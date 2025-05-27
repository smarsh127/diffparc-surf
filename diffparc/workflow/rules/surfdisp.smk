# to replicate earlier surface morphometry pipeline, we want to:
# 1. rigid reg target probseg to template probseg
# 2. run lddmm/greedy/ants to get mapping from template to target
# 3. apply transform to the template surface
# 4. calculate displacement based on target minus template
# 5. normalize by the average displacement vector in local neighborhood
#   - (was 10mm radius cube before, now we use 8mm fwhm surf smoothing - can tune this hyperparameter)
# 6. calculate inward/outward displacement using dot product with surface normal


rule rigid_shape_reg:
    """ rigidly register subj to template shape using moment-invariants """
    input:
        template=lambda wildcards: os.path.join(
            workflow.basedir, "..", config["seeds"][wildcards.seed]["template_probseg"]
        ),
        target=get_subject_seed_probseg,
    params:
        general_opts="-d 3",
        rigid_opts="-m SSD -moments 2 -det 1",
    output:
        xfm_ras=bids(
            root=root,
            suffix="xfm.txt",
            hemi="{hemi}",
            from_=config["template"],
            to="subj",
            desc="rigid",
            type_="ras",
            label="{seed}",
            datatype="warps",
            **subj_wildcards
        ),
        xfm_itk=bids(
            root=root,
            suffix="xfm.txt",
            hemi="{hemi}",
            from_=config["template"],
            to="subj",
            desc="rigid",
            type_="itk",
            label="{seed}",
            datatype="warps",
            **subj_wildcards
        ),
        warped_target=bids(
            root=root,
            **subj_wildcards,
            desc="rigid",
            hemi="{hemi}",
            label="{seed}",
            space=config["template"],
            datatype="morph",
            suffix="probseg.nii.gz"
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    threads: 8
    shell:
        "greedy -threads {threads} {params.general_opts} {params.rigid_opts}  -i {input.template} {input.target} -o {output.xfm_ras} &&  "
        "greedy -threads {threads} {params.general_opts} -rf {input.template} -rm {input.target} {output.warped_target}  -r {output.xfm_ras} && "
        "c3d_affine_tool {output.xfm_ras} -oitk {output.xfm_itk}"


rule fluid_shape_reg:
    """ fluid registration to get the displacements """
    input:
        template=lambda wildcards: os.path.join(
            workflow.basedir, "..", config["seeds"][wildcards.seed]["template_probseg"]
        ),
        target=bids(
            root=root,
            **subj_wildcards,
            desc="rigid",
            hemi="{hemi}",
            label="{seed}",
            space=config["template"],
            datatype="morph",
            suffix="probseg.nii.gz"
        ),
    output:
        warp=bids(
            root=root,
            datatype="warps",
            suffix="warp.nii.gz",
            hemi="{hemi}",
            from_="subject",
            to=config["template"],
            label="{seed}",
            **subj_wildcards
        ),
        warped=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            desc="fluid",
            label="{seed}",
            space=config["template"],
            datatype="morph",
            suffix="probseg.nii.gz"
        ),
    threads: 8
    resources:
        mem_mb=16000,  # right now these are on the high-end -- could implement benchmark rules to do this at some point..
        time=60,  # 1 hrs
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        "greedy -d 3 -threads {threads} -m SSD -i {input.template} {input.target} -o {output.warp} -n 100x50x10 && "
        "greedy -d 3 -threads {threads} -rf {input.template} -rm {input.target} {output.warped} -r {output.warp} "


rule convert_warpfield:
    input:
        warp=bids(
            root=root,
            datatype="warps",
            suffix="warp.nii.gz",
            hemi="{hemi}",
            from_="subject",
            to=config["template"],
            label="{seed}",
            **subj_wildcards
        ),
    output:
        warp=bids(
            root=root,
            datatype="warps",
            hemi="{hemi}",
            suffix="surfwarp.nii.gz",
            to_="subject",
            from_=config["template"],
            label="{seed}",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "wb_command -convert-warpfield -from-itk {input} -to-world {output}"


rule deform_surface:
    input:
        surf_gii=rules.set_surface_structure.output,
        warp=bids(
            root=root,
            datatype="warps",
            hemi="{hemi}",
            suffix="surfwarp.nii.gz",
            to_="subject",
            from_=config["template"],
            label="{seed}",
            **subj_wildcards
        ),
    output:
        surf_warped=bids(
            root=root,
            hemi="{hemi}",
            **subj_wildcards,
            desc="fluid",
            datatype="morph",
            suffix="{seed}.surf.gii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "wb_command -surface-apply-warpfield {input.surf_gii} {input.warp} {output.surf_warped}"


rule compute_displacement_metrics:
    input:
        ref_surf=rules.set_surface_structure.output,
        comp_surf=bids(
            root=root,
            hemi="{hemi}",
            **subj_wildcards,
            desc="fluid",
            datatype="morph",
            suffix="{seed}.surf.gii"
        ),
    output:
        scalar=temp(
            bids(
                root=root,
                hemi="{hemi}",
                **subj_wildcards,
                datatype="morph",
                label="{seed}",
                suffix="scalardisp.shape.gii"
            )
        ),
        vector=temp(
            bids(
                root=root,
                hemi="{hemi}",
                **subj_wildcards,
                datatype="morph",
                label="{seed}",
                suffix="vectordisp.shape.gii"
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "wb_command  -surface-to-surface-3d-distance {input.comp_surf} {input.ref_surf} {output.scalar} -vectors {output.vector}"


rule calc_template_surf_normals:
    input:
        ref_surf=rules.set_surface_structure.output,
    output:
        normals=get_template_prefix(
            root=root, subj_wildcards=subj_wildcards, template=config["template"]
        )
        + "_hemi-{hemi}_label-{seed}_normals.shape.gii",
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "wb_command -surface-normals {input} {output}"


rule smooth_displacement_vectors:
    input:
        ref_surf=rules.set_surface_structure.output,
        vector=bids(
            root=root,
            hemi="{hemi}",
            **subj_wildcards,
            datatype="morph",
            label="{seed}",
            suffix="vectordisp.shape.gii"
        ),
    params:
        fwhm=config["surfdisp_normalization_fwhm_mm"],
    output:
        vector=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            desc="smoothed",
            datatype="morph",
            label="{seed}",
            suffix="vectordisp.shape.gii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "wb_command -metric-smoothing {input.ref_surf} {input.vector} {params.fwhm} {output.vector} -fwhm"


rule normalize_displacement_by_smoothed:
    """subtracts the displacement by the smoothed displacement to provide a spatially-local
    displacement estimate"""
    input:
        vector=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            datatype="morph",
            label="{seed}",
            suffix="vectordisp.shape.gii"
        ),
        smoothed=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            desc="smoothed",
            datatype="morph",
            label="{seed}",
            suffix="vectordisp.shape.gii"
        ),
    output:
        normalized=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            desc="normalized",
            datatype="morph",
            label="{seed}",
            suffix="vectordisp.shape.gii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "wb_command -metric-math 'vec - vecsmooth' {output.normalized}"
        " -var vec {input.vector} "
        " -var vecsmooth {input.smoothed}  "


rule calc_inout_displacement:
    """uses dot-product with surface normal to get inward/outward displacement"""
    input:
        vec=bids(
            root=root,
            hemi="{hemi}",
            **subj_wildcards,
            desc="normalized",
            datatype="morph",
            label="{seed}",
            suffix="vectordisp.shape.gii"
        ),
        norm=rules.calc_template_surf_normals.output,
    output:
        inout=bids(
            root=root,
            **subj_wildcards,
            hemi="{hemi}",
            datatype="surf",
            label="{seed}",
            suffix="inout.shape.gii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    shell:
        "wb_command -metric-math '(v1*n1)+(v2*n2)+(v3*n3)' {output.inout}"
        " -var v1 {input.vec} -column 1  "
        " -var v2 {input.vec} -column 2  "
        " -var v3 {input.vec} -column 3  "
        " -var n1 {input.norm} -column 1  "
        " -var n2 {input.norm} -column 2  "
        " -var n3 {input.norm} -column 3  "
