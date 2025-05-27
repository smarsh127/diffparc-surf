# if we skip dwi preproc, then we are loading dwi preproc in T1w space:

if config["in_prepdwi_dir"]:

    rule import_prepdwi_dwi:
        input:
            dwi_nii=bids(
                root=os.path.join(config["in_prepdwi_dir"] + "prepdwi"),
                suffix="dwi_space-T1w_preproc.nii.gz",
                datatype="dwi",
                **subj_wildcards
            ),
            dwi_bval=bids(
                root=os.path.join(config["in_prepdwi_dir"] + "prepdwi"),
                suffix="dwi_space-T1w_preproc.bval",
                datatype="dwi",
                **subj_wildcards
            ),
            dwi_bvec=bids(
                root=os.path.join(config["in_prepdwi_dir"] + "prepdwi"),
                suffix="dwi_space-T1w_preproc.bvec",
                datatype="dwi",
                **subj_wildcards
            ),
            brainmask_nii=bids(
                root=os.path.join(config["in_prepdwi_dir"] + "prepdwi"),
                suffix="dwi_space-T1w_brainmask.nii.gz",
                datatype="dwi",
                **subj_wildcards
            ),
        output:
            dwi=expand(
                bids(
                    root=root,
                    suffix="dwi{ext}",
                    desc="preproc",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                ),
                ext=[".nii.gz", ".bval", ".bvec"],
                allow_missing=True,
            ),
            brainmask=bids(
                root=root,
                suffix="mask.nii.gz",
                desc="brain",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
        group:
            "subj"
        run:
            for in_file, out_file in zip(input, output):
                shell("cp -v {in_file} {out_file}")




elif config["in_snakedwi_dir"]:

    rule import_snakedwi_dwi:
        input:
            dwi_nii=bids(
                root=config["in_snakedwi_dir"],
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                desc="eddy",
                suffix="dwi.nii.gz",
                datatype="dwi",
                **subj_wildcards
            ),
            dwi_bval=bids(
                root=config["in_snakedwi_dir"],
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                desc="eddy",
                suffix="dwi.bval",
                datatype="dwi",
                **subj_wildcards
            ),
            dwi_bvec=bids(
                root=config["in_snakedwi_dir"],
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                desc="eddy",
                suffix="dwi.bvec",
                datatype="dwi",
                **subj_wildcards
            ),
            brainmask_nii=bids(
                root=config["in_snakedwi_dir"],
                suffix="mask.nii.gz",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                desc="brain",
                datatype="dwi",
                **subj_wildcards
            ),
        output:
            dwi=expand(
                bids(
                    root=root,
                    suffix="dwi{ext}",
                    desc="preproc",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                ),
                ext=[".nii.gz", ".bval", ".bvec"],
                allow_missing=True,
            ),
            brainmask=bids(
                root=root,
                suffix="mask.nii.gz",
                desc="brain",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
        group:
            "subj"
        run:
            for in_file, out_file in zip(input, output):
                shell("cp -v {in_file} {out_file}")




else:

    rule resample_dwi_to_t1w:
        input:
            ref=bids(
                root=root,
                suffix="avgb0.nii.gz",
                space="T1w",
                desc="preproc",
                proc="crop",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
            dwi=bids(
                root=root,
                suffix="dwi.nii.gz",
                desc="preproc",
                datatype="dwi",
                **subj_wildcards
            ),
            xfm_itk=bids(
                root=root,
                suffix="xfm.txt",
                from_="dwi",
                to="T1w",
                type_="itk",
                datatype="dwi",
                **subj_wildcards
            ),
        params:
            interpolation="Linear",
        output:
            dwi=bids(
                root=root,
                suffix="dwi.nii.gz",
                desc="preproc",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
        container:
            config["singularity"]["diffparc"]
        resources:
            mem_mb=32000,  #-- this is going to be dependent on size of image.. 
        group:
            "subj"
        shell:
            "antsApplyTransforms -d 3 --input-image-type 3 --input {input.dwi} --reference-image {input.ref} --transform {input.xfm_itk} --interpolation {params.interpolation} --output {output.dwi} --verbose "

    rule resample_brainmask_to_t1w:
        input:
            ref=bids(
                root=root,
                suffix="avgb0.nii.gz",
                space="T1w",
                desc="preproc",
                proc="crop",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
            brainmask=get_dwi_mask(),
            xfm_itk=bids(
                root=root,
                suffix="xfm.txt",
                from_="dwi",
                to="T1w",
                type_="itk",
                datatype="dwi",
                **subj_wildcards
            ),
        params:
            interpolation="NearestNeighbor",
        output:
            brainmask=bids(
                root=root,
                suffix="mask.nii.gz",
                desc="brain",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
        container:
            config["singularity"]["diffparc"]
        resources:
            mem_mb=32000,  #-- this is going to be dependent on size of image.. 
        group:
            "subj"
        shell:
            "antsApplyTransforms -d 3 --input-image-type 0 --input {input.brainmask} --reference-image {input.ref} --transform {input.xfm_itk} --interpolation {params.interpolation} --output {output.brainmask} --verbose"

    rule rotate_bvecs_to_t1w:
        input:
            bvecs=bids(
                root=root,
                suffix="dwi.bvec",
                desc="preproc",
                datatype="dwi",
                **subj_wildcards
            ),
            xfm_fsl=bids(
                root=root,
                suffix="xfm.txt",
                from_="dwi",
                to="T1w",
                type_="fsl",
                datatype="dwi",
                **subj_wildcards
            ),
            bvals=bids(
                root=root,
                suffix="dwi.bval",
                desc="preproc",
                datatype="dwi",
                **subj_wildcards
            ),
        params:
            script=os.path.join(workflow.basedir, "scripts/rotate_bvecs.sh"),
        output:
            bvecs=bids(
                root=root,
                suffix="dwi.bvec",
                desc="preproc",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
            bvals=bids(
                root=root,
                suffix="dwi.bval",
                desc="preproc",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            ),
        container:
            config["singularity"]["diffparc"]  #fsl
        group:
            "subj"
        shell:
            "chmod a+x {params.script} && "
            "{params.script} {input.bvecs} {input.xfm_fsl} {output.bvecs} && "
            "cp -v {input.bvals} {output.bvals}"


rule reg_dwi_to_t1:
    input:
        t1w=bids(
            root=root,
            suffix="T1w.nii.gz",
            desc="preproc",
            datatype="anat",
            **subj_wildcards
        ),
        avgb0=bids(
            root=root,
            suffix="b0.nii.gz",
            desc="preproc",
            datatype="dwi",
            **subj_wildcards
        ),
    params:
        general_opts="-d 3",
        rigid_opts="-m NMI -a -dof 6 -ia-identity",
    output:
        warped_avgb0=bids(
            root=root,
            suffix="avgb0.nii.gz",
            space="T1w",
            desc="preproc",
            datatype="dwi",
            **subj_wildcards
        ),
        xfm_ras=bids(
            root=root,
            suffix="xfm.txt",
            from_="dwi",
            to="T1w",
            type_="ras",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    log:
        bids(root="logs", suffix="reg_b0_to_t1.txt", datatype="dwi", **subj_wildcards),
    threads: 8
    shell:
        "greedy -threads {threads} {params.general_opts} {params.rigid_opts}  -i {input.t1w} {input.avgb0} -o {output.xfm_ras}  &> {log}  && "
        "greedy -threads {threads} {params.general_opts} -rf {input.t1w} -rm {input.avgb0} {output.warped_avgb0}  -r {output.xfm_ras} &>> {log}"


rule qc_reg_dwi_t1:
    input:
        ref=bids(
            root=root,
            suffix="avgb0.nii.gz",
            space="T1w",
            desc="preproc",
            datatype="dwi",
            **subj_wildcards
        ),
        flo=bids(
            root=root,
            suffix="T1w.nii.gz",
            desc="preproc",
            datatype="anat",
            **subj_wildcards
        ),
    output:
        png=report(
            bids(
                root="qc", suffix="reg.png", **subj_wildcards, from_="dwiref", to="T1w"
            ),
            caption="../report/reg_dwi_t1.rst",
            category="B0 T1w registration",
        ),
        html=bids(
            root="qc", suffix="reg.html", from_="dwiref", to="T1w", **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    script:
        "../scripts/vis_regqc.py"


rule convert_xfm_ras2itk:
    input:
        xfm_ras=bids(
            root=root,
            suffix="xfm.txt",
            from_="dwi",
            to="T1w",
            type_="ras",
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        xfm_itk=bids(
            root=root,
            suffix="xfm.txt",
            from_="dwi",
            to="T1w",
            type_="itk",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        "c3d_affine_tool {input.xfm_ras}  -oitk {output.xfm_itk}"


rule convert_xfm_ras2fsl:
    input:
        t1w=bids(
            root=root,
            suffix="T1w.nii.gz",
            desc="preproc",
            datatype="anat",
            **subj_wildcards
        ),
        avgb0=bids(
            root=root,
            suffix="b0.nii.gz",
            desc="preproc",
            datatype="dwi",
            **subj_wildcards
        ),
        xfm_ras=bids(
            root=root,
            suffix="xfm.txt",
            from_="dwi",
            to="T1w",
            type_="ras",
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        xfm_fsl=bids(
            root=root,
            suffix="xfm.txt",
            from_="dwi",
            to="T1w",
            type_="fsl",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        "c3d_affine_tool {input.xfm_ras} -ref {input.t1w} -src {input.avgb0} -ras2fsl -o {output.xfm_fsl}"
# tight crop around b0 after rotating into T1w space



rule create_cropped_ref:
    input:
        warped_avgb0=bids(
            root=root,
            suffix="avgb0.nii.gz",
            space="T1w",
            desc="preproc",
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        cropped_avgb0=bids(
            root=root,
            suffix="avgb0.nii.gz",
            space="T1w",
            desc="preproc",
            proc="crop",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        "c3d {input} -trim 0vox {output}"


# for later resampling..
rule write_nii_resolution_to_txt:
    input:
        "{prefix}.nii.gz",
    output:
        "{prefix}.resolution_mm.txt",
    group:
        "subj"
    container:
        config["singularity"]["diffparc"]
    script:
        "../scripts/write_nii_resolution_to_txt.py"


# rules for creating reference image for each resampling scheme -- only the rules that are required will be run..
rule create_cropped_ref_t1_resolution:
    input:
        cropped_avgb0=bids(
            root=root,
            suffix="avgb0.nii.gz",
            space="T1w",
            desc="preproc",
            proc="crop",
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        avgb0_crop_resample=bids(
            root=root,
            suffix="avgb0.nii.gz",
            space="T1w",
            desc="preproc",
            proc="crop",
            res="T1w",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    shell:
        "cp {input} {output}"


rule create_cropped_ref_dwi_resolution:
    input:
        cropped=bids(
            root=root,
            suffix="avgb0.nii.gz",
            space="T1w",
            desc="preproc",
            proc="crop",
            datatype="dwi",
            **subj_wildcards
        ),
        res_txt_orig=bids(
            root=root,
            suffix="b0.resolution_mm.txt",
            desc="preproc",
            datatype="dwi",
            **subj_wildcards
        ),
    output:
        resampled=bids(
            root=root,
            suffix="avgb0.nii.gz",
            space="T1w",
            desc="preproc",
            proc="crop",
            res="orig",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        "c3d {input.cropped} -resample-mm `cat {input.res_txt_orig}` {output}"


rule create_cropped_ref_custom_resolution:
    input:
        cropped=bids(
            root=root,
            suffix="avgb0.nii.gz",
            space="T1w",
            desc="preproc",
            proc="crop",
            datatype="dwi",
            **subj_wildcards
        ),
    params:
        resolution="x".join(
            [str(vox) for vox in config["resample_dwi"]["resample_mm"]]
        )
        + "mm",
    output:
        resampled=bids(
            root=root,
            suffix="avgb0.nii.gz",
            space="T1w",
            desc="preproc",
            proc="crop",
            res="custom",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    shell:
        "c3d {input} -resample-mm {params.resolution} {output}"


# dti fitting on dwi in t1w space


rule dtifit_resampled_t1w:
    input:
        dwi=bids(
            root=root,
            suffix="dwi.nii.gz",
            desc="preproc",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
        bvals=bids(
            root=root,
            suffix="dwi.bval",
            desc="preproc",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
        bvecs=bids(
            root=root,
            suffix="dwi.bvec",
            desc="preproc",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
        brainmask=bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **subj_wildcards
        ),
    params:
        out_basename=lambda wildcards, output: os.path.join(output.out_folder, "dti"),
    output:
        out_folder=directory(
            bids(
                root=root,
                suffix="dtifit",
                desc="preproc",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **subj_wildcards
            )
        ),
        out_fa=os.path.join(
            directory(
                bids(
                    root=root,
                    suffix="dtifit",
                    desc="preproc",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **subj_wildcards
                )
            ),
            "dti_FA.nii.gz",
        ),
    container:
        config["singularity"]["diffparc"]  #fsl
    group:
        "subj"
    shell:
        "mkdir -p {output.out_folder} && "
        "dtifit --data={input.dwi} --bvecs={input.bvecs} --bvals={input.bvals} --mask={input.brainmask} --out={params.out_basename}"
