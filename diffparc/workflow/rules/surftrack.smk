rule convert_rigid_to_world:
    input:
        xfm_itk=bids(
            root="work",
            suffix="xfm.txt",
            hemi="{hemi}",
            from_="{template}",
            to="subj",
            desc="rigid",
            type_="itk",
            label="{seed}",
            datatype="morph",
            **config["subj_wildcards"]
        ),
    output:
        rigid_world=bids(
            root="work",
            suffix="xfm.txt",
            hemi="{hemi}",
            from_="{template}",
            to="subj",
            desc="rigid",
            type_="world",
            label="{seed}",
            datatype="surftrack",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -convert-affine -from-itk {input} -to-world {output} -inverse"


rule transform_template_surf_to_t1:
    """ transforms the template surface to the subject T1 """
    input:
        surf_gii=bids(
            root="work",
            hemi="{hemi}",
            **config["subj_wildcards"],
            desc="fluid",
            from_="{template}",
            datatype="morph",
            suffix="{seed}.surf.gii"
        ),
        rigid_world=bids(
            root="work",
            suffix="xfm.txt",
            hemi="{hemi}",
            from_="{template}",
            to="subj",
            desc="rigid",
            type_="world",
            label="{seed}",
            datatype="surftrack",
            **config["subj_wildcards"]
        ),
    output:
        surf_warped=bids(
            root="work",
            **config["subj_wildcards"],
            space="individual",
            hemi="{hemi}",
            from_="{template}",
            datatype="surftrack",
            suffix="{seed}.surf.gii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -surface-apply-affine {input.surf_gii} {input.rigid_world} {output.surf_warped}"


# for seeding, create a csv with vertex coords
rule create_surf_seed_csv:
    input:
        surf=bids(
            root="work",
            **config["subj_wildcards"],
            hemi="{hemi}",
            space="individual",
            from_="{template}",
            datatype="surftrack",
            suffix="{seed}.surf.gii"
        ),
    output:
        csv=bids(
            root="work",
            **config["subj_wildcards"],
            hemi="{hemi}",
            space="individual",
            from_="{template}",
            datatype="surftrack",
            label="{seed}",
            suffix="seeds.csv"
        ),
    group:
        "subj"
    script:
        "../scripts/surf_to_seed_csv.py"


rule track_from_vertices:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
    input:
        wm_fod=get_fod_for_tracking,
        dwi=bids(
            root="results",
            datatype="dwi",
            suffix="dwi.mif",
            **config["subj_wildcards"],
        ),
        mask=bids(
            root="results",
            datatype="dwi",
            suffix="mask.mif",
            **config["subj_wildcards"],
        ),
        csv=bids(
            root="work",
            **config["subj_wildcards"],
            space="individual",
            hemi="{hemi}",
            from_=config["template"],
            datatype="surftrack",
            label="{seed}",
            suffix="seeds.csv"
        ),
    params:
        radius="0.5",
        seedspervertex="{seedspervertex}",
    output:
        tck_dir=temp(
            directory(
                bids(
                    root="work",
                    datatype="surftrack",
                    hemi="{hemi}",
                    label="{seed}",
                    seedspervertex="{seedspervertex}",
                    suffix="vertextracts",
                    **config["subj_wildcards"],
                )
            )
        ),
    threads: 32
    resources:
        mem_mb=128000,
        time=1440,
    group:
        "subj"
    container:
        config["singularity"]["diffparc_deps"]
    shell:
        "mkdir -p {output.tck_dir} && "
        "parallel --bar --link --jobs {threads} "
        "tckgen -quiet -nthreads 0  -algorithm iFOD2 -mask {input.mask} "
        " {input.wm_fod} {output.tck_dir}/vertex_{{1}}.tck "
        " -seed_sphere {{2}},{params.radius} -seeds {params.seedspervertex} "
        " :::  `seq -w $(cat {input.csv} | wc -l)` ::: `cat {input.csv}` "

def get_dseg_targets(wildcards):

    if config['use_synthseg']:
        return bids(
            root="results",
            **config["subj_wildcards"],
            space="individual",
            desc="{targets}",
            from_='synthseg',
            datatype="anat",
            suffix="dseg.mif"
        ),

    else:
        return bids(
            root="results",
            **config["subj_wildcards"],
            space="individual",
            desc="{targets}",
            from_=config["template"],
            datatype="anat",
            suffix="dseg.mif"
        ),



rule connectivity_from_vertices:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
    input:
        tck_dir=bids(
            root="work",
            datatype="surftrack",
            hemi="{hemi}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="vertextracts",
            **config["subj_wildcards"],
        ),
        targets=get_dseg_targets
    output:
        conn_dir=temp(
            directory(
                bids(
                    root="work",
                    datatype="surftrack",
                    hemi="{hemi}",
                    desc="{targets}",
                    label="{seed}",
                    seedspervertex="{seedspervertex}",
                    suffix="vertexconn",
                    **config["subj_wildcards"],
                )
            )
        ),
    threads: 32
    resources:
        mem_mb=128000,
        time=1440,
    group:
        "subj"
    container:
        config["singularity"]["diffparc_deps"]
    shell:
        "mkdir -p {output.conn_dir} && "
        "parallel --eta --jobs {threads} "
        "tck2connectome -nthreads 0 -quiet {input.tck_dir}/vertex_{{1}}.tck {input.targets} {output.conn_dir}/conn_{{1}}.csv -vector"
        " ::: `ls {input.tck_dir} | grep -Po '(?<=vertex_)[0-9]+'`"


rule gen_vertex_conn_csv:
    input:
        conn_dir=bids(
            root="work",
            datatype="surftrack",
            desc="{targets}",
            hemi="{hemi}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="vertexconn",
            **config["subj_wildcards"],
        ),
    params:
        header_line=lambda wildcards: ",".join(
            config["targets"][wildcards.targets]["labels"]
        ),
    output:
        conn_csv=bids(
            root="work",
            datatype="surftrack",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="conn.csv",
            **config["subj_wildcards"],
        ),
    group:
        "subj"
    script:
        "../scripts/gather_csv_files.py"


rule conn_csv_to_metric:
    input:
        csv=bids(
            root="work",
            datatype="surftrack",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="conn.csv",
            **config["subj_wildcards"],
        ),
    output:
        gii_metric=temp(
            bids(
                root="work",
                datatype="surftrack",
                hemi="{hemi}",
                desc="{targets}",
                label="{seed}",
                seedspervertex="{seedspervertex}",
                suffix="nostructconn.shape.gii",
                **config["subj_wildcards"],
            )
        ),
    group:
        "subj"
    script:
        "../scripts/conn_csv_to_gifti_metric.py"


rule set_structure_conn_metric:
    input:
        gii_metric=bids(
            root="work",
            datatype="surftrack",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="nostructconn.shape.gii",
            **config["subj_wildcards"],
        ),
    params:
        structure=lambda wildcards: config["hemi_to_structure"][wildcards.hemi],
    output:
        gii_metric=bids(
            root="work",
            datatype="surftrack",
            hemi="{hemi}",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="conn.shape.gii",
            **config["subj_wildcards"],
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "cp {input} {output} && "
        "wb_command -set-structure {output} {params.structure}"


rule create_cifti_conn_dscalar:
    input:
        left_metric=bids(
            root="work",
            datatype="surftrack",
            hemi="L",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="conn.shape.gii",
            **config["subj_wildcards"],
        ),
        right_metric=bids(
            root="work",
            datatype="surftrack",
            hemi="R",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="conn.shape.gii",
            **config["subj_wildcards"],
        ),
    output:
        cifti_dscalar=bids(
            root="work",
            datatype="surftrack",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="conn.dscalar.nii",
            **config["subj_wildcards"],
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-create-dense-scalar {output} -left-metric {input.left_metric} -right-metric {input.right_metric}"


rule create_cifti_conn_dscalar_maxprob:
    input:
        cifti_dscalar=bids(
            root="work",
            datatype="surftrack",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="conn.dscalar.nii",
            **config["subj_wildcards"],
        ),
    output:
        cifti_dscalar=bids(
            root="work",
            datatype="surftrack",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="maxprob.dscalar.nii",
            **config["subj_wildcards"],
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-reduce {input} INDEXMAX {output}"


# need to then convert that into a label, then can use parcellate
rule create_cifti_maxprob_dlabel:
    input:
        cifti_dscalar=bids(
            root="work",
            datatype="surftrack",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="maxprob.dscalar.nii",
            **config["subj_wildcards"],
        ),
        label_list_txt=lambda wildcards: os.path.join(
            workflow.basedir,
            "..",
            config["targets"][wildcards.targets]["label_list_txt"],
        ),
    output:
        cifti_dlabel=bids(
            root="work",
            datatype="surftrack",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="maxprob.dlabel.nii",
            **config["subj_wildcards"],
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-label-import {input.cifti_dscalar} {input.label_list_txt} {output.cifti_dlabel}"


rule parcellate_inout_displacement:
    input:
        cifti_dscalar=bids(
            root="work",
            **config["subj_wildcards"],
            desc="inout",
            from_="{template}",
            datatype="morph",
            label="{seed}",
            suffix="surfdisp.dscalar.nii"
        ),
        cifti_dlabel=bids(
            root="work",
            datatype="surftrack",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="maxprob.dlabel.nii",
            **config["subj_wildcards"],
        ),
    output:
        cifti_pscalar=bids(
            root="work",
            **config["subj_wildcards"],
            desc="inout",
            from_="{template}",
            datatype="morph",
            label="{seed}",
            parcel="{targets}",
            seedspervertex="{seedspervertex}",
            suffix="surfdisp.pscalar.nii"
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -cifti-parcellate {input.cifti_dscalar} {input.cifti_dlabel} COLUMN "
        " {output}"
