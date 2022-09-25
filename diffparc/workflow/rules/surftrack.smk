rule convert_warpfield_template_to_indiv:
    input:
        warp=bids(
            root="work",
            datatype="anat",
            suffix="warp.nii.gz",
            from_="subject",
            to="{template}",
            **config["subj_wildcards"]
        ),
    output:
        warp=bids(
            root="work",
            datatype="surftrack",
            suffix="surfwarp.nii.gz",
            to_="subject",
            from_="{template}",
            **config["subj_wildcards"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -convert-warpfield -from-itk {input} -to-world {output}"


rule convert_affine_to_world:
    input:
        affine_itk=bids(
            root="work",
            datatype="anat",
            suffix="affine.txt",
            from_="subject",
            to="{template}",
            desc="itk",
            **config["subj_wildcards"]
        ),
    output:
        affine_world=bids(
            root="work",
            datatype="anat",
            suffix="affine.txt",
            from_="subject",
            to="{template}",
            desc="world",
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
        surf_gii="results/tpl-{template}/tpl-{template}_hemi-{hemi}_{seed}.surf.gii",
        warp=bids(
            root="work",
            datatype="surftrack",
            suffix="surfwarp.nii.gz",
            to_="subject",
            from_="{template}",
            **config["subj_wildcards"]
        ),
        affine_world=bids(
            root="work",
            datatype="anat",
            suffix="affine.txt",
            from_="subject",
            to="{template}",
            desc="world",
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
    shadow:
        "minimal"
    group:
        "subj"
    container:
        config["singularity"]["autotop"]
    shell:
        "wb_command -surface-apply-warpfield {input.surf_gii} {input.warp} transformed_with_warpfield.surf.gii && "
        "wb_command -surface-apply-affine transformed_with_warpfield.surf.gii {input.affine_world} {output.surf_warped}"


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
        targets=bids(
            root="results",
            **config["subj_wildcards"],
            space="individual",
            desc="{targets}",
            from_=config["template"],
            datatype="anat",
            suffix="dseg.mif"
        ),
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
    shell:
        "cp {input} {output} && "
        "wb_command -set-structure {output} OTHER"
