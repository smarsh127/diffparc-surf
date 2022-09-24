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
        'wb_command -convert-affine -from-itk {input} -to-world {output} -inverse'

rule transform_template_surf_to_t1:
    """ transforms the template surface to the subject T1 """
    input:
        surf_gii="results/tpl-{template}/tpl-{template}_{seed}.surf.gii",
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
            space='individual',
            from_="{template}",
            datatype="surftrack",
            suffix="{seed}.surf.gii"
        ),
    shadow: 'minimal'
    group:
        "subj"
    container:
        config["singularity"]["autotop"]

    shell:
        'wb_command -surface-apply-warpfield {input.surf_gii} {input.warp} transformed_with_warpfield.surf.gii && '
        'wb_command -surface-apply-affine transformed_with_warpfield.surf.gii {input.affine_world} {output.surf_warped}' 



