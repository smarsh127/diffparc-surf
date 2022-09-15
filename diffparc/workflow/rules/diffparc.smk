rule create_upsampled_cropped_seed_ref:
    input:
        ref=bids(
            root="results",
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **config["subj_wildcards"]
        ),
    params:
        resample_res=config["resample_seed_res"],
    output:
        ref=bids(
            root="results",
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res="upsampled",
            datatype="dwi",
            **config["subj_wildcards"]
        ),
    shell:
        "c3d {input} -resample-mm {params.resample_res}  -o {output}"


rule transform_seed_to_subject:
    input:
        seed=lambda wildcards: os.path.join(
            workflow.basedir, "..", config["seeds"][wildcards.seed]["template_probseg"]
        ),
        ref=bids(
            root="results",
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res="upsampled",
            datatype="dwi",
            **config["subj_wildcards"]
        ),
        inv_warp=bids(
            root="work",
            datatype="anat",
            suffix="invwarp.nii.gz",
            from_="subject",
            to="{template}",
            **config["subj_wildcards"]
        ),
        affine_xfm_itk=bids(
            root="work",
            datatype="anat",
            suffix="affine.txt",
            from_="subject",
            to="{template}",
            desc="itk",
            **config["subj_wildcards"]
        ),
    output:
        seed=bids(
            root="results",
            **config["subj_wildcards"],
            space="individual",
            label="{seed}",
            from_="{template}",
            datatype="anat",
            suffix="probseg.nii.gz"
        ),
    envmodules:
        "ants",
    container:
        config["singularity"]["ants"]
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
            space="individual",
            label="{seed}",
            from_="{template}",
            suffix="transformseedtosubject.log"
        ),
    group:
        "participant1"
    threads: 8
    shell:
        "antsApplyTransforms -d 3 --interpolation Linear -i {input.seed} -o {output} -r {input.ref} -t [{input.affine_xfm_itk},1] {input.inv_warp}  &> {log}"


rule transform_targets_to_subject:
    input:
        targets=lambda wildcards: os.path.join(
            workflow.basedir,
            "..",
            config["targets"][wildcards.targets]["template_dseg"],
        ),
        ref=bids(
            root="results",
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **config["subj_wildcards"]
        ),
        inv_warp=bids(
            root="work",
            datatype="anat",
            suffix="invwarp.nii.gz",
            from_="subject",
            to="{template}",
            **config["subj_wildcards"]
        ),
        affine_xfm_itk=bids(
            root="work",
            datatype="anat",
            suffix="affine.txt",
            from_="subject",
            to="{template}",
            desc="itk",
            **config["subj_wildcards"]
        ),
    output:
        targets=bids(
            root="results",
            **config["subj_wildcards"],
            space="individual",
            desc="{targets}",
            from_="{template}",
            datatype="anat",
            suffix="dseg.nii.gz"
        ),
    envmodules:
        "ants",
    container:
        config["singularity"]["ants"]
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
            desc="{targets}",
            from_="{template}",
            suffix="transformtargetstosubjects.log"
        ),
    group:
        "participant1"
    threads: 8
    shell:
        "antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.targets} -o {output} -r {input.ref} -t [{input.affine_xfm_itk},1] {input.inv_warp}  &> {log}"


rule binarize_trim_subject_seed:
    input:
        seed_res=bids(
            root="results",
            **config["subj_wildcards"],
            space="individual",
            label="{seed}",
            from_="{template}",
            datatype="anat",
            suffix="probseg.nii.gz"
        ),
    params:
        threshold=lambda wildcards: config["seeds"][wildcards.seed]["probseg_threshold"],
    output:
        seed_thr=bids(
            root="results",
            **config["subj_wildcards"],
            space="individual",
            label="{seed}",
            from_="{template}",
            datatype="anat",
            suffix="mask.nii.gz"
        ),
    container:
        config["singularity"]["prepdwi"]
    log:
        bids(
            root="logs",
            **config["subj_wildcards"],
            label="{seed}",
            from_="{template}",
            suffix="binarizetrimsubjectseed.log"
        ),
    container:
        config["singularity"]["prepdwi"]
    group:
        "participant1"
    shell:
        "c3d {input.seed_res} -threshold 0.5 inf 1 0 -trim 0vox -type uchar -o {output} &> {log}"


def get_fod_for_tracking(wildcards):
    if config["fod_algorithm"] == "csd":
        return (
            bids(
                root="results",
                datatype="dwi",
                alg="csd",
                desc="wm",
                suffix="fod.mif",
                **config["subj_wildcards"],
            ),
        )
    elif config["fod_algorithm"] == "msmt_csd":
        return (
            bids(
                root="results",
                datatype="dwi",
                desc="wmnorm",
                alg="msmt",
                suffix="fod.mif",
                **config["subj_wildcards"],
            ),
        )


rule create_voxel_seed_images:
    input:
        seed=bids(
            root="results",
            **config["subj_wildcards"],
            space="individual",
            label="{seed}",
            from_=config["template"],
            datatype='anat',
            suffix="mask.nii.gz"
        ),
    output:
        voxseeds_dir=temp(
            directory(
                bids(
                    root="results",
                    **config["subj_wildcards"],
                    space="individual",
                    label="{seed}",
                    from_=config["template"],
                    datatype="dwi",
                    suffix="voxseeds"
                )
            )
        ),
    script:
        "../scripts/create_voxel_seed_images.py"


rule track_from_voxels:
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
        vox_seeds_dir=bids(
            root="results",
            **config["subj_wildcards"],
            space="individual",
            label="{seed}",
            datatype="dwi",
            from_=config["template"],
            suffix="voxseeds"
        ),
    params:
        seeds_per_voxel="{seedpervox}",
    output:
        tck_dir=temp(directory(
            bids(
                root="results",
                datatype="dwi",
                label="{seed}",
                seedpervox="{seedpervox}",
                suffix="voxtracts",
                **config["subj_wildcards"],
            ))
        ),
    threads: 32
    resources:
        mem_mb=128000,
        time=1440,
    group:
        "subj2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mkdir -p {output.tck_dir} && "
        "parallel --bar --jobs {threads} "
        "tckgen -quiet -nthreads 0  -algorithm iFOD2 -mask {input.mask} "
        " {input.wm_fod} {output.tck_dir}/vox_{{1}}.tck "
        " -seed_random_per_voxel {input.vox_seeds_dir}/seed_{{1}}.nii {params.seeds_per_voxel} "
        " ::: `ls {input.vox_seeds_dir} | grep -Po '(?<=seed_)[0-9]+'`"


rule connectivity_from_voxels:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
    input:
        tck_dir=bids(
            root="results",
            datatype="dwi",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="voxtracts",
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
                    root="results",
                    datatype="dwi",
                    desc="{targets}",
                    label="{seed}",
                    seedpervox="{seedpervox}",
                    suffix="voxconn",
                    **config["subj_wildcards"],
                )
            )
        ),
    threads: 32
    resources:
        mem_mb=128000,
        time=1440,
    group:
        "subj2"
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mkdir -p {output.conn_dir} && "
        "parallel --eta --jobs {threads} "
        "tck2connectome -nthreads 0 -quiet {input.tck_dir}/vox_{{1}}.tck {input.targets} {output.conn_dir}/conn_{{1}}.csv -vector"
        " ::: `ls {input.tck_dir} | grep -Po '(?<=vox_)[0-9]+'`"


rule dseg_nii2mif:
    input:
        "{file}_dseg.nii.gz",
    output:
        "{file}_dseg.mif",
    container:
        config["singularity"]["mrtrix"]
    shell:
        "mrconvert {input} {output} -nthreads {threads}"


rule gen_conn_csv:
    input:
        conn_dir=bids(
            root="results",
            datatype="dwi",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="voxconn",
            **config["subj_wildcards"],
        ),
    params:
        header_line=lambda wildcards: ",".join(
            config["targets"][wildcards.targets]["labels"]
        ),
    output:
        conn_csv=bids(
            root="results",
            datatype="dwi",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="conn.csv",
            **config["subj_wildcards"],
        ),
    group:
        "subj2"
    container:
        config["singularity"]["mrtrix"]
    script:
        "../scripts/gather_csv_files.py"


#    shell:
#        'echo {params.header_line} > {output.conn_csv} && '
#        'for f in `ls {input.conn_dir}/*.csv`; do tail -n 1 $f; done >> {output.conn_csv}'


rule conn_csv_to_image:
    input:
        conn_csv=bids(
            root="results",
            datatype="dwi",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="conn.csv",
            **config["subj_wildcards"],
        ),
        seed_nii=bids(
            root="results",
            **config["subj_wildcards"],
            space="individual",
            label="{seed}",
            from_=config["template"],
            datatype="anat",
            suffix="mask.nii.gz"
        ),
    output:
        conn_nii=bids(
            root="results",
            datatype="dwi",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="conn.nii.gz",
            **config["subj_wildcards"],
        ),
    script:
        "../scripts/conn_csv_to_image.py"


rule transform_conn_to_template:
    input:
        conn_nii=bids(
            root="results",
            datatype="dwi",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="conn.nii.gz",
            **config["subj_wildcards"],
        ),
        warp=bids(
            root="work",
            datatype="anat",
            suffix="warp.nii.gz",
            from_="subject",
            to="{template}",
            **config["subj_wildcards"]
        ),
        affine_xfm_itk=bids(
            root="work",
            datatype="anat",
            suffix="affine.txt",
            from_="subject",
            to="{template}",
            desc="itk",
            **config["subj_wildcards"]
        ),
        ref=os.path.join(workflow.basedir, "..", config["template_t1w"]),
    output:
        conn_nii=bids(
            root="results",
            datatype="dwi",
            space="{template}",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="conn.nii.gz",
            **config["subj_wildcards"],
        ),
    container:
        config["singularity"]["ants"]
    threads: 8
    resources:
        mem_mb=8000,
    log:
        bids(
            root="logs",
            space="{template}",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="transformconntotemplate.log",
            **config["subj_wildcards"],
        ),
    group:
        "participant1"
    shell:
        #using nearestneighbor to avoid bluring with background -- background set as -1
        "antsApplyTransforms -d 3 -e 3  --interpolation NearestNeighbor -i {input.conn_nii}  -o {output.conn_nii}  -r {input.ref} -t {input.warp} -t {input.affine_xfm_itk} &> {log} "


rule maxprob_conn:
    """ generate maxprob connectivity, adding outside striatum, and inside striatum (at a particular "streamline count" threshold) to identify unlabelled regions """
    input:
        conn_nii=bids(
            root="results",
            datatype="dwi",
            space="{template}",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
            suffix="conn.nii.gz",
            **config["subj_wildcards"],
        ),
    output:
        conn_nii=bids(
            root="results",
            datatype="dwi",
            space="{template}",
            desc="{targets}",
            label="{seed}",
            seedpervox="{seedpervox}",
            segtype="maxprob",
            suffix="dseg.nii.gz",
            **config["subj_wildcards"],
        ),
    shell:
        "c4d {input} -slice w 0:-1 -vote -o {output} "


# isosurf template probseg
rule gen_template_surface:
    input:
        nii=lambda wildcards: config["seeds"][wildcards.seed]["template_probseg"],
    params:
        threshold=lambda wildcards: config["seeds"][wildcards.seed]["probseg_threshold"],
    output:
        vtk="results/tpl-{template}/tpl-{template}_label-{seed}_isosurf.vtk",
    script:
        "../scripts/gen_isosurface.py"


"""
rule save_connmap_template_npz:
    input:
        mask = bids(root='results',template='{template}',label='{seed}',suffix='mask.nii.gz'),
        probtrack_dir = bids(root='results',**config['subj_wildcards'],label='{seed}',space='{template}',suffix='probtrack')
    params:
        connmap_3d = lambda wildcards, input: expand(bids(root=input.probtrack_dir,include_subject_dir=False,subject=wildcards.subject,prefix='seeds_to',label='{target}',suffix='mask.nii.gz'), target=targets),
    output:
        connmap_npz = bids(root='results',**config['subj_wildcards'],label='{seed}',space='{template}',suffix='connMap.npz')
    log: 'logs/save_connmap_to_template_npz/sub-{subject}_{seed}_{template}.log'
    group: 'participant1'
    conda: '../envs/sklearn.yml'
    script: '../scripts/save_connmap_template_npz.py'

#space-{template}
rule gather_connmap_group:
    input:
        connmap_npz = expand(bids(root='results',**config['subj_wildcards'],label='{seed}',space='{template}',suffix='connMap.npz'), subject=subjects,allow_missing=True)
    output:
        connmap_group_npz = bids(root='results',template='{template}',desc='concat',label='{seed}',from_='group',suffix='connMap.npz')
    log: 'logs/gather_connmap_group/{seed}_{template}.log'
    conda: '../envs/sklearn.yml'
    group: 'group1'
    script: '../scripts/gather_connmap_group.py'
     

#space-{template},  dseg
rule spectral_clustering:
    input:
        connmap_group_npz = bids(root='results',template='{template}',desc='concat',label='{seed}',from_='group',suffix='connMap.npz')
    params:
        max_k = config['max_k']
    output:
        cluster_k = expand(bids(root='results',template='{template}',label='{seed}',from_='group',method='spectralcosine',k='{k}',suffix='dseg.nii.gz'),k=range(2,config['max_k']+1),allow_missing=True)
    log: 'logs/spectral_clustering/{seed}_{template}.log'
    conda: '../envs/sklearn.yml'
    group: 'group1'
    script: '../scripts/spectral_clustering.py'
        
 
  
"""
