
rule transform_seed_to_subject:
    input: 
        seed = lambda wildcards: config['seeds'][wildcards.seed]['template_probseg'],
        ref= bids(
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
        seed = bids(root='results',subject='{subject}',space='individual',label='{seed}',from_='{template}',suffix='probseg.nii.gz'),
    envmodules: 'ants'
    container: config['singularity']['ants']
    log: 'logs/transform_to_subject/{template}_sub-{subject}_{seed}.log'
    group: 'participant1'
    threads: 8
    shell:
        'antsApplyTransforms -d 3 --interpolation Linear -i {input.seed} -o {output} -r {input.ref} -t [{input.affine_xfm_itk},1] {input.inv_warp}  &> {log}'


rule transform_targets_to_subject:
    input: 
        targets = lambda wildcards: config['targets'][wildcards.targets]['template_dseg'],
        ref= bids(
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
        targets = bids(root='results',subject='{subject}',space='individual',desc='{targets}',from_='{template}',suffix='dseg.nii.gz'),
    envmodules: 'ants'
    container: config['singularity']['ants']
    log: 'logs/transform_to_subject/{template}_sub-{subject}_{targets}.log'
    group: 'participant1'
    threads: 8
    shell:
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.targets} -o {output} -r {input.ref} -t [{input.affine_xfm_itk},1] {input.inv_warp}  &> {log}'




rule binarize_subject_seed:
    input: 
        seed_res = bids(root='results',subject='{subject}',space='individual',label='{seed}',from_='{template}',suffix='probseg.nii.gz')
    params:
        threshold = lambda wildcards: config['seeds'][wildcards.seed]['probseg_threshold']
    output: 
        seed_thr = bids(root='results',subject='{subject}',space='individual',label='{seed}',from_='{template}',suffix='mask.nii.gz')
    container: config['singularity']['prepdwi']
    log: 'logs/binarize_subject_seed/{template}_sub-{subject}_{seed}.log'
    container: config['singularity']['prepdwi']
    group: 'participant1'
    shell:
        'c3d {input.seed_res} -threshold 0.5 inf 1 0 -o {output} &> {log}'
        


rule split_targets:
    input: 
        targets = bids(root='results',subject='{subject}',space='individual',desc='{targets}',from_='{template}',suffix='dseg.nii.gz'),
    params:
        target_nums = lambda wildcards: [str(i+1) for i in range(len(config['targets'][wildcards.targets]['labels']))],
        target_seg = lambda wildcards, output: expand('{target_seg_dir}/sub-{subject}_label-{target}_mask.nii.gz',target_seg_dir=output.target_seg_dir,subject=wildcards.subject,target=config['targets'][wildcards.targets]['labels'])
    output:
        target_seg_dir = directory(bids(root='results',subject='{subject}',desc='{targets}',from_='{template}',suffix='targets'))
    container: config['singularity']['prepdwi']
    log: 'logs/split_targets/sub-{subject}_{targets}_{template}.log'
    threads: 32 
    group: 'participant1'
    shell:  #TODO: could do this in c3d with less effort.. 
        'mkdir -p {output} && parallel  --jobs {threads} fslmaths {input.targets} -thr {{1}} -uthr {{1}} -bin {{2}} &> {log} ::: {params.target_nums} :::+ {params.target_seg}'


rule gen_targets_txt:
    input:
        target_seg_dir = bids(root='results',subject='{subject}',desc='{targets}',from_='{template}',suffix='targets')
    params:
        target_seg = lambda wildcards, input: expand('{target_seg_dir}/sub-{subject}_label-{target}_mask.nii.gz',target_seg_dir=input.target_seg_dir,subject=wildcards.subject,target=config['targets'][wildcards.targets]['labels'])
    output:
        target_txt = bids(root='results',subject='{subject}',suffix='targets.txt',desc='{targets}',from_='{template}')
    log: 
        bids(root='logs',subject='{subject}',suffix='targets.txt',desc='{targets}',from_='{template}')
    group: 'participant1'
    run:
        f = open(output.target_txt,'w')
        for s in params.target_seg:
            f.write(f'{s}\n')
        f.close()

rule run_probtrack:
    input:
        seed = bids(root='results',subject='{subject}',space='individual',label='{seed}',from_='{template}',suffix='mask.nii.gz'),
        target_txt = rules.gen_targets_txt.output,
        mask = bids(
            root="results",
            suffix="mask.nii.gz",
            desc="brain",
            space="T1w",
            res=config["resample_dwi"]["resample_scheme"],
            datatype="dwi",
            **config["subj_wildcards"]
        ),
        target_seg_dir = bids(root='results',subject='{subject}',desc='{targets}',from_='{template}',suffix='targets'),
        bedpost_dir=
            bids(
                root="results",
                desc="eddy",
                suffix="diffusion.bedpostX",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **config["subj_wildcards"]
            ),
    params:
        bedpost_merged = lambda wildcards, input: os.path.join(input.bedpost_dir,'merged'),
        probtrack_opts = config['probtrack']['opts'],
        out_target_seg = lambda wildcards, output: expand(bids(root=output.probtrack_dir,include_subject_dir=False,prefix='seeds_to',label='{target}',suffix='mask.nii.gz'), target=config['targets'][wildcards.targets]['labels']),
        nsamples = config['probtrack']['nsamples'],
#        container = config['singularity']['fsl_603']
        
    container: config['singularity']['fsl_603']
    output:
        probtrack_dir = directory(bids(root='results',subject='{subject}',label='{seed}',from_='{template}',desc='{targets}',suffix='probtrack'))
    threads: 8
    resources: 
        mem_mb = 8000, 
        time = 30, #30 mins
        gpus = 1 #1 gpu
    log: 'logs/run_probtrack/sub-{subject}_{seed}_{template}_{targets}.log'
    group: 'participant1'
    shell:
        #'mkdir -p {output.probtrack_dir} && singularity exec -e --nv {params.container} probtrackx2_gpu --samples={params.bedpost_merged}  --mask={input.mask} --seed={input.seed} ' 
        'mkdir -p {output.probtrack_dir} &&  probtrackx2 --samples={params.bedpost_merged}  --mask={input.mask} --seed={input.seed} ' 
        '--targetmasks={input.target_txt} --seedref={input.seed} --nsamples={params.nsamples} '
        '--dir={output.probtrack_dir} {params.probtrack_opts} -V 2  &> {log}'



rule track_from_seed:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
    input:
        wm_fod=bids(
            root="results",
            datatype='response',
            desc='normalized',
            suffix='wm_fod.mif',
            **config['subj_wildcards'],
        ),
        dwi=bids(
            root="results",
            datatype='dwi',
            suffix='dwi.mif',
            **config['subj_wildcards'],
        ),
        mask=bids(
            root="results",
            datatype='dwi',
            suffix='mask.mif',
            **config['subj_wildcards'],
        ),
        seed = bids(root='results',subject='{subject}',space='individual',label='{seed}',from_=config['template'],suffix='mask.nii.gz'),
    params:
        streamlines="",
#        seed_strategy=lambda wildcards,input: f'-seed_image {input.seed}'
        seed_strategy=lambda wildcards,input: f'-seed_random_per_voxel {input.seed} 50'

    output:
        tck=bids(
            root="results",
            datatype='tractography',
            label='{seed}',
            suffix='tractography.tck',
            **config['subj_wildcards'],
        ),
        seed_locs=bids(
            root="results",
            datatype='tractography',
            label='{seed}',
            suffix='seedlocs.txt',
            **config['subj_wildcards'],
        )

    threads: 32
    resources:
        mem_mb=128000,
	time=1440
    group: "subj2"
    container:
        config['singularity']['mrtrix']
    shell:
        'tckgen -nthreads {threads} -algorithm iFOD2 -mask {input.mask} '
        ' {input.wm_fod} {output.tck} '
        ' -seed_grid_per_voxel {input.seed} 1 '
        ' -output_seeds {output.seed_locs} '

rule create_voxel_seed_images:
    input:
        seed = bids(root='results',subject='{subject}',space='individual',label='{seed}',from_=config['template'],suffix='mask.nii.gz'),
    output:
        voxseeds_dir = directory(bids(root='results',subject='{subject}',space='individual',label='{seed}',from_=config['template'],suffix='voxseeds'))
    script: '../scripts/create_voxel_seed_images.py'

rule track_from_voxels:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
    input:
        wm_fod=bids(
            root="results",
            datatype='response',
            desc='normalized',
            suffix='wm_fod.mif',
            **config['subj_wildcards'],
        ),
        dwi=bids(
            root="results",
            datatype='dwi',
            suffix='dwi.mif',
            **config['subj_wildcards'],
        ),
        mask=bids(
            root="results",
            datatype='dwi',
            suffix='mask.mif',
            **config['subj_wildcards'],
        ),
        vox_seeds_dir = bids(root='results',subject='{subject}',space='individual',label='{seed}',from_=config['template'],suffix='voxseeds')
    params:
        seeds_per_voxel=config['seeds_per_voxel']
    output:
        tck_dir=directory(bids(
            root="results",
            datatype='tractography',
            label='{seed}',
            suffix='voxtracts',
            **config['subj_wildcards'],
        )),

    threads: 32
    resources:
        mem_mb=128000,
	time=1440
    group: "subj2"
    container:
        config['singularity']['mrtrix']
    shell:
        'mkdir -p {output.tck_dir} && '
        'parallel --jobs {threads} '
        'tckgen -quiet -nthreads 0  -algorithm iFOD2 -mask {input.mask} '
        ' {input.wm_fod} {output.tck_dir}/vox_{{1}}.tck '
        ' -seed_random_per_voxel {input.vox_seeds_dir}/seed_{{1}}.nii {params.seeds_per_voxel} '
        " ::: `ls {input.vox_seeds_dir} | grep -Po '(?<=seed_)[0-9]+'`"

rule connectivity_from_voxels:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
    input:
        tck_dir=bids(
            root="results",
            datatype='tractography',
            label='{seed}',
            suffix='voxtracts',
            **config['subj_wildcards'],
        ),
        targets = bids(root='results',subject='{subject}',space='individual',desc='{targets}',from_=config['template'],suffix='dseg.nii.gz'),
    output:
         conn_dir=directory(bids(
            root="results",
            datatype='tractography',
            desc='{targets}',
            label='{seed}',
            suffix='voxconn',
            **config['subj_wildcards'],
        )),

   
    threads: 32
    resources:
        mem_mb=128000,
	time=1440
    group: "subj2"
    container:
        config['singularity']['mrtrix']
    shell:
        'mkdir -p {output.conn_dir} && '
        'parallel --jobs {threads} '
        'tck2connectome -nthreads 0 -quiet {input.tck_dir}/vox_{{1}}.tck {input.targets} {output.conn_dir}/conn_{{1}}.csv -vector'
        " ::: `ls {input.tck_dir} | grep -Po '(?<=vox_)[0-9]+'`"


rule gen_conn_csv:
    # Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670
    input:
         conn_dir=bids(
            root="results",
            datatype='tractography',
            desc='{targets}',
            label='{seed}',
            suffix='voxconn',
            **config['subj_wildcards'],
        ),
    params:
        header_line=lambda wildcards: ','.join(config['targets'][wildcards.targets]['labels'])
    output:
         conn_csv=bids(
            root="results",
            datatype='tractography',
            desc='{targets}',
            label='{seed}',
            suffix='conn.csv',
            **config['subj_wildcards'],
        ),
    group: "subj2"
    container:
        config['singularity']['mrtrix']
    shell:
        'echo {params.header_line} > {output.conn_csv} && '
        'for f in `ls {input.conn_dir}/*.csv`; do tail -n 1 $f; done >> {output.conn_csv}'


rule conn_csv_to_image:
    input:
         conn_csv=bids(
            root="results",
            datatype='tractography',
            desc='{targets}',
            label='{seed}',
            suffix='conn.csv',
            **config['subj_wildcards'],
         ),
         seed_nii = bids(root='results',subject='{subject}',space='individual',label='{seed}',from_=config['template'],suffix='mask.nii.gz'),
    output:
          conn_nii=bids(
            root="results",
            datatype='tractography',
            desc='{targets}',
            label='{seed}',
            suffix='conn.nii.gz',
            **config['subj_wildcards'],
        ),
    script: '../scripts/conn_csv_to_image.py'



rule transform_conn_to_template:
    input:
        conn_nii=bids(
            root="results",
            datatype='tractography',
            desc='{targets}',
            label='{seed}',
            suffix='conn.nii.gz',
            **config['subj_wildcards'],
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
        ref = config['template_t1w']
    output:
        conn_nii=bids(
            root="results",
            datatype='tractography',
            space='{template}',
            desc='{targets}',
            label='{seed}',
            suffix='conn.nii.gz',
            **config['subj_wildcards'],
        ),
    container: config['singularity']['ants']
    threads: 8
    resources:
        mem_mb = 8000
    log: 'logs/transform_conn_to_template/sub-{subject}_{seed}_{template}_{targets}.log'
    group: 'participant1'
    shell:
        #using nearestneighbor to avoid bluring with background -- background set as -1
        'antsApplyTransforms -d 3 -e 3  --interpolation NearestNeighbor -i {input.conn_nii}  -o {output.conn_nii}  -r {input.ref} -t {input.warp} -t {input.affine_xfm_itk} &> {log} '


rule maxprob_conn:
    """ generate maxprob connectivity, adding outside striatum, and inside striatum (at a particular "streamline count" threshold) to identify unlabelled regions """
    input:
        conn_nii=bids(
            root="results",
            datatype='tractography',
            space='{template}',
            desc='{targets}',
            label='{seed}',
            suffix='conn.nii.gz',
            **config['subj_wildcards'],
        ),
    output:
        conn_nii=bids(
            root="results",
            datatype='tractography',
            space='{template}',
            desc='{targets}',
            label='{seed}',
            segtype='maxprob',
            suffix='dseg.nii.gz',
            **config['subj_wildcards'],
        ),
    shell:
        'c4d {input} -slice w 0:-1 -vote -o {output} '


"""
rule save_connmap_template_npz:
    input:
        mask = bids(root='results',template='{template}',label='{seed}',suffix='mask.nii.gz'),
        probtrack_dir = bids(root='results',subject='{subject}',label='{seed}',space='{template}',suffix='probtrack')
    params:
        connmap_3d = lambda wildcards, input: expand(bids(root=input.probtrack_dir,include_subject_dir=False,subject=wildcards.subject,prefix='seeds_to',label='{target}',suffix='mask.nii.gz'), target=targets),
    output:
        connmap_npz = bids(root='results',subject='{subject}',label='{seed}',space='{template}',suffix='connMap.npz')
    log: 'logs/save_connmap_to_template_npz/sub-{subject}_{seed}_{template}.log'
    group: 'participant1'
    conda: '../envs/sklearn.yml'
    script: '../scripts/save_connmap_template_npz.py'

#space-{template}
rule gather_connmap_group:
    input:
        connmap_npz = expand(bids(root='results',subject='{subject}',label='{seed}',space='{template}',suffix='connMap.npz'), subject=subjects,allow_missing=True)
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
