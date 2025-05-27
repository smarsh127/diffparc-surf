def get_dscalar_nii(wildcards):
    metric = wildcards.metric
    if metric == "surfarea" or metric == "inout" or metric == "surfarearatio":
        dscalar = bids(
            root=root,
            datatype="surf",
            label="{seed}",
            suffix="{metric}.dscalar.nii",
            **subj_wildcards,
        )
    elif (
        metric[:6] == "bundle"
    ):  # bundleFA, bundleMD, bundle...FA (normalized versions)
        dscalar = bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            label="{seed}",
            suffix="{metric}.dscalar.nii",
            **subj_wildcards,
        )
    elif (
        metric[:4] == "surf"
    ):  # this captures surfFA, surfMD, surf...FA (normalized versions) ...
        dscalar = bids(
            root=root,
            datatype="surf",
            label="{seed}",
            suffix="{metric}.dscalar.nii",
            **subj_wildcards,
        )

    return dscalar.format(**wildcards)


def get_dlabel_nii(wildcards):
    if config["use_template_parcellation"]:
        return bids(
            root=os.path.join(workflow.basedir, "..", "resources", "tpl-ctrlavg"),
            prefix="tpl-ctrlavg",
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            reduce="median",
            suffix="maxprob.dlabel.nii",
        ).format(**wildcards)

    else:
        return bids(
            root=root,
            datatype="surf",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="maxprob.dlabel.nii",
            **subj_wildcards,
        ).format(**wildcards)


rule write_surf_metrics_csv:
    """ for backwards compatiblity with old diffparc - 
    separate file for each metric, using identical column names (parcels), 
    and index column as "subj", formatted as sub-{subject}_ses-{session} """
    input:
        dscalar=get_dscalar_nii,
        dlabel=get_dlabel_nii,
    params:
        index_col_value=bids(
            **subj_wildcards, include_subject_dir=False, include_session_dir=False
        ),
        index_col_name="subj",
    output:
        csv=bids(
            root=root,
            datatype="tabular",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="{metric}.csv",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    script:
        "../scripts/write_surf_metrics.py"


def get_maxprob_dsegs(wildcards):
    out_dict = {}
    if wildcards.metric == "vol":
        for hemi in config["hemispheres"]:
            out_dict[f"maxprob_{hemi}"] = bids(
                root=root,
                datatype="anat",
                hemi=hemi,
                desc="{targets}",
                label="{seed}",
                seedspervoxel="{seedspervoxel}",
                method="{method}",
                segtype="maxprob",
                suffix="dseg.nii.gz",
                **subj_wildcards,
            ).format(**wildcards)
    elif wildcards.metric == "volmni":
        for hemi in config["hemispheres"]:
            out_dict[f"maxprob_{hemi}"] = bids(
                root=root,
                datatype="anat",
                hemi=hemi,
                space=config["template"],
                warp="linear",
                desc="{targets}",
                label="{seed}",
                seedspervoxel="{seedspervoxel}",
                method="{method}",
                segtype="maxprob",
                suffix="dseg.nii.gz",
                **subj_wildcards,
            ).format(**wildcards)
    return out_dict


def get_aux_dseg_for_metrics(wildcards):
    if wildcards.metric == "vol":
        return bids(
            root=root,
            datatype="anat",
            desc="{dseg_method}",
            suffix="dseg.nii.gz",
            **subj_wildcards,
        ).format(**wildcards)
    elif wildcards.metric == "volmni":
        return bids(
            root=root,
            datatype="anat",
            desc="{dseg_method}",
            space=config["template"],
            warp="linear",
            suffix="dseg.nii.gz",
            **subj_wildcards,
        ).format(**wildcards)


rule write_dseg_vol_metrics_csv:
    """ for backwards compatiblity with old diffparc - 
    separate file for each metric, using identical column names (parcels), 
    and index column as "subj", formatted as sub-{subject}_ses-{session} """
    input:
        dseg=get_aux_dseg_for_metrics,
        labels_tsv=lambda wildcards: os.path.join(
            workflow.basedir,
            "..",
            config["aux_dseg"][wildcards.dseg_method]["label_tsv"],
        ),
    params:
        index_col_value=bids(
            **subj_wildcards, include_subject_dir=False, include_session_dir=False
        ),
        index_col_name="subj",
    output:
        csv=bids(
            root=root,
            datatype="tabular",
            method="{dseg_method}",
            suffix="{metric,vol|volmni}.csv",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    script:
        "../scripts/write_aux_dseg_vol_metrics.py"


rule write_dseg_dti_metrics_csv:
    """ calcs mean dti metric in each dseg label """
    input:
        dseg=bids(
            root=root,
            datatype="anat",
            desc="{dseg_method}",
            suffix="dseg.nii.gz",
            **subj_wildcards
        ),
        metric=bids(
            root=root,
            datatype="dwi",
            resliced="{dseg_method}",
            suffix="{metric}.nii.gz",
            **subj_wildcards,
        ),
        labels_tsv=lambda wildcards: os.path.join(
            workflow.basedir,
            "..",
            config["aux_dseg"][wildcards.dseg_method]["label_tsv"],
        ),
    params:
        index_col_value=bids(
            **subj_wildcards, include_subject_dir=False, include_session_dir=False
        ),
        index_col_name="subj",
    output:
        csv=bids(
            root=root,
            datatype="tabular",
            method="{dseg_method}",
            suffix="{metric}.csv",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    script:
        "../scripts/write_aux_dseg_dti_metrics.py"


rule write_vol_metrics_csv:
    """ for backwards compatiblity with old diffparc - 
    separate file for each metric, using identical column names (parcels), 
    and index column as "subj", formatted as sub-{subject}_ses-{session} """
    input:
        unpack(get_maxprob_dsegs),
    params:
        index_col_value=bids(
            **subj_wildcards, include_subject_dir=False, include_session_dir=False
        ),
        index_col_name="subj",
        parc_list=lambda wildcards: config["targets"][wildcards.targets]["labels"],
    output:
        csv=bids(
            root=root,
            datatype="tabular",
            desc="{targets}",
            label="{seed}",
            seedspervoxel="{seedspervoxel}",
            method="{method}",
            suffix="{metric,vol|volmni}.csv",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    script:
        "../scripts/write_vol_metrics.py"


rule write_indepconn_metric_csv:
    """ this reads in the raw connectivity
    before normalization, and calculates the
    mean number of streamlines to each target
    from non-zero values"""
    input:
        csv_left=bids(
            root=root,
            datatype="tracts",
            hemi="L",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="conn.csv",
            **subj_wildcards,
        ),
        csv_right=bids(
            root=root,
            datatype="tracts",
            hemi="R",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="conn.csv",
            **subj_wildcards,
        ),
    params:
        target_labels=lambda wildcards: config["targets"][wildcards.targets]["labels"],
        index_col_value=bids(
            **subj_wildcards, include_subject_dir=False, include_session_dir=False
        ),
        index_col_name="subj",
    output:
        csv=bids(
            root=root,
            datatype="tabular",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            method="{method}",
            suffix="indepconn.csv",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    script:
        "../scripts/write_indepconn_metric.py"


rule write_surf_volumes_csv:
    input:
        vol_txts=expand(
            bids(
                root=root,
                **subj_wildcards,
                hemi="{hemi}",
                datatype="surf",
                suffix="{seed}.surfvolume.txt"
            ),
            hemi=config["hemispheres"],
            seed=config["seeds"].keys(),
            allow_missing=True,
        ),
    params:
        col_names=expand(
            "{hemi}_{seed}",
            hemi=config["hemispheres"],
            seed=config["seeds"].keys(),
            allow_missing=True,
        ),
        index_col_value=bids(
            **subj_wildcards, include_subject_dir=False, include_session_dir=False
        ),
        index_col_name="subj",
    output:
        csv=bids(
            root=root,
            datatype="tabular",
            method="{method,mrtrix|fsl}",
            suffix="surfvol.csv",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    script:
        "../scripts/write_surface_volume_metric.py"


rule write_surf_volumes_mni_csv:
    """ surf volume but after linear xfm to mni space"""
    input:
        vol_txts=expand(
            bids(
                root=root,
                **subj_wildcards,
                hemi="{hemi}",
                datatype="surf",
                space=config["template"],
                warp="linear",
                suffix="{seed}.surfvolume.txt"
            ),
            hemi=config["hemispheres"],
            seed=config["seeds"].keys(),
            allow_missing=True,
        ),
    params:
        col_names=expand(
            "{hemi}_{seed}",
            hemi=config["hemispheres"],
            seed=config["seeds"].keys(),
            allow_missing=True,
        ),
        index_col_value=bids(
            **subj_wildcards, include_subject_dir=False, include_session_dir=False
        ),
        index_col_name="subj",
    output:
        csv=bids(
            root=root,
            datatype="tabular",
            method="{method,mrtrix|fsl}",
            suffix="surfvolmni.csv",
            **subj_wildcards,
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "subj"
    script:
        "../scripts/write_surface_volume_metric.py"


rule concat_subj_csv:
    input:
        csvs=expand(
            bids(
                root=root,
                datatype="tabular",
                desc="{targets}",
                label="{seed}",
                seedspervertex="{seedspervertex}",
                suffix="{suffix}.csv",
                **subj_wildcards
            ),
            zip,
            **subj_zip_list,
            allow_missing=True
        ),
        #loop over subjects and sessions 
    output:
        csv=bids(
            root=root,
            subject="group",
            datatype="tabular",
            desc="{targets}",
            label="{seed}",
            seedspervertex="{seedspervertex}",
            suffix="{suffix}.csv",
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "agg"
    script:
        "../scripts/concat_csv.py"
