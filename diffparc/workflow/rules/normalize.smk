# rules for normalizing dti metrics to help deal with site effects

# uses a tissue segmentation to sample the normalization factors from WM only,
# to potentially avoid confounds with the amount of CSF etc..


rule zscore_norm:
    input:
        metric=bids(
            root=root, datatype="dwi", suffix="{metric}.nii.gz", **subj_wildcards
        ),
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="tissue",
            method="synthseg",
            suffix="dseg.nii.gz"
        ),
    params:
        mask_label=config["vbm"]["tissue_lut"]["WM"],
    output:
        metric=bids(
            root=root, datatype="dwi", suffix="znorm{metric}.nii.gz", **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/normalize_by_zscore.py"


rule perc_norm:
    input:
        metric=bids(
            root=root, datatype="dwi", suffix="{metric}.nii.gz", **subj_wildcards
        ),
        dseg=bids(
            root=root,
            datatype="anat",
            **subj_wildcards,
            desc="tissue",
            method="synthseg",
            suffix="dseg.nii.gz"
        ),
    params:
        mask_label=config["vbm"]["tissue_lut"]["WM"],
        lower_perc=0.02,
        upper_perc=0.98,
    output:
        metric=bids(
            root=root, datatype="dwi", suffix="pnorm{metric}.nii.gz", **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["pythondeps"]
    script:
        "../scripts/normalize_by_percentile.py"
