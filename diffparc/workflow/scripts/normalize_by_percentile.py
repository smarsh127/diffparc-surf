import nibabel as nib
import numpy as np

nib_metric = nib.load(snakemake.input.metric)
metric_vol = nib_metric.get_fdata()

dseg_vol = nib.load(snakemake.input.dseg).get_fdata()

lower_perc = np.percentile(
    metric_vol[dseg_vol == snakemake.params.dseg_label], snakemake.params.lower_perc
)
upper_perc = np.percentile(
    metric_vol[dseg_vol == snakemake.params.dseg_label], snakemake.params.upper_perc
)

metric_vol = metric_vol - lower_perc / (
    upper_perc - lower_perc
)  # rescale from lower-upper to 0-1


out_nib = nib.Nifti1Image(
    metric_vol, affine=nib_metric.affine, header=nib_metric.header
)
out_nib.to_filename(snakemake.output.metric)
