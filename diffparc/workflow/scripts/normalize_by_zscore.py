import nibabel as nib
import numpy as np

nib_metric = nib.load(snakemake.input.metric)
metric_vol = nib_metric.get_fdata()

dseg_vol = nib.load(snakemake.input.dseg).get_fdata()


mean = np.mean(metric_vol[dseg_vol == snakemake.params.dseg_label])
std = np.std(metric_vol[dseg_vol == snakemake.params.dseg_label])

metric_vol = (metric_vol - mean) / std

out_nib = nib.Nifti1Image(
    metric_vol, affine=nib_metric.affine, header=nib_metric.header
)
out_nib.to_filename(snakemake.output.metric)
