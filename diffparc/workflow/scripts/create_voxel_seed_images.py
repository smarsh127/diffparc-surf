import nibabel as nib
import numpy as np
import os

seed_nib  = nib.load(snakemake.input.seed)
seed_vol = seed_nib.get_fdata()
os.mkdir(snakemake.output.voxseeds_dir)


for i,vox in enumerate(np.argwhere(seed_vol)):
    out_vol = np.zeros(seed_vol.shape,dtype='uint8')
    out_vol[vox[0],vox[1],vox[2]] = 1
    out_nii = nib.Nifti1Image(out_vol, affine=seed_nib.affine)
    out_nii.to_filename(os.path.join(snakemake.output.voxseeds_dir,f'seed_{i:05d}.nii'))
