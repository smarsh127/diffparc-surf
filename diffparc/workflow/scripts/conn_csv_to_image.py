import nibabel as nib
import numpy as np


seed_nii = snakemake.input.seed_nii
seed_nib = nib.load(seed_nii)
seed_vol = seed_nib.get_fdata()

mask = seed_vol>0

conn = np.loadtxt(snakemake.input.conn_csv,skiprows=1,delimiter=',')

num_targets = conn.shape[1]
conn_vol = np.zeros((seed_vol.shape[0],seed_vol.shape[1],seed_vol.shape[2],num_targets+1))

#we add a background channel too
bg_vol = np.zeros(seed_vol.shape)
bg_vol[~mask] = 1
conn_vol[:,:,:,0] = bg_vol

for i in range(num_targets):
    conn_single = np.zeros(seed_vol.shape)
    conn_single[mask] = conn[:,i]
    conn_vol[:,:,:,i+1] = conn_single


nii_conn_vol = nib.Nifti1Image(conn_vol.astype('int16'),affine=seed_nib.affine)

nii_conn_vol.to_filename(snakemake.output.conn_nii)

