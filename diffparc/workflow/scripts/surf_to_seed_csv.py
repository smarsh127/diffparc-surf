import nibabel as nib
import numpy as np


surf_verts = nib.load(snakemake.input.surf).agg_data()[0]
np.savetxt(snakemake.output.csv,surf_verts,delimiter=',',fmt='%0.4f')

