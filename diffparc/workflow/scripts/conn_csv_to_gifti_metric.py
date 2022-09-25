import nibabel as nib
import numpy as np

conn = np.loadtxt(snakemake.input.csv,delimiter=',',skiprows=1)
conn_darray = nib.gifti.GiftiDataArray(
    data=conn,
    intent='NIFTI_INTENT_TIMESERIES',
    datatype='NIFTI_TYPE_INT32'
    )
gifti = nib.GiftiImage()
gifti.add_gifti_data_array(conn_darray)

gifti.to_filename(snakemake.output.gii_metric)




