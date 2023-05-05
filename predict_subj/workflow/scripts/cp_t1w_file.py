import os
import re
import shutil

def copy_best_t1w_match(directory, output_path):
    t1w_patterns = [
        r't1w',
        r't1_weighted',
        r't1-w',
        r't1_mpr',
        r'mprage',
        r'sag_mpr',
        r'sagittal_mpr',
        r'sagittal_t1',
        r'sag_t1',
        r't1_sag',
        r't1_sagittal',
        r'spgr',
        r'fse_t1',
        r't1_fse',
        r't1_se',
        r't1_3d',
        r't1_volumetric',
    ]

    t1w_regex = re.compile('|'.join(t1w_patterns), re.IGNORECASE)

    best_t1w_file = None
    for filename in os.listdir(directory):
        if filename.lower().endswith('.nii') or filename.lower().endswith('.nii.gz'):
            if t1w_regex.search(filename):
                best_t1w_file = filename
                break

    if best_t1w_file is not None:
        shutil.copy(os.path.join(directory, best_t1w_file), output_path)
        print(f"Copied {best_t1w_file} to {output_path}")
    else:
        print("No T1-weighted MRI file found.")

# Specify the directory containing the NIfTI files
nifti_directory = snakemake.input.nifti_dir
output_path = snakemake.output.t1

# Call the function to copy the best match
copy_best_t1w_match(nifti_directory, output_path)

