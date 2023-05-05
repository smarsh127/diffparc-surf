import os
import tarfile
import zipfile
import shutil

def is_tarfile(file_path):
    return tarfile.is_tarfile(file_path)

def is_zipfile(file_path):
    return zipfile.is_zipfile(file_path)

def flatten_and_extract(archive_file_path, output_folder):
    temp_extract_folder = os.path.join(output_folder, 'temp_extract')
    os.makedirs(temp_extract_folder, exist_ok=True)

    if is_tarfile(archive_file_path):
        with tarfile.open(archive_file_path, 'r') as tar:
            tar.extractall(temp_extract_folder)
    elif is_zipfile(archive_file_path):
        with zipfile.ZipFile(archive_file_path, 'r') as zip_ref:
            zip_ref.extractall(temp_extract_folder)
    else:
        raise ValueError(f"Unsupported archive format for file: {archive_file_path}")

    for root, dirs, files in os.walk(temp_extract_folder):
        for file in files:
            file_path = os.path.join(root, file)
            shutil.move(file_path, os.path.join(output_folder, file))

    shutil.rmtree(temp_extract_folder)

archive_file_path = snakemake.input.archive
output_folder = snakemake.output.dicom_dir
flatten_and_extract(archive_file_path, output_folder)

