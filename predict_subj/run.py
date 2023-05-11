#!/usr/bin/env python3

import argparse
import os
import yaml
from snakemake import snakemake

def run_app(args):
    config = {
        'in_archive': args.in_archive,
        's3_bucket_models': args.s3_bucket_models,
        's3_model_zip': args.s3_model_zip,
        'out_csv': args.out_csv,
        'age': args.age,
        'sex': args.sex,
    }

    snakemake(os.path.join(os.path.dirname(__file__), 'workflow','Snakefile'),
                workdir=os.path.dirname(__file__), 
                config=config,
                force_incomplete=True,
                )

def main():
    parser = argparse.ArgumentParser(description='End-to-end diffparc prediction for a single subject, from dicom to classification output')
    parser.add_argument('--in-archive', required=True, help='Input archive file (tar, tgz, zip, tar.gz) containing dicom files for a single subject')
    parser.add_argument('--s3-model-zip', required=False, help='Path to input model zip inside S3 bucket')
    parser.add_argument('--s3-bucket-models', required=False, help='Name of S3 bucket containing photonai models')
    parser.add_argument('--age', required=True, help='Age of subject')
    parser.add_argument('--sex', required=True, help='Sex of subject')
    parser.add_argument('--out-csv', required=True, help='Output csv file')

    args = parser.parse_args()
    run_app(args)


if __name__ == '__main__':
    main()


