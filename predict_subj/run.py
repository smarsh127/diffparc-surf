#!/usr/bin/env python3

import argparse
import os
import yaml
from snakemake import snakemake

def run_app(args):
    config = {
        'in_archive': args.in_archive,
        'in_model': args.in_model,
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
    parser.add_argument('--in-model', required=True, help='Input photonai model directory')
    parser.add_argument('--age', required=True, help='Age of subject')
    parser.add_argument('--sex', required=True, help='Sex of subject')
    parser.add_argument('--out-csv', required=True, help='Output csv file')

    args = parser.parse_args()
    run_app(args)


if __name__ == '__main__':
    main()


