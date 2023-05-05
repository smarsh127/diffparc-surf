#!/usr/bin/env python3

import argparse
import os
import yaml
from snakemake import snakemake

def run_app(args):
    config = {
        'in_archive': args.in_archive,
        'in_model': args.in_model,
        'out_csv': args.out_csv
    }

    snakemake(os.path.join(os.path.dirname(__file__), 'workflow','Snakefile'),
                workdir=os.path.dirname(__file__), 
                config=config,
                force_incomplete=True,
                )

def main():
    parser = argparse.ArgumentParser(description='Run Snakemake workflow with specified parameters.')
    parser.add_argument('--in-archive', required=True, help='Input archive file')
    parser.add_argument('--in-model', required=True, help='Input model file')
    parser.add_argument('--out-csv', required=True, help='Output CSV file')

    args = parser.parse_args()
    run_app(args)


if __name__ == '__main__':
    main()


