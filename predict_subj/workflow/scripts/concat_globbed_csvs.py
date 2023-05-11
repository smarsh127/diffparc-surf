import os
from glob import glob

import pandas as pd

# glob for the csvs we want in the folder
csv_glob = snakemake.params.feature_csv_globs
csvs = glob(
    os.path.join(snakemake.input.extract_dir, csv_glob), recursive=True
)
# then concatenate them
pd.concat([pd.read_csv(csv) for csv in csvs]).to_csv(
    snakemake.output.csv, index=False
)
