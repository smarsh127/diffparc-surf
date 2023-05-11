import pandas as pd

# read the features table
df = pd.read_csv(snakemake.input.csv)

for key,value in snakemake.params.metadata.items():
    df[key] = value


# save it
df.to_csv(snakemake.output.csv, index=False)
