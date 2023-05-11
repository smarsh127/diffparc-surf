import pandas as pd

df = pd.read_csv(snakemake.input.csv)

# prepend kind to all but 1st column
df.columns = [df.columns[0]] + [
    snakemake.wildcards.kind + "_" + col for col in df.columns[1:]
]

df.to_csv(snakemake.output.csv, index=False)
