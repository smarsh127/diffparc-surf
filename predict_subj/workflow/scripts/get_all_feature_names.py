from snakemake.shell import shell
import pandas as pd

shell("mkdir -p {snakemake.output.features}/touch")
df = pd.read_csv(snakemake.input.csv)
for feat in df.columns[1:-2]:  # leave out first col, and last 2 cols
    shell("touch {snakemake.output.features}/touch/{feat}")

