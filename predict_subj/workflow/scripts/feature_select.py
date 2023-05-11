import pandas as pd
from snakemake.io import glob_wildcards


(feature_names,) = glob_wildcards(
    snakemake.input.features + "/touch/{feature}"
)

selected_features = feature_names

feature_df = pd.DataFrame(selected_features, columns=["features"])
feature_df.to_csv(snakemake.output.selected_features, index=False)
