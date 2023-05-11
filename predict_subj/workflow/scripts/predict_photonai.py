from glob import glob

import pandas as pd
from photonai import Hyperpipe

# get permutation importances posthoc
pipe = Hyperpipe.load_optimum_pipe(
    glob(f"{snakemake.input.project_folder}/*/best_model.photonai")[0]
)

# get list of features
feature_df = pd.read_csv(snakemake.input.features)
feature_names = feature_df["features"].values

# load test data
test_df = pd.read_csv(snakemake.input.test_csv)

X = test_df[feature_names]


predictions = pipe.predict(X)

out_df = test_df[["subj", "age", "sex"]].copy()
out_df["predicted_group"] = predictions

out_df.to_csv(snakemake.output.predictions_csv, index=False)
