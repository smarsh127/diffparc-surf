from functools import reduce
from os.path import join

import pandas as pd

extract_ids = snakemake.params.extract_ids
subj_sess_list = list(extract_ids["patterns"].keys())

df_list = []

for csv in snakemake.input.csvs:

    df = (
        pd.read_table(csv, sep=",")
        # harmonize the subj column name:
        .rename(columns=extract_ids["column_rename"])
    )

    df_list.append(df)

# merge the columns (using the subj (sub-{subject}_ses-{session}) column as the index
df = reduce(lambda df1, df2: pd.merge(df1, df2, on="subj"), df_list)

# add subject column by pattern match from subj column:
df = df.assign(
    subject=lambda df_: df_.subj.str.extract(
        extract_ids["patterns"]["subject"]
    )
)

# optionally add session similarly to subject
if "session" in extract_ids["patterns"].keys():
    df = df.assign(
        session=lambda df_: df_.subj.str.extract(
            extract_ids["patterns"]["session"]
        )
    )


# apply custom query (e.g. to select the session):
if "query_filter" in extract_ids.keys():
    df = df.query(extract_ids["query_filter"])

# remove bad subjects
df = df[~df["subject"].isin(snakemake.params.bad_subjects)]


# save to csv
df.to_csv(snakemake.output.csv, index=False)
