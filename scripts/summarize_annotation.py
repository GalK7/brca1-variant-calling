import pandas as pd

df = pd.read_csv("results/annotation/variants_table.tsv", sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "gt", "gnomad_af", "clnsig", "clnrevstat", "ann"])

df["first_ann"] = df["ann"].str.split(",").str[0]
fields = df["first_ann"].str.split("|", expand=True)
df["effect"] = fields[1]
df["impact"] = fields[2]
df["gene"] = fields[3]
df["hgvs_c"] = fields[9]
df["gnomad_af"] = pd.to_numeric(df["gnomad_af"], errors="coerce")

print(df.shape)
print(df.head())