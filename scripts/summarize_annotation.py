import pandas as pd

df = pd.read_csv("results/annotation/variants_table.tsv", sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "gt", "gnomad_af", "clnsig", "clnrevstat", "ann"])

df = df[df["alt"] != "*"]

df["first_ann"] = df["ann"].str.split(",").str[0]
fields = df["first_ann"].str.split("|", expand=True)
df["effect"] = fields[1]
df["impact"] = fields[2]
df["gene"] = fields[3]
df["hgvs_c"] = fields[9]
df["gnomad_af"] = pd.to_numeric(df["gnomad_af"], errors="coerce")

df["ba1"] = df["gnomad_af"] > 0.05

# print(f"Effects: {df["effect"].value_counts()}\n")
# print(f"ClinVar classification: {df["clnsig"].value_counts()}\n")
# print(f"Variants passing BA1 (gnomAD AF > 5%):{df["ba1"].sum()}")

review = df[(~df["ba1"]) & (df["impact"].isin(["HIGH","MODERATE"])| df["clnsig"].str.contains("Pathogenic|Uncertain", na=False))]
# print(len(review))
# print(review[["pos", "gene", "hgvs_c", "impact", "gnomad_af", "clnsig"]])

conflict = df[df["ba1"] & df["clnsig"].str.contains("Uncertain", na=False)]
# print(len(conflict))
# print(conflict[["pos", "gene", "hgvs_c", "impact", "gnomad_af", "clnsig"]])

summary = df[["chrom", "pos", "ref", "alt", "gt", "gene", "effect", "impact", "hgvs_c", "gnomad_af", "ba1", "clnsig", "clnrevstat"]]
summary.to_csv("results/annotation/variant_summary.tsv", sep="\t", index=False)