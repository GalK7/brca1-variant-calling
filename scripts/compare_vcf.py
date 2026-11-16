import gzip

with gzip.open("data/truth_brca1.vcf.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        print(fields[0:5])
        break