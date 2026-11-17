import gzip

def read_vcf(path):
    variants = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            var_details = (fields[0], int(fields[1]), fields[3], fields[4])
            variants.append(var_details)
    return variants

truth_variants = read_vcf('data/truth_brca1.vcf.gz')
bcf_variants = read_vcf('results/variants/na12878.bcftools.vcf.gz')
gatk_variants = read_vcf('results/variants/na12878.gatk.vcf.gz')

print(len(truth_variants), len(bcf_variants), len(gatk_variants))