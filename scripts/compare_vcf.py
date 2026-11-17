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


set_truth_variants = set(truth_variants)
set_bcf_variants = set(bcf_variants)

tp = set_truth_variants&set_bcf_variants
fn = set_truth_variants-set_bcf_variants
fp = set_bcf_variants-set_truth_variants

recall = len(tp)/(len(tp)+len(fn))
precision = len(tp)/(len(tp)+len(fp))

print(f"tp: {len(tp)}, fp: {len(fp)}, fn: {len(fn)} \nPercision: {round(precision,3)}, Recall: {round(recall,3)}" )