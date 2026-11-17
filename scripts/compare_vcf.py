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

def read_bed(path):
    intervals = []
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            region_start = int(fields[1])
            region_end = int(fields[2])
            if chrom == "chr17" and region_start < 43140000 and region_end > 43030000:
                intervals.append((chrom, region_start, region_end))
    return intervals

def in_regions(variant, intervals):
    for chrom, start, end in intervals:
        if variant[0] == chrom and start < variant[1] <= end:
            return True
    return False

truth_variants = read_vcf('data/truth_brca1.vcf.gz')
bcf_variants = read_vcf('results/variants/na12878.bcftools.vcf.gz')
gatk_variants = read_vcf('results/variants/na12878.gatk.vcf.gz')


# set_truth_variants = set(truth_variants)
# set_bcf_variants = set(bcf_variants)

# tp = set_truth_variants&set_bcf_variants
# fn = set_truth_variants-set_bcf_variants
# fp = set_bcf_variants-set_truth_variants

# recall = len(tp)/(len(tp)+len(fn))
# precision = len(tp)/(len(tp)+len(fp))

# print(f"tp: {len(tp)}, fp: {len(fp)}, fn: {len(fn)} \nPercision: {round(precision,3)}, Recall: {round(recall,3)}" )


truth_bed = read_bed("data/HG001_benchmark.bed")
# print(len(truth_bed), truth_bed[0], truth_bed[-1])

print(in_regions(truth_variants[0], truth_bed))