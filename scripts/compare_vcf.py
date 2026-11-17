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

def filter_variants(variants, intervals):
    kept = set()
    for variant in variants:
        if in_regions(variant, intervals):
            kept.add(variant)
    return kept

truth_variants = read_vcf('data/truth_brca1.vcf.gz')
bcf_variants = read_vcf('results/variants/na12878.bcftools.vcf.gz')
gatk_variants = read_vcf('results/variants/na12878.gatk.vcf.gz')
truth_bed = read_bed("data/HG001_benchmark.bed")

truth_filtered = filter_variants(truth_variants, truth_bed)
bcf_filtered = filter_variants(bcf_variants, truth_bed)

tp = truth_filtered&bcf_filtered
fn = truth_filtered-bcf_filtered
fp = bcf_filtered-truth_filtered
recall = len(tp)/(len(tp)+len(fn))
precision = len(tp)/(len(tp)+len(fp))

print(f"# truth variants before filtering: {len(truth_variants)}\t # truth variants after filtering: {len(truth_filtered)}, \n# bcf variants before filtering: {len(bcf_variants)}\t # bcf variants after filtering: {len(bcf_filtered)} ")
print(f"tp: {len(tp)}, fp: {len(fp)}, fn: {len(fn)} \nPrecision: {round(precision,3)}, Recall: {round(recall,3)}" )