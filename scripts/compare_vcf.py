import gzip

def read_vcf(path):
    variants = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            var_details = (fields[0], int(fields[1]), fields[3].upper(), fields[4].upper())
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

def compare(truth, calls):
    tp = truth&calls
    fn = truth-calls
    fp = calls-truth
    recall = len(tp)/(len(tp)+len(fn))
    precision = len(tp)/(len(tp)+len(fp))
    return tp, fp, fn, precision, recall

truth_variants = read_vcf('data/truth_brca1.norm.vcf.gz')
bcf_variants = read_vcf('results/variants/na12878.bcftools.norm.vcf.gz')
gatk_variants = read_vcf('results/variants/na12878.gatk.norm.vcf.gz')
truth_bed = read_bed("data/HG001_benchmark.bed")

truth_filtered = filter_variants(truth_variants, truth_bed)
bcf_filtered = filter_variants(bcf_variants, truth_bed)
gatk_filtered = filter_variants(gatk_variants, truth_bed)

tools = [("bcftools", bcf_filtered), ("GATK", gatk_filtered)]

print(f"# truth variants before filtering: {len(truth_variants)}\t # truth variants after filtering: {len(truth_filtered)}, \n# bcf variants before filtering: {len(bcf_variants)}\t # bcf variants after filtering: {len(bcf_filtered)}, \n# gatk variants before filtering: {len(gatk_variants)}\t # gatk variants after filtering: {len(gatk_filtered)}")

for tool_name, tool_calls in tools:
    tp, fp, fn, precision, recall = compare(truth_filtered, tool_calls)
    print(f"\n{tool_name} results:")
    print(f"tp: {len(tp)}, fp: {len(fp)}, fn: {len(fn)} \nPrecision: {round(precision,3)}, Recall: {round(recall,3)}" )
    print("FP:")
    for var in sorted(fp):
        print(var)
    print("FN:")
    for var in sorted(fn):
        print(var)