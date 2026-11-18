# Germline variant calling on NA12878 (BRCA1), benchmarked against GIAB

A reproducible germline variant-calling pipeline, from paired-end FASTQ to VCF, run with two variant callers (bcftools and GATK HaplotypeCaller) and evaluated against the Genome in a Bottle (GIAB) truth set for NA12878. The GATK calls are then annotated with SnpEff, gnomAD allele frequencies and ClinVar classifications.

The pipeline is written in Snakemake, and the benchmarking and summary scripts are written in plain Python.

## Data

| | Source |
|---|---|
| Sample | NA12878 (GIAB HG001), Illumina HiSeq 2×148 bp |
| Region | chr17:43,030,000–43,140,000 (BRCA1 ± ~15 kb) |
| Reads | Streamed from the GIAB 300x BAM and randomly subsampled to 10% (~30x), then converted back to FASTQ |
| Reference | GRCh38 chr17 (UCSC hg38) |
| Truth set | GIAB HG001 v4.2.1 benchmark VCF and confident-region BED |
| Gene annotation | SnpEff database `GRCh38.mane.1.5.refseq` (MANE Select transcripts) |
| Population frequencies | gnomAD v4.1 genomes |
| Clinical classifications | ClinVar GRCh38 VCF (latest weekly release) |

All downloads and preparation steps are in [`scripts/prepare_data.sh`](scripts/prepare_data.sh). The SnpEff database is downloaded by the pipeline itself.

## Pipeline

```
FASTQ ─ fastp ─ BWA-MEM ─ sort/index ─ MarkDuplicates ─┬─ bcftools mpileup/call ─┐
                                                       └─ GATK HaplotypeCaller ──┴─ bcftools norm ─ compare_vcf.py

GATK normalized VCF ─ SnpEff ─ gnomAD AF ─ ClinVar ─ summarize_annotation.py
```

| Step | Tool |
|---|---|
| Read QC and adapter/quality trimming | fastp |
| Alignment | BWA-MEM, samtools sort/index |
| Duplicate marking | GATK MarkDuplicates |
| Variant calling | bcftools mpileup + call; GATK HaplotypeCaller |
| Normalization (left-align, trim, split multi-allelic sites) | bcftools norm |
| Benchmarking | [`scripts/compare_vcf.py`](scripts/compare_vcf.py) |
| Functional annotation | SnpEff |
| Population frequency and ClinVar annotation | bcftools annotate |
| Annotation summary | [`scripts/summarize_annotation.py`](scripts/summarize_annotation.py) |

The workflow is defined in [`workflow/Snakefile`](workflow/Snakefile). SnpEff runs in its own conda environment ([`workflow/envs/snpeff.yaml`](workflow/envs/snpeff.yaml)), because it needs a newer Java version than GATK.

### Benchmarking

`compare_vcf.py` reads the normalized VCFs and the GIAB BED file, keeps only variants inside GIAB confident regions (for both the truth set and the calls), and compares them by chromosome, position, REF and ALT. It reports true positives (TP), false positives (FP), false negatives (FN), precision and recall for each caller, and lists every FP and FN.

### Annotation

`summarize_annotation.py` takes the first (highest-impact) SnpEff annotation of each variant and applies two filters:
- **BA1**: gnomAD allele frequency above 5%. Under the ACMG/AMP guidelines this is stand-alone evidence for a benign classification.
- **Review**: variants that do not pass BA1 and either have HIGH or MODERATE predicted impact, or are classified in ClinVar as pathogenic or of uncertain significance (VUS).

It also lists ClinVar VUS that pass BA1, since a common variant classified as VUS is a conflict worth noting. The script writes a summary table and a text report.

## Results

### Benchmarking

170 truth variants fall inside the confident regions.

| Caller | TP | FP | FN | Precision | Recall |
|---|---|---|---|---|---|
| bcftools | 168 | 0 | 2 | 1.000 | 0.988 |
| GATK HaplotypeCaller | 169 | 0 | 1 | 1.000 | 0.994 |

Missed variants:
- **chr17:43109041**, a 16 bp insertion in a short tandem repeat: missed by bcftools, found by GATK. HaplotypeCaller reassembles local haplotypes, which helps with longer indels.
- **chr17:43038123**, a multi-allelic site in an AC repeat (two different insertion alleles): both callers found one allele and missed the other.

### Why normalization matters

The first comparison inside the confident regions, without normalization, gave bcftools 25 FP and 26 FN (precision 0.85, recall 0.845). Looking at them showed that almost every FP had an FN at the same position: the same variant, written differently by the caller and by GIAB. 49 of the 51 "errors" were representation differences:
- indels written with the whole repeat instead of the minimal form, or not left-aligned;
- lower-case reference bases copied from the soft-masked UCSC reference;
- multi-allelic sites written on one line.

Normalizing all VCFs with `bcftools norm` and comparing alleles in upper case removed them, leaving the 2 real misses above. GATK output needed almost no realignment, while bcftools had 52 realigned records.

### Annotation

236 GATK variants were annotated (after removing spanning-deletion `*` alleles).
- Most are non-coding: 173 intronic, 25 intergenic, 27 upstream/downstream and 4 in 3' UTRs. 4 are missense and 3 are synonymous.
- 230 of 236 pass BA1. The 4 missense variants are all common (gnomAD AF 0.30–0.49).
- ClinVar: 176 benign or likely benign, 1 VUS, 59 not in ClinVar.
- **0 variants were flagged for review.** No rare protein-altering variants and no pathogenic variants were found, which is the expected result for this sample.
- The single VUS, BRCA1 c.\*872_\*873delAA (3' UTR, one ClinVar submitter), has a gnomAD AF of 0.28 and passes BA1, so it is most likely benign.

## How to run

```bash
conda env create -f environment.yml
conda activate varcall
bash scripts/prepare_data.sh    # download and prepare input data (once)
snakemake -c 4 --use-conda      # run the pipeline
cat results/benchmark.txt                         # benchmark report
cat results/annotation/annotation_report.txt      # annotation report
```

On the first run, Snakemake creates the SnpEff environment and downloads the SnpEff database (several hundred MB).

## Repository structure

```
workflow/Snakefile                Snakemake pipeline
workflow/envs/snpeff.yaml         conda environment for SnpEff
scripts/prepare_data.sh           data download and preparation
scripts/compare_vcf.py            benchmarking against the truth set
scripts/summarize_annotation.py   annotation summary and filters
environment.yml                   main conda environment
```

`data/` and `results/` are created locally and are not tracked in Git.

## Limitations

- One sample, one ~110 kb region, ~30x coverage: a learning-scale benchmark, not a full validation.
- Reads were extracted from a BAM that was already aligned to this region, so the alignment rate (100%) is not representative of whole-genome data.
- Variants are matched by allele only. Genotypes are not compared, and the comparison is not haplotype-aware like dedicated tools such as hap.py or RTG vcfeval.
- After 10% subsampling almost no PCR duplicates remain, so MarkDuplicates has little effect here.
- Only the GATK calls are annotated, and only the first SnpEff annotation of each variant is used.
- BA1 is the only ACMG/AMP criterion applied; this is a frequency filter, not a full variant classification.
- The ClinVar download is the latest release, so classifications may change between runs.