#!/usr/bin/env bash
# Prepare input data for the pipeline: reference, reads and truth set.
# Run once from the project root:  bash scripts/prepare_data.sh
set -euo pipefail

mkdir -p data
cd data

# --- Reference: chr17 of GRCh38 (UCSC hg38) and the indexes the tools need ---
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz
gunzip chr17.fa.gz
samtools faidx chr17.fa                    # .fai index (samtools, bcftools, GATK)
bwa index chr17.fa                         # BWA alignment index
gatk CreateSequenceDictionary -R chr17.fa  # .dict file (GATK)

# --- Reads: NA12878 (GIAB HG001) paired-end reads over BRCA1 ---
# Stream only the BRCA1 region (+-~15 kb) from the GIAB 300x Illumina BAM.
# -s 42.1 keeps a random 10% of reads (seed 42), giving ~30x coverage.
BAM="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/NHGRI_Illumina300X_novoalign_bams/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam"
samtools view -b -s 42.1 "$BAM" chr17:43030000-43140000 > na12878_brca1_subset.bam
# Sort by read name so mates are adjacent, then convert back to paired FASTQ.
samtools sort -n -o na12878_brca1_namesorted.bam na12878_brca1_subset.bam
samtools fastq -1 na12878_R1.fastq.gz -2 na12878_R2.fastq.gz -0 /dev/null -s /dev/null -n na12878_brca1_namesorted.bam

# --- Truth set: GIAB HG001 v4.2.1 benchmark variants and confident regions ---
TRUTH="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38"
bcftools view -r chr17:43030000-43140000 "$TRUTH/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" -Oz -o truth_brca1.vcf.gz
bcftools index -t truth_brca1.vcf.gz
curl -s -o HG001_benchmark.bed "$TRUTH/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"