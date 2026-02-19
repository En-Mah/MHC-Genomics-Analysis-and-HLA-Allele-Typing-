#!/usr/bin/env bash
set -euo pipefail

mkdir -p project/{counts,genes,seq,csv,logs}

for bam in project/bam/*.bam; do
  sample=$(basename "$bam" .bam)
  echo "=== Processing $sample ==="

  featureCounts -p -B \
    -a project/gtf/gencode.v37lift37.annotation.gtf \
    -o "project/counts/${sample}.counts.txt" \
    -g gene_name -t exon \
    "$bam" > "project/logs/${sample}.featureCounts.log" 2>&1

  awk 'BEGIN{FS="\t"} $1 !~ /^#/ && NR>2 { if($NF>0) print $1 }' \
    "project/counts/${sample}.counts.txt" \
    > "project/genes/${sample}.genes.txt"

  ./make_gene_table.py \
    project/gtf/gencode.v37lift37.annotation.gtf \
    "project/genes/${sample}.genes.txt" \
    "project/genes/${sample}.genes.bed" \
    "project/genes/${sample}.genes.tsv" \
    > "project/logs/${sample}.make_gene_table.log" 2>&1

  bedtools getfasta \
    -fi project/ref/ucsc.hg19.fasta \
    -bed "project/genes/${sample}.genes.bed" \
    -s -name \
    -fo "project/seq/${sample}.genes.fa" \
    2> "project/logs/${sample}.bedtools_getfasta.log"

  ./merge_sequence.py \
    "project/genes/${sample}.genes.tsv" \
    "project/seq/${sample}.genes.fa" \
    "project/csv/${sample}.task1.csv" \
    > "project/logs/${sample}.merge.log" 2>&1
done

echo "All done. CSV count:"
ls -1 project/csv/*.task1.csv 2>/dev/null | wc -l
