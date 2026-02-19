samtools view -H MOT36308.bam

samtools idxstats MOT36308.bam

gunzip gencode.v37lift37.annotation.gtf.gz

head -n 20 gencode.v37lift37.annotation.gtf

mkdir -p project/{bam,gtf,ref,counts,genes,seq,csv,logs}

samtools index project/bam/MOT36308.bam

awk '$3=="gene"' project/gtf/gencode.v37lift37.annotation.gtf > project/gtf/genes_only.gtf

samtools faidx project/ref/ucsc.hg19.fasta


for bam in project/bam/*.bam; do
  sample=$(basename "$bam" .bam)

  featureCounts -a project/gtf/gencode.v37lift37.annotation.gtf \
    -o "project/counts/${sample}.counts.txt" -g gene_name -t exon "$bam"

  awk 'BEGIN{FS="\t"} $1 !~ /^#/ && NR>2 { if($NF>0) print $1 }' \
    "project/counts/${sample}.counts.txt" > "project/genes/${sample}.genes.txt"

  python3 make_gene_table.py project/gtf/gencode.v37lift37.annotation.gtf \
    "project/genes/${sample}.genes.txt" \
    "project/genes/${sample}.genes.bed" \
    "project/genes/${sample}.genes.tsv"

  bedtools getfasta -fi project/ref/ucsc.hg19.fasta \
    -bed "project/genes/${sample}.genes.bed" -s -name \
    -fo "project/seq/${sample}.genes.fa"

  python3 merge_sequence.py \
    "project/genes/${sample}.genes.tsv" \
    "project/seq/${sample}.genes.fa" \
    "project/csv/${sample}.task1.csv"
done



ls -lh project/csv/*.task1.csv | head

brew install subread

featureCounts -v

conda install -c bioconda subread

featureCounts -v

bedtools --version

ls project/gtf/gencode.v37lift37.annotation.gtf
ls project/bam/MOT36308.bam

ls -lh project/counts/

conda install -c bioconda subread

featureCounts -v


featureCounts \
  -p \
  -a project/gtf/gencode.v37lift37.annotation.gtf \
  -o project/counts/MOT36308.counts.txt \
  -g gene_name \
  -t exon \
  project/bam/MOT36308.bam

ls -lh project/counts/MOT36308.counts.txt
head -n 5 project/counts/MOT36308.counts.txt

mkdir -p project/genes

awk 'BEGIN{FS="\t"} $1 !~ /^#/ && NR>2 { if($NF>0) print $1 }' project/counts/MOT36308.counts.txt \
  > project/genes/MOT36308.genes.txt

wc -l project/genes/MOT36308.genes.txt
head project/genes/MOT36308.genes.txt

paste -sd'|' project/genes/MOT36308.genes.txt > project/genes/MOT36308.genes.regex

wc -l project/genes/MOT36308.genes.tsv
head project/genes/MOT36308.genes.tsv

awk 'BEGIN{OFS="\t"} {print $2, $3-1, $4, $1, 0, $5}' project/genes/MOT36308.genes.tsv \
> project/genes/MOT36308.genes.bed

head project/genes/MOT36308.genes.bed

samtools faidx project/ref/ucsc.hg19.fasta

bedtools getfasta -fi project/ref/ucsc.hg19.fasta \
  -bed project/genes/MOT36308.genes.bed \
  -s -name \
  -fo project/seq/MOT36308.genes.fa

wc -l project/genes/MOT36308.genes.tsv
head project/genes/MOT36308.genes.tsv

ls -lh project/gtf/gencode.v37lift37.annotation.gtf
head -n 3 project/gtf/gencode.v37lift37.annotation.gtf

grep -m 1 'gene_name "' project/gtf/gencode.v37lift37.annotation.gtf

grep -m 5 'gene_name "TFAP2A"' project/gtf/gencode.v37lift37.annotation.gtf

mkdir -p project/genes
./make_gene_table.py \
  project/gtf/gencode.v37lift37.annotation.gtf \
  project/genes/MOT36308.genes.txt \
  project/genes/MOT36308.genes.bed \
  project/genes/MOT36308.genes.tsv

wc -l project/genes/MOT36308.genes.tsv
head project/genes/MOT36308.genes.tsv

grep -m 1 'gene_name "' project/gtf/gencode.v37lift37.annotation.gtf

./make_gene_table.py project/gtf/gencode.v37lift37.annotation.gtf project/genes/MOT36308.genes.txt project/genes/MOT36308.genes.bed project/genes/MOT36308.genes.tsv

head project/genes/MOT36308.genes.tsv
head project/genes/MOT36308.genes.bed

bedtools --version

conda install -c bioconda bedtools

ls -lh project/ref/ucsc.hg19.fasta

samtools faidx project/ref/ucsc.hg19.fasta

mkdir -p project/seq

bedtools getfasta \
  -fi project/ref/ucsc.hg19.fasta \
  -bed project/genes/MOT36308.genes.bed \
  -s -name \
  -fo project/seq/MOT36308.genes.fa

head -n 4 project/seq/MOT36308.genes.fa

mkdir -p project/csv
./merge_sequence.py \
  project/genes/MOT36308.genes.tsv \
  project/seq/MOT36308.genes.fa \
  project/csv/MOT36308.task1.csv

head project/csv/MOT36308.task1.csv

awk -F',' 'NR>1 && $6=="" {c++} END{print "empty sequences:", c+0}' project/csv/MOT36308.task1.csv

bedtools --version

ls -lh project/ref/ucsc.hg19.fasta

mkdir -p project/ref
cd project/ref

curl -L -o chromFa.tar.gz https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/chromFa.tar.gz

tar -xzf chromFa.tar.gz

cat chr*.fa > ucsc.hg19.fasta

cd ../../

grep -m 3 '^>' project/ref/ucsc.hg19.fasta

samtools faidx project/ref/ucsc.hg19.fasta

ls -lh project/ref/ucsc.hg19.fasta.fai

mkdir -p project/seq

bedtools getfasta \
  -fi project/ref/ucsc.hg19.fasta \
  -bed project/genes/MOT36308.genes.bed \
  -s -name \
  -fo project/seq/MOT36308.genes.fa

grep -m 3 '^>' project/ref/ucsc.hg19.fasta

head -n 4 project/seq/MOT36308.genes.fa

samtools faidx project/ref/ucsc.hg19.fasta

ls -lh project/ref/ucsc.hg19.fasta.fai

mkdir -p project/seq

bedtools getfasta \
  -fi project/ref/ucsc.hg19.fasta \
  -bed project/genes/MOT36308.genes.bed \
  -s -name \
  -fo project/seq/MOT36308.genes.fa

head -n 6 project/seq/MOT36308.genes.fa

mkdir -p project/csv
./merge_sequence.py \
  project/genes/MOT36308.genes.tsv \
  project/seq/MOT36308.genes.fa \
  project/csv/MOT36308.task1.csv

head project/csv/MOT36308.task1.csv

awk -F',' 'NR>1 && $6=="" {c++} END{print "empty sequences:", c+0}' project/csv/MOT36308.task1.csv

head -n 6 project/seq/MOT36308.genes.fa

./merge_sequence.py project/genes/MOT36308.genes.tsv project/seq/MOT36308.genes.fa project/csv/MOT36308.task1.csv

grep -m 10 '^>' project/seq/MOT36308.genes.fa

bedtools getfasta \
  -fi project/ref/ucsc.hg19.fasta \
  -bed project/genes/MOT36308.genes.bed \
  -s -name \
  -fo project/seq/MOT36308.genes.fa \
  2> project/logs/bedtools_getfasta_MOT36308.log

wc -l project/logs/bedtools_getfasta_MOT36308.log
head project/logs/bedtools_getfasta_MOT36308.log

bedtools getfasta \
  -fi project/ref/ucsc.hg19.fasta \
  -bed project/genes/MOT36308.genes.bed \
  -s -nameOnly \
  -fo project/seq/MOT36308.genes.fa

grep -m 5 '^>' project/seq/MOT36308.genes.fa

./merge_sequence.py \
  project/genes/MOT36308.genes.tsv \
  project/seq/MOT36308.genes.fa \
  project/csv/MOT36308.task1.csv

head -n 3 project/csv/MOT36308.task1.csv

bedtools getfasta \
  -fi project/ref/ucsc.hg19.fasta \
  -bed project/genes/MOT36308.genes.bed \
  -s -nameOnly \
  -fo project/seq/MOT36308.genes.fa

./merge_sequence.py project/genes/MOT36308.genes.tsv project/seq/MOT36308.genes.fa project/csv/MOT36308.task1.csv

head -n 3 project/csv/MOT36308.task1.csv

wc -l project/csv/MOT36308.task1.csv

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



ls project/csv/*.task1.csv | wc -l



cat > run_task1.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

./run_task1.sh

bash -n run_task1.sh

ls project/csv/*.task1.csv | wc -l

ls project/bam/*.bam | wc -l
ls project/csv/*.task1.csv | wc -l

