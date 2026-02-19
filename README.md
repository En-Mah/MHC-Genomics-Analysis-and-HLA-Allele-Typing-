# MHC-Genomics-Analysis-and-HLA-Allele-Typing

**Course:** Bioinformatics (Fall 2025)
**Project Date:** January 2026

---

## 1) Project Overview

This project implements a complete bioinformatics pipeline for analyzing sequencing data from the **Major Histocompatibility Complex (MHC)** region. The MHC region is located on **chromosome 6** and contains many immune-related genes, including **HLA genes**, which are essential for antigen presentation and immune response.

The MHC region is biologically important but technically challenging because it is:

* highly polymorphic
* gene-dense
* rich in homologous sequences
* difficult to analyze with short-read sequencing

The main purpose of this project is to go from raw alignment files (**BAM**) to interpretable results by completing three sequential tasks:

1. **Gene-level annotation and sequence extraction**
2. **Mapping quality (MAPQ) quality control**
3. **MHC-specific gene extraction and allele typing via alignment** 

---

## 2) Inputs and Dataset

### Primary Input Files

The project dataset contains multiple aligned sequencing samples:

* **20 BAM files** (paired-end sequencing alignments)

### Required Reference Resources

To interpret BAM alignments correctly, two reference files were required:

* **GTF gene annotation file** (must match BAM genome build)
* **Reference genome FASTA** (for extracting gene sequences)

Because the BAM files were aligned to **UCSC hg19**, the annotation and reference genome also had to be hg19-compatible. 

---

## 3) Tools and Software Used

This pipeline was executed using standard command-line genomics tools and Python scripting:

### Command-line tools

* **samtools**

  * Inspect BAM headers
  * Extract MAPQ values
  * Index reference FASTA

* **featureCounts** (Subread package)

  * Count fragments per gene (paired-end mode)

* **bedtools**

  * Extract gene sequences from hg19 FASTA using BED coordinates

### Programming

* **Python 3**

  * Parse GTF
  * Filter gene sets
  * Merge tables
  * Convert CSV → FASTA
  * Implement global alignment and allele typing

---

## 4) Task 1 — Gene-Level Annotation from BAM using GTF + FASTA

### 4.1 Objective

The goal of Task 1 was to annotate each BAM file with gene information and produce a structured table of genes covered by sequencing reads.

A gene was considered **covered/present** if:

> gene is covered ⇐⇒ gene fragment count > 0

The required output was **one CSV per sample**, containing gene metadata and reference sequences. 

---

### 4.2 Key Methodology

Task 1 required building a reproducible pipeline with the following steps.

---

### Step 1 — Verify BAM Reference Build

Before using any annotation, I checked the BAM header to confirm the genome build:

```bash
samtools view -H <sample>.bam
```

The BAM header included a reference path containing:

`/genomics/opt/RefGenomes/hg19/ucsc.hg19.fasta`

This confirmed that all BAM files were aligned to **UCSC hg19**, so the annotation and FASTA must match hg19 coordinate space. 

---

### Step 2 — Choose the Correct Gene Annotation (GTF)

Because the reference was hg19, I selected an hg19-compatible GTF:

**GENCODE v37lift37**:

`gencode.v37lift37.annotation.gtf`

This ensured gene coordinates aligned correctly with BAM coordinates and chromosome naming conventions (chr6, not 6). 

---

### Step 3 — Count Gene-Level Fragments with featureCounts

The BAM files are paired-end, so gene counting must be done in paired-end mode (fragment counting).

Command used:

```bash
featureCounts -p -B \
-a project/gtf/gencode.v37lift37.annotation.gtf \
-o project/counts/<sample>.counts.txt \
-g gene_name \
-t exon \
project/bam/<sample>.bam
```

Key settings:

* `-p`: paired-end fragment counting
* `-B`: require properly paired fragments
* `-t exon`: count exons and aggregate to gene-level
* `-g gene_name`: group features by gene name

This produced a gene-level count table for each sample. 

---

### Step 4 — Extract Genes with Count > 0

To get only covered genes:

```bash
awk 'BEGIN{FS="\t"} $1 !~ /^#/ && NR>2 { if($NF>0) print $1 }' \
project/counts/<sample>.counts.txt \
> project/genes/<sample>.genes.txt
```

Example result:

* Sample **MOT36308** contained **734 genes** with count > 0. 

---

### Step 5 — Extract Gene Coordinates from GTF (Python)

Counts alone are not enough. The CSV must contain:

* chromosome
* start
* end
* strand

A Python script was implemented to:

1. Load `<sample>.genes.txt`
2. Parse the GTF
3. Keep only `feature == "gene"`
4. Extract gene_name from attributes
5. Output a coordinate table (TSV)
6. Output a BED file for sequence extraction

Coordinate conversion:

* GTF uses **1-based inclusive**
* BED uses **0-based half-open**

So:

* BED start = GTF start − 1
* BED end = GTF end

This step was essential for correct sequence extraction. 

---

### Step 6 — Build and Index the hg19 FASTA

I downloaded UCSC hg19 chromosome FASTA files and concatenated them:

```bash
cat chr*.fa > ucsc.hg19.fasta
```

Then indexed the FASTA:

```bash
samtools faidx project/ref/ucsc.hg19.fasta
```

This enabled fast extraction of sequences by coordinate. 

---

### Step 7 — Extract Gene Sequences with bedtools

Using the BED intervals from Step 5:

```bash
bedtools getfasta \
-fi project/ref/ucsc.hg19.fasta \
-bed project/genes/<sample>.genes.bed \
-s -name \
-fo project/seq/<sample>.genes.fa
```

Important flags:

* `-s`: strand-aware extraction
* `-name`: use gene names in FASTA headers

---

### Step 8 — Merge Coordinates + Sequences into Final CSV

The final deliverable is one CSV per sample:

`project/csv/<sample>.task1.csv`

A key technical issue occurred:

bedtools FASTA headers looked like:

`>ABCF1::chr6:30539169-30564956(+)`

But the gene table contained only:

`ABCF1`

So the merge initially failed and returned:

> Missing sequences for 734 genes.

Fix:
The merge script was updated to normalize FASTA headers by splitting on `::` and keeping only the gene name.

After the fix, the pipeline produced:

> Wrote 734 rows. Missing sequences for 0 genes. 

---

### 4.3 Final Outputs (Task 1)

For each BAM file:

* `project/csv/<sample>.task1.csv`

And intermediate files for reproducibility:

* `project/counts/<sample>.counts.txt`
* `project/genes/<sample>.genes.txt`
* `project/genes/<sample>.genes.tsv`
* `project/genes/<sample>.genes.bed`
* `project/seq/<sample>.genes.fa`

The pipeline was successfully executed for **20 BAM files**, producing **20 final CSV outputs**. 

---

## 5) Task 2 — Mapping Quality (MAPQ) Analysis

### 5.1 Objective

The goal of Task 2 was to assess whether read alignments are reliable enough for downstream MHC interpretation.

This step is especially important because the MHC region often produces ambiguous alignments due to homology and repeats.

MAPQ (mapping quality) was used as the main metric:

* MAPQ ≥ 30 → confidently mapped reads
* MAPQ = 0 → ambiguous multi-mapped reads 

---

### 5.2 Methods

For each BAM file:

1. Extract MAPQ values using samtools

2. Generate plots for each sample:

   * histogram (counts)
   * histogram (normalized)
   * KDE (density)
   * CDF (cumulative distribution)

3. Across all samples:

   * bar plot: % reads with MAPQ ≥ 30
   * bar plot: % reads with MAPQ = 0

These plots were used to identify outlier samples and evaluate alignment confidence. 

---

### 5.3 Results Summary

Across the dataset:

* Most samples showed strong mapping quality
* Typically **~80–85%** of reads had MAPQ ≥ 30
* Some samples had elevated MAPQ = 0 fractions, indicating more ambiguity

A representative example sample (MOT36308) showed a **bimodal distribution**:

* a strong high-MAPQ peak
* a smaller MAPQ=0 peak

This pattern matches known characteristics of MHC sequencing. 

---

### 5.4 Interpretation

Conclusion from Task 2:

* The dataset is generally high-quality and suitable for downstream analysis.
* Samples with higher MAPQ=0 require caution because ambiguous mappings can bias:

  * coverage depth estimates
  * gene-level comparisons
  * allele typing results 

---

## 6) Task 3 — MHC Extraction and Allele Typing

### 6.1 Objective

Task 3 focused on MHC-specific analysis by:

1. Extracting MHC genes from Task 1 outputs
2. Converting extracted genes to FASTA
3. Performing global alignment to reference HLA alleles
4. Selecting the best allele for each gene based on:

   * alignment score
   * percent identity 

---

### 6.2 Inputs

Task 3 used:

* Task 1 per-sample gene tables:

  * `project/csv/MOT*.task1.csv`
* Reference HLA allele sequences:

  * downloaded from the IPD-IMGT/HLA database mirror
  * combined into:

    * `project/ref/hla_alleles.fa` 

---

### 6.3 Step 1 — Extract MHC Genes

Each Task 1 table was filtered to keep only MHC-related genes (Class I, II, III).

Output:

`project/task3/step1_mhc_tables/MOTxxxx.mhc_step1.csv`

Technical issue:
Some sequence fields exceeded Python’s default CSV field size limit.

Fix:
Increased CSV field size limit in the Step 1 script to allow robust parsing of long sequences. 

---

### 6.4 Step 2 — Convert MHC Tables to FASTA

Each filtered MHC table was converted into FASTA.

Example header format:

`>HLA-A|class=I|chr6:...`

Output:

`project/task3/step2_mhc_fasta/MOTxxxx.mhc_step2.fa`

Sanity check:
The number of FASTA records (`>`) was verified to match the number of extracted MHC genes. 

---

### 6.5 Step 3 — Pairwise Alignment and Allele Typing

#### Reference database preparation

Reference HLA FASTA files were assembled from locus-specific IMGT/HLA nucleotide FASTAs and concatenated into one file.

A key challenge:
IMGT/HLA headers often include internal IDs like `HLA00001`, meaning allele names are not always the first token.

Fix:
A robust parser was implemented that searches the full header line for allele patterns like:

* `A*01:01:01:01`
* `DRB1*15:01:01` 

---

#### Alignment method

For each sample gene sequence:

* perform global alignment against allele sequences for the same locus
* use Needleman–Wunsch dynamic programming
* scoring:

  * match = +2
  * mismatch = −1
  * gap = −2

Best allele selection:

1. maximum alignment score
2. tie-breaker: maximum percent identity 

---

### 6.6 Output Format

Final allele typing table per sample:

`project/task3/step3_alignment/MOTxxxx.mhc_allele_typing.real.csv`

Columns:

* `sample_id`
* `gene_name`
* `mhc_class`
* `best_allele`
* `read_support`
* `mean_alignment_score`
* `percent_identity` 

---

### 6.7 Results Summary

Task 3 was executed successfully for all 20 samples.

Outputs generated:

* 20 MHC filtered tables (Step 1)
* 20 MHC FASTA files (Step 2)
* allele typing results per sample (Step 3)

Observations:

* HLA Class I and II genes usually produced valid allele matches
* Many Class III genes produced empty allele results (not present in the reference database) 

---

## 7) Final Combined Interpretation (What the Whole Project Shows)

### 7.1 Sequencing and alignment quality is sufficient

Most samples contain high MAPQ values, suggesting reliable alignments.

### 7.2 Gene annotation confirms broad coverage

Each sample contains hundreds of genes with nonzero coverage (example: 734 genes in MOT36308).

### 7.3 MHC coverage is uneven but expected

Coverage differences across HLA genes reflect known MHC complexity and mapping difficulty.

### 7.4 Allele typing is possible but limited

High-confidence allele calls were obtained for some loci (example: HLA-A), but many loci remain difficult to type using short-read data alone. 

---

## 8) Limitations

This project is scientifically meaningful, but several limitations apply:

1. **Short-read sequencing limits HLA typing resolution**

   * many alleles are too similar to distinguish reliably

2. **MAPQ is aligner-dependent**

   * MAPQ values are not absolute truth, only confidence estimates

3. **Reference-based gene sequences are not sample haplotypes**

   * Task 1 sequences come from hg19 reference, not individual variants

4. **MHC multi-mapping affects coverage interpretation**

   * ambiguous reads can bias gene-level depth and allele calls 

---

## 9) Deliverables Produced

### Task 1 Deliverables

* `project/csv/<sample>.task1.csv` (20 files)

### Task 2 Deliverables

* MAPQ distribution plots (per sample)
* summary bar plots across samples
* written QC interpretation

### Task 3 Deliverables

* MHC-only tables (20 files)
* MHC FASTA outputs (20 files)
* allele typing outputs per sample

---

## 10) Conclusion

This project successfully built and executed a full MHC-focused genomics pipeline:

* BAM files were correctly annotated using hg19-compatible GTF + FASTA
* mapping quality was evaluated to validate alignment reliability
* MHC genes were extracted and typed using pairwise global alignment against IMGT/HLA alleles

The results demonstrate that:

* most samples have strong alignment quality
* MHC gene coverage is present but uneven (expected)
* allele typing is feasible for some loci but remains challenging overall

This workflow provides a reproducible foundation for deeper MHC analysis, including improved allele typing methods, variant calling, or long-read sequencing validation.
