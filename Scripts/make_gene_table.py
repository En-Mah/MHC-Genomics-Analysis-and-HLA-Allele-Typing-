#!/usr/bin/env python3
import re, sys, csv

gtf_path = sys.argv[1]
genes_list_path = sys.argv[2]
out_bed_path = sys.argv[3]
out_tsv_path = sys.argv[4]

genes = set(line.strip() for line in open(genes_list_path) if line.strip())
if not genes:
    raise SystemExit("Gene list is empty.")

def get_attr(attr_str, key):
    m = re.search(rf'{key} "([^"]+)"', attr_str)
    return m.group(1) if m else None

# gene_name -> (gene_name, chrom, start, end, strand)
uniq = {}

with open(gtf_path) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attrs = parts
        if feature != "gene":
            continue

        gene_name = get_attr(attrs, "gene_name") or get_attr(attrs, "gene_id")
        if not gene_name:
            continue
        if gene_name not in genes:
            continue

        uniq[gene_name] = (gene_name, chrom, int(start), int(end), strand)

# TSV (1-based)
with open(out_tsv_path, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["gene_name","chromosome","start","end","strand"])
    for gene_name in sorted(uniq):
        w.writerow(list(uniq[gene_name]))

# BED (0-based)
with open(out_bed_path, "w") as out:
    for gene_name in sorted(uniq):
        _, chrom, start, end, strand = uniq[gene_name]
        out.write(f"{chrom}\t{start-1}\t{end}\t{gene_name}\t0\t{strand}\n")

print(f"Found {len(uniq)} genes in GTF (out of {len(genes)} input genes).")
