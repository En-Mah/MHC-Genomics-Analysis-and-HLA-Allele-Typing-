#!/usr/bin/env python3
import sys, csv

tsv_path = sys.argv[1]
fasta_path = sys.argv[2]
out_csv = sys.argv[3]

# FASTA -> dict(gene_name -> sequence)
seq = {}
name = None
chunks = []

def normalize(header: str) -> str:
    return header.split("::")[0]

with open(fasta_path) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                seq[name] = "".join(chunks)
            raw = line[1:].split()[0]
            name = normalize(raw)
            chunks = []
        else:
            chunks.append(line)
    if name is not None:
        seq[name] = "".join(chunks)

with open(tsv_path) as inp, open(out_csv, "w", newline="") as out:
    r = csv.DictReader(inp, delimiter="\t")
    fieldnames = ["gene_name","chromosome","start","end","strand","sequence"]
    w = csv.DictWriter(out, fieldnames=fieldnames)
    w.writeheader()

    total = 0
    missing = 0
    for row in r:
        total += 1
        g = row["gene_name"]
        s = seq.get(g, "")
        if s == "":
            missing += 1
        w.writerow({
            "gene_name": g,
            "chromosome": row["chromosome"],
            "start": row["start"],
            "end": row["end"],
            "strand": row["strand"],
            "sequence": s
        })

print(f"Wrote {total} rows. Missing sequences for {missing} genes.")
