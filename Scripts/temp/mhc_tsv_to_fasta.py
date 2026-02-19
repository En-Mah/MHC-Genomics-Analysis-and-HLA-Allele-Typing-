#!/usr/bin/env python3
import csv, glob, os, sys, textwrap

# allow very large fields (sequences)
max_int = sys.maxsize
while True:
    try:
        csv.field_size_limit(max_int)
        break
    except OverflowError:
        max_int = int(max_int / 10)

IN_GLOB = "project/task3/mhc_tables/*.mhc_genes.tsv"
OUT_DIR = "project/task3/mhc_fasta"
os.makedirs(OUT_DIR, exist_ok=True)

def wrap_seq(s, width=80):
    return "\n".join(textwrap.wrap(s, width=width))

def main():
    files = sorted(glob.glob(IN_GLOB))
    if not files:
        print(f"ERROR: No MHC TSV files found at {IN_GLOB}", file=sys.stderr)
        sys.exit(1)

    for path in files:
        sample = os.path.basename(path).replace(".mhc_genes.tsv", "")
        out_fa = os.path.join(OUT_DIR, f"{sample}.mhc_genes.fa")

        n = 0
        with open(path, newline="") as f, open(out_fa, "w") as out:
            r = csv.DictReader(f, delimiter="\t")
            required = {"gene_name","mhc_class","chromosome","start","end","strand","sequence"}
            if not required.issubset(set(r.fieldnames or [])):
                raise RuntimeError(f"{path} missing required columns. Found: {r.fieldnames}")

            for row in r:
                gene = row["gene_name"].strip()
                mhc_class = row["mhc_class"].strip()
                chrom = row["chromosome"].strip()
                start = row["start"].strip()
                end = row["end"].strip()
                strand = row["strand"].strip()
                seq = (row["sequence"] or "").strip().upper()

                if not gene or not seq:
                    continue

                header = f">{gene}|class={mhc_class}|{chrom}:{start}-{end}({strand})"
                out.write(header + "\n")
                out.write(wrap_seq(seq, 80) + "\n")
                n += 1

        print(f"[OK] {sample}: wrote {n} FASTA records -> {out_fa}")

if __name__ == "__main__":
    main()
