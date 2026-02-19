#!/usr/bin/env python3
import csv, glob, os, sys

# --- Fix: allow very large CSV fields (gene sequences can be huge) ---
max_int = sys.maxsize
while True:
    try:
        csv.field_size_limit(max_int)
        break
    except OverflowError:
        max_int = int(max_int / 10)

TASK1_GLOB = "project/csv/*.task1.csv"
MHC_LIST   = "project/task3/lists/mhc_gene_list.tsv"
OUT_DIR    = "project/task3/mhc_tables"

os.makedirs(OUT_DIR, exist_ok=True)

def load_mhc_map(path):
    mhc = {}
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            g = row["gene_name"].strip()
            if not g or g.startswith("#"):
                continue
            mhc[g] = row["mhc_class"].strip()
    return mhc

def main():
    mhc_map = load_mhc_map(MHC_LIST)
    if not mhc_map:
        print("ERROR: MHC gene list is empty. Check mhc_gene_list.tsv", file=sys.stderr)
        sys.exit(1)

    files = sorted(glob.glob(TASK1_GLOB))
    if not files:
        print(f"ERROR: No Task1 CSVs found at {TASK1_GLOB}", file=sys.stderr)
        sys.exit(1)

    for path in files:
        sample = os.path.basename(path).replace(".task1.csv", "")
        out_path = os.path.join(OUT_DIR, f"{sample}.mhc_genes.tsv")

        with open(path, newline="") as f_in, open(out_path, "w", newline="") as f_out:
            r = csv.DictReader(f_in)

            if not r.fieldnames or "gene_name" not in r.fieldnames:
                raise RuntimeError(f"{path} does not look like a Task1 CSV (missing gene_name).")

            fieldnames = ["gene_name", "mhc_class"] + [c for c in r.fieldnames if c != "gene_name"]
            w = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t")
            w.writeheader()

            kept = 0
            for row in r:
                g = (row.get("gene_name") or "").strip()
                if g in mhc_map:
                    row_out = {"gene_name": g, "mhc_class": mhc_map[g]}
                    for k in r.fieldnames:
                        if k == "gene_name":
                            continue
                        row_out[k] = row.get(k, "")
                    w.writerow(row_out)
                    kept += 1

        print(f"[OK] {sample}: wrote {kept} MHC genes -> {out_path}")

if __name__ == "__main__":
    main()
