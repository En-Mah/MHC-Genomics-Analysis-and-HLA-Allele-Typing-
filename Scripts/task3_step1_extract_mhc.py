#!/usr/bin/env python3

from __future__ import annotations
import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

try:
    csv.field_size_limit(1024 * 1024 * 1024)  # 1GB
except Exception:
    pass


MHC_CLASS_I = {
    "HLA-A", "HLA-B", "HLA-C",
    "HLA-E", "HLA-F", "HLA-G",
    "HLA-H", "HLA-J", "HLA-K", "HLA-L",
}

MHC_CLASS_II = {
    "HLA-DRA",
    "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
    "HLA-DQA1", "HLA-DQA2",
    "HLA-DQB1", "HLA-DQB2",
    "HLA-DPA1", "HLA-DPB1",
    "HLA-DMA", "HLA-DMB",
    "HLA-DOA", "HLA-DOB",
    "HLA-DCA", "HLA-DCB", 
}

MHC_CLASS_III = {
    "C2", "C4A", "C4B", "CFB",
    "TNF", "LTA", "LTB",
    "HSPA1A", "HSPA1B", "HSPA1L",
    "CYP21A2",
    "SKIV2L", "NELFE", "TNXB",
    "NOTCH4", "BATF2",
}


def normalize_gene_name(name: str) -> str:
    """Normalize gene symbol for matching."""
    n = (name or "").strip()
    if "::" in n:
        n = n.split("::", 1)[0].strip()
    return n


def load_custom_mhc_list(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with path.open("r", newline="") as f:
        sample = f.read(2048)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        reader = csv.DictReader(f, dialect=dialect)
        if reader.fieldnames is None:
            raise ValueError("Custom MHC list has no header.")

        fields = {h.strip().lower(): h for h in reader.fieldnames}
        if "gene_name" not in fields or "mhc_class" not in fields:
            raise ValueError("Custom MHC list must contain columns: gene_name, mhc_class")

        gene_col = fields["gene_name"]
        class_col = fields["mhc_class"]

        for row in reader:
            g = normalize_gene_name(row.get(gene_col, ""))
            c = (row.get(class_col, "") or "").strip().upper()
            if not g:
                continue
            if c in {"1"}:
                c = "I"
            if c in {"2"}:
                c = "II"
            if c in {"3"}:
                c = "III"
            if c not in {"I", "II", "III"}:
                continue
            mapping[g] = c
    return mapping


def built_in_mapping() -> Dict[str, str]:
    m: Dict[str, str] = {}
    for g in MHC_CLASS_I:
        m[g] = "I"
    for g in MHC_CLASS_II:
        m[g] = "II"
    for g in MHC_CLASS_III:
        m[g] = "III"
    return m


def detect_input_files(input_path: Path) -> List[Path]:
    if input_path.is_file():
        return [input_path]
    if input_path.is_dir():
        files = sorted(input_path.glob("*.task1.csv"))
        return files
    raise FileNotFoundError(f"Input path not found: {input_path}")


def infer_sample_name(csv_path: Path) -> str:
    name = csv_path.name
    if name.endswith(".task1.csv"):
        return name[:-len(".task1.csv")]
    return csv_path.stem


def read_csv_header(path: Path) -> List[str]:
    with path.open("r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, None)
    if not header:
        raise ValueError(f"Empty CSV or missing header: {path}")
    return header


def extract_mhc_rows(
    csv_in: Path,
    csv_out: Path,
    mhc_map: Dict[str, str],
    gene_column_candidates: Tuple[str, ...] = ("gene_name", "gene", "Gene", "GENE", "geneSymbol"),
) -> Tuple[int, int]:
    """
    Returns: (total_rows, mhc_rows)
    """
    with csv_in.open("r", newline="") as fin:
        reader = csv.DictReader(fin)
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in {csv_in}")

        fieldnames_lower = {h.lower(): h for h in reader.fieldnames}
        gene_col: Optional[str] = None
        for cand in gene_column_candidates:
            if cand.lower() in fieldnames_lower:
                gene_col = fieldnames_lower[cand.lower()]
                break
        if gene_col is None:
            raise ValueError(
                f"Could not find gene column in {csv_in}. "
                f"Expected one of: {gene_column_candidates}. Found: {reader.fieldnames}"
            )

        out_fields = list(reader.fieldnames)
        if "mhc_class" not in out_fields:
            out_fields.append("mhc_class")

        csv_out.parent.mkdir(parents=True, exist_ok=True)
        with csv_out.open("w", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=out_fields)
            writer.writeheader()

            total = 0
            kept = 0
            for row in reader:
                total += 1
                g = normalize_gene_name(row.get(gene_col, ""))
                if not g:
                    continue

                mhc_class = mhc_map.get(g)
                if mhc_class is None:
                    if g.startswith("HLA-"):
                        if g.startswith("HLA-D"):
                            mhc_class = "II"
                        else:
                            mhc_class = "I"

                if mhc_class is None:
                    continue

                row["mhc_class"] = mhc_class
                writer.writerow(row)
                kept += 1

    return total, kept


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Task3 Step1: Extract MHC genes (Class I/II/III) from Task1 CSV tables."
    )
    ap.add_argument(
        "-i", "--input",
        required=True,
    )
    ap.add_argument(
        "-o", "--output_dir",
        default="project/task3/step1_mhc_tables",
    )
    ap.add_argument(
        "--mhc_list",
        default=None,
    )
    ap.add_argument(
        "--log",
        default="project/task3/logs/task3_step1_extract_mhc.log",
    )

    args = ap.parse_args()

    in_path = Path(args.input)
    out_dir = Path(args.output_dir)
    log_path = Path(args.log)

    if args.mhc_list:
        mhc_map = load_custom_mhc_list(Path(args.mhc_list))
    else:
        mhc_map = built_in_mapping()

    files = detect_input_files(in_path)
    if not files:
        print(f"No input files found under: {in_path}", file=sys.stderr)
        return 2

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        log.write("sample\tinput_csv\ttotal_rows\tmhc_rows\toutput_csv\n")

        for f in files:
            sample = infer_sample_name(f)
            out_csv = out_dir / f"{sample}.mhc_step1.csv"
            total, kept = extract_mhc_rows(f, out_csv, mhc_map)
            log.write(f"{sample}\t{f}\t{total}\t{kept}\t{out_csv}\n")

            print(f"[OK] {sample}: kept {kept} / {total} rows -> {out_csv}")

    print(f"\nWrote log: {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
