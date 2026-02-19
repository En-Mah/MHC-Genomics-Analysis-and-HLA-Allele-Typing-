#!/usr/bin/env python3

from __future__ import annotations
import argparse
import csv
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple


try:
    csv.field_size_limit(1024 * 1024 * 1024)  # 1GB
except Exception:
    pass


def detect_input_files(input_path: Path) -> List[Path]:
    if input_path.is_file():
        return [input_path]
    if input_path.is_dir():
        return sorted(input_path.glob("*.mhc_step1.csv"))
    raise FileNotFoundError(f"Input path not found: {input_path}")


def infer_sample_name(step1_csv: Path) -> str:
    name = step1_csv.name
    if name.endswith(".mhc_step1.csv"):
        return name[:-len(".mhc_step1.csv")]
    return step1_csv.stem


def wrap_fasta(seq: str, width: int = 60) -> str:
    seq = re.sub(r"\s+", "", seq)
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


def csv_to_fasta(csv_in: Path, fasta_out: Path) -> Tuple[int, int]:
    """
    (rows_read, records_written)
    """
    with csv_in.open("r", newline="") as fin:
        reader = csv.DictReader(fin)
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in {csv_in}")

        required = ["gene_name", "sequence"]
        for r in required:
            if r not in reader.fieldnames:
                raise ValueError(f"{csv_in} missing required column: {r}. Found: {reader.fieldnames}")

        has_chr = "chromosome" in reader.fieldnames
        has_start = "start" in reader.fieldnames
        has_end = "end" in reader.fieldnames
        has_strand = "strand" in reader.fieldnames
        has_class = "mhc_class" in reader.fieldnames

        fasta_out.parent.mkdir(parents=True, exist_ok=True)
        rows = 0
        recs = 0

        with fasta_out.open("w") as fout:
            for row in reader:
                rows += 1
                gene = (row.get("gene_name") or "").strip()
                seq = (row.get("sequence") or "").strip()
                if not gene or not seq:
                    continue

                mhc_class = (row.get("mhc_class") or "").strip() if has_class else ""
                chrom = (row.get("chromosome") or "").strip() if has_chr else ""
                start = (row.get("start") or "").strip() if has_start else ""
                end = (row.get("end") or "").strip() if has_end else ""
                strand = (row.get("strand") or "").strip() if has_strand else ""

                loc = ""
                if chrom and start and end:
                    loc = f"{chrom}:{start}-{end}"
                    if strand:
                        loc += f"({strand})"

                header_parts = [gene]
                if mhc_class:
                    header_parts.append(f"class={mhc_class}")
                if loc:
                    header_parts.append(loc)

                header = ">" + "|".join(header_parts)
                fout.write(header + "\n")
                fout.write(wrap_fasta(seq) + "\n")

                recs += 1

        return rows, recs


def main() -> int:
    ap = argparse.ArgumentParser(description="Task3 Step2: Convert MHC Step1 CSV to FASTA.")
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output_dir", default="project/task3/step2_mhc_fasta")
    ap.add_argument("--log", default="project/task3/logs/task3_step2_csv_to_fasta.log")

    args = ap.parse_args()

    in_path = Path(args.input)
    out_dir = Path(args.output_dir)
    log_path = Path(args.log)

    files = detect_input_files(in_path)
    if not files:
        print(f"No Step1 CSV files found under: {in_path}", file=sys.stderr)
        return 2

    log_path.parent.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    with log_path.open("w") as log:
        log.write("sample\tinput_csv\trows_read\trecords_written\toutput_fasta\n")
        for f in files:
            sample = infer_sample_name(f)
            out_fa = out_dir / f"{sample}.mhc_step2.fa"
            rows, recs = csv_to_fasta(f, out_fa)
            log.write(f"{sample}\t{f}\t{rows}\t{recs}\t{out_fa}\n")
            print(f"[OK] {sample}: wrote {recs} FASTA records -> {out_fa}")

    print(f"\nWrote log: {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
