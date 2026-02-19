#!/usr/bin/env python3

from __future__ import annotations
import argparse
import hashlib
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


FASTA_HEADER_RE = re.compile(r"^>(?P<gene>[^|]+)(?:\|class=(?P<class>[^|]+))?(?:\|(?P<loc>.+))?$")


def detect_input_files(input_path: Path) -> List[Path]:
    if input_path.is_file():
        return [input_path]
    if input_path.is_dir():
        return sorted(input_path.glob("*.mhc_step2.fa"))
    raise FileNotFoundError(f"Input path not found: {input_path}")


def infer_sample_name(fa: Path) -> str:
    name = fa.name
    if name.endswith(".mhc_step2.fa"):
        return name[:-len(".mhc_step2.fa")]
    return fa.stem


def parse_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    header = None
    seq_chunks: List[str] = []
    with path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def gc_pct(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return 100.0 * gc / len(seq)


def sha1_sig(seq: str) -> str:
    return hashlib.sha1(seq.encode("utf-8")).hexdigest()


def extract_hla_records(fa_in: Path) -> List[Dict[str, str]]:
    out: List[Dict[str, str]] = []
    for header, seq in parse_fasta(fa_in):
        m = FASTA_HEADER_RE.match(">" + header)
        gene = header.split("|", 1)[0].strip() if m is None else (m.group("gene") or "").strip()
        mhc_class = ""
        if m is not None and m.group("class"):
            mhc_class = m.group("class").strip()

        if not gene.startswith("HLA-"):
            continue

        seq = re.sub(r"\s+", "", seq)
        rec = {
            "gene_name": gene,
            "mhc_class": mhc_class,
            "length": str(len(seq)),
            "gc_pct": f"{gc_pct(seq):.2f}",
            "signature_sha1": sha1_sig(seq),
        }
        out.append(rec)
    return out


def write_csv(rows: List[Dict[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["gene_name", "mhc_class", "length", "gc_pct", "signature_sha1"]
    with out_path.open("w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(r.get(c, "") for c in cols) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description="Task3 Step3: Simple HLA typing-like summary from MHC FASTA.")
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output_dir", default="project/task3/step3_typing")
    ap.add_argument("--write_all_summary", action="store_true")
    args = ap.parse_args()

    in_path = Path(args.input)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = detect_input_files(in_path)
    if not fasta_files:
        print(f"No input FASTA found under: {in_path}", file=sys.stderr)
        return 2

    all_rows: List[Dict[str, str]] = []

    for fa in fasta_files:
        sample = infer_sample_name(fa)
        rows = extract_hla_records(fa)
        out_csv = out_dir / f"{sample}.hla_typing.csv"
        write_csv(rows, out_csv)
        print(f"[OK] {sample}: HLA records={len(rows)} -> {out_csv}")

        if args.write_all_summary:
            for r in rows:
                rr = dict(r)
                rr["sample"] = sample
                all_rows.append(rr)

    if args.write_all_summary:
        summary_path = out_dir / "ALL.hla_typing_summary.csv"
        cols = ["sample", "gene_name", "mhc_class", "length", "gc_pct", "signature_sha1"]
        with summary_path.open("w") as f:
            f.write(",".join(cols) + "\n")
            for r in all_rows:
                f.write(",".join(r.get(c, "") for c in cols) + "\n")
        print(f"[OK] Wrote summary -> {summary_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
