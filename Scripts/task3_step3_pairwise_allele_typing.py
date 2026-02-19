#!/usr/bin/env python3
"""
Task 3 - Step 3: Pairwise alignment & allele typing (IMGT/HLA reference)
Output format matches the provided example CSV:

sample_id,gene_name,mhc_class,best_allele,read_support,mean_alignment_score,percent_identity

- We align sample gene sequences (from Step2 FASTA) to reference allele sequences.
- For each sample gene (typically HLA-*), we pick the best matching allele by:
    1) alignment score (global Needleman-Wunsch)
    2) percent identity as tie-breaker
- For genes without reference alleles: best_allele empty, support=0, score=0, pid=0.
"""

from __future__ import annotations
import argparse
import csv
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Optional
import hashlib


# ---------------- FASTA parsing ----------------
def parse_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    header = None
    seq_chunks: List[str] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def clean_seq(seq: str) -> str:
    return re.sub(r"[^ACGTNacgtn]", "", seq).upper()


# ---------------- sample header parsing ----------------
# Sample FASTA header: >HLA-A|class=I|chr... OR >C2|class=III|...
SAMPLE_HDR_RE = re.compile(r"^(?P<gene>[^|]+)(?:\|class=(?P<class>[^|]+))?.*$")

def parse_sample_header(h: str) -> Tuple[str, str]:
    m = SAMPLE_HDR_RE.match(h or "")
    gene = (m.group("gene") if m else (h.split("|", 1)[0] if h else "")).strip()
    mhc_class = (m.group("class") if (m and m.group("class")) else "").strip()
    return gene, mhc_class


def normalize_class_label(mhc_class: str) -> str:
    c = (mhc_class or "").strip().upper()
    if c == "I":
        return "Class I"
    if c == "II":
        return "Class II"
    if c == "III":
        return "Class III"
    return ""


# ---------------- reference header parsing (robust for IMGT/HLA) ----------------
def parse_ref_header(h: str) -> Tuple[str, str]:
    """
    Robust parser for IMGT/HLA FASTA headers.
    Many headers look like:
      HLA00001 A*01:01:01:01 ...
      HLA12345 DRB1*15:01:01 ...
    So we search the whole header for a locus*allele pattern.

    Returns:
      gene: HLA-A, HLA-DRB1, ...
      allele_id: HLA-A*01:01:01:01, ...
    """
    if not h:
        return "NA", "NA"

    loci = r"(A|B|C|E|F|G|DRA|DRB[0-9]+|DQA[0-9]+|DQB[0-9]+|DPA[0-9]+|DPB[0-9]+|DMA|DMB|DOA|DOB)"
    m = re.search(rf"(?<![A-Z0-9])(?P<gene>{loci})\*(?P<al>[0-9:]+)", h)
    if not m:
        m2 = re.search(rf"(?P<gene>HLA-{loci})\*(?P<al>[0-9:]+)", h)
        if not m2:
            return "NA", "NA"
        gene_part = m2.group("gene").replace("HLA-", "")
        al = m2.group("al")
    else:
        gene_part = m.group("gene")
        al = m.group("al")

    gene = f"HLA-{gene_part}"
    allele_id = f"{gene}*{al}"
    return gene, allele_id


# ---------------- Alignment (global Needleman-Wunsch) ----------------
def nw_global_align_score_and_identity(
    a: str, b: str,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -2
) -> Tuple[int, float]:
    a = clean_seq(a)
    b = clean_seq(b)
    n = len(a)
    m = len(b)
    if n == 0 or m == 0:
        return (0, 0.0)

    prev = [0] * (m + 1)
    curr = [0] * (m + 1)
    ptr = bytearray((n + 1) * (m + 1))  # 0 diag, 1 up, 2 left

    for j in range(1, m + 1):
        prev[j] = prev[j - 1] + gap
        ptr[j] = 2
    for i in range(1, n + 1):
        curr[0] = prev[0] + gap
        ptr[i * (m + 1)] = 1
        ai = a[i - 1]
        for j in range(1, m + 1):
            bj = b[j - 1]
            s_diag = prev[j - 1] + (match if ai == bj else mismatch)
            s_up = prev[j] + gap
            s_left = curr[j - 1] + gap

            best = s_diag
            direction = 0
            if s_up > best:
                best = s_up
                direction = 1
            if s_left > best:
                best = s_left
                direction = 2

            curr[j] = best
            ptr[i * (m + 1) + j] = direction
        prev, curr = curr, prev

    score = prev[m]

    # Traceback for identity
    i, j = n, m
    matches = 0
    aligned = 0
    while i > 0 or j > 0:
        direction = ptr[i * (m + 1) + j]
        if direction == 0:
            aligned += 1
            if a[i - 1] == b[j - 1]:
                matches += 1
            i -= 1
            j -= 1
        elif direction == 1:
            aligned += 1
            i -= 1
        else:
            aligned += 1
            j -= 1

    pid = (100.0 * matches / aligned) if aligned else 0.0
    return score, pid


# ---------------- Main typing ----------------
def infer_sample_id(sample_fasta: Path) -> str:
    name = sample_fasta.name
    if name.endswith(".mhc_step2.fa"):
        return name[:-len(".mhc_step2.fa")]
    return sample_fasta.stem


def main() -> int:
    ap = argparse.ArgumentParser(description="Task3 Step3: Pairwise alignment & allele typing (IMGT/HLA).")
    ap.add_argument("--sample_fasta", required=True)
    ap.add_argument("--ref_fasta", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--log", default=None)

    ap.add_argument("--match", type=int, default=2)
    ap.add_argument("--mismatch", type=int, default=-1)
    ap.add_argument("--gap", type=int, default=-2)

    args = ap.parse_args()

    sample_fa = Path(args.sample_fasta)
    ref_fa = Path(args.ref_fasta)
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    sample_id = infer_sample_id(sample_fa)

    # Load sample genes (we will output one row per gene record in sample FASTA)
    sample_records: List[Tuple[str, str, str]] = []  # (gene, class_label, seq)
    sample_genes_set = set()

    for h, s in parse_fasta(sample_fa):
        gene, mhc_class = parse_sample_header(h)
        class_label = normalize_class_label(mhc_class)
        seq = clean_seq(s)
        if not gene or not seq:
            continue
        sample_records.append((gene, class_label, seq))
        sample_genes_set.add(gene)

    # Load reference alleles, but ONLY for genes present in sample to save memory/time
    ref_by_gene: Dict[str, List[Tuple[str, str]]] = {}
    total_ref = 0
    genes_seen = set()

    for h, s in parse_fasta(ref_fa):
        gene, allele_id = parse_ref_header(h)
        if gene == "NA":
            continue
        if gene not in sample_genes_set:
            continue
        seq = clean_seq(s)
        if not seq:
            continue
        ref_by_gene.setdefault(gene, []).append((allele_id, seq))
        total_ref += 1
        genes_seen.add(gene)

    # Align and produce output rows
    out_rows: List[Dict[str, str]] = []

    # We can have multiple records per gene; we aggregate per gene_name:
    # - read_support = number of records for that gene in sample FASTA
    # - mean_alignment_score = mean of best scores across records
    # - percent_identity = mean of pids across records
    # - best_allele = allele from the record with highest score
    per_gene_acc: Dict[str, Dict[str, object]] = {}

    for gene, class_label, seq in sample_records:
        ref_list = ref_by_gene.get(gene, [])
        if not ref_list:
            acc = per_gene_acc.setdefault(gene, {
                "sample_id": sample_id,
                "gene_name": gene,
                "mhc_class": class_label,
                "best_allele": "",
                "read_support": 0,
                "scores": [],
                "pids": [],
                "best_score": None,
                "best_pid": None,
                "best_allele_id": "",
            })
            acc["read_support"] = int(acc["read_support"]) + 0  # keep 0 for no-ref genes
            continue

        best_score = None
        best_pid = None
        best_allele = ""

        for allele_id, ref_seq in ref_list:
            score, pid = nw_global_align_score_and_identity(
                seq, ref_seq,
                match=args.match, mismatch=args.mismatch, gap=args.gap
            )
            if (best_score is None) or (score > best_score) or (score == best_score and pid > (best_pid or -1)):
                best_score = score
                best_pid = pid
                best_allele = allele_id

        acc = per_gene_acc.setdefault(gene, {
            "sample_id": sample_id,
            "gene_name": gene,
            "mhc_class": class_label,
            "best_allele": "",
            "read_support": 0,
            "scores": [],
            "pids": [],
            "best_score": None,
            "best_pid": None,
            "best_allele_id": "",
        })
        acc["read_support"] = int(acc["read_support"]) + 1
        acc["scores"].append(best_score or 0)
        acc["pids"].append(best_pid or 0.0)

        # track global best allele per gene
        if acc["best_score"] is None or (best_score is not None and best_score > acc["best_score"]):
            acc["best_score"] = best_score
            acc["best_pid"] = best_pid
            acc["best_allele_id"] = best_allele

    # Build final rows in sample gene order (as they appeared)
    seen = set()
    for gene, class_label, _ in sample_records:
        if gene in seen:
            continue
        seen.add(gene)
        acc = per_gene_acc.get(gene)
        if not acc:
            # shouldn't happen
            out_rows.append({
                "sample_id": sample_id,
                "gene_name": gene,
                "mhc_class": class_label,
                "best_allele": "",
                "read_support": "0",
                "mean_alignment_score": "0",
                "percent_identity": "0",
            })
            continue

        rs = int(acc["read_support"])
        if rs == 0:
            out_rows.append({
                "sample_id": sample_id,
                "gene_name": gene,
                "mhc_class": acc["mhc_class"],
                "best_allele": "",
                "read_support": "0",
                "mean_alignment_score": "0",
                "percent_identity": "0",
            })
        else:
            scores = acc["scores"]
            pids = acc["pids"]
            mean_score = sum(scores) / len(scores) if scores else 0
            mean_pid = sum(pids) / len(pids) if pids else 0.0
            out_rows.append({
                "sample_id": sample_id,
                "gene_name": gene,
                "mhc_class": acc["mhc_class"],
                "best_allele": acc["best_allele_id"],
                "read_support": str(rs),
                "mean_alignment_score": f"{mean_score:.1f}",
                "percent_identity": f"{mean_pid:.6f}",
            })

    # Write CSV
    cols = ["sample_id", "gene_name", "mhc_class", "best_allele", "read_support", "mean_alignment_score", "percent_identity"]
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    if args.log:
        lp = Path(args.log)
        lp.parent.mkdir(parents=True, exist_ok=True)
        lp.write_text(
            f"sample_id={sample_id}\n"
            f"sample_records={len(sample_records)}\n"
            f"ref_loaded={total_ref}\n"
            f"ref_genes_matched={sorted(genes_seen)}\n"
        )

    print(f"[OK] Wrote: {out_csv}  (genes={len(out_rows)})")
    print(f"[INFO] Reference alleles loaded (filtered to sample genes): {total_ref}  (genes_matched={len(genes_seen)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
