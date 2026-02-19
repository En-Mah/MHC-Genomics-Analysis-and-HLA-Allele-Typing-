#!/usr/bin/env python3
import os
import sys
import glob
import zlib
import warnings
from collections import defaultdict, Counter

from Bio import SeqIO
from Bio import BiopythonDeprecationWarning
from Bio import pairwise2

SAMPLE_FASTA_GLOB = "project/task3/mhc_fasta/*.mhc_genes.fa"
REF_FASTA = "project/task3/ref/hla_gen.fasta"
OUT_DIR = "project/task3/allele_calls"
os.makedirs(OUT_DIR, exist_ok=True)

# Silence Biopython deprecation warnings (pairwise2 is deprecated but works reliably)
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)

def make_sketch(seq: str, k: int = 17, step: int = 25, max_items: int = 4000):
    """Fast k-mer sketch for prefiltering."""
    seq = seq.upper()
    n = len(seq)
    out = set()
    if n < k:
        return out
    for i in range(0, n - k + 1, step):
        out.add(zlib.crc32(seq[i:i+k].encode("ascii", "ignore")))
        if len(out) >= max_items:
            break
    return out

def gene_token_from_sample(gene_name: str) -> str:
    """Convert sample gene name like 'HLA-A' -> 'A'."""
    gene_name = gene_name.strip()
    return gene_name[4:] if gene_name.startswith("HLA-") else gene_name

def parse_sample_header(header: str):
    """
    Sample FASTA header example:
      >HLA-A|class=I|chr6:...
    Returns (gene_name, mhc_class).
    """
    parts = header.split("|")
    gene = parts[0].replace(">", "").strip()
    mhc_class = "NA"
    for p in parts:
        if p.startswith("class="):
            mhc_class = p.split("=", 1)[1].strip()
    return gene, mhc_class

def collect_needed_tokens():
    """Scan sample FASTAs to learn which HLA tokens we must support."""
    needed = set()
    for fa in glob.glob(SAMPLE_FASTA_GLOB):
        for rec in SeqIO.parse(fa, "fasta"):
            gene, _ = parse_sample_header(rec.description)
            if gene.startswith("HLA-"):
                needed.add(gene_token_from_sample(gene))
    return needed

def parse_token_from_ref_header(desc: str):
    """
    IMGT/HLA hla_gen.fasta header typically:
      >HLA:HLA00001 A*01:01:01:01 3503 bp
      >HLA:HLAxxxxx DRB1*15:01:01:01 ...
    Pick first whitespace token containing '*' and ':'.
    """
    for t in desc.split():
        if "*" in t and ":" in t:
            return t.strip()
    return None

def load_ref_by_token_manual(ref_fa: str, needed_tokens: set):
    """
    Manual FASTA parsing (robust, matches what grep sees).
    Store only tokens in needed_tokens -> faster, less RAM.
    Returns:
      ref[token] = list of (allele_name, seq, sketch)
      counts[token] = count
    """
    ref = defaultdict(list)
    counts = Counter()

    header = None
    token = None
    allele_name = None
    seq_chunks = []

    def flush_current():
        nonlocal token, allele_name, seq_chunks
        if token and token in needed_tokens and seq_chunks:
            seq = "".join(seq_chunks).upper()
            ref[token].append((allele_name, seq, make_sketch(seq)))
            counts[token] += 1
        seq_chunks = []

    with open(ref_fa, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # allow comment lines before the first record
            if header is None and line[0] in (";", "#", "!"):
                continue

            if line.startswith(">"):
                if header is not None:
                    flush_current()

                header = line[1:]
                token = None
                allele_name = None
                seq_chunks = []

                allele_tok = parse_token_from_ref_header(header)
                if allele_tok and "*" in allele_tok:
                    gene_part, rest = allele_tok.split("*", 1)
                    token = gene_part.replace("HLA-", "").strip()
                    allele_name = f"HLA-{token}*{rest.strip()}"

            else:
                if token and token in needed_tokens:
                    seq_chunks.append(line)

    if header is not None:
        flush_current()

    return ref, counts

def percent_identity_from_pairwise2(alnA: str, alnB: str) -> float:
    """Percent identity ignoring gap positions."""
    matches = 0
    aligned = 0
    for a, b in zip(alnA, alnB):
        if a == "-" or b == "-":
            continue
        aligned += 1
        if a == b:
            matches += 1
    return 0.0 if aligned == 0 else 100.0 * matches / aligned

def main():
    if not os.path.exists(REF_FASTA):
        print(f"ERROR: Missing reference FASTA: {REF_FASTA}", file=sys.stderr)
        sys.exit(1)

    needed = collect_needed_tokens()
    if not needed:
        print(f"ERROR: No needed HLA tokens found from sample FASTAs at {SAMPLE_FASTA_GLOB}", file=sys.stderr)
        sys.exit(1)

    ref, counts = load_ref_by_token_manual(REF_FASTA, needed)

    # Debug summary
    for t in ["A","B","C","DRB1","DQB1","DQA1","DPB1","DPA1","DRA","E","F","G"]:
        if t in needed:
            print(f"[REF] {t}: {counts.get(t,0)} alleles (needed)", file=sys.stderr)

    # Scoring scheme (same idea as before)
    match = 2
    mismatch = -1
    gap_open = -5
    gap_extend = -1

    sample_files = sorted(glob.glob(SAMPLE_FASTA_GLOB))
    if not sample_files:
        print(f"ERROR: No sample FASTAs found at {SAMPLE_FASTA_GLOB}", file=sys.stderr)
        sys.exit(1)

    for fa in sample_files:
        sample = os.path.basename(fa).replace(".mhc_genes.fa", "")
        out_tsv = os.path.join(OUT_DIR, f"{sample}.alleles.tsv")

        with open(out_tsv, "w") as out:
            out.write("gene_name\tmhc_class\tbest_matching_allele\talignment_score\tpercent_identity\tref_candidates\tprefilter_used\n")

            for rec in SeqIO.parse(fa, "fasta"):
                gene_name, mhc_class = parse_sample_header(rec.description)
                tok = gene_token_from_sample(gene_name)
                seq = str(rec.seq).upper()

                if tok not in ref:
                    out.write(f"{gene_name}\t{mhc_class}\tNA\tNA\tNA\t0\t0\n")
                    continue

                candidates = ref[tok]
                n_candidates = len(candidates)

                # Prefilter top 50 by sketch overlap
                sample_sk = make_sketch(seq)
                scored = []
                for allele_name, refseq, refsk in candidates:
                    scored.append((len(sample_sk & refsk), allele_name, refseq))
                scored.sort(reverse=True, key=lambda x: x[0])
                top = scored[:50]
                prefilter_used = len(top)

                # 1) score-only over top candidates (fast, no alignment explosion)
                best = None  # (score, allele_name, refseq)
                for overlap, allele_name, refseq in top:
                    score = pairwise2.align.globalms(
                        seq, refseq,
                        match, mismatch,
                        gap_open, gap_extend,
                        score_only=True
                    )
                    if best is None or score > best[0]:
                        best = (score, allele_name, refseq)

                if best is None:
                    out.write(f"{gene_name}\t{mhc_class}\tNA\tNA\tNA\t{n_candidates}\t{prefilter_used}\n")
                    continue

                best_score, best_allele, best_refseq = best

                # 2) compute percent identity only for the best candidate (one alignment only)
                aln = pairwise2.align.globalms(
                    seq, best_refseq,
                    match, mismatch,
                    gap_open, gap_extend,
                    one_alignment_only=True
                )[0]
                pid = percent_identity_from_pairwise2(aln.seqA, aln.seqB)

                out.write(f"{gene_name}\t{mhc_class}\t{best_allele}\t{best_score:.2f}\t{pid:.2f}\t{n_candidates}\t{prefilter_used}\n")

        print(f"[OK] {sample} -> {out_tsv}")

if __name__ == "__main__":
    main()
