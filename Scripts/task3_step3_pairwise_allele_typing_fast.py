#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

def parse_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    h=None; seq=[]
    with path.open() as f:
        for line in f:
            line=line.strip()
            if not line: 
                continue
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(seq)
                h=line[1:].strip()
                seq=[]
            else:
                seq.append(line)
    if h is not None:
        yield h, "".join(seq)

def clean_seq(s: str) -> str:
    return re.sub(r"[^ACGTNacgtn]", "", s).upper()

SAMPLE_HDR_RE = re.compile(r"^(?P<gene>[^|]+)(?:\|class=(?P<class>[^|]+))?.*$")
def parse_sample_header(h: str) -> Tuple[str, str]:
    m = SAMPLE_HDR_RE.match(h or "")
    gene = (m.group("gene") if m else (h.split("|",1)[0] if h else "")).strip()
    c = (m.group("class") if (m and m.group("class")) else "").strip().upper()
    if c == "I":  cls="Class I"
    elif c == "II": cls="Class II"
    elif c == "III": cls="Class III"
    else: cls=""
    return gene, cls

def parse_ref_header_imgt(h: str) -> Tuple[str, str]:
    # Find allele pattern like A*01:01:01:01 or DRB1*15:01:01 in full header
    loci = r"(A|B|C|E|F|G|DRA|DRB[0-9]+|DQA[0-9]+|DQB[0-9]+|DPA[0-9]+|DPB[0-9]+|DMA|DMB|DOA|DOB)"
    m = re.search(rf"(?<![A-Z0-9])(?P<gene>{loci})\*(?P<al>[0-9:]+)", h or "")
    if not m:
        return "NA","NA"
    gene = f"HLA-{m.group('gene')}"
    allele_id = f"{gene}*{m.group('al')}"
    return gene, allele_id

def nw_score_pid(a: str, b: str, match=2, mismatch=-1, gap=-2) -> Tuple[int, float]:
    a=clean_seq(a); b=clean_seq(b)
    n=len(a); m=len(b)
    if n==0 or m==0: return 0,0.0
    prev=[0]*(m+1); curr=[0]*(m+1)
    ptr=bytearray((n+1)*(m+1))
    for j in range(1,m+1):
        prev[j]=prev[j-1]+gap
        ptr[j]=2
    for i in range(1,n+1):
        curr[0]=prev[0]+gap
        ptr[i*(m+1)]=1
        ai=a[i-1]
        for j in range(1,m+1):
            bj=b[j-1]
            sd=prev[j-1]+(match if ai==bj else mismatch)
            su=prev[j]+gap
            sl=curr[j-1]+gap
            best=sd; d=0
            if su>best: best=su; d=1
            if sl>best: best=sl; d=2
            curr[j]=best
            ptr[i*(m+1)+j]=d
        prev,curr=curr,prev
    score=prev[m]
    i,j=n,m; matches=0; aligned=0
    while i>0 or j>0:
        d=ptr[i*(m+1)+j]
        if d==0:
            aligned+=1
            if a[i-1]==b[j-1]: matches+=1
            i-=1; j-=1
        elif d==1:
            aligned+=1; i-=1
        else:
            aligned+=1; j-=1
    pid=100.0*matches/aligned if aligned else 0.0
    return score,pid

def infer_sample_id(p: Path) -> str:
    n=p.name
    return n[:-len(".mhc_step2.fa")] if n.endswith(".mhc_step2.fa") else p.stem

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sample_fasta", required=True)
    ap.add_argument("--ref_dir", required=True, help="Directory containing locus fasta like A_nuc.fasta, DRB1_nuc.fasta ...")
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--log", default=None)
    ap.add_argument("--match", type=int, default=2)
    ap.add_argument("--mismatch", type=int, default=-1)
    ap.add_argument("--gap", type=int, default=-2)
    args=ap.parse_args()

    sample_fa=Path(args.sample_fasta)
    ref_dir=Path(args.ref_dir)
    out_csv=Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    sample_id=infer_sample_id(sample_fa)

    # Map sample gene -> locus file name
    def gene_to_file(gene: str) -> Path | None:
        if not gene.startswith("HLA-"): return None
        locus=gene.replace("HLA-","")
        return ref_dir / f"{locus}_nuc.fasta"

    # Read sample records
    sample_records=[]
    for h,s in parse_fasta(sample_fa):
        gene,cls=parse_sample_header(h)
        seq=clean_seq(s)
        if gene and seq:
            sample_records.append((gene,cls,seq))

    # Group sample records by gene
    by_gene: Dict[str, List[Tuple[str,str]]] = {}
    for gene,cls,seq in sample_records:
        by_gene.setdefault(gene, []).append((cls,seq))

    # Cache reference alleles per gene
    ref_cache: Dict[str, List[Tuple[str,str]]] = {}

    out_rows=[]
    for gene, recs in by_gene.items():
        cls = recs[0][0] if recs else ""
        ref_path = gene_to_file(gene)
        if ref_path is None or not ref_path.exists():
            out_rows.append({
                "sample_id": sample_id,
                "gene_name": gene,
                "mhc_class": cls,
                "best_allele": "",
                "read_support": "0",
                "mean_alignment_score": "0",
                "percent_identity": "0",
            })
            continue

        if gene not in ref_cache:
            alleles=[]
            for h,s in parse_fasta(ref_path):
                g,aid=parse_ref_header_imgt(h)
                if g != gene:  # safety
                    continue
                seq=clean_seq(s)
                if seq:
                    alleles.append((aid,seq))
            ref_cache[gene]=alleles

        ref_list=ref_cache[gene]
        if not ref_list:
            out_rows.append({
                "sample_id": sample_id,
                "gene_name": gene,
                "mhc_class": cls,
                "best_allele": "",
                "read_support": "0",
                "mean_alignment_score": "0",
                "percent_identity": "0",
            })
            continue

        scores=[]; pids=[]
        best_allele=""; best_score=None; best_pid=None

        # Usually 1 record per gene, but keep general
        for _, seq in recs:
            L=len(seq); lo=int(L*0.90); hi=int(L*1.10)
            filtered=[(aid,rs) for aid,rs in ref_list if lo <= len(rs) <= hi]
            if not filtered:
                filtered=ref_list

            local_best_score=None; local_best_pid=None; local_best_allele=""
            for aid, rseq in filtered:
                sc,pid = nw_score_pid(seq, rseq, match=args.match, mismatch=args.mismatch, gap=args.gap)
                if local_best_score is None or sc>local_best_score or (sc==local_best_score and pid>(local_best_pid or -1)):
                    local_best_score=sc; local_best_pid=pid; local_best_allele=aid
            scores.append(local_best_score or 0)
            pids.append(local_best_pid or 0.0)
            if best_score is None or (local_best_score is not None and local_best_score > best_score):
                best_score=local_best_score; best_pid=local_best_pid; best_allele=local_best_allele

        mean_score = sum(scores)/len(scores) if scores else 0
        mean_pid = sum(pids)/len(pids) if pids else 0.0
        out_rows.append({
            "sample_id": sample_id,
            "gene_name": gene,
            "mhc_class": cls,
            "best_allele": best_allele,
            "read_support": str(len(recs)),
            "mean_alignment_score": f"{mean_score:.1f}",
            "percent_identity": f"{mean_pid:.6f}",
        })

    cols=["sample_id","gene_name","mhc_class","best_allele","read_support","mean_alignment_score","percent_identity"]
    with out_csv.open("w", newline="") as f:
        w=csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    if args.log:
        Path(args.log).write_text(f"sample_id={sample_id}\nGenes={len(out_rows)}\n")

    print(f"[OK] Wrote: {out_csv} (genes={len(out_rows)})")

if __name__=="__main__":
    raise SystemExit(main())
