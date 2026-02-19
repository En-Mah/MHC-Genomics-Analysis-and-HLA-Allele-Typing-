#!/usr/bin/env python3
import os, glob, subprocess, csv
import numpy as np
import matplotlib.pyplot as plt

BAM_DIR = "project/bam"
MAPQ_DIR = "project/mapq"
PLOTS_DIR = "project/task2_plots"
QC_DIR = "project/qc"

os.makedirs(MAPQ_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(QC_DIR, exist_ok=True)

BINS = np.arange(-0.5, 61.5, 1.0)   # MAPQ usually 0..60; 1-unit bins
BIN_CENTERS = (BINS[:-1] + BINS[1:]) / 2

def extract_mapq_from_bam(bam_path: str, out_txt: str) -> None:
    # Write MAPQ column (5th field) for every read
    cmd = ["samtools", "view", bam_path]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.stdout is None:
        raise RuntimeError("samtools stdout unavailable")

    with open(out_txt, "w") as f:
        for line in p.stdout:
            parts = line.split("\t")
            if len(parts) >= 5:
                f.write(parts[4] + "\n")

    _, err = p.communicate()
    if p.returncode != 0:
        raise RuntimeError(f"samtools view failed: {bam_path}\n{err}")

def load_mapq(txt_path: str) -> np.ndarray:
    vals = []
    with open(txt_path) as f:
        for line in f:
            line = line.strip()
            if line:
                vals.append(int(line))
    return np.asarray(vals, dtype=int)

def kde_gaussian(x: np.ndarray, grid: np.ndarray) -> np.ndarray:
    x = x.astype(float)
    n = len(x)
    if n < 2:
        return np.zeros_like(grid, dtype=float)
    std = np.std(x, ddof=1)
    bw = 1.06 * std * (n ** (-1/5)) if std > 0 else 1.0
    bw = max(bw, 0.5)  

    diffs = (grid[:, None] - x[None, :]) / bw
    dens = np.mean(np.exp(-0.5 * diffs**2) / (bw * np.sqrt(2*np.pi)), axis=1)
    return dens

def plot_hist_counts(x: np.ndarray, sample: str, out_png: str) -> None:
    plt.figure()
    plt.hist(x, bins=BINS)
    plt.xlabel("Mapping Quality (MAPQ)")
    plt.ylabel("Number of Reads")
    plt.title(f"MAPQ Histogram (Counts) - {sample}")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

def plot_hist_density(x: np.ndarray, sample: str, out_png: str) -> None:
    plt.figure()
    plt.hist(x, bins=BINS, density=True)
    plt.xlabel("Mapping Quality (MAPQ)")
    plt.ylabel("Probability Density")
    plt.title(f"MAPQ Histogram (Normalized) - {sample}")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

def plot_kde(x: np.ndarray, sample: str, out_png: str) -> None:
    plt.figure()
    grid = np.linspace(0, 60, 601)
    dens = kde_gaussian(x, grid)
    plt.plot(grid, dens)
    plt.xlabel("Mapping Quality (MAPQ)")
    plt.ylabel("Density (KDE)")
    plt.title(f"MAPQ Density (KDE) - {sample}")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

def plot_cdf(x: np.ndarray, sample: str, out_png: str) -> None:
    plt.figure()
    xs = np.sort(x)
    ys = np.arange(1, len(xs) + 1) / len(xs)
    plt.plot(xs, ys)
    plt.xlabel("Mapping Quality (MAPQ)")
    plt.ylabel("Cumulative Fraction of Reads")
    plt.title(f"MAPQ CDF - {sample}")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

def main():
    bam_files = sorted(glob.glob(os.path.join(BAM_DIR, "*.bam")))
    if not bam_files:
        raise SystemExit(f"No BAM files found in {BAM_DIR}")

    all_samples = []
    all_mapq = []
    summary_rows = []

    for bam in bam_files:
        sample = os.path.basename(bam).replace(".bam", "")
        all_samples.append(sample)

        mapq_txt = os.path.join(MAPQ_DIR, f"{sample}.mapq.txt")
        if not os.path.exists(mapq_txt):
            extract_mapq_from_bam(bam, mapq_txt)

        x = load_mapq(mapq_txt)
        all_mapq.append(x)

        plot_hist_counts(x, sample, os.path.join(PLOTS_DIR, f"{sample}.hist_counts.png"))
        plot_hist_density(x, sample, os.path.join(PLOTS_DIR, f"{sample}.hist_norm.png"))
        plot_kde(x, sample, os.path.join(PLOTS_DIR, f"{sample}.kde.png"))
        plot_cdf(x, sample, os.path.join(PLOTS_DIR, f"{sample}.cdf.png"))

        n = len(x)
        pct0 = (np.sum(x == 0) / n * 100) if n else np.nan
        pct30 = (np.sum(x >= 30) / n * 100) if n else np.nan
        summary_rows.append({
            "sample": sample,
            "n_reads": n,
            "mean_mapq": float(np.mean(x)) if n else np.nan,
            "median_mapq": float(np.median(x)) if n else np.nan,
            "pct_mapq0": float(pct0),
            "pct_mapq_ge_30": float(pct30),
        })

    summary_csv = os.path.join(QC_DIR, "task2_mapq_summary.csv")
    with open(summary_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        w.writerows(summary_rows)

    # 1) Boxplot: MAPQ per sample
    plt.figure(figsize=(max(10, len(all_samples) * 0.25), 6))
    plt.boxplot(all_mapq, labels=all_samples, showfliers=False)
    plt.xticks(rotation=90)
    plt.xlabel("Sample")
    plt.ylabel("Mapping Quality (MAPQ)")
    plt.title("MAPQ Distribution Across Samples (Boxplot)")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "ALL.boxplot_mapq.png"), dpi=200)
    plt.close()

    # 2) Bar plot: % MAPQ=0
    plt.figure(figsize=(max(10, len(all_samples) * 0.25), 5))
    pct0s = [r["pct_mapq0"] for r in summary_rows]
    plt.bar(all_samples, pct0s)
    plt.xticks(rotation=90)
    plt.xlabel("Sample")
    plt.ylabel("% Reads with MAPQ = 0")
    plt.title("Fraction of Ambiguously Mapped Reads (MAPQ=0) per Sample")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "ALL.bar_pct_mapq0.png"), dpi=200)
    plt.close()

    # 3) Bar plot: % MAPQ>=30
    plt.figure(figsize=(max(10, len(all_samples) * 0.25), 5))
    pct30s = [r["pct_mapq_ge_30"] for r in summary_rows]
    plt.bar(all_samples, pct30s)
    plt.xticks(rotation=90)
    plt.xlabel("Sample")
    plt.ylabel("% Reads with MAPQ ≥ 30")
    plt.title("Fraction of Confidently Mapped Reads (MAPQ≥30) per Sample")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "ALL.bar_pct_mapq_ge30.png"), dpi=200)
    plt.close()

    # 4) Heatmap
    mat = []
    for x in all_mapq:
        counts, _ = np.histogram(x, bins=BINS)
        mat.append(counts)
    mat = np.asarray(mat)

    plt.figure(figsize=(12, max(5, len(all_samples) * 0.18)))
    plt.imshow(mat, aspect="auto")
    plt.yticks(np.arange(len(all_samples)), all_samples)
    plt.xticks(np.arange(0, 61, 5), [str(i) for i in range(0, 61, 5)])
    plt.xlabel("MAPQ Bin")
    plt.ylabel("Sample")
    plt.title("MAPQ Binned Counts Across Samples (Heatmap)")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, "ALL.heatmap_binned_counts.png"), dpi=200)
    plt.close()

    print("Done.")
    print(f"Per-sample plots + global plots saved to: {PLOTS_DIR}")
    print(f"Summary CSV saved to: {summary_csv}")

if __name__ == "__main__":
    main()
