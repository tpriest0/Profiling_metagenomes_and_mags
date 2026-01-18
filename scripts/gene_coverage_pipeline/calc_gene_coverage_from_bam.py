#!/usr/bin/env python3
import argparse
import subprocess
import sys
from collections import Counter, defaultdict
import pysam
import pandas as pd
import numpy as np


def median_nonzero_from_hist(depth_counts: Counter, covered_positions: int) -> float:
    """Median depth across covered bases only (depth > 0)."""
    if covered_positions <= 0:
        return 0.0

    k = (covered_positions - 1) // 2  # 0-based index
    running = 0
    for d in sorted(depth_counts):
        running += depth_counts[d]
        if running > k:
            return float(d)
    return 0.0


def percentile_from_hist(depth_counts: Counter, n: int, p: float) -> float:
    """Return the pth percentile from a depth histogram (covered bases only)."""
    if n <= 0:
        return 0.0

    target = int(np.ceil((p / 100.0) * n))
    running = 0
    for d in sorted(depth_counts):
        running += depth_counts[d]
        if running >= target:
            return float(d)
    return float(max(depth_counts)) if depth_counts else 0.0


def winsorized_mean_from_hist(depth_counts: Counter, length: int, covered_positions: int, p: float) -> float:
    """
    Winsorized mean depth across ALL bases (zeros included).
    Depths above the pth percentile (computed on covered bases) are capped.
    """
    if length <= 0:
        return np.nan
    if covered_positions <= 0:
        return 0.0

    cap = percentile_from_hist(depth_counts, covered_positions, p)

    clipped_sum = 0.0
    for d, c in depth_counts.items():
        clipped_sum += min(d, cap) * c

    return clipped_sum / length


def main():
    #####
    # Command-line arguments
    #####
    parser = argparse.ArgumentParser(
        description="Fast per-gene coverage stats using samtools depth streaming."
    )
    parser.add_argument("-i", "--in-bam", required=True, help="Input BAM (indexed)")
    parser.add_argument("-o", "--out-coverage", required=True, help="Output TSV")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Threads for samtools depth")
    parser.add_argument("-r", "--read-length", type=int, default=150, help="Read length (reported only)")
    parser.add_argument("--hcov-threshold", type=float, default=0.5,
                        help="Horizontal coverage threshold (fraction, default=0.5)")
    parser.add_argument("--winsor-p", type=float, default=95.0,
                        help="Percentile for winsorised mean (default=95)")
    args = parser.parse_args()

    #####
    # Load references and lengths from BAM
    #####
    with pysam.AlignmentFile(args.in_bam, "rb") as bam:
        refs = list(bam.references)
        lengths = dict(zip(bam.references, bam.lengths))

    #####
    # Aggregators
    #####
    covered_positions = defaultdict(int)
    total_depth = defaultdict(int)
    depth_hists = defaultdict(Counter)

    #####
    # Stream depth from samtools
    #####
    cmd = ["samtools", "depth", "-d", "0", "-a", args.in_bam]

    try:
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, bufsize=1
        )
    except FileNotFoundError:
        print("[ERROR] samtools not found in PATH.", file=sys.stderr)
        sys.exit(2)

    for line in proc.stdout:
        ref, _pos, depth_str = line.strip().split("\t")
        d = int(depth_str)
        if d > 0:
            covered_positions[ref] += 1
            total_depth[ref] += d
            depth_hists[ref][d] += 1

    stderr = proc.stderr.read()
    if proc.wait() != 0:
        print(stderr, file=sys.stderr)
        raise RuntimeError("samtools depth failed")

    #####
    # Build output table
    #####
    rows = []
    for ref in refs:
        L = lengths.get(ref, 0)
        cov_pos = covered_positions.get(ref, 0)
        tdepth = total_depth.get(ref, 0)
        hist = depth_hists.get(ref, Counter())

        horiz_frac = cov_pos / L if L > 0 else np.nan
        horiz_pct = horiz_frac * 100.0 if L > 0 else np.nan
        scale = 1000.0 / L if L > 0 else np.nan

        # Defaults
        mean_p95 = 0.0
        median_nz = 0.0

        if L > 0 and horiz_frac >= args.hcov_threshold:
            mean_p95 = winsorized_mean_from_hist(hist, L, cov_pos, args.winsor_p)
            median_nz = median_nonzero_from_hist(hist, cov_pos)

        rows.append({
            "Gene_name": ref,
            "Gene_length": L,
            "Horizontal_coverage": horiz_pct,
            "Total_bases_aligned": tdepth,
            "Bases_per_kbp": tdepth * scale if L > 0 else np.nan,

            # Requested output names
            "Mean_coverage_p95": mean_p95,
            "Mean_coverage_p95_per_kbp": mean_p95 * scale if L > 0 else np.nan,
            "Median_coverage": median_nz,
            "Median_coverage_per_kbp": median_nz * scale if L > 0 else np.nan,
        })

    pd.DataFrame(rows).to_csv(args.out_coverage, sep="\t", index=False)
    print(
        f"[INFO] Wrote coverage table: {args.out_coverage} ({len(rows)} genes)",
        file=sys.stderr
    )


if __name__ == "__main__":
    main()