#!/usr/bin/env python3
"""
Estimate coverage modes from a read-to-assembly depth profile.

For a dikaryon we expect roughly two modes: haplotype-specific sequence at ~1x the
per-haplotype depth, and collapsed/shared sequence at ~2x it. With 60.4 Gb of reads over a
a dikaryon they land near 1x and 2x the expected per-haplotype depth, which the caller
derives from --read_yield_bases and --genome_size.

Those numbers are EXPECTATIONS, not thresholds. The modes are estimated empirically from the
data and only then compared against expectation, so an unexpected result shows up as a
disagreement to investigate rather than being silently forced into the expected shape.

A substantial third mode well below the haplotype peak is the signature worth investigating
as contaminant or under-represented sequence.
"""

import argparse
import gzip
import json
import sys
from collections import Counter


def load_depths(path, max_depth):
    """Histogram of binned per-base depths."""
    hist = Counter()
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                d = int(parts[2])
            except ValueError:
                continue
            if 0 < d <= max_depth:
                hist[d] += 1
    return hist


def smooth(hist, max_depth, window=5):
    """Simple moving average, so local jitter doesn't register as a peak."""
    counts = [hist.get(d, 0) for d in range(max_depth + 1)]
    out = []
    for i in range(len(counts)):
        lo, hi = max(0, i - window), min(len(counts), i + window + 1)
        out.append(sum(counts[lo:hi]) / (hi - lo))
    return out


def find_peaks(sm, min_depth, prominence_frac=0.15):
    """Local maxima that are a meaningful fraction of the global maximum."""
    peaks = []
    peak_max = max(sm[min_depth:]) if len(sm) > min_depth else 0
    if peak_max <= 0:
        return peaks
    for d in range(min_depth + 1, len(sm) - 1):
        if sm[d] >= sm[d - 1] and sm[d] > sm[d + 1] and sm[d] >= prominence_frac * peak_max:
            # Collapse near-adjacent detections into one peak.
            if peaks and d - peaks[-1][0] < 5:
                if sm[d] > peaks[-1][1]:
                    peaks[-1] = (d, sm[d])
            else:
                peaks.append((d, sm[d]))
    return peaks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depth", required=True, help="binned depth tsv(.gz): contig, pos, depth")
    ap.add_argument("--coverage", help="samtools coverage output, for per-contig breadth")
    ap.add_argument("--assembly-id", required=True)
    ap.add_argument("--expected-haplotype-cov", type=float, default=None,
                    help="expected per-haplotype depth; omit to report modes without "
                         "scoring them against an expectation")
    ap.add_argument("--max-depth", type=int, default=500)
    ap.add_argument("--outprefix", default="coverage_summary")
    ap.add_argument("--histogram", help="write the depth histogram as TSV (depth, "
                                        "positions, pct_of_positions)")
    args = ap.parse_args()

    hist = load_depths(args.depth, args.max_depth)
    if not hist:
        print("WARNING: no depth data parsed", file=sys.stderr)
        result = {"assembly_id": args.assembly_id, "error": "no depth data"}
        with open(f"{args.outprefix}.json", "w") as fh:
            json.dump(result, fh, indent=2)
        return

    total_positions = sum(hist.values())
    mean_depth = sum(d * n for d, n in hist.items()) / total_positions

    sm = smooth(hist, args.max_depth)
    peaks = find_peaks(sm, min_depth=3)
    peaks.sort(key=lambda p: p[1], reverse=True)
    top = sorted(d for d, _ in peaks[:4])

    exp_hap = args.expected_haplotype_cov
    exp_col = exp_hap * 2 if exp_hap else None

    interpretation = []
    hap_peak = None

    # No expectation supplied: report what was measured and stop there. A default figure
    # would be worse than nothing -- it reads as a real comparison while being a guess.
    if exp_hap is None:
        interpretation.append(
            "Detected modes: %s. No expected depth available (the sample's yield is "
            "normally measured by READ_QC; --qc_only skips it, in which case set "
            "--read_yield_bases), so these are reported unscored."
            % (", ".join("%sx" % d for d in top) if top else "none")
        )

    for d in top if exp_hap else []:
        if 0.6 * exp_hap <= d <= 1.4 * exp_hap:
            hap_peak = d
            interpretation.append(
                f"{d}x is consistent with haplotype-specific sequence (expected ~{exp_hap:.0f}x)"
            )
        elif 0.6 * exp_col <= d <= 1.4 * exp_col:
            interpretation.append(
                f"{d}x is consistent with collapsed/shared sequence (expected ~{exp_col:.0f}x)"
            )
        elif d < 0.5 * exp_hap:
            interpretation.append(
                f"{d}x is well below the expected haplotype depth — investigate as "
                f"contaminant or under-represented sequence"
            )
        else:
            interpretation.append(f"{d}x does not match either expectation")

    if not peaks:
        interpretation.append("No clear coverage modes detected.")
    elif exp_col and hap_peak and not any(0.6 * exp_col <= d <= 1.4 * exp_col for d in top):
        interpretation.append(
            "No collapsed-sequence mode detected, which is what a well-separated "
            "(non-collapsed) dikaryotic assembly should look like."
        )

    result = {
        "assembly_id": args.assembly_id,
        "mean_depth": round(mean_depth, 2),
        "positions_sampled": total_positions,
        "detected_modes": top,
        "expected_haplotype_cov": exp_hap,
        "expected_collapsed_cov": exp_col,
        "interpretation": interpretation,
    }

    if args.coverage:
        try:
            with open(args.coverage) as fh:
                rows = [l.split("\t") for l in fh if not l.startswith("#")]
            rows = [r for r in rows if len(r) >= 6]
            if rows:
                breadth = [float(r[5]) for r in rows]
                result["n_contigs"] = len(rows)
                result["mean_breadth_pct"] = round(sum(breadth) / len(breadth), 2)
        except (OSError, ValueError, IndexError) as e:
            print(f"WARNING: could not parse coverage file: {e}", file=sys.stderr)

    with open(f"{args.outprefix}.json", "w") as fh:
        json.dump(result, fh, indent=2)

    # The histogram itself, not just the modes extracted from it. A summary says "modes
    # at 57x"; the curve shows whether that is a clean unimodal peak or a shoulder on something
    # else, which is the difference between a retained dikaryon and a collapsed one.
    if args.histogram:
        with open(args.histogram, "w") as fh:
            fh.write("depth\tpositions\tpct_of_positions\n")
            for d in range(1, args.max_depth + 1):
                n = hist.get(d, 0)
                fh.write("%d\t%d\t%.6f\n" % (d, n, 100.0 * n / total_positions))

    with open(f"{args.outprefix}.txt", "w") as fh:
        fh.write(f"Assembly: {args.assembly_id}\n")
        fh.write(f"Mean depth: {result['mean_depth']}x\n")
        fh.write(f"Detected coverage modes: {', '.join(str(d) + 'x' for d in top) or 'none'}\n")
        if exp_hap:
            fh.write(f"Expected: ~{exp_hap:.0f}x haplotype-specific, "
                     f"~{exp_col:.0f}x collapsed\n\n")
        else:
            fh.write("Expected: not computed (no --read_yield_bases)\n\n")
        for line in interpretation:
            fh.write(f"  - {line}\n")

    print(f"Coverage modes for {args.assembly_id}: {top} (mean {result['mean_depth']}x)")


if __name__ == "__main__":
    main()
