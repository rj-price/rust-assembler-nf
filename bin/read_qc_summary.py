#!/usr/bin/env python3
"""
Combine per-run seqkit stats into a single read QC summary.

The point of this script is the *per-run comparison*: the two runs have very different
character (run1 ~11 kb / Q40, run2 ~27.5 kb / Q29) and that difference is a QC observation
to be surfaced, never a trigger to filter (DECISION D2).

Outputs:
  read_summary.tsv        one row per run, plus a combined row
  read_summary.json       same, machine-readable for downstream reporting
  read_summary_mqc.tsv    MultiQC custom-content table
"""

import argparse
import csv
import json
import sys
from pathlib import Path


def parse_seqkit(path):
    """seqkit stats --all --tabular -> dict of the fields we care about."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty seqkit stats file: {path}")
    return rows[0]


def num(row, *keys, default=0.0):
    """seqkit column names drift between versions, so accept several spellings."""
    for k in keys:
        if k in row and row[k] not in ("", "N/A", None):
            try:
                return float(str(row[k]).replace(",", ""))
            except ValueError:
                continue
    return default


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", nargs="+", required=True,
                    help="run_id:path pairs of seqkit stats TSVs")
    ap.add_argument("--outprefix", default="read_summary")
    args = ap.parse_args()

    records = []
    for item in args.stats:
        run_id, _, path = item.partition(":")
        if not path:
            sys.exit(f"expected run_id:path, got '{item}'")
        row = parse_seqkit(path)

        records.append({
            "run": run_id,
            "num_reads": int(num(row, "num_seqs")),
            "total_bases": int(num(row, "sum_len")),
            "total_gb": round(num(row, "sum_len") / 1e9, 2),
            "min_len": int(num(row, "min_len")),
            "mean_len": round(num(row, "avg_len"), 1),
            "median_len": round(num(row, "Q2"), 1),
            "max_len": int(num(row, "max_len")),
            "n50": int(num(row, "N50")),
            "avg_qual": round(num(row, "AvgQual"), 2),
            "pct_q20": round(num(row, "Q20(%)"), 2),
            "pct_q30": round(num(row, "Q30(%)"), 2),
            "gc_pct": round(num(row, "GC(%)"), 2),
        })

    records.sort(key=lambda r: r["run"])

    # Combined row. Length stats are weighted by yield; N50 cannot be pooled from summaries,
    # so it is left blank rather than fabricated from a wrong assumption.
    total_reads = sum(r["num_reads"] for r in records)
    total_bases = sum(r["total_bases"] for r in records)
    combined = {
        "run": "combined",
        "num_reads": total_reads,
        "total_bases": total_bases,
        "total_gb": round(total_bases / 1e9, 2),
        "min_len": min((r["min_len"] for r in records), default=0),
        "mean_len": round(total_bases / total_reads, 1) if total_reads else 0,
        "median_len": "",
        "max_len": max((r["max_len"] for r in records), default=0),
        "n50": "",
        "avg_qual": round(
            sum(r["avg_qual"] * r["total_bases"] for r in records) / total_bases, 2
        ) if total_bases else 0,
        "pct_q20": round(
            sum(r["pct_q20"] * r["total_bases"] for r in records) / total_bases, 2
        ) if total_bases else 0,
        "pct_q30": round(
            sum(r["pct_q30"] * r["total_bases"] for r in records) / total_bases, 2
        ) if total_bases else 0,
        "gc_pct": round(
            sum(r["gc_pct"] * r["total_bases"] for r in records) / total_bases, 2
        ) if total_bases else 0,
    }

    all_rows = records + [combined]
    fields = list(all_rows[0].keys())

    with open(f"{args.outprefix}.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(all_rows)

    with open(f"{args.outprefix}.json", "w") as fh:
        json.dump({"runs": records, "combined": combined}, fh, indent=2)

    # Flag the run-to-run divergence as an observation for the report.
    notes = []
    if len(records) > 1:
        q30 = [r["pct_q30"] for r in records]
        if max(q30) - min(q30) > 25:
            notes.append(
                "Runs differ markedly in base quality (Q30 range "
                f"{min(q30):.1f}%-{max(q30):.1f}%). This is recorded as an observation; "
                "no filtering is applied by default."
            )
        lens = [r["mean_len"] for r in records]
        if max(lens) > 1.5 * min(lens):
            notes.append(
                f"Runs differ markedly in read length (mean {min(lens):.0f}-{max(lens):.0f} bp). "
                "The longer reads are potentially valuable for repeat and haplotype resolution."
            )
    if notes:
        with open(f"{args.outprefix}_notes.txt", "w") as fh:
            fh.write("\n".join(notes) + "\n")
        for n in notes:
            print(f"NOTE: {n}", file=sys.stderr)

    with open(f"{args.outprefix}_mqc.tsv", "w") as fh:
        fh.write("# id: read_summary\n")
        fh.write("# section_name: 'HiFi read summary'\n")
        fh.write("# description: 'Per-run and combined HiFi read statistics. "
                 "Quality differences between runs are observations, not filters.'\n")
        fh.write("# format: 'tsv'\n")
        fh.write("# plot_type: 'table'\n")
        fh.write("\t".join(fields) + "\n")
        for r in all_rows:
            fh.write("\t".join(str(r[f]) for f in fields) + "\n")

    print(f"Wrote {args.outprefix}.tsv / .json / _mqc.tsv ({len(all_rows)} rows)")


if __name__ == "__main__":
    main()
