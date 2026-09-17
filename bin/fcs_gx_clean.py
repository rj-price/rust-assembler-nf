#!/usr/bin/env python3
"""
Act on an FCS-GX contamination report: write a cleaned assembly, and say exactly what was
removed.

Every other contamination step in this pipeline reports and stops. This one is the single
place that acts, so it is deliberately conservative and fully auditable:

  * the input assembly is never modified — a NEW FASTA is written (DECISION D6);
  * EXCLUDE contigs are dropped whole;
  * TRIM and FIX regions are cut out of their contig, and the surviving flanks are kept as
    separate records suffixed `:cleaned-<n>`, rather than being joined across the excised
    region (joining would invent a junction that no read supports);
  * REVIEW and INFO are left alone. They are FCS-GX saying "look at this", not "this is
    contamination", and a dikaryon has enough genuinely odd sequence that acting on them
    automatically would remove real biology;
  * a removal manifest is written alongside, so any contig that disappeared can be traced to
    the report line and taxon that caused it.

The FCS-GX action report is TSV with a `#` header line and the columns:

    seq_id  start_pos  end_pos  seq_len  action  div  agg_cont_cov  top_tax_name

Coordinates in that report are 1-based inclusive.
"""

import argparse
import gzip
import io
import json
import os
import sys

# Acted on. Anything else in the report is recorded and left in place.
DROP_ACTIONS = ("EXCLUDE",)
CUT_ACTIONS = ("TRIM", "FIX")


def open_text(path):
    """Open a FASTA whether or not it is gzipped."""
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r")


def read_fasta(path):
    """Yield (name, description_line, sequence). Kept streaming: assemblies run to gigabases."""
    name, header, chunks = None, None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, header, "".join(chunks)
                header = line[1:]
                name = header.split()[0] if header else ""
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield name, header, "".join(chunks)


def parse_report(path):
    """
    Return {seq_id: [ {action, start, end, taxon}, ... ]} from an FCS-GX action report.

    Tolerant of a missing or extra column: FCS-GX's report has gained columns between
    releases, and the first five are the ones that matter here.
    """
    per_seq = {}
    if not path or not os.path.exists(path):
        return per_seq

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 5:
                continue
            seq_id, start, end, _seq_len, action = f[0], f[1], f[2], f[3], f[4].strip().upper()
            taxon = f[7] if len(f) > 7 else ""
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            per_seq.setdefault(seq_id, []).append(
                {"action": action, "start": start_i, "end": end_i, "taxon": taxon}
            )
    return per_seq


def merge_intervals(intervals):
    """Merge overlapping 1-based inclusive (start, end) pairs."""
    if not intervals:
        return []
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def surviving_segments(seq_len, cut_intervals):
    """Complement of the cut intervals within 1..seq_len, as 1-based inclusive pairs."""
    segments = []
    cursor = 1
    for start, end in cut_intervals:
        if start > cursor:
            segments.append((cursor, start - 1))
        cursor = max(cursor, end + 1)
    if cursor <= seq_len:
        segments.append((cursor, seq_len))
    return segments


def write_wrapped(fh, seq, width=60):
    for i in range(0, len(seq), width):
        fh.write(seq[i : i + width] + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True, help="assembly FASTA (.fa/.fasta, optionally .gz)")
    ap.add_argument("--report", help="FCS-GX action report TSV; absent means nothing to do")
    ap.add_argument("--assembly-id", required=True)
    ap.add_argument("--out-fasta", required=True, help="cleaned FASTA (gzipped)")
    ap.add_argument("--out-removed", help="FASTA of the removed sequence, gzipped")
    ap.add_argument("--out-manifest", required=True, help="TSV manifest of every action taken")
    ap.add_argument("--out-json", required=True, help="machine-readable summary")
    ap.add_argument(
        "--min-segment",
        type=int,
        default=1000,
        help="drop a surviving flank shorter than this after a trim (default 1000)",
    )
    args = ap.parse_args()

    per_seq = parse_report(args.report)

    stats = {
        "assembly_id": args.assembly_id,
        "contigs_in": 0,
        "contigs_out": 0,
        "contigs_excluded": 0,
        "contigs_trimmed": 0,
        "bases_in": 0,
        "bases_out": 0,
        "bases_removed": 0,
        "actions_seen": {},
        "report_present": bool(per_seq),
    }

    manifest = open(args.out_manifest, "w")
    manifest.write("seq_id\taction\tstart\tend\tbases\ttop_tax_name\toutcome\n")

    out = gzip.open(args.out_fasta, "wt")
    removed = gzip.open(args.out_removed, "wt") if args.out_removed else None

    for name, header, seq in read_fasta(args.fasta):
        seq_len = len(seq)
        stats["contigs_in"] += 1
        stats["bases_in"] += seq_len

        hits = per_seq.get(name, [])
        for h in hits:
            stats["actions_seen"][h["action"]] = stats["actions_seen"].get(h["action"], 0) + 1

        drop = [h for h in hits if h["action"] in DROP_ACTIONS]
        cuts = [h for h in hits if h["action"] in CUT_ACTIONS]

        # Whole-contig exclusion wins: if FCS-GX wants the contig gone, trimming it is moot.
        if drop:
            stats["contigs_excluded"] += 1
            stats["bases_removed"] += seq_len
            for h in drop:
                manifest.write(
                    "%s\t%s\t%d\t%d\t%d\t%s\t%s\n"
                    % (name, h["action"], h["start"], h["end"], seq_len, h["taxon"], "contig_dropped")
                )
            if removed:
                removed.write(">%s\n" % header)
                write_wrapped(removed, seq)
            continue

        if not cuts:
            stats["contigs_out"] += 1
            stats["bases_out"] += seq_len
            out.write(">%s\n" % header)
            write_wrapped(out, seq)
            continue

        # Trim/fix: excise the flagged regions and keep the flanks as separate records.
        stats["contigs_trimmed"] += 1
        merged = merge_intervals([[h["start"], h["end"]] for h in cuts])
        for h in cuts:
            manifest.write(
                "%s\t%s\t%d\t%d\t%d\t%s\t%s\n"
                % (
                    name,
                    h["action"],
                    h["start"],
                    h["end"],
                    h["end"] - h["start"] + 1,
                    h["taxon"],
                    "region_excised",
                )
            )
            if removed:
                removed.write(">%s:%d-%d %s\n" % (name, h["start"], h["end"], h["action"]))
                write_wrapped(removed, seq[h["start"] - 1 : h["end"]])

        segments = surviving_segments(seq_len, merged)
        kept = [s for s in segments if (s[1] - s[0] + 1) >= args.min_segment]
        dropped_short = [s for s in segments if s not in kept]

        for start, end in dropped_short:
            manifest.write(
                "%s\t%s\t%d\t%d\t%d\t%s\t%s\n"
                % (name, "SHORT_FLANK", start, end, end - start + 1, "", "flank_dropped")
            )

        stats["bases_removed"] += seq_len - sum(e - s + 1 for s, e in kept)

        # A single surviving segment spanning the whole contig keeps the original name, so
        # that a contig which lost only a terminal sliver stays traceable across the summary.
        for idx, (start, end) in enumerate(kept, start=1):
            sub = seq[start - 1 : end]
            new_name = name if len(kept) == 1 and start == 1 else "%s:cleaned-%d" % (name, idx)
            out.write(">%s len=%d src=%s:%d-%d\n" % (new_name, len(sub), name, start, end))
            write_wrapped(out, sub)
            stats["contigs_out"] += 1
            stats["bases_out"] += len(sub)

    out.close()
    if removed:
        removed.close()
    manifest.close()

    stats["pct_bases_removed"] = (
        round(100.0 * stats["bases_removed"] / stats["bases_in"], 4) if stats["bases_in"] else 0.0
    )

    with open(args.out_json, "w") as fh:
        json.dump(stats, fh, indent=2, sort_keys=True)
        fh.write("\n")

    sys.stderr.write(
        "FCS-GX clean %s: %d/%d contigs kept, %d excluded, %d trimmed, %.4f%% of bases removed\n"
        % (
            args.assembly_id,
            stats["contigs_out"],
            stats["contigs_in"],
            stats["contigs_excluded"],
            stats["contigs_trimmed"],
            stats["pct_bases_removed"],
        )
    )

    if not stats["report_present"]:
        sys.stderr.write(
            "NOTE: no FCS-GX findings for this assembly; the cleaned FASTA is a copy of the input.\n"
        )


if __name__ == "__main__":
    main()
