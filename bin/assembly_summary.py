#!/usr/bin/env python3
"""
Build the machine-readable assembly comparison.

Two modes:
  record  — distil one assembly's QC outputs into a single JSON record
  merge   — combine those records into assembly_summary.tsv / .json

Design notes:
  * Assemblies are NEVER ranked by N50 (DECISION D7). N50 is reported, and that is all.
  * `size_flag` scores assembly size against the expected 0.8-1.0 Gb dikaryon. A ~450 Mb
    result is flagged `collapsed` — the specific failure this project exists to avoid — and
    an oversized one suggests retained contamination or duplication artefacts.
  * The BUSCO breakdown is kept intact. High duplication is EXPECTED in a correctly phased
    dikaryotic assembly, so it is reported without judgement rather than penalised.
"""

import argparse
import csv
import json
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------- parsers

def parse_gfastats(path):
    """
    gfastats emits 'Key: value' lines.

    Keys are matched case-insensitively: gfastats capitalises them ('Total scaffold length',
    'Scaffold N50') and the capitalisation has varied between versions. A silent mismatch
    here would leave assembly_size unset and every size_flag stuck on 'unknown', so the
    lookup must not depend on getting the case exactly right.
    """
    out = {}
    if not path or not Path(path).exists():
        return out
    key_map = {
        "total scaffold length": "assembly_size",
        "total contig length": "total_contig_length",
        "# scaffolds": "n_scaffolds",
        "# contigs": "n_contigs",
        "scaffold n50": "n50",
        "contig n50": "contig_n50",
        "scaffold l50": "l50",
        "contig l50": "contig_l50",
        "largest scaffold": "largest_contig",
        "gc content %": "gc_pct",
        "# gaps in scaffolds": "n_gaps",
    }
    with open(path) as fh:
        for line in fh:
            if ":" not in line:
                continue
            k, _, v = line.partition(":")
            k, v = k.strip().lower(), v.strip()
            if k in key_map:
                try:
                    out[key_map[k]] = float(v) if "." in v else int(v)
                except ValueError:
                    out[key_map[k]] = v

    # An assembly with no scaffold-length line but a contig-length one still has a size.
    if "assembly_size" not in out and "total_contig_length" in out:
        out["assembly_size"] = out["total_contig_length"]
    return out


def parse_busco(path):
    """Retain the full breakdown, not a single percentage."""
    out = {}
    if not path or not Path(path).exists():
        return out
    try:
        with open(path) as fh:
            data = json.load(fh)
    except (ValueError, OSError) as e:
        # Warn rather than return an empty dict silently. A BUSCO breakdown that vanishes
        # without trace is exactly the kind of gap this project cannot afford: duplication
        # is the primary evidence that both nuclei were retained.
        print("WARNING: could not parse BUSCO JSON %s: %s" % (path, e), file=sys.stderr)
        return out

    results = data.get("results", data)
    field_map = {
        "Complete percentage": "busco_complete_pct",
        "Single copy percentage": "busco_single_pct",
        "Multi copy percentage": "busco_duplicated_pct",
        "Fragmented percentage": "busco_fragmented_pct",
        "Missing percentage": "busco_missing_pct",
        "Complete BUSCOs": "busco_complete_n",
        "Single copy BUSCOs": "busco_single_n",
        "Multi copy BUSCOs": "busco_duplicated_n",
        "Fragmented BUSCOs": "busco_fragmented_n",
        "Missing BUSCOs": "busco_missing_n",
        "n_markers": "busco_n_markers",
    }
    for src, dest in field_map.items():
        if src in results:
            out[dest] = results[src]
    if "lineage_dataset" in data:
        ld = data["lineage_dataset"]
        out["busco_lineage"] = ld.get("name") if isinstance(ld, dict) else ld
    return out


def parse_qv(path):
    """Merqury .qv is tab-separated; the QV value is column 4."""
    if not path or not Path(path).exists():
        return {}
    try:
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 4:
                    try:
                        return {"qv": round(float(f[3]), 2)}
                    except ValueError:
                        continue
    except OSError:
        pass
    return {}


def parse_telomeres(path):
    """
    Telomere placement, from telomere_scan.py.

    `telomere_interstitial_arrays` is the one to watch: a telomere in the middle of a contig
    is a mis-join. Zero is the expected and healthy value, so a non-zero count belongs in the
    comparison table rather than buried in a per-assembly file.
    """
    if not path or not Path(path).exists():
        return {}
    try:
        with open(path) as fh:
            d = json.load(fh)
    except (json.JSONDecodeError, OSError):
        return {}
    keep = (
        "telomere_capped_ends",
        "telomere_pct_ends",
        "telomere_t2t_contigs",
        "telomere_interstitial_arrays",
        "telomere_interstitial_contigs",
        "telomere_longest_array_units",
        "telomere_motif",
    )
    return {k: d[k] for k in keep if k in d}


def parse_coverage(path):
    if not path or not Path(path).exists():
        return {}
    try:
        with open(path) as fh:
            d = json.load(fh)
    except (json.JSONDecodeError, OSError):
        return {}
    return {
        "mean_depth": d.get("mean_depth"),
        "coverage_modes": ",".join(str(m) for m in d.get("detected_modes", [])),
        "coverage_interpretation": " | ".join(d.get("interpretation", [])),
    }


def _load_json(path):
    if not path or not Path(path).exists():
        return {}
    try:
        with open(path) as fh:
            return json.load(fh)
    except (json.JSONDecodeError, OSError):
        return {}


def parse_cleaning(fcs_path, mito_path):
    """
    What the QC columns were measured on, and what cleaning removed to get there.

    Every other metric in the record describes the assembly as handed on -- after FCS-GX,
    and after the mito screen when it ran. That is a different assembly from the raw one
    (M. larici-populina hifiasm hap1: 136.03 Mb raw, 104.40 Mb cleaned), so the row has to
    say which it is, and keep the raw size so the removal itself stays visible: a large
    fraction removed is evidence about the sample, not just a housekeeping step.
    """
    fcs = _load_json(fcs_path)
    mito = _load_json(mito_path)
    out = {}
    if not fcs and not mito:
        out["qc_input"] = "raw"
        return out

    out["qc_input"] = "fcs_gx_cleaned_mito_free" if mito else "fcs_gx_cleaned"
    raw = fcs.get("bases_in") or mito.get("bases_in")
    if raw:
        out["raw_assembly_size"] = raw
    if fcs:
        out["contaminant_bases_removed"] = fcs.get("bases_removed")
        out["contaminant_contigs_removed"] = fcs.get("contigs_excluded")
    if mito:
        out["mito_bases_removed"] = mito.get("mito_bases")
        out["mito_contigs_removed"] = mito.get("mito_contigs")
    removed = (fcs.get("bases_removed") or 0) + (mito.get("mito_bases") or 0)
    if raw:
        out["cleaning_removed_pct"] = round(100.0 * removed / raw, 2)
    return out


# ---------------------------------------------------------------------------- scoring

def size_flag(assembly_size, meta, lo, hi):
    """
    Score assembly size against the dikaryon expectation.

    Haplotype-resolved assemblies (hap1/hap2) hold ONE nucleus, so they are judged against
    half the dikaryon range. Primary assemblies may legitimately span from one collapsed
    haploid copy up to both haplotypes retained.
    """
    if not assembly_size:
        return "unknown"

    atype = meta.get("assembly_type", "")
    if "hap1" in atype or "hap2" in atype:
        lo, hi = lo / 2, hi / 2

    if assembly_size < 0.75 * lo:
        return "collapsed"
    if assembly_size > 1.15 * hi:
        return "oversized"
    return "expected"


# ---------------------------------------------------------------------------- commands

def cmd_record(args):
    meta = json.loads(args.meta)
    rec = {
        "assembly_id": meta.get("id"),
        "sample": meta.get("sample"),
        "assembler": meta.get("assembler"),
        "assembly_type": meta.get("assembly_type"),
        "readset": meta.get("readset"),
    }
    rec.update(parse_gfastats(args.gfastats))
    rec.update(parse_busco(args.busco))
    rec.update(parse_qv(args.qv))
    rec.update(parse_coverage(args.coverage))
    rec.update(parse_telomeres(args.telomeres))
    rec.update(parse_cleaning(args.fcs_clean, args.mito))

    rec["size_flag"] = size_flag(
        rec.get("assembly_size"), meta, args.dikaryon_min, args.dikaryon_max
    )

    with open(args.out, "w") as fh:
        json.dump(rec, fh, indent=2)
    print(f"Wrote record for {rec['assembly_id']} (size_flag={rec['size_flag']})")


# ---------------------------------------------------------------------------- scaffolds

def scaffold_lengths(path):
    """
    Sequence lengths from a FASTA, descending.

    gfastats reports '# scaffolds', which for M. larici-populina hap1 is 55 — but only 18 of
    those are chromosome-scale, and 18 is the number the paper quotes and the number this
    branch exists to produce. That distribution is not in gfastats' output at any verbosity,
    so it is measured here from the sequence itself.
    """
    lengths, cur, seen = [], 0, False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seen:
                    lengths.append(cur)
                cur, seen = 0, True
            else:
                cur += len(line.strip())
    if seen:
        lengths.append(cur)
    lengths.sort(reverse=True)
    return lengths


def l_count(lengths, frac):
    """How many of the longest sequences it takes to cover `frac` of the total."""
    total = sum(lengths)
    if not total:
        return None
    target, run = total * frac, 0
    for i, n in enumerate(lengths, start=1):
        run += n
        if run >= target:
            return i
    return len(lengths)


def cmd_scaffold_record(args):
    meta = json.loads(args.meta)
    gf = parse_gfastats(args.gfastats)
    lengths = scaffold_lengths(args.fasta)
    total = sum(lengths)

    big = [n for n in lengths if n >= args.chromosome_min_length]
    span = sum(big)

    rec = {
        "assembly_id": meta.get("id"),
        # The id carries a '_scaffolds' suffix so it can never be mistaken for the contig-level
        # row of the same assembly. Keep the unsuffixed name too, so the two can be joined.
        "source_assembly": re.sub(r"_scaffolds$", "", meta.get("id") or ""),
        "sample": meta.get("sample"),
        "assembler": meta.get("assembler"),
        "assembly_type": meta.get("assembly_type"),
        "readset": meta.get("readset"),
        "scaffolder": args.scaffolder,
        "n_scaffolds": gf.get("n_scaffolds", len(lengths)),
        "scaffold_total_length": gf.get("assembly_size", total),
        "scaffold_n50": gf.get("n50"),
        "scaffold_l50": gf.get("l50"),
        "largest_scaffold": gf.get("largest_contig"),
        # The before/after contrast. Scaffolding does not change these: identical contig N50
        # with a larger scaffold N50 is the signature of joins and nothing else.
        "n_contigs": gf.get("n_contigs"),
        "contig_n50": gf.get("contig_n50"),
        "n_gaps": gf.get("n_gaps"),
        "chromosome_min_length": args.chromosome_min_length,
        "n_chromosome_scale": len(big),
        "chromosome_scale_span": span,
        "chromosome_scale_pct": round(100.0 * span / total, 2) if total else None,
        # Threshold-free companion to the count above, so the table still says something
        # useful when the 1 Mb default is wrong for the species.
        "l90": l_count(lengths, 0.90),
    }

    with open(args.out, "w") as fh:
        json.dump(rec, fh, indent=2)
    print("Wrote scaffold record for %s (%s chromosome-scale of %s scaffolds)" % (
        rec["assembly_id"], rec["n_chromosome_scale"], rec["n_scaffolds"]))


SCAFFOLD_ORDER = [
    "assembly_id", "source_assembly", "sample", "assembler", "assembly_type", "readset",
    "scaffolder", "n_scaffolds", "scaffold_total_length", "n_chromosome_scale",
    "chromosome_scale_span", "chromosome_scale_pct", "chromosome_min_length",
    "scaffold_n50", "scaffold_l50", "l90", "largest_scaffold",
    "n_contigs", "contig_n50", "n_gaps",
]


def cmd_scaffold_merge(args):
    records = []
    for p in args.records:
        try:
            with open(p) as fh:
                records.append(json.load(fh))
        except (json.JSONDecodeError, OSError) as e:
            print("WARNING: skipping unreadable scaffold record %s: %s" % (p, e),
                  file=sys.stderr)

    if not records:
        sys.exit("No scaffold records to merge")

    records.sort(key=lambda r: str(r.get("assembly_id") or ""))

    keys = [k for k in SCAFFOLD_ORDER if any(k in r for r in records)]
    fields = keys + sorted({k for r in records for k in r} - set(keys))

    with open("scaffold_summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in records:
            w.writerow({k: r.get(k, "") for k in fields})

    with open("scaffold_summary.json", "w") as fh:
        json.dump({"scaffolds": records}, fh, indent=2)

    mqc_cols = [c for c in [
        "assembly_id", "scaffolder", "n_scaffolds", "n_chromosome_scale",
        "chromosome_scale_pct", "scaffold_n50", "contig_n50", "n_gaps",
        "scaffold_total_length",
    ] if c in fields]

    with open("scaffold_summary_mqc.tsv", "w") as fh:
        fh.write("# id: scaffold_summary\n")
        fh.write("# section_name: 'Hi-C scaffolding'\n")
        fh.write("# description: 'Scaffolded PHASED haplotypes. These rows are NOT candidate "
                 "assemblies and are absent from the assembly comparison on purpose: a "
                 "scaffold N50 is raised by joins alone and is not comparable with the "
                 "contig N50 of an unscaffolded assembly. Compare scaffold_n50 against "
                 "contig_n50 in the SAME row to see what scaffolding actually did. "
                 "n_chromosome_scale counts scaffolds at or above chromosome_min_length.'\n")
        fh.write("# format: 'tsv'\n")
        fh.write("# plot_type: 'table'\n")
        fh.write("\t".join(mqc_cols) + "\n")
        for r in records:
            fh.write("\t".join(str(r.get(c, "")) for c in mqc_cols) + "\n")

    print("Wrote scaffold_summary.tsv / .json (%d scaffolded assemblies)" % len(records))


# Column order: identity, then size/contiguity, then completeness, then quality.
PREFERRED_ORDER = [
    "assembly_id", "sample", "assembler", "assembly_type", "readset",
    "qc_input", "assembly_size", "size_flag", "n_contigs", "n_scaffolds",
    "raw_assembly_size", "cleaning_removed_pct",
    "contaminant_bases_removed", "contaminant_contigs_removed",
    "mito_bases_removed", "mito_contigs_removed",
    "n50", "contig_n50", "l50", "largest_contig", "gc_pct",
    "busco_lineage", "busco_complete_pct", "busco_single_pct",
    "busco_duplicated_pct", "busco_fragmented_pct", "busco_missing_pct",
    "busco_complete_n", "busco_single_n", "busco_duplicated_n",
    "busco_fragmented_n", "busco_missing_n", "busco_n_markers",
    "qv", "mean_depth", "coverage_modes", "coverage_interpretation",
    "telomere_capped_ends", "telomere_pct_ends", "telomere_t2t_contigs",
    "telomere_interstitial_arrays", "telomere_interstitial_contigs",
    "telomere_longest_array_units", "telomere_motif",
]


def cmd_merge(args):
    records = []
    for p in args.records:
        try:
            with open(p) as fh:
                records.append(json.load(fh))
        except (json.JSONDecodeError, OSError) as e:
            print(f"WARNING: skipping unreadable record {p}: {e}", file=sys.stderr)

    if not records:
        sys.exit("No assembly records to merge")

    # Sort for a stable, readable table — explicitly NOT by N50 or any quality metric,
    # so the file never implies a ranking.
    records.sort(key=lambda r: (
        str(r.get("assembler") or ""),
        str(r.get("readset") or ""),
        str(r.get("assembly_type") or ""),
    ))

    keys = [k for k in PREFERRED_ORDER if any(k in r for r in records)]
    extra = sorted({k for r in records for k in r} - set(keys))
    fields = keys + extra

    with open("assembly_summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in records:
            w.writerow({k: r.get(k, "") for k in fields})

    with open("assembly_summary.json", "w") as fh:
        json.dump({"assemblies": records}, fh, indent=2)

    # MultiQC custom table, trimmed to the columns a human actually compares.
    mqc_cols = [c for c in [
        "assembly_id", "assembler", "assembly_type", "readset", "qc_input", "assembly_size",
        "size_flag", "raw_assembly_size", "cleaning_removed_pct", "n_contigs", "n50", "busco_complete_pct", "busco_duplicated_pct",
        "qv", "mean_depth", "coverage_modes", "telomere_capped_ends",
        "telomere_t2t_contigs", "telomere_interstitial_arrays",
    ] if c in fields]

    with open("assembly_summary_mqc.tsv", "w") as fh:
        fh.write("# id: assembly_summary\n")
        fh.write("# section_name: 'Assembly comparison'\n")
        fh.write("# description: 'Candidate assemblies with evidence. Assemblies are NOT "
                 "ranked by N50. Every metric describes the assembly named in qc_input -- "
                 "the CLEANED one whenever FCS-GX cleaning ran -- and raw_assembly_size / "
                 "cleaning_removed_pct show what cleaning took out. size_flag scores size "
                 "against the expected 0.8-1.0 Gb dikaryon; high BUSCO duplication is expected "
                 "in a correctly phased dikaryotic assembly.'\n")
        fh.write("# format: 'tsv'\n")
        fh.write("# plot_type: 'table'\n")
        fh.write("\t".join(mqc_cols) + "\n")
        for r in records:
            fh.write("\t".join(str(r.get(c, "")) for c in mqc_cols) + "\n")

    mis_joined = [r for r in records if r.get("telomere_interstitial_arrays")]
    if mis_joined:
        print("\nAssemblies with telomeres INSIDE contigs (possible mis-joins):",
              file=sys.stderr)
        for r in mis_joined:
            print("  %s: %s interstitial array(s) on %s contig(s)" % (
                r.get("assembly_id"), r.get("telomere_interstitial_arrays"),
                r.get("telomere_interstitial_contigs")), file=sys.stderr)

    flagged = [r for r in records if r.get("size_flag") in ("collapsed", "oversized")]
    if flagged:
        print("\nAssemblies flagged on size:", file=sys.stderr)
        for r in flagged:
            print(f"  {r.get('assembly_id')}: {r.get('size_flag')} "
                  f"({r.get('assembly_size')} bp)", file=sys.stderr)

    print(f"Wrote assembly_summary.tsv / .json ({len(records)} assemblies)")


def main():
    ap = argparse.ArgumentParser()
    # NOTE: no `required=True` here. This script runs under whatever python3 the executing
    # environment provides, and the cluster's `nextflow` conda env is python 3.6, where
    # add_subparsers() does not accept `required`. It raises TypeError at import time, which
    # is how every ASSEMBLY_RECORD task failed in run 32466998. Check the result instead.
    sub = ap.add_subparsers(dest="cmd")

    r = sub.add_parser("record")
    r.add_argument("--meta", required=True, help="assembly meta as JSON")
    r.add_argument("--gfastats")
    r.add_argument("--busco")
    r.add_argument("--qv")
    r.add_argument("--coverage")
    r.add_argument("--telomeres")
    r.add_argument("--fcs-clean", help="fcs_gx_clean.json; its presence means QC ran on the cleaned assembly")
    r.add_argument("--mito", help="mito_screen.json; its presence means QC ran on the mito-free assembly")
    r.add_argument("--dikaryon-min", type=float, default=800e6)
    r.add_argument("--dikaryon-max", type=float, default=1000e6)
    r.add_argument("--out", required=True)
    r.set_defaults(func=cmd_record)

    m = sub.add_parser("merge")
    m.add_argument("records", nargs="+")
    m.set_defaults(func=cmd_merge)

    sr = sub.add_parser("scaffold-record")
    sr.add_argument("--meta", required=True, help="assembly meta as JSON")
    sr.add_argument("--fasta", required=True, help="scaffolded FASTA")
    sr.add_argument("--gfastats", required=True)
    sr.add_argument("--scaffolder", default="")
    sr.add_argument("--chromosome-min-length", type=int, default=1000000)
    sr.add_argument("--out", required=True)
    sr.set_defaults(func=cmd_scaffold_record)

    sm = sub.add_parser("scaffold-merge")
    sm.add_argument("records", nargs="+")
    sm.set_defaults(func=cmd_scaffold_merge)

    args = ap.parse_args()
    if not getattr(args, "cmd", None):
        ap.error("a subcommand is required (record | merge)")
    args.func(args)


if __name__ == "__main__":
    main()
