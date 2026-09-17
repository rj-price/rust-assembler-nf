#!/usr/bin/env python3
"""
Locate telomeric repeat arrays in an assembly, at contig ends and internally.

Two things are reported, and the second is the point of the exercise:

  * how many contig ends carry a telomere, and how many contigs are capped at BOTH ends
    (telomere-to-telomere) — a contiguity measure that, unlike N50, means something
    biological;
  * every INTERSTITIAL array, i.e. a telomere sitting in the middle of a contig. That is the
    signature of a mis-join, where a chromosome end has been stitched into the interior.
    Flye in particular is often said to do this, so the pipeline measures it rather than
    repeating the claim.

Distance from the nearest contig end is computed AFTER locating each array anywhere in the
sequence, so the terminal/interstitial split is a measurement rather than a definition.

The motif is a parameter. TTAGGG is canonical for fungi including the rusts, but it is not
universal, and a wrong motif yields a confidently empty result. The array length distribution
in the output is the check: real telomeres appear as tens to hundreds of clean tandem units,
whereas a wrong motif gives only scattered 3-4 unit hits at chance frequency.
"""

import argparse
import gzip
import json
import re
import sys

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path):
    """Stream (name, sequence). Handles plain or gzipped input."""
    opener = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def find_arrays(seq, motif, min_units):
    """
    Yield (start, end, units, strand) for every run of >= min_units tandem copies.

    Both strands are searched, since a telomere at the 5' end of a contig appears as the
    reverse complement. The work is handed to `re` rather than a Python loop: these are
    gigabase assemblies, and a per-position scan in the interpreter is minutes-to-hours
    slower than the same scan in C.
    """
    for unit, strand in ((motif, "+"), (revcomp(motif), "-")):
        pattern = re.compile("(?:%s){%d,}" % (unit, min_units))
        for m in pattern.finditer(seq):
            yield m.start(), m.end(), (m.end() - m.start()) // len(unit), strand


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", help="assembly FASTA, optionally gzipped")
    ap.add_argument("--assembly-id", required=True)
    ap.add_argument("--motif", default="TTAGGG",
                    help="telomeric repeat unit (default: %(default)s, canonical for fungi)")
    ap.add_argument("--min-units", type=int, default=3,
                    help="minimum tandem copies to call an array (default: %(default)s)")
    ap.add_argument("--interstitial-min-units", type=int, default=20,
                    help="stricter threshold for calling an INTERNAL array real, since short "
                         "G-rich runs occur by chance in the interior (default: %(default)s)")
    ap.add_argument("--end-window", type=int, default=1000,
                    help="an array within this many bp of a contig end is terminal "
                         "(default: %(default)s)")
    ap.add_argument("--arrays", help="TSV of every array retained")
    ap.add_argument("--out", required=True, help="JSON summary for the assembly record")
    args = ap.parse_args()

    motif = args.motif.upper()
    if set(motif) - set("ACGT") or not motif:
        sys.exit("--motif must be a non-empty ACGT string, got %r" % args.motif)

    arrays_fh = open(args.arrays, "w") if args.arrays else None
    if arrays_fh:
        arrays_fh.write("assembly_id\tcontig\tcontig_len\tstart\tend\tunits\tstrand"
                        "\tdist_from_end\tposition\n")

    n_contigs = n_bp = 0
    capped_left = capped_right = t2t = 0
    interstitial = interstitial_bp = 0
    interstitial_contigs = set()
    longest_units = 0

    for name, seq in read_fasta(args.fasta):
        n_contigs += 1
        length = len(seq)
        n_bp += length
        seq = seq.upper()
        has_left = has_right = False

        for start, end, units, strand in find_arrays(seq, motif, args.min_units):
            dist = min(start, length - end)
            if start < args.end_window:
                position = "start"
                has_left = True
            elif length - end < args.end_window:
                position = "end"
                has_right = True
            elif units >= args.interstitial_min_units:
                position = "interstitial"
                interstitial += 1
                interstitial_bp += end - start
                interstitial_contigs.add(name)
            else:
                continue        # short interior hit: chance occurrence, not a telomere
            longest_units = max(longest_units, units)
            if arrays_fh:
                arrays_fh.write("%s\t%s\t%d\t%d\t%d\t%d\t%s\t%d\t%s\n" % (
                    args.assembly_id, name, length, start, end, units, strand,
                    dist, position))

        capped_left += has_left
        capped_right += has_right
        t2t += has_left and has_right

    if arrays_fh:
        arrays_fh.close()

    ends = 2 * n_contigs
    capped = capped_left + capped_right
    summary = {
        "assembly_id": args.assembly_id,
        "telomere_motif": motif,
        "n_contigs": n_contigs,
        "assembly_bp": n_bp,
        "telomere_capped_ends": capped,
        "telomere_pct_ends": round(100.0 * capped / ends, 2) if ends else None,
        "telomere_t2t_contigs": t2t,
        "telomere_longest_array_units": longest_units,
        "telomere_interstitial_arrays": interstitial,
        "telomere_interstitial_contigs": len(interstitial_contigs),
        "telomere_interstitial_bp": interstitial_bp,
    }

    with open(args.out, "w") as fh:
        json.dump(summary, fh, indent=2)

    print("%s: %d/%d contig ends capped (%s%%), %d T2T contigs, %d interstitial arrays" % (
        args.assembly_id, capped, ends, summary["telomere_pct_ends"], t2t, interstitial))

    # A confidently empty result usually means the wrong motif, not a telomere-free genome.
    if capped == 0 and n_contigs:
        print("WARNING: no telomeres found at any contig end. Check --telomere_motif: %s "
              "may not be the repeat for this species." % motif, file=sys.stderr)


if __name__ == "__main__":
    main()
