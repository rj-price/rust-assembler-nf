#!/usr/bin/env python3
"""
Call mitochondrial contigs in an assembly, from alignment against a reference mitogenome.

Why this exists
---------------
Duplessis et al. (2026) removed the mitochondrion from their Melampsora haplotypes as a
separate step. This pipeline did not, so every haplotype it has produced carries a mito
contig inside a "nuclear" assembly -- it inflates the assembly size, it is the reason a
single contig can sit at 50x the nuclear depth, and it is one of the contigs in the short
tail that made our contig counts higher than the published ones.

FCS-GX does not catch it. FCS-GX asks "does this sequence belong to a DIFFERENT organism",
and a mitochondrion belongs to this one. It is the right sequence in the wrong assembly,
which is a different question and needs a different tool.

Why BLASTn and not minimap2
---------------------------
The first version of this used minimap2 -x asm20. Measured against the RefSeq mitochondrion
database on the real mlp98AG31 haplotypes, it called 31 contigs / 3.41 Mb of hap1
mitochondrial -- against a true mitogenome of ~47 kb on one contig. The false calls were
tandem-repeat arrays chaining onto unrelated mitogenomes at roughly 14% identity, which
minimap2 reports as one long alignment block and does not filter, because it has no identity
floor. BLASTn with DUST masking at >=90% identity, the paper's own parameters, separates the
real mitochondrion (aligned fraction 0.35) from the worst false positive (0.04) with room to
spare. Low-complexity masking is the load-bearing half: mitogenomes and fungal repeat arrays
are both AT-rich, so without it the search matches composition rather than homology.

Calling rule
------------
A contig is called MITO when the fraction of ITS OWN length covered by alignments to the
reference mitogenome reaches --min-aligned-frac (default 0.2), and it is no longer than
--max-length (default 500 kb).

Both conditions matter, and the second is not redundant:

  * The fraction test, not a raw alignment length or a best-hit identity, is what separates
    a real mitochondrion from a NUMT. Nuclear-mitochondrial insertions are common and can
    align over many kilobases, but they sit inside a chromosome-scale contig, so the
    fraction of that contig which is mitochondrial stays tiny. A raw-length threshold would
    delete a whole chromosome because 8 kb of it is an old insertion.

  * The length ceiling is the backstop for the case where the fraction test is wrong --
    a badly fragmented assembly where a 2 Mb contig somehow reaches 50%, or a reference
    mitogenome that contains nuclear contamination of its own.

Coverage intervals are merged before the fraction is computed, so a contig that aligns to
the same reference region three times over is not counted three times -- mitogenomes are
circular and repetitive, and overlapping alignments are the norm rather than the exception.

Everything with ANY alignment is written to the TSV with its call, so NUMTs and near-miss
contigs are visible and auditable. Only the MITO rows are acted on. The removed sequence is
written out rather than deleted, and the input assembly is never modified -- the same three
safeguards FCS_GX_CLEAN uses.
"""

import argparse
import gzip
import json
import sys


def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def read_fasta(path):
    """Yield (name, description_line, sequence). Names are split on whitespace."""
    name, header, chunks = None, None, []
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, header, "".join(chunks)
                header = line.rstrip("\n")
                name = header[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, header, "".join(chunks)


def merge_intervals(intervals):
    """Merge overlapping [start, end) intervals and return the total length covered."""
    if not intervals:
        return 0
    intervals.sort()
    total = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
    total += cur_end - cur_start
    return total


def parse_blast_hits(path):
    """
    Return {query_name: [query_length, [(start, end), ...]]}.

    Reads BLAST tabular output written as
        6 qseqid qlen qstart qend pident length sseqid

    BLAST coordinates are 1-based inclusive and, on a minus-strand hit, qstart > qend. Both
    are normalised here to 0-based half-open ascending intervals so merge_intervals sees one
    consistent convention -- getting this wrong would produce negative interval lengths and
    silently deflate every aligned fraction.

    Hits are kept regardless of which reference they came from. The nearest mitogenome in the
    database may be several genera away, a reference may itself be in pieces, and a circular
    molecule linearised at an arbitrary point produces hits to both ends of the query.
    """
    hits = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            qname, qlen = f[0], int(f[1])
            qstart, qend = int(f[2]), int(f[3])
            if qstart > qend:
                qstart, qend = qend, qstart
            entry = hits.setdefault(qname, [qlen, []])
            entry[1].append((qstart - 1, qend))
    return hits


def gc_fraction(seq):
    if not seq:
        return 0.0
    gc = sum(seq.count(base) for base in "GCgc")
    at = sum(seq.count(base) for base in "ATat")
    return gc / (gc + at) if (gc + at) else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="assembly to screen")
    ap.add_argument("--hits", required=True,
                    help="BLASTn tabular: 6 qseqid qlen qstart qend pident length sseqid")
    ap.add_argument("--assembly-id", required=True)
    # 0.2, not 0.5, and measured rather than guessed: on mlp98AG31 the mitochondrion and
    # its double-length concatemer both sit at 0.349 while the worst false positive reaches
    # 0.041. The threshold belongs in that gap, nearer the noise than the signal because the
    # reference may be a distant relative and a real mitochondrion can align over less of
    # itself than this one did.
    ap.add_argument("--min-aligned-frac", type=float, default=0.2)
    ap.add_argument("--max-length", type=int, default=500000)
    ap.add_argument("--out-fasta", required=True, help="assembly with mito contigs removed")
    ap.add_argument("--out-mito", required=True, help="the removed contigs")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args()

    hits = parse_blast_hits(args.hits)

    calls = {}
    for qname, (qlen, intervals) in hits.items():
        covered = merge_intervals(intervals)
        frac = covered / qlen if qlen else 0.0
        is_mito = frac >= args.min_aligned_frac and qlen <= args.max_length
        # Say WHY, so a surprising call can be argued with rather than just re-run.
        if is_mito:
            reason = "aligned_fraction>=%.2f and length<=%d" % (args.min_aligned_frac,
                                                               args.max_length)
        elif frac >= args.min_aligned_frac:
            reason = "aligned fraction high but contig longer than --max-length"
        else:
            reason = "partial alignment only (NUMT, shared repeat or distant homology)"
        calls[qname] = {
            "length": qlen,
            "aligned_bp": covered,
            "aligned_fraction": round(frac, 4),
            "call": "MITO" if is_mito else "KEEP",
            "reason": reason,
        }

    kept_n = kept_bp = mito_n = mito_bp = 0
    total_n = total_bp = 0

    with open_maybe_gzip(args.out_fasta, "wt") as keep_fh, \
         open_maybe_gzip(args.out_mito, "wt") as mito_fh:
        for name, header, seq in read_fasta(args.fasta):
            total_n += 1
            total_bp += len(seq)
            call = calls.get(name)
            if call is not None:
                call["gc"] = round(gc_fraction(seq), 4)
            target = mito_fh if (call and call["call"] == "MITO") else keep_fh
            if call and call["call"] == "MITO":
                mito_n += 1
                mito_bp += len(seq)
            else:
                kept_n += 1
                kept_bp += len(seq)
            target.write(header + "\n")
            for i in range(0, len(seq), 60):
                target.write(seq[i:i + 60] + "\n")

    with open(args.out_tsv, "w") as fh:
        fh.write("assembly_id\tcontig\tlength\taligned_bp\taligned_fraction\tgc\tcall\treason\n")
        for name in sorted(calls, key=lambda n: -calls[n]["aligned_fraction"]):
            c = calls[name]
            fh.write("%s\t%s\t%d\t%d\t%.4f\t%s\t%s\t%s\n" % (
                args.assembly_id, name, c["length"], c["aligned_bp"],
                c["aligned_fraction"],
                ("%.4f" % c["gc"]) if "gc" in c else "NA",
                c["call"], c["reason"]))

    summary = {
        "assembly_id": args.assembly_id,
        "min_aligned_frac": args.min_aligned_frac,
        "max_length": args.max_length,
        "contigs_in": total_n,
        "bases_in": total_bp,
        "contigs_kept": kept_n,
        "bases_kept": kept_bp,
        "mito_contigs": mito_n,
        "mito_bases": mito_bp,
        "mito_percent_of_assembly": round(100.0 * mito_bp / total_bp, 4) if total_bp else 0.0,
        "contigs_with_any_alignment": len(calls),
        "calls": calls,
    }
    with open(args.out_json, "w") as fh:
        json.dump(summary, fh, indent=2)

    print("%s: %d contig(s), %d bp called mitochondrial; %d contig(s) had partial "
          "alignments and were kept." % (args.assembly_id, mito_n, mito_bp,
                                         len(calls) - mito_n), file=sys.stderr)


if __name__ == "__main__":
    main()
