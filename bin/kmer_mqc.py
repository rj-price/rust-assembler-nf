#!/usr/bin/env python3
"""
Turn a GenomeScope2 summary into a MultiQC table and a machine-readable JSON.

GenomeScope2's own summary.txt is a fixed-width human-readable block that no MultiQC module
parses, so its genome size, heterozygosity and repeat content never reached the report even
though they are the numbers every assembly decision here is scored against. This lifts them
out into both forms.

CAVEAT, carried into the table itself rather than left in a config comment: GenomeScope2
models a DIPLOID. A dikaryon can break that fit, so these figures are a cross-check on the
size estimate, not the source of truth. Smudgeplot and the observed coverage modes are the
other two checks.

That caveat was not enough. On M. larici-populina 98AG31 GenomeScope2 reported a haploid
length of 48 Mb with no indication anything was wrong, and the pipeline carried that number
for weeks; the assembled haplotypes came out at 102-104 Mb, a factor of two out. The fit had
failed, and everything needed to see that was already on disk:

    Model Fit          67.73% - 80.09%     (a converged fit is >95%)
    Homozygous (aa)    0% - 100%           unconstrained
    Heterozygous (ab)  0% - 100%           unconstrained
    length   3.300e+07 +/- 3.602e+07   p = 0.360   not distinguishable from zero
    kmercov  2.081e+02                 vs smudgeplot 1n = 104.6

The mechanism is peak misassignment: GenomeScope2 locked onto the homozygous AABB peak at
multiplicity ~210 and called it 1n, so every length came out halved. A dikaryon is exactly
the case that provokes this, which makes it this pipeline's problem rather than a curiosity.

So this script now WARNS. It does not correct anything and does not fail the run -- the
estimate may still be the best available, and it is the user who has to judge that. It makes
the failure visible instead of letting a confident wrong number travel downstream.
"""

import argparse
import json
import re
import sys

# GenomeScope2 summary lines look like:
#   Genome Haploid Length     464,283,624 bp     526,048,914 bp
#   Heterozygosity            1.31%              1.42%
# i.e. a label followed by a min and a max estimate.
# The heterozygosity row is labelled "Heterozygous (ab)" for a diploid model and
# "Heterozygosity" in some builds, so accept either — matching only one of them yielded a
# confident "NA" in the report for the single number this pipeline cares about most.
FIELDS = [
    ("genome_haploid_length", r"Genome Haploid Length"),
    ("genome_repeat_length", r"Genome Repeat Length"),
    ("genome_unique_length", r"Genome Unique Length"),
    ("heterozygosity", r"(?:Heterozygous \(ab\)|Heterozygosity)"),
    ("homozygosity", r"Homozygous \(aa\)"),
    ("model_fit", r"Model Fit"),
    ("read_error_rate", r"Read Error Rate"),
]


def parse_number(token):
    token = token.replace(",", "").replace("%", "").replace("bp", "").strip()
    try:
        return float(token)
    except ValueError:
        return None


def parse_summary(path):
    """Return {field: {'min': x, 'max': y}} for the fields we report."""
    text = open(path).read()
    out = {}
    for key, label in FIELDS:
        m = re.search(r"^%s\s+(.*)$" % label, text, re.MULTILINE)
        if not m:
            continue
        # Split the value part into the min/max pair; units are dropped by parse_number.
        tokens = re.findall(r"[\d,.]+\s*(?:bp|%)?", m.group(1))
        values = [v for v in (parse_number(t) for t in tokens) if v is not None]
        if not values:
            continue
        out[key] = {"min": values[0], "max": values[-1]}
    return out


def parse_model_kmercov(path):
    """
    Pull `kmercov` out of GenomeScope2's model.txt (the nls fit summary).

    The row is a fixed-width R coefficient table:
        kmercov 2.081e+02  2.294e-01 907.273   <2e-16 ***
    Only the estimate (first column after the name) is wanted.
    """
    try:
        text = open(path).read()
    except OSError:
        return None
    m = re.search(r"^kmercov\s+([\d.eE+-]+)", text, re.MULTILINE)
    return float(m.group(1)) if m else None


def parse_smudgeplot_1n(path):
    """
    Pull the 1n coverage smudgeplot settled on from its verbose summary.

    The line reads:
        1n coverage used in smudgeplot (one of the three above):\t104.6
    which is smudgeplot's own choice between the user, subset and highest-peak estimates --
    the right one to compare against, rather than re-deriving a preference here.
    """
    try:
        text = open(path).read()
    except OSError:
        return None
    m = re.search(r"^1n coverage used in smudgeplot[^:]*:\s*([\d.]+)", text, re.MULTILINE)
    return float(m.group(1)) if m else None


# A converged GenomeScope2 fit sits well above 90%; the failed 98AG31 fit topped out at 80.1%.
MODEL_FIT_FLOOR = 90.0

# How close to an exact 2x (or 0.5x) the kmercov/1n ratio has to be before it is called peak
# misassignment rather than ordinary disagreement. 98AG31 was 208.1/104.6 = 1.99.
PLOIDY_RATIO_TOLERANCE = 0.15


def check_fit(parsed, kmercov, smudge_1n):
    """
    Return a list of human-readable warnings about the fit. Empty means nothing detected --
    NOT that the fit is good, since only a few specific failures are testable here.
    """
    warnings = []

    fit = parsed.get("model_fit")
    if fit and fit.get("max") is not None and fit["max"] < MODEL_FIT_FLOOR:
        warnings.append(
            "model fit only %.1f%% (a converged fit is >%.0f%%) - treat the genome size as "
            "unreliable" % (fit["max"], MODEL_FIT_FLOOR)
        )

    # An unconstrained parameter is reported as the full 0-100% range. It means the optimiser
    # could not pin it down at all, which is a clearer failure signal than the fit percentage.
    for key, label in (("heterozygosity", "heterozygosity"), ("homozygosity", "homozygosity")):
        entry = parsed.get(key)
        if entry and entry.get("min") == 0.0 and entry.get("max") == 100.0:
            warnings.append(
                "%s is unconstrained (0-100%%) - the model did not converge" % label
            )

    # The one that actually catches a halved genome. GenomeScope2 fitting the homozygous peak
    # as 1n is the classic dikaryon/high-heterozygosity failure, and smudgeplot is the
    # independent estimate that exposes it.
    if kmercov and smudge_1n:
        ratio = kmercov / smudge_1n
        for factor, effect in ((2.0, "HALVED"), (0.5, "DOUBLED")):
            if abs(ratio - factor) <= PLOIDY_RATIO_TOLERANCE:
                warnings.append(
                    "GenomeScope kmercov %.1f is %.2fx smudgeplot's 1n coverage %.1f - it "
                    "has fitted the wrong peak, so the haploid length is likely %s"
                    % (kmercov, ratio, smudge_1n, effect)
                )
                break
        else:
            if ratio < 0.8 or ratio > 1.25:
                warnings.append(
                    "GenomeScope kmercov %.1f and smudgeplot 1n coverage %.1f disagree "
                    "(%.2fx) - cross-check the genome size before relying on it"
                    % (kmercov, smudge_1n, ratio)
                )

    return warnings


def fmt_bp(value):
    if value is None:
        return "NA"
    return "%.1f Mb" % (value / 1e6)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", required=True, help="GenomeScope2 <prefix>_summary.txt")
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-mqc", required=True)
    ap.add_argument("--model", help="GenomeScope2 <prefix>_model.txt, for kmercov")
    ap.add_argument("--smudgeplot-summary",
                    help="smudgeplot <prefix>_verbose_summary.txt, for the 1n coverage estimate")
    args = ap.parse_args()

    parsed = parse_summary(args.summary)
    if not parsed:
        sys.stderr.write(
            "WARNING: no recognised fields in %s — GenomeScope2's summary layout may have "
            "changed. Writing an empty table rather than failing the run.\n" % args.summary
        )

    # Both cross-check inputs are optional: smudgeplot can be off, and model.txt is an
    # optional GenomeScope2 output. A missing file means that check is skipped, not that the
    # fit passed it -- which is why the JSON records what was actually available.
    kmercov = parse_model_kmercov(args.model) if args.model and not args.model.startswith("NO_FILE") else None
    smudge_1n = (parse_smudgeplot_1n(args.smudgeplot_summary)
                 if args.smudgeplot_summary and "NO_FILE" not in args.smudgeplot_summary else None)

    warnings = check_fit(parsed, kmercov, smudge_1n)
    for w in warnings:
        sys.stderr.write("WARNING [%s]: %s\n" % (args.sample_id, w))

    record = {
        "sample": args.sample_id,
        "genomescope2": parsed,
        "kmercov": kmercov,
        "smudgeplot_1n_coverage": smudge_1n,
        "fit_warnings": warnings,
    }
    with open(args.out_json, "w") as fh:
        json.dump(record, fh, indent=2, sort_keys=True)
        fh.write("\n")

    def rng(key, formatter=lambda v: "%.2f" % v if v is not None else "NA"):
        entry = parsed.get(key)
        if not entry:
            return "NA"
        lo, hi = entry.get("min"), entry.get("max")
        if lo == hi:
            return formatter(lo)
        return "%s - %s" % (formatter(lo), formatter(hi))

    with open(args.out_mqc, "w") as fh:
        fh.write("# id: genomescope2\n")
        fh.write("# section_name: 'k-mer genome model (GenomeScope2)'\n")
        fh.write(
            "# description: 'Reference-free genome characterisation from the meryl k-mer "
            "spectrum, computed BEFORE assembly. GenomeScope2 fits a DIPLOID model, which a "
            "dikaryon can break — read these as a cross-check on the assumed genome size, "
            "alongside smudgeplot and the observed coverage modes, not as the source of "
            "truth. The Fit warnings column reports specific, testable failures -- a poor "
            "model fit, an unconstrained parameter, or a kmercov that is a clean multiple of "
            "smudgeplot\'s 1n coverage, which means the wrong peak was fitted and the "
            "genome size is out by that factor. \"none detected\" is not a clean bill of "
            "health: only those few failures are checked.'\n"
        )
        fh.write("# format: 'tsv'\n")
        fh.write("# plot_type: 'table'\n")
        fh.write(
            "Sample\tHaploid length\tUnique length\tRepeat length\tHeterozygosity (%)\t"
            "Model fit (%)\tRead error rate (%)\tFit warnings\n"
        )
        # The warnings share the row rather than sitting in a separate section, because the
        # number they qualify is in this row. A reader who sees "48.0 Mb" needs to see
        # "the model did not converge" in the same glance, not two tables further down.
        fh.write(
            "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
            % (
                args.sample_id,
                rng("genome_haploid_length", fmt_bp),
                rng("genome_unique_length", fmt_bp),
                rng("genome_repeat_length", fmt_bp),
                rng("heterozygosity"),
                rng("model_fit"),
                rng("read_error_rate"),
                "; ".join(warnings) if warnings else "none detected",
            )
        )


if __name__ == "__main__":
    main()
