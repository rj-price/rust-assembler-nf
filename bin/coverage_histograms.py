#!/usr/bin/env python3
"""
Combine per-assembly depth histograms into one overlaid MultiQC line plot.

One curve per candidate, on a shared axis, because the comparison is the point: a retained
dikaryon shows a single clean mode at the per-nucleus depth, a collapsed one shows its mass
at twice that, and an assembly carrying spurious duplicate sequence shows a shoulder at half
it. Those three shapes are obvious side by side and easy to miss in eleven separate files.

Curves are normalised to percent of positions so assemblies of different sizes are
comparable -- otherwise the largest assembly simply plots highest everywhere.
"""

import argparse
import os
import sys


def read_hist(path):
    """(depth -> pct) from a coverage_summary.py --histogram TSV."""
    out = {}
    with open(path) as fh:
        header = fh.readline()
        if "depth" not in header:
            fh.seek(0)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                out[int(f[0])] = float(f[2])
            except ValueError:
                continue
    return out


def read_expected_map(path):
    """(assembly_id -> expected per-nucleus depth) from a two-column TSV."""
    out = {}
    if not path:
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            try:
                out[f[0]] = float(f[1])
            except ValueError:
                continue
    return out


def svg_small_multiples(series, max_depth, expected=None, expected_map=None, cols=3):
    """
    One panel per assembly, shared axes, hand-written SVG.

    Small multiples rather than one overlaid plot on purpose: a categorical palette runs out
    at eight hues, and a run can easily produce eleven candidates. Panels also make the
    shape comparison — one clean peak, versus a low-depth shoulder — easier than eleven
    overlapping curves. The interactive overlay lives in the MultiQC report instead.

    No plotting library: this runs in a stock python container, and adding matplotlib to it
    would cost more than the plot is worth.
    """
    names = sorted(series)
    rows = (len(names) + cols - 1) // cols
    pw, ph, gx, gy = 250, 130, 26, 34          # panel box and gaps
    ml, mt = 42, 26                            # margins for axis labels and title
    W = ml + cols * pw + (cols - 1) * gx + 14
    H = mt + rows * ph + (rows - 1) * gy + 20

    ymax = max(max(h.values()) for h in series.values() if h) or 1.0
    out = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" '
           'height="%d" font-family="sans-serif">' % (W, H, W, H),
           '<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H),
           '<text x="8" y="16" font-size="13" fill="#111">Coverage depth distribution '
           '(%% of assembled positions)</text>']

    for idx, name in enumerate(names):
        r, c = divmod(idx, cols)
        x0 = ml + c * (pw + gx)
        y0 = mt + r * (ph + gy)
        h = series[name]
        out.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" '
                   'stroke="#dddddd"/>' % (x0, y0, pw, ph))
        # Per-assembly expectation where one is known, falling back to the global figure.
        # Yield differs per sample, so one line across every panel would be wrong for all
        # but one of them in a multi-sample run, which is why this is not a single global param.
        exp = (expected_map or {}).get(name, expected)
        if exp:
            xe = x0 + min(exp, max_depth) / max_depth * pw
            out.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#bbbbbb" '
                       'stroke-dasharray="3 3"/>' % (xe, y0, xe, y0 + ph))
        pts = []
        for d in range(1, max_depth + 1):
            x = x0 + d / max_depth * pw
            y = y0 + ph - (h.get(d, 0.0) / ymax) * ph
            pts.append("%s%.1f %.1f" % ("M" if d == 1 else "L", x, y))
        out.append('<path d="%s" fill="none" stroke="#1baf7a" stroke-width="1.4"/>'
                   % " ".join(pts))
        out.append('<text x="%d" y="%d" font-size="9" fill="#333">%s</text>'
                   % (x0, y0 - 4, name))
        out.append('<text x="%d" y="%d" font-size="8" fill="#777" text-anchor="end">%.1f%%</text>'
                   % (x0 - 3, y0 + 7, ymax))
        out.append('<text x="%d" y="%d" font-size="8" fill="#777" text-anchor="end">0</text>'
                   % (x0 - 3, y0 + ph))
        out.append('<text x="%d" y="%d" font-size="8" fill="#777">0</text>'
                   % (x0, y0 + ph + 11))
        out.append('<text x="%d" y="%d" font-size="8" fill="#777" text-anchor="end">%d'
                   '&#215;</text>' % (x0 + pw, y0 + ph + 11, max_depth))
    out.append("</svg>")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("histograms", nargs="+")
    ap.add_argument("--max-depth", type=int, default=250,
                    help="x-axis limit; the informative range is around the expected depth, "
                         "and a long flat repeat tail squashes it (default: %(default)s)")
    ap.add_argument("--out", default="coverage_histogram_mqc.tsv")
    ap.add_argument("--svg", default="coverage_histograms.svg",
                    help="small-multiples plot, one panel per assembly")
    ap.add_argument("--expected-depth", type=float, default=None,
                    help="expected per-nucleus depth, drawn as a reference line")
    ap.add_argument("--expected-depth-map", default=None,
                    help="TSV of assembly_id<TAB>expected depth, one row per assembly; "
                         "takes precedence over --expected-depth per panel")
    args = ap.parse_args()

    series = {}
    for p in args.histograms:
        name = os.path.basename(p).replace(".depth_histogram.tsv", "")
        h = read_hist(p)
        if h:
            series[name] = h
        else:
            print("WARNING: no data in %s" % p, file=sys.stderr)

    if not series:
        sys.exit("No usable depth histograms")

    names = sorted(series)
    with open(args.out, "w") as fh:
        fh.write("# id: coverage_histogram\n")
        fh.write("# section_name: 'Coverage depth distribution'\n")
        fh.write("# description: 'Per-base read depth, one curve per assembly, normalised to "
                 "percent of positions. A retained dikaryon shows a single mode at the "
                 "per-nucleus depth; mass at twice that indicates collapse, and a shoulder at "
                 "half it indicates spurious duplicate sequence.'\n")
        fh.write("# plot_type: 'linegraph'\n")
        fh.write("# pconfig:\n")
        fh.write("#     id: 'coverage_histogram_plot'\n")
        fh.write("#     title: 'Coverage depth distribution'\n")
        fh.write("#     xlab: 'Depth (x)'\n")
        fh.write("#     ylab: '% of assembled positions'\n")
        fh.write("Depth\t" + "\t".join(names) + "\n")
        for d in range(1, args.max_depth + 1):
            fh.write("%d\t%s\n" % (d, "\t".join("%.6f" % series[n].get(d, 0.0)
                                                for n in names)))

    with open(args.svg, "w") as fh:
        fh.write(svg_small_multiples(series, args.max_depth, args.expected_depth,
                                     read_expected_map(args.expected_depth_map)))

    print("Wrote %s and %s (%d assemblies, depth 1-%d)"
          % (args.out, args.svg, len(names), args.max_depth))


if __name__ == "__main__":
    main()
