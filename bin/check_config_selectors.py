#!/usr/bin/env python3
"""Pre-flight guards against traps that fail SILENTLY or only on real data.

CHECK 1 -- the same process selector in more than one Nextflow config file.

Nextflow does NOT deep-merge `withName:` / `withLabel:` blocks across included config files.
If the same selector string appears in two files, the later file's block REPLACES the earlier
one wholesale, silently, keeping none of its settings.

That is not a hypothetical. In run 32466998, conf/apptainer.config's
    withName: 'FCS_GX' { container = ... }
wiped conf/base.config's
    withName: 'FCS_GX' { cpus = 16; memory = 700.GB; time = 48.h }
so FCS-GX fell back to its generic label's 220 GB, thrashed against a 498 GB index and hit the
wall twice. Six other processes were silently mis-sized the same way, and conf/test.config was
deleting FLYE's and BUSCO's containers out of the smoke test whose job is to prove containers
work.

CHECK 2 -- a bin/ script a module invokes that is not executable, or has no shebang.

Nextflow stages bin/ onto PATH but does NOT chmod it, so a script committed without its
executable bit dies with exit 126 ("Permission denied") the first time the process runs.

A stub run cannot catch this: stub blocks touch their outputs instead of calling the script,
so the whole test suite passes and the failure waits for real data. In the 2026-08-26 full
run, KMER_SUMMARY died at 2 h 31 m -- after MERYL_COUNT had chewed through all 60.29 Gb --
because bin/kmer_mqc.py was mode 644. bin/fcs_gx_clean.py was queued to fail the same way
several days later, once the assemblies existed.

CHECK 3 -- a lone backslash escape inside a process script block.

A Nextflow script block is a Groovy string BEFORE it is a shell script, so Groovy consumes
`\\n` and `\\t` first. Inside a single-quoted awk program a Groovy-eaten `\\n` arrives as a real
newline, which is the syntax error "Unexpected end of string". This applies to COMMENTS in the
block too -- Nextflow does not strip them -- so prose describing an escape sequence breaks the
script just as code does. Both happened: the awk in KRAKEN2, and then the comment written to
warn about the awk, each costing an hour-plus of finished classification.

CHECK 4 -- a process directive set from a bare `params.x` rather than a closure.

Nextflow evaluates a config file as it parses it. conf/base.config is included by
nextflow.config BEFORE any profile, so every param a site profile supplies is still null at
that point and `container = params.fcs_gx_sif` bakes in `container = null` for good. The task
then runs on the bare host with no container and no warning of any kind.

On 2026-09-01 that meant FCS-GX did not run on any of the five assemblies: /app/bin/run_gx
exists only inside the image, run_gx exited 2, a `|| echo WARNING` swallowed it, and both
report outputs were declared `optional`, so the task exited 0 and FCS_GX_CLEAN passed each
assembly through untouched while reporting "no findings". Wrapping the value in a closure --
`container = { params.fcs_gx_sif }` -- defers it to task submission, after profiles resolve.

Usage:  python3 bin/check_config_selectors.py [conf_dir]
"""
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

SELECTOR = re.compile(r"""with(?:Name|Label)\s*:\s*(?:'([^']+)'|"([^"]+)"|([A-Za-z_]\w*))""")


def strip_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)   # block comments
    return re.sub(r"//[^\n]*", "", text)                 # line comments


def check_selectors(conf_dir):
    where = defaultdict(set)
    for f in sorted(conf_dir.glob("*.config")):
        for m in SELECTOR.finditer(strip_comments(f.read_text())):
            where[next(g for g in m.groups() if g)].add(f.name)

    dupes = {s: fs for s, fs in where.items() if len(fs) > 1}
    if dupes:
        print("ERROR: process selector defined in more than one config file.", file=sys.stderr)
        print("Nextflow REPLACES rather than merges these, so settings are lost silently.\n",
              file=sys.stderr)
        for sel in sorted(dupes):
            print("  %-16s <-- %s" % (sel, ", ".join(sorted(dupes[sel]))), file=sys.stderr)
        print("\nFix: declare the container and the resources together in ONE block.",
              file=sys.stderr)
        return 1

    print("OK: %d selectors, none defined in more than one config file." % len(where))
    return 0


def check_bin_scripts(repo):
    """Every bin/ script a module actually invokes must be executable and have a shebang."""
    bin_dir = repo / "bin"
    if not bin_dir.is_dir():
        print("OK: no bin/ directory to check.")
        return 0

    # Only scripts a process actually calls matter. A helper run by hand (this file,
    # make_test_data.sh) is never staged onto a task's PATH, so its mode is irrelevant --
    # flagging it would be noise that trains you to ignore the check.
    nf_text = "\n".join(
        f.read_text() for f in sorted(repo.glob("**/*.nf"))
        if ".nextflow" not in f.parts and "work" not in f.parts
    )

    problems = []
    checked = 0
    for script in sorted(bin_dir.iterdir()):
        if not script.is_file():
            continue
        # Word-boundary match: kmer_mqc.py must not be satisfied by some other_kmer_mqc.py.
        if not re.search(r"(?<![\w./-])%s\b" % re.escape(script.name), nf_text):
            continue
        checked += 1
        why = []
        if not os.access(script, os.X_OK):
            why.append("not executable (mode %04o)" % (script.stat().st_mode & 0o7777))
        with open(script, "rb") as fh:
            if not fh.read(2) == b"#!":
                why.append("no shebang")
        if why:
            problems.append((script.name, "; ".join(why)))

    if problems:
        print("\nERROR: bin/ script invoked by a module is not runnable.", file=sys.stderr)
        print("Nextflow stages bin/ onto PATH but does not chmod it -- these die with exit 126",
              file=sys.stderr)
        print("on real data, and a stub run cannot catch it (stubs never call the script).\n",
              file=sys.stderr)
        for name, why in problems:
            print("  %-24s <-- %s" % (name, why), file=sys.stderr)
        print("\nFix: chmod +x bin/<script>  (and commit the mode: git update-index --chmod=+x)",
              file=sys.stderr)
        return 1

    print("OK: %d bin/ scripts invoked by modules, all executable with a shebang." % checked)
    return 0


def check_script_escapes(repo):
    """Flag a single backslash escape inside a process script block.

    Only `\\n`, `\\t` and `\\r` are reported. `\\$` is deliberately ignored: Groovy renders it as a
    plain `$`, which is what the author wanted and is harmless in both code and comments. The
    three flagged here are the ones that inject a real newline or tab and so change the
    STRUCTURE of the emitted script.
    """
    escape = re.compile(r"(?<!\\)\\(?!\\)[ntr]")
    problems = []
    for nf in sorted(repo.glob("**/*.nf")):
        if ".nextflow" in nf.parts or "work" in nf.parts:
            continue
        text = nf.read_text()
        # Only blocks that become a SHELL SCRIPT. A triple-quoted string elsewhere (the
        # workflow.onComplete message in main.nf, say) is pure Groovy, where a lone \\n is a
        # newline on purpose and entirely correct.
        for blk in re.finditer(r'(?:script|shell|stub)\s*:\s*(?:.*?\n)??\s*"""(.*?)"""',
                               text, re.S):
            first = text[: blk.start(1)].count("\n") + 1
            for offset, line in enumerate(blk.group(1).split("\n")):
                if escape.search(line):
                    kind = "comment" if line.strip().startswith("#") else "code"
                    problems.append((nf.relative_to(repo), first + offset, kind, line.strip()))

    if problems:
        print("\nERROR: single backslash escape inside a process script block.", file=sys.stderr)
        print("Groovy consumes it before bash or awk ever sees it, injecting a real newline or",
              file=sys.stderr)
        print("tab and changing the script's structure. Comments are NOT exempt.\n",
              file=sys.stderr)
        for path, line, kind, text in problems:
            print("  %s:%d (%s)" % (path, line, kind), file=sys.stderr)
            print("      %s" % text[:96], file=sys.stderr)
        print("\nFix: double it (\\\\n), or in prose spell it out in words.", file=sys.stderr)
        return 1

    print("OK: no un-doubled backslash escapes in process script blocks.")
    return 0


# Directives whose value is resolved once at config-parse time, so a bare `params.x` in a file
# included before the site profile silently becomes null. Nextflow accepts a closure for all
# of them and evaluates it per task instead.
LAZY_DIRECTIVES = ("container", "queue", "clusterOptions", "containerOptions", "beforeScript")
# Not anchored to line start: conf/apptainer.config writes whole selector blocks inline,
# as `withName: 'X' { container = ... }`.
DIRECTIVE_ASSIGN = re.compile(
    r"(?<![\w.])(%s)\s*=\s*(?P<val>\S[^\n]*)" % "|".join(LAZY_DIRECTIVES)
)


def check_lazy_directives(conf_dir):
    """CHECK 4: a parse-time directive set from params without a closure."""
    bad = []
    for path in sorted(conf_dir.glob("*.config")):
        text = strip_comments(path.read_text())
        for m in DIRECTIVE_ASSIGN.finditer(text):
            # Trim the trailing brace of an inline selector block: `{ container = X }`.
            val = m.group("val").strip().rstrip("}").strip()
            # A closure is fine -- it is evaluated per task, long after profiles resolve.
            if val.startswith("{") or "params." not in val:
                continue
            line = text[: m.start()].count("\n") + 1
            bad.append((path.name, line, m.group(1), val))

    if bad:
        print("FAIL: process directive read from params without a closure.\n")
        for name, line, directive, val in bad:
            print("  %s:%d  %s = %s" % (name, line, directive, val))
        print(
            "\nThis file may be parsed before the profile that sets the param, in which case the\n"
            "value is null and stays null -- for `container` that means the task runs on the bare\n"
            "host, silently. Wrap it so it is evaluated per task instead:\n"
            "    container = { params.fcs_gx_sif }"
        )
        return 1

    print("OK: no parse-time process directives read params without a closure.")
    return 0


def main():
    repo = Path(__file__).resolve().parent.parent
    conf_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else repo / "conf"
    if not conf_dir.is_dir():
        sys.exit("No such directory: %s" % conf_dir)

    # Run EVERY check before returning, rather than short-circuiting on the first failure:
    # finding one problem, fixing it, and then discovering the second on the next run is
    # exactly the slow loop these guards exist to prevent.
    rc = check_selectors(conf_dir)
    rc |= check_bin_scripts(repo)
    rc |= check_script_escapes(repo)
    rc |= check_lazy_directives(conf_dir)
    return rc


if __name__ == "__main__":
    sys.exit(main())
