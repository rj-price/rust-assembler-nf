process SMUDGEPLOT {
    tag "${meta.id}"
    label 'process_himem'

    publishDir "${params.outdir}/kmer/smudgeplot", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(kmers)

    output:
    tuple val(meta), path("*.png")        , emit: plots,   optional: true
    tuple val(meta), path("*_summary.txt"), emit: summary, optional: true
    path "versions.yml"                   , emit: versions

    when:
    params.run_smudgeplot

    script:
    // Step 3 of 3. Ploidy sanity check that does NOT assume a diploid, which is exactly why
    // it is here: GenomeScope2's diploid model may fit a dikaryon poorly.
    //
    // The `meryl print` that used to sit in the middle of this script now lives in
    // MERYL_PRINT, because meryl is not in the smudgeplot container — the 2026-08-26 run died
    // here with "meryl: command not found" after 3 h of queueing.
    """
    # The subcommand is verified before it is run. A blanket `|| skip` here previously turned
    # "hetkmers is not a valid task name" -- a wrong-container bug -- into the reassuring
    # message "found no usable k-mer pairs", and the run reported a clean skip for days. Only a
    # genuine absence of pairs may be skipped; anything else must fail loudly.
    if ! smudgeplot.py --help 2>&1 | grep -qw hetkmers; then
        echo "ERROR: this smudgeplot has no 'hetkmers' task. Version 0.4.x renamed it to" >&2
        echo "       'hetmers' and reads a FastK .ktab, not a k-mer dump. Pin 0.2.5 -- see" >&2
        echo "       conf/apptainer.config." >&2
        smudgeplot.py --version >&2 || true
        exit 1
    fi

    # Do NOT blanket-skip on failure. hetkmers loads the whole dump into memory and is a
    # realistic OOM candidate on a dikaryon (9.1 GB of k-mers here), and an OOM swallowed as
    # "no usable pairs" is both a wrong answer and a defeat of the retry ladder, which would
    # otherwise double the memory and succeed. A kill signal must propagate; only an ordinary
    # non-zero exit is treated as "this data has no smudge".
    set +e
    smudgeplot.py hetkmers -o ${meta.id} < ${kmers}
    RC=\$?
    set -e

    if [ "\$RC" -ge 128 ]; then
        echo "ERROR: hetkmers was killed by signal \$(( RC - 128 )) (exit \$RC)." >&2
        echo "       Almost certainly out of memory: it holds the entire k-mer dump in RAM." >&2
        echo "       Failing so the retry ladder can re-run it with more." >&2
        exit \$RC
    fi

    if [ "\$RC" -ne 0 ]; then
        echo "hetkmers exited \$RC; treating as no usable k-mer pairs" >&2
        touch ${meta.id}_skipped_summary.txt
    fi

    if [ -f "${meta.id}_coverages.tsv" ]; then
        smudgeplot.py plot -o ${meta.id} ${meta.id}_coverages.tsv
    elif [ ! -f "${meta.id}_skipped_summary.txt" ]; then
        # hetkmers returned 0 but produced nothing: a real "no pairs" result, recorded as such
        # rather than left as a silently missing output.
        echo "hetkmers produced no coverages file; no smudge to plot" >&2
        touch ${meta.id}_skipped_summary.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smudgeplot: \$(smudgeplot.py --version 2>&1 | sed 's/smudgeplot //' | head -1)
    END_VERSIONS
    """

    stub:
    """
    # _verbose_summary.txt, matching what smudgeplot actually writes. The stub used to emit
    # <id>_summary.txt, which collides with GenomeScope2's summary of the same name once both
    # are staged into KMER_SUMMARY -- a collision that cannot happen in a real run and so hid
    # behind the stub until KMER_SUMMARY started consuming both.
    touch ${meta.id}_verbose_summary.txt
    echo '"${task.process}": {smudgeplot: stub}' > versions.yml
    """
}
