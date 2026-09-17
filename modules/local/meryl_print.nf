process MERYL_PRINT {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(meryl_db), path(cutoffs)

    output:
    tuple val(meta), path("${meta.id}.kmers_in_range.tsv"), emit: kmers
    path "versions.yml"                                   , emit: versions

    when:
    params.run_smudgeplot

    script:
    // Step 2 of 3: the meryl half of smudgeplot's workflow, which cannot share a container
    // with the smudgeplot half (see SMUDGEPLOT_CUTOFF). Dumps only the k-mers between the
    // coverage cutoffs — the full database is far too large to print.
    """
    L=\$(cut -f1 ${cutoffs})
    U=\$(cut -f2 ${cutoffs})

    meryl print \\
        [ less-than \${U} [ greater-than \${L} ${meryl_db} ] ] \\
        > ${meta.id}.kmers_in_range.tsv

    if [ ! -s ${meta.id}.kmers_in_range.tsv ]; then
        echo "ERROR: no k-mers between the smudgeplot cutoffs L=\$L U=\$U." >&2
        echo "       The coverage model and the k-mer database disagree; check the histogram." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$(meryl --version 2>&1 | sed 's/^meryl //' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.kmers_in_range.tsv
    echo '"${task.process}": {meryl: stub}' > versions.yml
    """
}
