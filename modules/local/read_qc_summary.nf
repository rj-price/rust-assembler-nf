process READ_QC_SUMMARY {
    tag "${sample}"
    label 'process_single'

    publishDir "${params.outdir}/read_qc", mode: params.publish_dir_mode

    input:
    tuple val(sample), val(run_ids), path(stats_files)

    output:
    path "read_summary.tsv"      , emit: tsv
    path "read_summary.json"     , emit: json
    path "read_summary_mqc.tsv"  , emit: mqc
    path "read_summary_notes.txt", emit: notes, optional: true
    path "versions.yml"          , emit: versions

    script:
    // Pair each run id with its stats file positionally; both lists come from the same
    // collected channel, so order is consistent.
    def pairs = [run_ids, stats_files]
        .transpose()
        .collect { rid, f -> "${rid}:${f}" }
        .join(' ')
    """
    read_qc_summary.py --stats ${pairs} --outprefix read_summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch read_summary.tsv read_summary.json read_summary_mqc.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
