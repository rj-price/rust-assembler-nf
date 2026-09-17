process SEQKIT_STATS {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/read_qc/seqkit", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.seqkit_stats.tsv"), emit: stats
    path "versions.yml",                                  emit: versions

    script:
    // -a gives N50, Q1/Q2/Q3, AvgQual etc. — the full picture, not just count/length.
    """
    seqkit stats \\
        --all \\
        --tabular \\
        --threads ${task.cpus} \\
        ${reads} \\
        > ${meta.id}.seqkit_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """

    stub:
    // A header and one row, not an empty file: the per-sample yield is read back out of
    // this TSV with splitCsv, and an empty stub output fails the whole run with
    // "Missing 'header' in CSV file" before any downstream logic is exercised.
    //
    // The escapes are DOUBLED. Singles worked only by accident: Groovy consumed them and
    // emitted real tabs and a real newline into the printf format string, which bash then
    // passed through unchanged. bin/check_config_selectors.py rejects that on sight, because
    // the same pattern in an awk program or a different quoting context silently restructures
    // the script instead of surviving it.
    """
    printf 'file\\tformat\\ttype\\tnum_seqs\\tsum_len\\tmin_len\\tavg_len\\tmax_len\\tQ1\\tQ2\\tQ3\\tsum_gap\\tN50\\tN50_num\\tQ20(%%)\\tQ30(%%)\\tAvgQual\\tGC(%%)\\n' \\
        > ${meta.id}.seqkit_stats.tsv
    printf '%s\\tFASTA\\tDNA\\t100000\\t1000000000\\t1000\\t10000.0\\t30000\\t6000.0\\t8000.0\\t10000.0\\t0\\t9000\\t3000\\t0.00\\t0.00\\t0.00\\t43.60\\n' \\
        stub.fasta.gz >> ${meta.id}.seqkit_stats.tsv
    echo '"${task.process}": {seqkit: stub}' > versions.yml
    """
}
