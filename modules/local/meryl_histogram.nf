process MERYL_HISTOGRAM {
    tag "${meta.id}"
    label 'process_single'

    publishDir "${params.outdir}/kmer/spectra", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(meryl_db)

    output:
    tuple val(meta), path("${meta.id}.hist"), emit: hist
    path "versions.yml"                     , emit: versions

    script:
    """
    meryl histogram \\
        threads=${task.cpus} \\
        ${meryl_db} \\
        > ${meta.id}.hist

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$(meryl --version 2>&1 | sed 's/meryl //' | head -1)
    END_VERSIONS
    """

    stub:
    """
    printf '1\\t1000\\n2\\t500\\n' > ${meta.id}.hist
    echo '"${task.process}": {meryl: stub}' > versions.yml
    """
}
