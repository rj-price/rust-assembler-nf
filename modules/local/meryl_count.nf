process MERYL_COUNT {
    tag "${meta.id}"
    label 'process_kmer'

    publishDir "${params.outdir}/kmer/meryl", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.k${params.kmer_size}.meryl"), emit: meryl_db
    path "versions.yml"                                           , emit: versions

    script:
    // meryl sizes its hash from the expected genome size; a bad guess costs memory, not
    // correctness. sized for a ~500 Mb haploid / ~1 Gb dikaryon rust.
    def mem_gb = task.memory ? (task.memory.toGiga() * 0.8) as int : 100
    """
    meryl count \\
        k=${params.kmer_size} \\
        threads=${task.cpus} \\
        memory=${mem_gb} \\
        output ${meta.id}.k${params.kmer_size}.meryl \\
        ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$(meryl --version 2>&1 | sed 's/meryl //' | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}.k${params.kmer_size}.meryl
    echo '"${task.process}": {meryl: stub}' > versions.yml
    """
}
