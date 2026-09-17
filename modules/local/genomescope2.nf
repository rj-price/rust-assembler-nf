process GENOMESCOPE2 {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/kmer/genomescope", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(hist)

    output:
    tuple val(meta), path("${meta.id}_summary.txt")     , emit: summary
    tuple val(meta), path("${meta.id}_model.txt")       , emit: model,  optional: true
    tuple val(meta), path("*.png")                      , emit: plots,  optional: true
    path "versions.yml"                                 , emit: versions

    script:
    // CAVEAT: GenomeScope2 models a DIPLOID. A dikaryon can break its fit, so treat its
    // genome-size estimate as a cross-check on --genome_size, not as the
    // source of truth. Smudgeplot and the observed coverage modes are the other checks.
    """
    genomescope2 \\
        --input ${hist} \\
        --output . \\
        --kmer_length ${params.kmer_size} \\
        --ploidy 2 \\
        --name_prefix ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomescope2: \$(genomescope2 --version 2>&1 | sed 's/GenomeScope //' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_summary.txt
    echo '"${task.process}": {genomescope2: stub}' > versions.yml
    """
}
