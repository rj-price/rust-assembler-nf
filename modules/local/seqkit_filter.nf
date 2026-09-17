process SEQKIT_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/derived/filtered_reads", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.filtered.fastq.gz"), emit: reads
    path "versions.yml"                                  , emit: versions

    when:
    params.min_read_length > 0 || params.min_read_q > 0

    script:
    // DECISION D6: this NEVER overwrites the input. The filtered dataset is a new,
    // separately-named file in derived/, leaving raw/ untouched.
    def len_arg = params.min_read_length > 0 ? "--min-len ${params.min_read_length}" : ''
    def q_arg   = params.min_read_q      > 0 ? "--min-qual ${params.min_read_q}"     : ''
    """
    seqkit seq \\
        ${len_arg} \\
        ${q_arg} \\
        --threads ${task.cpus} \\
        ${reads} \\
        | gzip -c > ${meta.id}.filtered.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.filtered.fastq.gz
    echo '"${task.process}": {seqkit: stub}' > versions.yml
    """
}
