process FLYE {
    tag "${meta.id}"
    label 'process_assembly'

    publishDir path: { "${params.outdir}/assembly/flye/${meta.id}" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.flye.fasta.gz")     , emit: fasta
    tuple val(meta), path("${meta.id}.assembly_graph.gfa"), emit: gfa, optional: true
    tuple val(meta), path("${meta.id}.flye.log")          , emit: log
    path "versions.yml"                                   , emit: versions

    when:
    params.run_flye

    script:
    // A fundamentally different assembly approach to hifiasm. The question this answers is
    // "do independent assemblers agree on the structure of the genome?", not "which has the
    // best N50" (DECISION D7).
    """
    flye \\
        --pacbio-hifi ${reads} \\
        --out-dir flye_out \\
        --threads ${task.cpus} \\
        --genome-size ${params.assembler_genome_size ?: params.genome_size} \\
        ${params.flye_extra}

    gzip -c flye_out/assembly.fasta > ${meta.id}.flye.fasta.gz
    cp flye_out/flye.log ${meta.id}.flye.log
    if [ -f flye_out/assembly_graph.gfa ]; then
        cp flye_out/assembly_graph.gfa ${meta.id}.assembly_graph.gfa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$(flye --version)
    END_VERSIONS
    """

    stub:
    """
    echo ">stub" | gzip -c > ${meta.id}.flye.fasta.gz
    touch ${meta.id}.flye.log
    echo '"${task.process}": {flye: stub}' > versions.yml
    """
}
