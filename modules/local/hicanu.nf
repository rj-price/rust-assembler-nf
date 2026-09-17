process HICANU {
    tag "${meta.id}"
    label 'process_assembly'

    publishDir "${params.outdir}/assembly/hicanu/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.hicanu.fasta.gz"), emit: fasta, optional: true
    tuple val(meta), path("${meta.id}.hicanu.log")     , emit: log,   optional: true
    path "versions.yml"                                , emit: versions, optional: true

    when:
    params.run_hicanu

    script:
    // KNOWN RISK: Canu self-submits to the grid by default, which fights the
    // Nextflow executor and produces jobs neither system is tracking. useGrid=false forces it
    // to run entirely inside this single allocation.
    //
    // It is also slow — potentially >1 week on 60 Gb — and its errorStrategy is 'ignore' in
    // base.config so a HiCanu timeout cannot sink the rest of the pipeline. hifiasm and Flye
    // already provide the independent-structure comparison if this one never lands.
    def mem_gb = task.memory.toGiga()
    """
    canu \\
        -p ${meta.id} \\
        -d canu_out \\
        genomeSize=${params.assembler_genome_size ?: params.genome_size} \\
        -pacbio-hifi ${reads} \\
        useGrid=false \\
        maxThreads=${task.cpus} \\
        maxMemory=${mem_gb} \\
        ${params.hicanu_extra} \\
        2>&1 | tee ${meta.id}.hicanu.log

    if [ -f canu_out/${meta.id}.contigs.fasta ]; then
        gzip -c canu_out/${meta.id}.contigs.fasta > ${meta.id}.hicanu.fasta.gz
    else
        echo "WARNING: HiCanu produced no contigs.fasta — check ${meta.id}.hicanu.log" >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        canu: \$(canu -version 2>&1 | head -1 | sed 's/Canu //')
    END_VERSIONS
    """

    stub:
    """
    echo ">stub" | gzip -c > ${meta.id}.hicanu.fasta.gz
    touch ${meta.id}.hicanu.log
    echo '"${task.process}": {canu: stub}' > versions.yml
    """
}
