process VERKKO {
    tag "${meta.id}"
    label 'process_assembly'

    publishDir "${params.outdir}/assembly/verkko/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.verkko.fasta.gz"), emit: fasta, optional: true
    tuple val(meta), path("${meta.id}*.gfa")           , emit: gfa,   optional: true
    path "versions.yml"                                , emit: versions, optional: true

    when:
    params.run_verkko

    script:
    // Off by default (DECISION D9). HiFi-only it mostly duplicates what hifiasm already tells
    // us; it earns its place with Hi-C, when NuclearPhaser can phase its primary.
    """
    # Verkko drives snakemake, which builds a SourceCache under \$HOME/.cache on startup.
    # \$HOME is read-only inside the container on the compute nodes, so that raised
    # "OSError: [Errno 30] Read-only file system: '/home/<user>/.cache'" and killed the task
    # before assembly began (job 34035432). Point HOME and the XDG cache at the task work
    # directory, which is writable and is discarded with the task.
    export HOME="\$PWD"
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p "\$XDG_CACHE_HOME"

    # --no-consensus-bam, for correctness rather than economy. Verkko's combineConsensus
    # rule builds a CRAM by globbing one pattern per consensus partition:
    #     Snakefiles/7-combineConsensus.sm:31
    #     expand("packages/part{nnnn}.bam*.bam", nnnn = ...part{xxxx}.cnspack...)
    # and hands the unexpanded patterns to `samtools merge`. utgcns only writes BAM records
    # for tigs it realigns, so a partition holding nothing but singletons legitimately
    # produces no BAM and exits 0. Its glob then matches nothing, bash passes the literal
    # `packages/part157.bam*.bam` through, samtools cannot open it, and snakemake's strict
    # mode kills the run after every other stage has succeeded. Job 34056916 died exactly
    # that way at 11 h 17 m, in 7-consensus, with 156 of 157 partitions fine and part157
    # logging "Processed 0 tigs and 782 singletons". This genome invites it: a long tail of
    # tiny contigs makes an all-singleton final partition likely.
    #
    # The flag sets withBAM=False, which skips that whole samtools block. Nothing is lost --
    # this pipeline consumes verkko_out/assembly.fasta and the GFAs, never the CRAM -- and
    # every other withBAM branch in verkko's driver concerns ONT or Hi-C read alignment,
    # neither of which we pass. It also halves verkko's layout memory request.
    verkko \\
        -d verkko_out \\
        --hifi ${reads} \\
        --threads ${task.cpus} \\
        --local-memory ${task.memory.toGiga()} \\
        --no-consensus-bam \\
        ${params.verkko_extra}

    if [ -f verkko_out/assembly.fasta ]; then
        gzip -c verkko_out/assembly.fasta > ${meta.id}.verkko.fasta.gz
    fi
    for g in verkko_out/*.gfa; do
        [ -e "\$g" ] && cp "\$g" ${meta.id}.\$(basename "\$g")
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verkko: \$(verkko --version 2>/dev/null | head -1)
    END_VERSIONS
    """

    stub:
    """
    echo ">stub" | gzip -c > ${meta.id}.verkko.fasta.gz
    echo '"${task.process}": {verkko: stub}' > versions.yml
    """
}
