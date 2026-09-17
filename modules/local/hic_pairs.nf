process HIC_PAIRS {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/phasing/hic_contacts", mode: params.publish_dir_mode,
        // The dedup stats are the record worth keeping; the pairs file is large and cheaply
        // reproducible from the alignment.
        saveAs: { fn -> fn.endsWith('.pairs.gz') ? null : fn }

    input:
    tuple val(meta), path(sam), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.id}.dedup.pairs.gz"), emit: pairs
    tuple val(meta), path("${meta.id}.hic_stats.txt") , emit: stats
    path "versions.yml"                               , emit: versions

    script:
    // Alignments into deduplicated ligation pairs.
    //
    // --walks-policy 5unique keeps the 5'-most unique alignment of a multi-fragment walk
    // rather than discarding the read outright, which matters on a library with many short
    // fragments.
    //
    // MAPQ is the one lever that changes what counts as a contact at all. NuclearPhaser's
    // author recommends MAPQ 10 with 100 kb bins, and MAPQ 30 with 20 kb bins as the
    // alternative — so both are parameters rather than hard-coded here, and they belong to
    // this process and HIC_MATRIX respectively.
    """
    zcat ${sam} \\
        | pairtools parse \\
            --min-mapq ${params.hic_min_mapq} \\
            --walks-policy 5unique \\
            --max-inter-align-gap 30 \\
            --chroms-path ${chrom_sizes} \\
            --nproc-in ${task.cpus} --nproc-out ${task.cpus} \\
        | pairtools sort --nproc ${task.cpus} --tmpdir=. \\
        | pairtools dedup \\
            --mark-dups \\
            --output-stats ${meta.id}.hic_stats.txt \\
            --output ${meta.id}.dedup.pairs.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version 2>&1 | sed 's/.*version //')
        hic_min_mapq: ${params.hic_min_mapq}
    END_VERSIONS
    """

    stub:
    """
    echo | gzip -c > ${meta.id}.dedup.pairs.gz
    touch ${meta.id}.hic_stats.txt
    echo '"${task.process}": {pairtools: stub}' > versions.yml
    """
}
