process HIC_MATRIX {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/phasing/hic_contacts", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(pairs), path(chrom_sizes)

    output:
    tuple val(meta), path("${meta.id}.contacts.ginteractions.tsv"), emit: contacts
    tuple val(meta), path("${meta.id}.cool")                      , emit: cool
    path "versions.yml"                                           , emit: versions

    script:
    // NuclearPhaser's documented route to a contact map is HiC-Pro followed by hicexplorer's
    // ginteractions conversion. HiC-Pro has no container anywhere, and it turns out neither
    // tool is needed: NuclearPhaser reads the map positionally as seven tab-separated columns
    //
    //     contig1  start1  end1  contig2  start2  end2  count
    //
    // which is exactly what `cooler dump --join` writes. That removes a whole unavailable
    // dependency from the path, and the .cool is kept as the reusable artefact — it is also
    // what you would load to look at the contact map yourself.
    """
    cooler cload pairs \\
        -c1 2 -p1 3 -c2 4 -p2 5 \\
        ${chrom_sizes}:${params.hic_matrix_resolution} \\
        ${pairs} \\
        ${meta.id}.cool

    cooler dump --join --table pixels ${meta.id}.cool \\
        > ${meta.id}.contacts.ginteractions.tsv

    if [ ! -s ${meta.id}.contacts.ginteractions.tsv ]; then
        echo "ERROR: empty Hi-C contact map for ${meta.id}." >&2
        echo "       Either the Hi-C reads do not belong to this assembly, or the MAPQ" >&2
        echo "       threshold (--hic_min_mapq ${params.hic_min_mapq}) removed everything." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(cooler --version 2>&1 | sed 's/cooler, version //')
        hic_matrix_resolution: ${params.hic_matrix_resolution}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.contacts.ginteractions.tsv ${meta.id}.cool
    echo '"${task.process}": {cooler: stub}' > versions.yml
    """
}
