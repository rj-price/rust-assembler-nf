process BUSCO {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/assembly_qc/busco", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(fasta), path(gfa)
    path  busco_db

    output:
    tuple val(meta), path("${meta.id}.busco.json"), emit: json
    tuple val(meta), path("${meta.id}.busco.txt") , emit: summary
    tuple val(meta), path("${meta.id}.busco_full_table.tsv"), emit: full_table, optional: true
    path "versions.yml"                           , emit: versions

    when:
    params.run_busco

    script:
    // DECISION D5: pucciniomycetes_odb12 is the correct clade for a rust. Lineages are
    // read from --busco_db, so --offline works. NOTE: odb12 lineages REQUIRE BUSCO v6 — v5 cannot
    // read them, which is why the container is pinned to 6.x.
    //
    // The full Complete/Single/Duplicated/Fragmented/Missing breakdown is retained rather
    // than collapsed to a "BUSCO %": in a correctly phased dikaryotic assembly, HIGH
    // DUPLICATION IS EXPECTED AND CORRECT, and near-zero duplication is a sign of collapse.
    """
    ASM=${fasta}
    if [[ "\$ASM" == *.gz ]]; then
        gunzip -c "\$ASM" > asm.fa
        ASM=asm.fa
    fi

    busco \\
        --in "\$ASM" \\
        --out ${meta.id}_busco \\
        --mode genome \\
        --lineage_dataset ${busco_db}/${params.busco_lineage} \\
        --cpu ${task.cpus} \\
        --offline \\
        --download_path ${busco_db}

    cp ${meta.id}_busco/short_summary.*.json ${meta.id}.busco.json
    cp ${meta.id}_busco/short_summary.*.txt  ${meta.id}.busco.txt

    # The per-gene table, not just the summary counts. NuclearPhaser reads it directly (BUSCO
    # id, status, contig) to find single-copy genes shared between candidate haplotypes, and
    # it is the only place the actual LOCATION of each duplicated BUSCO is recorded — which is
    # the evidence for whether duplication reflects two retained nuclei or one collapsed pair.
    cp ${meta.id}_busco/run_*/full_table.tsv ${meta.id}.busco_full_table.tsv || \
        echo "WARNING: no BUSCO full_table.tsv found for ${meta.id}" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$(busco --version 2>&1 | sed 's/BUSCO //')
        busco_lineage: ${params.busco_lineage}
    END_VERSIONS
    """

    stub:
    """
    echo '{"results": {"Complete percentage": 0}}' > ${meta.id}.busco.json
    touch ${meta.id}.busco.txt ${meta.id}.busco_full_table.tsv
    echo '"${task.process}": {busco: stub}' > versions.yml
    """
}
