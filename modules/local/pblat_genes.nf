process PBLAT_GENES {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/phasing/gene_mapping", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(fasta)
    path  genes

    output:
    tuple val(meta), path("${meta.id}.gene_mapping.psl"), emit: psl
    path "versions.yml"                                 , emit: versions

    script:
    // NuclearPhaser's documented gene-mapping step uses BioKanga blitz. BioKanga is
    // unmaintained and has no container anywhere (checked against the Galaxy depot: no
    // biokanga image at all), which would otherwise make this whole branch unbuildable on a
    // containers-first pipeline.
    //
    // It does not have to be BioKanga. NuclearPhaser reads the mapping file positionally in
    // GeneBinning.read_gene_mapping(): column 0 is the match score, column 9 the gene name
    // and column 13 the contig — which is plain PSL. `blitz` emits PSL, and so does pblat, a
    // maintained, containerised, threaded BLAT. The file NuclearPhaser sees is the same shape
    // either way.
    //
    // The header PSL normally carries would be read as data, so -noHead is required, not
    // cosmetic.
    """
    ASM=${fasta}
    if [[ "\$ASM" == *.gz ]]; then
        gunzip -c "\$ASM" > asm.fa
        ASM=asm.fa
    fi

    GENES=${genes}
    if [[ "\$GENES" == *.gz ]]; then
        gunzip -c "\$GENES" > genes.fa
        GENES=genes.fa
    fi

    pblat \\
        -threads=${task.cpus} \\
        -noHead \\
        "\$ASM" \\
        "\$GENES" \\
        ${meta.id}.gene_mapping.psl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pblat: \$(pblat 2>&1 | grep -io 'pblat[^ ]* v[^ ]*' | head -1 || echo 'unknown')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.gene_mapping.psl
    echo '"${task.process}": {pblat: stub}' > versions.yml
    """
}
