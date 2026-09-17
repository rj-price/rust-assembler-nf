process HIC_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(fasta)
    tuple val(hic_meta), path(hic_r1), path(hic_r2)

    output:
    tuple val(meta), path("${meta.id}.hic.sam.gz"), emit: sam
    tuple val(meta), path("${meta.id}.chrom.sizes"), emit: chrom_sizes
    path "versions.yml"                            , emit: versions

    script:
    // Hi-C reads aligned against ONE candidate assembly. This is per-assembly rather than
    // once for the run: contact coordinates are meaningless against a different set of contigs.
    //
    // -SP5M is not optional, and getting it wrong fails quietly rather than loudly:
    //   -S -P  stop bwa applying paired-end rescue and proper-pair logic, which assume a
    //          fragment-length distribution that Hi-C data does not have;
    //   -5     reports the 5'-most alignment as primary, which is what places a ligation
    //          junction;
    //   -M     marks split hits as secondary.
    // Without them the trans contacts degrade — and trans contacts are the entire phasing
    // signal, so the run would finish successfully and simply phase badly.
    //
    // Kept separate from pair-calling and matrix-building so that changing the matrix
    // resolution or the MAPQ threshold does not re-run the alignment, which is the expensive
    // part of this branch.
    """
    ASM=${fasta}
    if [[ "\$ASM" == *.gz ]]; then
        gunzip -c "\$ASM" > asm.fa
        ASM=asm.fa
    fi

    # Contig lengths, computed with awk because the bwa image carries no samtools.
    awk '/^>/ { if (name != "") print name "\\t" len; name = substr(\$1, 2); len = 0; next }
         { len += length(\$0) }
         END { if (name != "") print name "\\t" len }' "\$ASM" > ${meta.id}.chrom.sizes

    bwa index "\$ASM"

    bwa mem -SP5M -t ${task.cpus} "\$ASM" ${hic_r1} ${hic_r2} \\
        | gzip -c > ${meta.id}.hic.sam.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -i '^Version' | sed 's/Version: //')
    END_VERSIONS
    """

    stub:
    """
    echo | gzip -c > ${meta.id}.hic.sam.gz
    touch ${meta.id}.chrom.sizes
    echo '"${task.process}": {bwa: stub}' > versions.yml
    """
}
