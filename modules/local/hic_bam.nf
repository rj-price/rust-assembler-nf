process HIC_BAM {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(fasta), path(sam)

    output:
    tuple val(meta), path("${meta.id}.contigs.fa"),
                     path("${meta.id}.contigs.fa.fai"), emit: contigs
    tuple val(meta), path("${meta.id}.hic.bam")        , emit: bam
    path "versions.yml"                                , emit: versions

    script:
    // The bridge between HIC_ALIGN and a scaffolder. Both YaHS and HapHiC want a BAM and an
    // indexed, uncompressed contigs FASTA; HIC_ALIGN emits gzipped SAM because that is what
    // pairtools reads. Converting once here rather than inside each scaffolder means the two
    // can be compared on byte-identical input, which is the whole point of offering both.
    //
    // The FASTA is emitted alongside the BAM, decompressed and indexed, for a reason that
    // bites otherwise: a scaffolder reads the .fai to get contig lengths, and Nextflow stages
    // inputs as symlinks, so `samtools faidx` inside the scaffolder task would try to write
    // the index next to the ORIGINAL file in another task's work directory. Emitting the pair
    // together also guarantees the index cannot drift from the sequence it describes.
    //
    // -F 0x904 drops secondary (0x100), supplementary (0x800) and unmapped (0x4) records. No
    // MAPQ filter here on purpose: YaHS and HapHiC each apply their own, and baking one in
    // would silently override the setting the user thinks they are choosing.
    //
    // Sorted by read name, not coordinate. Both scaffolders read a Hi-C pair as two adjacent
    // records and neither will tell you if that assumption is broken -- they will simply find
    // fewer contacts and scaffold worse.
    """
    ASM=${fasta}
    if [[ "\$ASM" == *.gz ]]; then
        gunzip -c "\$ASM" > ${meta.id}.contigs.fa
    else
        cp -L "\$ASM" ${meta.id}.contigs.fa
    fi
    samtools faidx ${meta.id}.contigs.fa

    samtools view -h -F 0x904 -@ ${task.cpus} ${sam} \\
        | samtools sort -n -@ ${task.cpus} -m 1G -T sorttmp -o ${meta.id}.hic.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.contigs.fa ${meta.id}.contigs.fa.fai ${meta.id}.hic.bam
    echo '"${task.process}": {samtools: stub}' > versions.yml
    """
}
