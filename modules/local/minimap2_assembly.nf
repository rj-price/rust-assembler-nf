process MINIMAP2_ASSEMBLY {
    tag "${meta.id}"
    label 'process_high'

    publishDir path: { "${params.outdir}/assembly_qc/coverage/${meta.id}" }, mode: params.publish_dir_mode,
        // BAMs are enormous and reproducible; keep the derived depth/coverage tables.
        saveAs: { fn -> (fn.endsWith('.bam') || fn.endsWith('.bai')) ? null : fn }

    input:
    tuple val(meta), path(fasta), path(gfa), path(reads)

    output:
    tuple val(meta), path("${meta.id}.coverage.txt"), emit: coverage
    tuple val(meta), path("${meta.id}.depth.tsv.gz"), emit: depth
    path "versions.yml"                             , emit: versions

    script:
    """
    ASM=${fasta}
    if [[ "\$ASM" == *.gz ]]; then
        gunzip -c "\$ASM" > asm.fa
        ASM=asm.fa
    fi

    # samtools sort's -m is PER THREAD and defaults to 768 MB, so the sort buffer alone was
    # ~24 GB at 32 threads and grew silently with any change to task.cpus. That, plus the
    # minimap2 index, is the whole 52 GB peak this process was measured at. Bound it: sorting
    # is I/O-bound well before 32 threads, so half the cores at 1 GB caps the buffer at 16 GB
    # and makes the memory request predictable instead of a function of the CPU request.
    # An `x && y` list is deliberately avoided here: Nextflow runs the script under `set -e`,
    # so a false test would return non-zero and abort the task in the ordinary case.
    SORT_THREADS=\$(( ${task.cpus} / 2 ))
    if [ "\$SORT_THREADS" -lt 1 ]; then
        SORT_THREADS=1
    fi

    minimap2 -ax map-hifi -t ${task.cpus} "\$ASM" ${reads} \\
        | samtools sort -@ "\$SORT_THREADS" -m 1G -o ${meta.id}.bam -

    samtools index -@ ${task.cpus} ${meta.id}.bam

    # Per-contig breadth/depth summary.
    samtools coverage ${meta.id}.bam > ${meta.id}.coverage.txt

    # Per-base depth, binned to keep the file tractable on a ~1 Gb assembly. This is what
    # the coverage-mode estimation consumes.
    samtools depth -a ${meta.id}.bam \\
        | awk 'NR % 100 == 0 { print \$1"\\t"\$2"\\t"\$3 }' \\
        | gzip -c > ${meta.id}.depth.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.coverage.txt
    echo "" | gzip -c > ${meta.id}.depth.tsv.gz
    echo '"${task.process}": {minimap2: stub}' > versions.yml
    """
}
