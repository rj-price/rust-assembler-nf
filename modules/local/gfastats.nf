process GFASTATS {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/assembly_qc/gfastats", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(fasta), path(gfa)

    output:
    tuple val(meta), path("${meta.id}.gfastats.txt"), emit: stats
    path "versions.yml"                             , emit: versions

    script:
    // Basic contiguity/composition stats. NOTE (DECISION D7): these feed the report but
    // nothing downstream ranks assemblies by N50.
    """
    gfastats ${fasta} > ${meta.id}.gfastats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfastats: \$(gfastats --version 2>&1 | sed 's/gfastats //' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.gfastats.txt
    echo '"${task.process}": {gfastats: stub}' > versions.yml
    """
}
