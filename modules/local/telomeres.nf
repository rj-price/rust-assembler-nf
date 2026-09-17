process TELOMERES {
    tag "${meta.id}"
    label 'process_single'

    publishDir "${params.outdir}/assembly_qc/telomeres", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(fasta), path(gfa)

    output:
    tuple val(meta), path("${meta.id}.telomeres.json"), emit: json
    tuple val(meta), path("${meta.id}.telomeres.tsv") , emit: arrays
    path "versions.yml"                               , emit: versions

    when:
    params.run_telomeres

    script:
    // Telomere-capped contig ends are a contiguity measure with a biological meaning, unlike
    // N50 (DECISION D7): a genome has a known number of chromosome ends, so the count can be
    // read against an expectation rather than only compared between assemblies.
    //
    // The interstitial count is the other half, and the reason this runs on every candidate:
    // a telomere in the middle of a contig is a mis-join. Assemblers are routinely accused of
    // this, so the pipeline measures it per assembly instead of taking the reputation on
    // trust.
    """
    telomere_scan.py ${fasta} \\
        --assembly-id ${meta.id} \\
        --motif ${params.telomere_motif} \\
        --min-units ${params.telomere_min_units} \\
        --interstitial-min-units ${params.telomere_interstitial_min_units} \\
        --end-window ${params.telomere_end_window} \\
        --arrays ${meta.id}.telomeres.tsv \\
        --out ${meta.id}.telomeres.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{}' > ${meta.id}.telomeres.json
    touch ${meta.id}.telomeres.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
