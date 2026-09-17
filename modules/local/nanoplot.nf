process NANOPLOT {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/read_qc/nanoplot/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html")                       , emit: html
    tuple val(meta), path("*.png")                        , emit: png,  optional: true
    tuple val(meta), path("${meta.id}_NanoStats.txt")     , emit: stats
    path "versions.yml"                                   , emit: versions

    script:
    // --tsv_stats keeps the summary machine-readable for the cross-run comparison.
    // PacBio deliveries are sometimes FASTA (quality stripped), so pick the input flag off
    // the extension: --fastq on a FASTA is a hard failure, and the quality-derived panels
    // are simply absent rather than wrong.
    def input_flag = reads.toString().replaceAll(/\.gz$/, '') ==~ /.*\.(fa|fasta|fsa)$/ ? '--fasta' : '--fastq'
    """
    NanoPlot \\
        ${input_flag} ${reads} \\
        --threads ${task.cpus} \\
        --prefix ${meta.id}_ \\
        --tsv_stats \\
        --N50 \\
        --loglength \\
        --format png

    # NanoPlot's prefix handling varies between versions; normalise the name we depend on.
    if [ ! -f "${meta.id}_NanoStats.txt" ]; then
        cp \$(ls *NanoStats.txt | head -1) ${meta.id}_NanoStats.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(NanoPlot --version | sed 's/NanoPlot //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_report.html ${meta.id}_NanoStats.txt
    echo '"${task.process}": {nanoplot: stub}' > versions.yml
    """
}
