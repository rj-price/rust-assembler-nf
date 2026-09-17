process COVERAGE_HISTOGRAMS {
    label 'process_single'

    publishDir "${params.outdir}/assembly_qc/coverage", mode: params.publish_dir_mode

    input:
    path histograms
    path expected_depths

    output:
    path "coverage_histogram_mqc.tsv", emit: mqc
    path "coverage_histograms.svg"   , emit: plot
    path "versions.yml"              , emit: versions

    script:
    // Two outputs, deliberately: a small-multiples SVG that stands alone, and an overlaid
    // interactive line plot for MultiQC. The comparison between candidates is what carries
    // the information here, and it is invisible when the curves sit in separate files.
    //
    // Reference line per panel, not one line across the whole figure: the expectation is
    // yield/(2 x genome_size) and yield is per sample. ASSEMBLY_QC computes it for each
    // assembly and collects it here as a TSV; the global param remains the fallback for
    // single-sample projects and for qc_only, where READ_QC never runs.
    def expected = GenomeSize.expectedHaplotypeDepth(params.read_yield_bases, params.genome_size)
    def exp_arg  = expected ? "--expected-depth ${expected}" : ''
    def map_arg  = expected_depths.name.startsWith('NO_FILE') ? '' : "--expected-depth-map ${expected_depths}"
    """
    coverage_histograms.py ${histograms} \\
        --max-depth ${params.coverage_plot_max_depth} \\
        ${exp_arg} \\
        ${map_arg} \\
        --svg coverage_histograms.svg \\
        --out coverage_histogram_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch coverage_histogram_mqc.tsv coverage_histograms.svg
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
