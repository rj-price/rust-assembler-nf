process COVERAGE_SUMMARY {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/assembly_qc/coverage/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(depth), path(coverage), val(yield_bases)

    output:
    tuple val(meta), path("${meta.id}.coverage_summary.json"), emit: json
    tuple val(meta), path("${meta.id}.coverage_summary.txt") , emit: txt
    tuple val(meta), path("${meta.id}.depth_histogram.tsv")   , emit: histogram
    path "versions.yml"                                      , emit: versions

    script:
    // Expected per-haplotype depth, DERIVED so it tracks --genome_size. The previous code
    // said it was derived and then hard-coded 67.0, which is how every coverage report in
    // run 32466998 compared a measured 57x against a stale expectation and called it low.
    //
    // A dikaryon carries two nuclear genomes, so the sequenced target is 2 x genome_size and
    // per-haplotype depth is yield / (2 x genome_size). Collapsed (shared) sequence draws
    // reads from both nuclei and so sits at twice that.
    //
    // The yield is PER SAMPLE, taken from that sample's own SEQKIT_STATS and carried in on
    // the input tuple. It used to be params.read_yield_bases, a single global figure -- fine
    // for a one-sample project, wrong for a multi-sample one. The seven Pst isolates span
    // 1.33-2.72 Gb, i.e. 8.9-18.1x expected, so one number would have scored six of the
    // seven against another sample's expectation.
    //
    // Falls back to the global param when no per-sample yield is available (qc_only skips
    // READ_QC). Null when neither is known: the modes are still detected and reported, they
    // are simply not scored. Inventing a number here is how every coverage report in one run
    // compared a measured 57x against a stale 67x.
    def expected_cov = GenomeSize.expectedHaplotypeDepth(
        yield_bases ?: params.read_yield_bases, params.genome_size)
    def expected_arg = expected_cov ? "--expected-haplotype-cov ${expected_cov}" : ''
    """
    coverage_summary.py \\
        --depth ${depth} \\
        --coverage ${coverage} \\
        --assembly-id ${meta.id} \\
        ${expected_arg} \\
        --histogram ${meta.id}.depth_histogram.tsv \\
        --outprefix ${meta.id}.coverage_summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{}' > ${meta.id}.coverage_summary.json
    touch ${meta.id}.coverage_summary.txt ${meta.id}.depth_histogram.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
