process MULTIQC {
    label 'process_low'

    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    input:
    path multiqc_files
    path versions

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data"       , emit: data, optional: true
    path "versions.yml"       , emit: versions

    script:
    """
    multiqc \\
        --force \\
        --filename multiqc_report.html \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    echo '"${task.process}": {multiqc: stub}' > versions.yml
    """
}
