process ASSEMBLY_SUMMARY {
    label 'process_single'

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path records

    output:
    path "assembly_summary.tsv"     , emit: tsv
    path "assembly_summary.json"    , emit: json
    path "assembly_summary_mqc.tsv" , emit: mqc
    path "versions.yml"             , emit: versions

    script:
    """
    assembly_summary.py merge ${records}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch assembly_summary.tsv assembly_summary.json assembly_summary_mqc.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
