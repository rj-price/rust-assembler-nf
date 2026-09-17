process FCS_GX {
    tag "${meta.id}"
    label 'process_himem'

    publishDir "${params.outdir}/assembly_qc/fcs_gx/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(fasta), path(gfa)
    path  gx_db

    output:
    tuple val(meta), path("*.fcs_gx_report.txt")   , emit: report,  optional: true
    tuple val(meta), path("*.taxonomy.rpt")        , emit: taxonomy, optional: true
    path "versions.yml"                            , emit: versions

    when:
    params.run_fcs_gx

    script:
    // DECISION D4: FCS-GX rather than BlobToolKit — NCBI's own screen, with a single
    // database and image to provide. This is the assembly-level half of contamination screening;
    // contamination that slips past read classification is often much cleaner to see here.
    //
    // Like every other contamination step, this REPORTS. It does not clean anything.
    """
    ASM=${fasta}
    if [[ "\$ASM" == *.gz ]]; then
        gunzip -c "\$ASM" > asm.fa
        ASM=asm.fa
    fi

    # NO `|| echo WARNING` HERE. It was there until 2026-09-02, and it turned a total
    # failure into a silent success: params.fcs_gx_sif resolved to null, the task ran on the
    # bare host where /app/bin/run_gx does not exist, run_gx exited 2, the `||` swallowed it,
    # the task exited 0, and both report outputs are `optional` so Nextflow saw nothing
    # missing. FCS_GX_CLEAN then honestly reported "no FCS-GX findings" and copied each
    # input through unchanged. Five assemblies were published as contamination-screened
    # without a single base having been screened.
    #
    # A contamination screen that cannot say whether it ran is worse than no screen, because
    # the output looks the same either way. Fail instead.
    python3 /app/bin/run_gx \\
        --fasta "\$ASM" \\
        --gx-db ${gx_db} \\
        --tax-id ${params.fcs_gx_taxid} \\
        --out-dir .

    if ! ls *.fcs_gx_report.txt >/dev/null 2>&1; then
        echo "ERROR: run_gx exited 0 but wrote no *.fcs_gx_report.txt." >&2
        echo "       Nothing downstream can tell that apart from a clean assembly, so stop here." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcs_gx: \$(python3 /app/bin/run_gx --version 2>&1 | head -1 || echo 'unknown')
        fcs_gx_taxid: ${params.fcs_gx_taxid}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.fcs_gx_report.txt
    echo '"${task.process}": {fcs_gx: stub}' > versions.yml
    """
}
