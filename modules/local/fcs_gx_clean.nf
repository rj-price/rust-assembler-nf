process FCS_GX_CLEAN {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/derived/cleaned_assemblies", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(fasta), path(report)

    output:
    tuple val(meta), path("${meta.id}.cleaned.fa.gz")           , emit: fasta
    tuple val(meta), path("${meta.id}.contaminants.fa.gz")      , emit: contaminants
    tuple val(meta), path("${meta.id}.fcs_gx_clean.tsv")        , emit: manifest
    tuple val(meta), path("${meta.id}.fcs_gx_clean.json")       , emit: json
    path "versions.yml"                                         , emit: versions

    when:
    params.run_fcs_gx && params.fcs_gx_clean

    script:
    // The one stage in this pipeline that ACTS on a contamination call rather than reporting
    // it. Three things keep that safe:
    //
    //   * the source assembly is untouched — this writes a new FASTA to derived/ (D6);
    //   * only EXCLUDE and TRIM/FIX are acted on. REVIEW and INFO are recorded and left
    //     alone, because on a dikaryon "unusual" is not the same as "contaminant";
    //   * the removed sequence is kept, not deleted, so any call can be second-guessed.
    //
    // The cleaned assemblies are deliverables, not new QC candidates: feeding them back into
    // ASSEMBLY_QC would run FCS-GX on its own output. Compare cleaned against original using
    // the manifest and the removal percentages instead.
    //
    // A missing report is normal (FCS-GX found nothing, or was skipped for this candidate);
    // the script then emits a copy and says so.
    def report_arg = report.name.startsWith('NO_FILE') ? '' : "--report ${report}"
    """
    # Refuse to write over the staged input. Nextflow stages inputs as SYMLINKS, so if the
    # output name ever equals the input name, ">" follows the link and truncates the real
    # file in the publish directory. That is not hypothetical: it destroyed all 35 cleaned
    # FASTAs on 2026-09-08 when an id-normalisation change made the two names identical.
    # Cheap to check, and the failure it prevents is silent and irreversible.
    if [ "\$(basename ${fasta})" = "${meta.id}.cleaned.fa.gz" ]; then
        echo "ERROR: output ${meta.id}.cleaned.fa.gz would overwrite the input ${fasta}" >&2
        echo "       (input is a symlink to the source; writing would destroy it)" >&2
        exit 1
    fi

    fcs_gx_clean.py \\
        --fasta ${fasta} \\
        ${report_arg} \\
        --assembly-id ${meta.id} \\
        --out-fasta ${meta.id}.cleaned.fa.gz \\
        --out-removed ${meta.id}.contaminants.fa.gz \\
        --out-manifest ${meta.id}.fcs_gx_clean.tsv \\
        --out-json ${meta.id}.fcs_gx_clean.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    # Refuse to write over the staged input. Nextflow stages inputs as SYMLINKS, so if the
    # output name ever equals the input name, ">" follows the link and truncates the real
    # file in the publish directory. That is not hypothetical: it destroyed all 35 cleaned
    # FASTAs on 2026-09-08 when an id-normalisation change made the two names identical.
    # Cheap to check, and the failure it prevents is silent and irreversible.
    if [ "\$(basename ${fasta})" = "${meta.id}.cleaned.fa.gz" ]; then
        echo "ERROR: output ${meta.id}.cleaned.fa.gz would overwrite the input ${fasta}" >&2
        echo "       (input is a symlink to the source; writing would destroy it)" >&2
        exit 1
    fi

    echo | gzip -c > ${meta.id}.cleaned.fa.gz
    echo | gzip -c > ${meta.id}.contaminants.fa.gz
    touch ${meta.id}.fcs_gx_clean.tsv
    echo '{}' > ${meta.id}.fcs_gx_clean.json
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
