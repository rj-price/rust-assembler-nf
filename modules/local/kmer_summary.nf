process KMER_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    publishDir "${params.outdir}/kmer", mode: params.publish_dir_mode

    input:
    // stageAs into subdirectories, not for tidiness. GenomeScope2's summary and smudgeplot's
    // both end in _summary.txt, and in the stub profile both are literally
    // <id>_summary.txt, which fails the task before it starts:
    //     Process `KMER_ANALYSIS:KMER_SUMMARY` input file name collision -- There are
    //     multiple input files for each of the following file names: rust_all_summary.txt
    // A subdirectory per input makes the collision impossible regardless of what upstream
    // names its output, and `*` preserves the original basename so the NO_FILE placeholder
    // check below still sees the real name.
    tuple val(meta), path(summary), path(model, stageAs: 'gs_model/*'),
                                    path(smudge_summary, stageAs: 'smudge/*')

    output:
    tuple val(meta), path("${meta.id}.genomescope.json"), emit: json
    path "${meta.id}.genomescope2_mqc.tsv"              , emit: mqc
    path "versions.yml"                                 , emit: versions

    script:
    // GenomeScope2's estimates drive every size judgement downstream, but its own summary.txt
    // is a human-readable block that no MultiQC module reads — so before this the numbers
    // existed on disk and appeared nowhere in the report. The diploid-model caveat travels
    // with the table rather than living only in a config comment.
    //
    // The MQC filename is namespaced by sample. It used to be a bare genomescope2_mqc.tsv,
    // which works only for a single-sample run: with more than one readset the files collide
    // when MultiQC stages them (and overwrite each other in publishDir). MultiQC still merges
    // them into ONE table because the section identity comes from the `# id: genomescope2`
    // header inside the file, not from its name.
    //
    // model.txt and smudgeplot's verbose summary are here so the script can CHECK the fit
    // rather than only report it. GenomeScope2 gave 48 Mb for mlp98AG31 against 102-104 Mb
    // of assembled haplotype, with a failed fit that nothing surfaced; the cross-check that
    // catches it is kmercov against smudgeplot's independent 1n coverage. Both are optional
    // outputs upstream, so both arrive as placeholders when absent and their checks are
    // skipped -- see bin/kmer_mqc.py.
    def model_arg  = model.name.startsWith('NO_FILE') ? '' : "--model ${model}"
    def smudge_arg = smudge_summary.name.startsWith('NO_FILE') ? '' : "--smudgeplot-summary ${smudge_summary}"
    """
    kmer_mqc.py \\
        --summary ${summary} \\
        ${model_arg} \\
        ${smudge_arg} \\
        --sample-id ${meta.id} \\
        --out-json ${meta.id}.genomescope.json \\
        --out-mqc ${meta.id}.genomescope2_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{}' > ${meta.id}.genomescope.json
    touch ${meta.id}.genomescope2_mqc.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
