process GFA2FASTA {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/assembly/${meta.assembler}/${meta.sample}_${meta.readset}/fasta",
        mode: params.publish_dir_mode

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("${meta.id}.fa.gz"), emit: fasta
    path "versions.yml"                      , emit: versions

    script:
    // GFA segment lines carry the sequence in field 3. Records with '*' there carry no
    // sequence at all (hifiasm's .noseq graphs) — emitting them would produce a FASTA that
    // looks valid but contains nothing, so they are skipped and their absence is fatal
    // rather than silent.
    """
    awk '\$1 == "S" && \$3 != "*" { print ">"\$2"\\n"\$3 }' ${gfa} \\
        | fold -w 60 \\
        | gzip -c > ${meta.id}.fa.gz

    n_seqs=\$(zcat ${meta.id}.fa.gz | grep -c '^>' || true)
    if [ "\${n_seqs}" -eq 0 ]; then
        echo "ERROR: no sequences extracted from ${gfa} — is it a sequence-free (.noseq) graph?" >&2
        exit 1
    fi
    echo "Extracted \${n_seqs} sequences from ${gfa}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    echo ">stub" | gzip -c > ${meta.id}.fa.gz
    echo '"${task.process}": {awk: stub}' > versions.yml
    """
}
