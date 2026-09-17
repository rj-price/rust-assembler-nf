process MERQURY {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/assembly_qc/merqury/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(fasta), path(gfa), path(meryl_db)

    output:
    tuple val(meta), path("${meta.id}.qv")            , emit: qv,       optional: true
    tuple val(meta), path("${meta.id}.completeness.stats"), emit: completeness, optional: true
    tuple val(meta), path("*.png")                    , emit: plots,    optional: true
    path "versions.yml"                               , emit: versions

    when:
    params.run_merqury

    script:
    // Reference-free QV and k-mer completeness — the tool that distinguishes "highly
    // contiguous but erroneous" from "accurate and complete". Reuses the meryl database
    // already built during KMER_ANALYSIS rather than recounting 60 Gb of reads.
    """
    ASM=${fasta}
    if [[ "\$ASM" == *.gz ]]; then
        gunzip -c "\$ASM" > asm.fa
        ASM=asm.fa
    fi

    # merqury.sh sources its helpers from \$MERQURY and the biocontainer does not set it.
    # Unset, it looks for '/eval/spectra-cn.sh' and produces NO output while still exiting 0 —
    # a silent failure that leaves every QV blank.
    export MERQURY="\${MERQURY:-/usr/local/share/merqury}"
    if [ ! -x "\${MERQURY}/eval/spectra-cn.sh" ]; then
        echo "ERROR: MERQURY=\${MERQURY} does not contain eval/spectra-cn.sh" >&2
        exit 1
    fi

    merqury.sh ${meryl_db} "\$ASM" ${meta.id} || {
        echo "WARNING: merqury returned non-zero; check outputs" >&2
    }

    # merqury exits 0 even when it produced nothing, so verify the QV actually exists.
    if [ ! -s "${meta.id}.qv" ] && [ -z "\$(ls *.qv 2>/dev/null)" ]; then
        echo "ERROR: merqury produced no QV output — see logs/" >&2
        cat logs/* >&2 2>/dev/null || true
        exit 1
    fi

    # Normalise the output names we depend on downstream.
    [ -f ${meta.id}.qv ] || { f=\$(ls *.qv 2>/dev/null | head -1); [ -n "\$f" ] && cp "\$f" ${meta.id}.qv; }
    [ -f ${meta.id}.completeness.stats ] || { f=\$(ls *completeness.stats 2>/dev/null | head -1); [ -n "\$f" ] && cp "\$f" ${meta.id}.completeness.stats; }

    # NOTE: this container ships busybox grep, which has no -P. Keep version capture to
    # portable POSIX tools or it silently yields an empty version string.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merqury: \$(merqury.sh 2>&1 | sed -n 's/.*[Mm]erqury *v\\{0,1\\}\\([0-9][0-9.]*\\).*/\\1/p' | head -1 || true)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.qv ${meta.id}.completeness.stats
    echo '"${task.process}": {merqury: stub}' > versions.yml
    """
}
