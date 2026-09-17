process IPA {
    tag "${meta.id}"
    label 'process_assembly'

    publishDir "${params.outdir}/assembly/ipa/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.ipa.p_ctg.fasta.gz"), emit: fasta
    tuple val(meta), path("${meta.id}.ipa.a_ctg.fasta.gz"), emit: alternate, optional: true
    path "versions.yml"                                   , emit: versions

    when:
    params.run_ipa

    script:
    // PacBio's Improved Phased Assembler. A fourth independent view of the same reads, and
    // the only one here whose overlapper is Pancake -- Flye, HiCanu and hifiasm all take
    // different routes but none take this one.
    //
    // IPA is unmaintained (last upstream commit 2022-03-11) and every stock biocontainer of
    // it is broken; conf/apptainer.config points this process at a locally built image and
    // containers/ipa.def explains what had to be fixed. If the container is missing, the
    // run fails at validation in main.nf rather than here.
    def njobs    = task.cpus >= 4 ? 4 : 1
    def nthreads = Math.max(1, (task.cpus / njobs) as int)
    // `ipa local` parses --input-fn with argparse's default single-value action, so a second
    // filename after it is an unrecognised argument, not a second input: "ipa: error:
    // unrecognized arguments: ...". Its own help says to repeat the flag -- `-i fn1 -i fn2`.
    // Every other assembler here takes a bare list, so `--input-fn ${reads}` looked right and
    // survived the single-file mlp98AG31 sample; job 34173415's rust_all sample has two HiFi
    // runs and died in 2 s. Built here rather than in the script body so that the one-file
    // case renders the identical string and does not invalidate a cached nine-hour task.
    def input_fn = (reads instanceof List ? reads : [reads]).collect { "--input-fn ${it}" }.join(' ')
    """
    # IPA drives snakemake, which builds a SourceCache under \$HOME/.cache at startup, and
    # \$HOME is read-only inside the container on the compute nodes. Same failure as VERKKO
    # hit in job 34035432. Point HOME and the XDG cache at the task work directory.
    export HOME="\$PWD"
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p "\$XDG_CACHE_HOME" ipa_tmp

    # --tmp-dir defaults to /tmp, which on these nodes is small and shared. The overlap sort
    # spills there and would fill it; keep the spill in the task directory instead.
    ipa local \\
        ${input_fn} \\
        --nthreads ${nthreads} \\
        --njobs ${njobs} \\
        --run-dir RUN \\
        --tmp-dir "\$PWD/ipa_tmp" \\
        ${params.ipa_extra}

    # The stage directories are numbered, and the numbering shifts with --no-polish and
    # friends, so the final directory is matched by suffix rather than named outright.
    p_ctg=\$(ls RUN/*-final/final.p_ctg.fasta 2>/dev/null | head -1)
    a_ctg=\$(ls RUN/*-final/final.a_ctg.fasta 2>/dev/null | head -1)

    if [ -z "\$p_ctg" ]; then
        echo "IPA produced no final.p_ctg.fasta -- looked in RUN/*-final/" >&2
        ls -R RUN >&2 || true
        exit 1
    fi

    # IPA names its contigs with slashes -- `ctg/p/c/000000/0`, and `hap_ctg/p/c/000018/0`
    # for the alternate. No other assembler here does. BUSCO 6.1.0 refuses such a file
    # outright ("The character \"/\" is present in the fasta header ... which will crash
    # Reader"), which is how both IPA BUSCO tasks died in job 34071020 after IPA itself had
    # spent 8 h 30 m succeeding. A slash in a sequence name is trouble well beyond BUSCO --
    # it is a path separator, so anything that derives a filename from a sequence name is at
    # risk -- so clean it here, at the source, rather than teaching each consumer to cope.
    # Headers only: the sequence lines are untouched.
    sed '/^>/ s|/|_|g' "\$p_ctg" | gzip -c > ${meta.id}.ipa.p_ctg.fasta.gz
    if [ -n "\$a_ctg" ] && [ -s "\$a_ctg" ]; then
        sed '/^>/ s|/|_|g' "\$a_ctg" | gzip -c > ${meta.id}.ipa.a_ctg.fasta.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ipa: \$(ipa --version 2>/dev/null | head -1)
    END_VERSIONS
    """

    stub:
    """
    echo ">stub" | gzip -c > ${meta.id}.ipa.p_ctg.fasta.gz
    echo ">stub" | gzip -c > ${meta.id}.ipa.a_ctg.fasta.gz
    echo '"${task.process}": {ipa: stub}' > versions.yml
    """
}
