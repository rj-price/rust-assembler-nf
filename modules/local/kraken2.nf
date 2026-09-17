process KRAKEN2 {
    tag "${meta.id}"
    label 'process_himem'

    publishDir "${params.outdir}/contamination/kraken2", mode: params.publish_dir_mode,
        // The classified/unclassified FASTQs are large; keep reports by default.
        saveAs: { fn -> fn.endsWith('.fastq.gz') ? null : fn }

    input:
    tuple val(meta), path(reads)
    path  db

    output:
    tuple val(meta), path("${meta.id}.kraken2.report.txt"), emit: report
    tuple val(meta), path("${meta.id}.kraken2.summary.tsv"), emit: summary
    path "${meta.id}.kraken2.summary_mqc.tsv"              , emit: mqc
    path "versions.yml"                                    , emit: versions

    when:
    params.run_kraken2

    script:
    // DECISION: classification is EVIDENCE, not an action. Nothing is discarded here, and
    // "plant" reads in particular are not dropped on Kraken2's say-so alone — host mapping
    // and assembly-level screening are the corroborating lines of evidence.
    //
    // --confidence defaults to 0.1 rather than Kraken2's own 0.0. At 0.0 a single matching
    // k-mer anywhere in a read wins, and per-read classification is winner-takes-all, so a
    // 15-20 kb HiFi read accumulates enough spurious hits against an over-represented genome
    // to be called confidently wrong. Raising the floor narrows that without pretending the
    // result is authoritative: FCS-GX post-assembly is still the source of truth.
    """
    zcat -f ${reads} \\
        | kraken2 \\
            --db ${db} \\
            --threads ${task.cpus} \\
            --report ${meta.id}.kraken2.report.txt \\
            --use-names \\
            --confidence ${params.kraken2_confidence} \\
            --memory-mapping \\
            /dev/stdin \\
        > /dev/null

    # Roll the report up into the broad categories that actually matter here.
    #
    # Deliberately awk, not python. The kraken2 biocontainer has no python3, so the previous
    # heredoc failed with exit 127 AFTER the 100-minute classification had already succeeded —
    # throwing away real work on a formatting step. awk is present in every container here.
    #
    # Every backslash below is DOUBLED. This whole block is a Groovy string BEFORE it is a
    # shell script, so a lone backslash-n is consumed by Groovy and reaches awk as a real
    # newline -- inside a single-quoted awk program that is the syntax error "Unexpected end
    # of string", and it cost a 3 h 28 m classification in the 2026-08-26 run.
    #
    # That applies to COMMENTS TOO, which is not obvious and cost a second run: Nextflow does
    # not strip these lines, so an escape sequence written in prose here is interpreted just
    # as the code is. Never write a backslash in this block unless you mean it. Spell escape
    # sequences out in words instead.
    #
    # The classification is also no longer allowed to die on its summary -- see below.
    awk -F'\\t' '
        {
            pct = \$1; name = \$6
            gsub(/^[ \\t]+|[ \\t]+\$/, "", name)
            if (name == "unclassified") { unc = pct }
            if (name == "Basidiomycota" || name == "Pucciniomycotina" || name == "Fungi") {
                if (pct > fungal) fungal = pct
            }
            if (name == "Viridiplantae" || name == "Streptophyta") {
                if (pct > plant) plant = pct
            }
            if (name == "Bacteria") { if (pct > bact) bact = pct }
            if (name == "Archaea")  { if (pct > arch) arch = pct }
            if (name == "Viruses")  { if (pct > viru) viru = pct }
            # Animal is NOT an expected category for a rust, which is exactly why it is
            # reported. In the 2026-08-26 run 62.03% of reads landed in Chordata — a long-read
            # Kraken2 artefact, but the old buckets covered only fungi/plant/bacteria/archaea/
            # virus and so summed to 34% while claiming nothing was wrong. A category you do
            # not print is a category you cannot notice.
            if (name == "Metazoa") { if (pct > animal) animal = pct }
            # Root anchors the catch-all below. Kraken2 reports percentages of ALL reads, so
            # "Other classified" is what root accounts for that no named bucket claimed.
            if (name == "root") { if (pct > rootpct) rootpct = pct }
        }
        END {
            other = rootpct - fungal - plant - bact - arch - viru - animal
            if (other < 0) other = 0
            printf "category\\tpercent_reads\\n"
            printf "Rust/fungal\\t%.2f\\n", fungal
            printf "Plant\\t%.2f\\n",       plant
            printf "Animal\\t%.2f\\n",      animal
            printf "Bacteria\\t%.2f\\n",    bact
            printf "Archaea\\t%.2f\\n",     arch
            printf "Viruses\\t%.2f\\n",     viru
            printf "Other classified\\t%.2f\\n", other
            printf "Unclassified\\t%.2f\\n", unc
        }
    ' ${meta.id}.kraken2.report.txt > ${meta.id}.kraken2.summary.tsv || SUMMARY_FAILED=1

    # The classification is the expensive, irreplaceable output; the summary is a cosmetic
    # roll-up of a file already on disk. Three separate runs have now lost an hour or more of
    # finished classification to a bug in this roll-up (a python heredoc in a container with no
    # python; a Groovy-eaten backslash in the awk; a Groovy-eaten backslash in a COMMENT above
    # the awk). The asymmetry is absurd, so it is now structural rather than a promise to be
    # more careful: a broken summary degrades to a placeholder and the task still succeeds.
    #
    # The report itself is emitted either way, so nothing is silently lost -- and the
    # placeholder says plainly that it needs looking at rather than reporting zeroes, which
    # would read as "no contamination found".
    if [ "\${SUMMARY_FAILED:-0}" = "1" ] || [ ! -s ${meta.id}.kraken2.summary.tsv ]; then
        echo "WARNING: kraken2 summary step failed; the classification itself SUCCEEDED." >&2
        echo "         See ${meta.id}.kraken2.report.txt for the full result." >&2
        printf 'category\\tpercent_reads\\n'            > ${meta.id}.kraken2.summary.tsv
        printf 'SUMMARY FAILED - see report\\tNA\\n'    >> ${meta.id}.kraken2.summary.tsv
    fi

    # Wrap the same numbers with MultiQC headers. Without this the buckets exist only as a TSV
    # on disk: MultiQC's native kraken module read the raw report and rendered no taxon names
    # at all, so the 62% Metazoa signal from the 2026-08-26 run appeared nowhere in the report
    # anyone actually opens. A category you do not print is one you cannot
    # notice, and that applies to WHERE it is printed too.
    {
        echo "# id: kraken2_summary"
        echo "# section_name: 'Read classification (Kraken2)'"
        echo "# description: 'Broad categories from the Kraken2 report, as a percentage of ALL"
        echo "#   reads, summing to 100%. Read-level classification of LONG reads over-calls:"
        echo "#   a 15-20 kb read accumulates enough spurious k-mer hits against heavily"
        echo "#   represented genomes (human above all) to be called confidently wrong, and"
        echo "#   Kraken2 per read is winner-takes-all. Treat this as evidence about the READS."
        echo "#   FCS-GX post-assembly is the source of truth for contamination.'"
        echo "# format: 'tsv'"
        echo "# plot_type: 'bargraph'"
        echo "# pconfig:"
        echo "#    id: 'kraken2_summary_plot'"
        echo "#    ylab: 'percent of reads'"
        printf 'Sample\\t'
        cut -f1 ${meta.id}.kraken2.summary.tsv | tail -n +2 | paste -sd '\\t' -
        printf '${meta.id}\\t'
        cut -f2 ${meta.id}.kraken2.summary.tsv | tail -n +2 | paste -sd '\\t' -
    } > ${meta.id}.kraken2.summary_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version | head -1 | sed 's/Kraken version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.kraken2.report.txt ${meta.id}.kraken2.summary.tsv
    touch ${meta.id}.kraken2.summary_mqc.tsv
    echo '"${task.process}": {kraken2: stub}' > versions.yml
    """
}
