process MINIMAP2_HOST {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/contamination/host_mapping", mode: params.publish_dir_mode,
        saveAs: { fn -> fn.endsWith('.fastq.gz') ? null : fn }

    input:
    tuple val(meta), path(reads)
    path  reference

    output:
    tuple val(meta), path("${meta.id}.host_mapping.tsv")      , emit: summary
    path "${meta.id}.host_mapping_mqc.tsv"                    , emit: mqc
    tuple val(meta), path("${meta.id}.host_removed.fastq.gz") , emit: host_removed, optional: true
    tuple val(meta), path("${meta.id}.host_reads.fastq.gz")   , emit: host_reads,   optional: true
    path "versions.yml"                                       , emit: versions

    script:
    // Direct alignment evidence, which is far more informative than k-mer classification for
    // deciding how much of the data is really host. DECISION D6: reads are SPLIT into two new
    // files, never removed in place — "did removing host actually help?" stays answerable.
    def do_split = params.host_filtered_assembly ? 'true' : 'false'
    """
    minimap2 -ax map-hifi -t ${task.cpus} ${reference} ${reads} \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.host.bam -

    samtools index ${meta.id}.host.bam

    total=\$(samtools view -c ${meta.id}.host.bam)
    mapped=\$(samtools view -c -F 0x904 ${meta.id}.host.bam)
    unmapped=\$(samtools view -c -f 0x4 ${meta.id}.host.bam)
    mapped_bases=\$(samtools stats ${meta.id}.host.bam | awk -F'\\t' '/^SN\\tbases mapped \\(cigar\\)/{print \$3}')
    total_bases=\$(samtools stats ${meta.id}.host.bam | awk -F'\\t' '/^SN\\ttotal length/{print \$3}')

    # awk, NOT python3. The mulled minimap2+samtools image contains neither python3 nor
    # python, so the heredoc this replaces died with exit 127 in job 34007886 after minimap2
    # had already spent an hour mapping. awk ships in the image and does this arithmetic fine.
    #
    # Note the escaped \${mapped:-0}: the shell must expand these, not Groovy. Unescaped, the
    # previous version had Groovy consume them at parse time and render a literal 0, so every
    # count in the table would have been zero -- a clean-looking table of nothing, published
    # into MultiQC. Escape every shell variable in a script block.
    awk -v mapped="\${mapped:-0}" -v unmapped="\${unmapped:-0}" \\
        -v mb="\${mapped_bases:-0}" -v tb="\${total_bases:-0}" '
        BEGIN {
            OFS = "\\t"
            denom = mapped + unmapped
            print "metric", "value"
            print "reads_mapped_to_host", mapped
            print "reads_unmapped", unmapped
            printf "pct_reads_mapped_to_host%s%.2f\\n", OFS, (denom ? 100*mapped/denom : 0)
            print "bases_mapped_to_host", mb
            printf "pct_bases_mapped_to_host%s%.2f\\n", OFS, (tb ? 100*mb/tb : 0)
        }' > ${meta.id}.host_mapping.tsv

    # Wrap the same numbers with MultiQC headers. Without this the host mapping existed only
    # as a TSV on disk: computed, published, and absent from the report anyone opens. That is
    # the hidden-category failure exactly, and it matters more here than most, because these numbers are
    # the direct-alignment counterweight to Kraken2's "62% human" -- the one line of evidence
    # that says how much of the data really is the plant it grew on.
    {
        echo "# id: host_mapping"
        echo "# section_name: 'Host mapping (minimap2)'"
        echo "# description: 'Reads aligned against the host reference the sample was"
        echo "#   harvested from. Direct alignment evidence, and far harder to argue with than"
        echo "#   read-level k-mer classification. NOTHING is filtered on these numbers: the"
        echo "#   host-removed reads are only split into a separate dataset, and only when"
        echo "#   --host_filtered_assembly asks for it.'"
        echo "# format: 'tsv'"
        echo "# plot_type: 'table'"
        printf 'Sample\\t'
        cut -f1 ${meta.id}.host_mapping.tsv | tail -n +2 | paste -sd '\\t' -
        printf '${meta.id}\\t'
        cut -f2 ${meta.id}.host_mapping.tsv | tail -n +2 | paste -sd '\\t' -
    } > ${meta.id}.host_mapping_mqc.tsv

    if [ "${do_split}" = "true" ]; then
        samtools fastq -f 0x4 -@ ${task.cpus} ${meta.id}.host.bam | gzip -c > ${meta.id}.host_removed.fastq.gz
        samtools fastq -F 0x904 -@ ${task.cpus} ${meta.id}.host.bam | gzip -c > ${meta.id}.host_reads.fastq.gz
    fi

    rm -f ${meta.id}.host.bam ${meta.id}.host.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    // The split FASTQs are stubbed too when they are asked for: without them the
    // host_filtered readset silently never reaches ASSEMBLY, and a stub run would report
    // success while proving nothing about the branch it was meant to validate.
    def do_split = params.host_filtered_assembly ? 'true' : 'false'
    """
    touch ${meta.id}.host_mapping.tsv
    # Wrap the same numbers with MultiQC headers. Without this the host mapping existed only
    # as a TSV on disk: computed, published, and absent from the report anyone opens. That is
    # the hidden-category failure exactly, and it matters more here than most, because these numbers are
    # the direct-alignment counterweight to Kraken2's "62% human" -- the one line of evidence
    # that says how much of the data really is the plant it grew on.
    {
        echo "# id: host_mapping"
        echo "# section_name: 'Host mapping (minimap2)'"
        echo "# description: 'Reads aligned against the host reference the sample was"
        echo "#   harvested from. Direct alignment evidence, and far harder to argue with than"
        echo "#   read-level k-mer classification. NOTHING is filtered on these numbers: the"
        echo "#   host-removed reads are only split into a separate dataset, and only when"
        echo "#   --host_filtered_assembly asks for it.'"
        echo "# format: 'tsv'"
        echo "# plot_type: 'table'"
        printf 'Sample\\t'
        cut -f1 ${meta.id}.host_mapping.tsv | tail -n +2 | paste -sd '\\t' -
        printf '${meta.id}\\t'
        cut -f2 ${meta.id}.host_mapping.tsv | tail -n +2 | paste -sd '\\t' -
    } > ${meta.id}.host_mapping_mqc.tsv

    if [ "${do_split}" = "true" ]; then
        echo | gzip -c > ${meta.id}.host_removed.fastq.gz
        echo | gzip -c > ${meta.id}.host_reads.fastq.gz
    fi
    echo '"${task.process}": {minimap2: stub}' > versions.yml
    """
}
