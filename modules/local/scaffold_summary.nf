//
// SCAFFOLD_SUMMARY — the scaffolding result as a table, deliberately SEPARATE from the
// candidate comparison.
//
// WHY A SECOND TABLE rather than extra rows in assembly_summary.tsv. A scaffold N50 is raised
// by joins alone. For mlp98AG31 hap1 the contig N50 is unchanged at 4.96 Mb while the scaffold
// N50 is 5.60 Mb, off four joins totalling 800 bp of N. (4.96 Mb is the N50 of the CLEANED
// assembly, which is what scaffolding receives -- 59 contigs and 104.40 Mb, after FCS-GX
// dropped 232 bacterial contigs. assembly_summary.tsv quotes 4.59 Mb for the same haplotype
// because it reports the assembly BEFORE cleaning; the two are not the same sequence set.)
// Listing that beside unscaffolded
// candidates in one table invites exactly the N50 ranking DECISION D7 forbids, and
// subworkflows/local/scaffolding.nf states as a design rule that scaffolds never re-enter the
// comparison. So they stay apart, and every row here carries its own contig N50 — the only
// comparison that means anything is scaffold vs contig WITHIN a row.
//
// WHAT THIS CLOSES. Until now the number this whole branch exists to produce — 18
// chromosome-scale scaffolds, the figure Duplessis et al. (2026) quote — appeared in no
// summary anywhere. gfastats reports '# scaffolds: 55' for that same haplotype, counting 37
// small unplaced ones alongside the 18 that matter, and the published run had to be measured
// by hand from the FASTA to state the result.
//

process SCAFFOLD_RECORD {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(fasta), path(gfastats)

    output:
    path "${meta.id}.scaffold_record.json", emit: record
    path "versions.yml"                   , emit: versions

    script:
    def meta_json = groovy.json.JsonOutput.toJson([
        id           : meta.id,
        sample       : meta.sample,
        assembler    : meta.assembler,
        assembly_type: meta.assembly_type,
        readset      : meta.readset
    ])
    """
    assembly_summary.py scaffold-record \\
        --meta '${meta_json}' \\
        --fasta ${fasta} \\
        --gfastats ${gfastats} \\
        --scaffolder ${params.scaffolder} \\
        --chromosome-min-length ${params.scaffold_chromosome_min_length} \\
        --out ${meta.id}.scaffold_record.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"assembly_id": "${meta.id}"}' > ${meta.id}.scaffold_record.json
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

process SCAFFOLD_SUMMARY {
    label 'process_single'

    publishDir "${params.outdir}/scaffolding", mode: params.publish_dir_mode

    input:
    path records

    output:
    path "scaffold_summary.tsv"     , emit: tsv
    path "scaffold_summary.json"    , emit: json
    path "scaffold_summary_mqc.tsv" , emit: mqc
    path "versions.yml"             , emit: versions

    script:
    """
    assembly_summary.py scaffold-merge ${records}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch scaffold_summary.tsv scaffold_summary.json scaffold_summary_mqc.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
