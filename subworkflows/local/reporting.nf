//
// REPORTING — machine-readable summaries suitable for automated comparison.
//
// The deliverable of this pipeline is not a FASTA, it is the evidence needed to decide which
// assembly is biologically credible. assembly_summary.tsv/.json is that evidence in a form
// that can be diffed, scripted against, and carried into a project write-up.
//

include { ASSEMBLY_RECORD  } from '../../modules/local/assembly_record'
include { ASSEMBLY_SUMMARY } from '../../modules/local/assembly_summary'
include { COVERAGE_HISTOGRAMS } from '../../modules/local/coverage_histograms'
include { MULTIQC          } from '../../modules/local/multiqc'

workflow REPORTING {

    take:
    ch_for_summary   // [ meta, fasta, gfastats, busco, qv, coverage, telomeres, fcs_clean, mito ]
    ch_depth_hists   // per-assembly depth histogram TSVs
    ch_expected_depths // "assembly_id\texpected depth" rows, one per assembly
    ch_multiqc_files
    ch_versions

    main:
    // Nextflow's join(remainder: true) yields nulls for missing optional QC outputs.
    // Substitute a placeholder so the process still stages, and let the script drop them.
    //
    // One placeholder PER SLOT, not one shared placeholder. Nextflow stages inputs by
    // filename, so if two optional outputs were missing for the same assembly, a single
    // shared NO_FILE was staged twice and the task died with
    //     input file name collision -- multiple input files for NO_FILE
    // That is not an edge case here: whenever an assembly lacks both a QV and a coverage
    // summary (exactly what a node fault produces) the record for it cannot be built at all.
    def placeholders = ['gfastats', 'busco', 'qv', 'coverage', 'telomeres', 'fcs_clean', 'mito']
                           .collect { file("${projectDir}/assets/NO_FILE_${it}") }

    ch_records_in = ch_for_summary.map { items ->
        def meta  = items[0]
        def fasta = items[1]                    // never null
        def opts  = (2..8).collect { i -> items[i] ?: placeholders[i - 2] }
        [meta, fasta] + opts
    }

    ASSEMBLY_RECORD(ch_records_in)

    ASSEMBLY_SUMMARY(ASSEMBLY_RECORD.out.record.collect())

    // Overlaid depth curves for all candidates on one axis, each panel carrying its own
    // sample's expected depth. ifEmpty keeps the plot running when no expectation is known
    // (no yield measured and no --read_yield_bases): the curves are the point, the reference
    // line is the annotation.
    ch_expected_tsv = ch_expected_depths
        .collectFile(name: 'expected_depths.tsv', newLine: true, sort: true)
        .ifEmpty(file("${projectDir}/assets/NO_FILE_expected_depths"))

    COVERAGE_HISTOGRAMS(ch_depth_hists.collect(), ch_expected_tsv)

    MULTIQC(
        ch_multiqc_files
            .mix(ASSEMBLY_SUMMARY.out.mqc)
            .mix(COVERAGE_HISTOGRAMS.out.mqc)
            .collect(),
        ch_versions
    )

    emit:
    summary_tsv = ASSEMBLY_SUMMARY.out.tsv
    report      = MULTIQC.out.report
}
