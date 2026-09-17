//
// KMER_ANALYSIS — reference-free genome characterisation, BEFORE assembly.
//
// This runs first because hifiasm makes coverage-based decisions when separating homozygous
// from heterozygous sequence, and a bad inferred homozygous-coverage threshold produces an
// assembly that is substantially too large or too small. Knowing the k-mer spectrum up front
// lets us spot that before committing to a multi-day run.
//
// The meryl database is also reused later by Merqury for QV, so it is emitted, not discarded.
//

include { MERYL_COUNT     } from '../../modules/local/meryl_count'
include { MERYL_HISTOGRAM } from '../../modules/local/meryl_histogram'
include { GENOMESCOPE2    } from '../../modules/local/genomescope2'
include { SMUDGEPLOT_CUTOFF } from '../../modules/local/smudgeplot_cutoff'
include { MERYL_PRINT     } from '../../modules/local/meryl_print'
include { SMUDGEPLOT      } from '../../modules/local/smudgeplot'
include { KMER_SUMMARY    } from '../../modules/local/kmer_summary'

workflow KMER_ANALYSIS {

    take:
    ch_reads   // [ meta(id,sample,readset), [ fastqs ] ]

    main:
    ch_versions = Channel.empty()

    MERYL_COUNT(ch_reads)
    ch_versions = ch_versions.mix(MERYL_COUNT.out.versions)

    MERYL_HISTOGRAM(MERYL_COUNT.out.meryl_db)
    ch_versions = ch_versions.mix(MERYL_HISTOGRAM.out.versions)

    GENOMESCOPE2(MERYL_HISTOGRAM.out.hist)
    ch_versions = ch_versions.mix(GENOMESCOPE2.out.versions)

    // Ploidy structure check that does not assume a diploid model.
    //
    // Three processes for what smudgeplot documents as one script, because the middle step is
    // meryl's and the outer two are smudgeplot's, and no biocontainer holds both. Splitting on
    // the container boundary is the only honest option.
    SMUDGEPLOT_CUTOFF(MERYL_HISTOGRAM.out.hist)
    ch_versions = ch_versions.mix(SMUDGEPLOT_CUTOFF.out.versions.ifEmpty(null))

    MERYL_PRINT(MERYL_COUNT.out.meryl_db.join(SMUDGEPLOT_CUTOFF.out.cutoffs))
    ch_versions = ch_versions.mix(MERYL_PRINT.out.versions.ifEmpty(null))

    SMUDGEPLOT(MERYL_PRINT.out.kmers)
    ch_versions = ch_versions.mix(SMUDGEPLOT.out.versions.ifEmpty(null))

    // Lift GenomeScope2's estimates into the report. Without this the genome size,
    // heterozygosity and repeat content the whole pipeline is scored against appeared
    // nowhere except a text file in the output tree.
    //
    // The summary is also where the fit is CHECKED, not just reported, so model.txt and
    // smudgeplot's verbose summary are joined in. Both are optional outputs — GenomeScope2
    // does not always write a model, and smudgeplot can fail or be skipped — so join on
    // remainder and substitute a placeholder, which the module reads as "skip that check".
    ch_kmer_summary_in = GENOMESCOPE2.out.summary
        .join(GENOMESCOPE2.out.model, remainder: true)
        .join(SMUDGEPLOT.out.summary, remainder: true)
        .map { meta, summary, model, smudge ->
            tuple(meta,
                  summary,
                  model ?: file("${projectDir}/assets/NO_FILE_gs_model"),
                  smudge ?: file("${projectDir}/assets/NO_FILE_smudge_summary"))
        }
        .filter { meta, summary, model, smudge -> summary != null }

    KMER_SUMMARY(ch_kmer_summary_in)
    ch_versions = ch_versions.mix(KMER_SUMMARY.out.versions.ifEmpty(null))

    emit:
    meryl_db  = MERYL_COUNT.out.meryl_db      // reused by Merqury
    hist      = MERYL_HISTOGRAM.out.hist
    summary   = GENOMESCOPE2.out.summary
    summary_json = KMER_SUMMARY.out.json
    for_multiqc  = KMER_SUMMARY.out.mqc
    versions  = ch_versions
}
