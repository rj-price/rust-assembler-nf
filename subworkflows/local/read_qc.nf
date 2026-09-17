//
// READ_QC — characterise the HiFi reads, per run and combined.
//
// DECISION D2: no filtering by default. Runs can differ sharply in quality and read length
// (94.9% vs ~36% of bases >=Q30, 11 kb vs 27.5 kb mean, on the dataset this was built for).
// That difference is reported as an observation. Filtering happens only if the user explicitly asks, and even then it
// produces a SEPARATE dataset — the raw FASTQs are never modified (DECISION D6).
//

include { SEQKIT_STATS    } from '../../modules/local/seqkit_stats'
include { NANOPLOT        } from '../../modules/local/nanoplot'
include { SEQKIT_FILTER   } from '../../modules/local/seqkit_filter'
include { READ_QC_SUMMARY } from '../../modules/local/read_qc_summary'

workflow READ_QC {

    take:
    ch_reads_per_run   // [ meta(id,sample,run), fastq ]

    main:
    ch_versions = Channel.empty()

    SEQKIT_STATS(ch_reads_per_run)
    ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions.first())

    NANOPLOT(ch_reads_per_run)
    ch_versions = ch_versions.mix(NANOPLOT.out.versions.first())

    // Cross-run comparison table. Run ids and stats files must stay in matching order, so
    // they are sorted together and passed as one tuple rather than as three channels.
    ch_collected = SEQKIT_STATS.out.stats
        .map { meta, stats -> tuple(meta.sample, meta.run, stats) }
        .toSortedList { a, b -> a[1] <=> b[1] }   // deterministic ordering by run id
        .filter { rows -> rows.size() > 0 }
        .map { rows ->
            tuple(rows[0][0], rows.collect { it[1] }, rows.collect { it[2] })
        }

    READ_QC_SUMMARY(ch_collected)
    ch_versions = ch_versions.mix(READ_QC_SUMMARY.out.versions)

    // Per-sample sequenced yield, summed over that sample's runs. This is what makes the
    // coverage expectation per sample rather than a single global --read_yield_bases: the
    // number is already measured here, so asking the user to restate it as one figure for a
    // multi-sample run could only ever be right for one of them.
    ch_sample_yield = SEQKIT_STATS.out.stats
        .splitCsv(sep: '\t', header: true, elem: 1)
        .map { meta, row -> tuple(meta.sample, (row['sum_len'] as BigInteger)) }
        .groupTuple()
        .map { sample, sums -> tuple(sample, sums.sum()) }

    // Optional filtered dataset — inert unless --min_read_length / --min_read_q are set.
    SEQKIT_FILTER(ch_reads_per_run)
    ch_versions = ch_versions.mix(SEQKIT_FILTER.out.versions.first().ifEmpty(null))

    // Pool the per-run filtered FASTQs into one `filtered` readset, mirroring how the raw
    // reads are pooled into `all`. Filtering is per-run because quality differs sharply
    // between runs, but the assembly question is about the dataset as a whole.
    ch_filtered_combined = SEQKIT_FILTER.out.reads
        .map { meta, fq -> tuple(meta.sample, fq) }
        .groupTuple()
        .map { sample, fqs ->
            tuple([ id: "${sample}_filtered".toString(), sample: sample, readset: 'filtered' ], fqs)
        }

    emit:
    sample_yield  = ch_sample_yield          // [ sample, total sequenced bases ]
    summary_tsv   = READ_QC_SUMMARY.out.tsv
    summary_json  = READ_QC_SUMMARY.out.json
    filtered      = SEQKIT_FILTER.out.reads
    filtered_combined = ch_filtered_combined
    for_multiqc   = NANOPLOT.out.stats.map { meta, f -> f }
                        .mix(READ_QC_SUMMARY.out.mqc)
    versions      = ch_versions
}
