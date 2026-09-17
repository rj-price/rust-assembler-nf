//
// CONTAMINATION (read level).
//
// Mandatory stage, not optional QC — contamination is assumed possible from the outset.
// Two complementary lines of evidence at read level, and a third (FCS-GX) after assembly,
// because contamination that survives read classification is often much cleaner to spot
// once contigs exist.
//
// Nothing is discarded automatically. Kraken2 alone is not grounds for dropping "plant"
// reads; host mapping either corroborates it or does not.
//

include { KRAKEN2       } from '../../modules/local/kraken2'
include { MINIMAP2_HOST } from '../../modules/local/minimap2_host'

workflow CONTAMINATION {

    take:
    ch_reads   // [ meta(id,sample,readset), [ fastqs ] ]

    main:
    ch_versions = Channel.empty()

    // A. k-mer/taxonomic classification.
    // The DB path is nullable (site-specific, no portable default). KRAKEN2's own `when:`
    // gates execution, but the input is still evaluated, and file(null) is a hard error —
    // so stage a per-slot placeholder when it is unset. Per-slot, not shared: one shared
    // NO_FILE staged twice into the same process collides on the staging name.
    ch_kraken_db = params.kraken_db
        ? file(params.kraken_db, checkIfExists: false)
        : file("${projectDir}/assets/NO_FILE_kraken_db")

    KRAKEN2(ch_reads, ch_kraken_db)
    ch_versions = ch_versions.mix(KRAKEN2.out.versions.ifEmpty(null))

    // B. Host-reference mapping. Inert until a reference is supplied. For a rust this is the
    //    plant it was harvested from, e.g. Phaseolus vulgaris for a bean rust.
    ch_host_ref = params.host_reference
        ? Channel.value(file(params.host_reference, checkIfExists: true))
        : Channel.empty()

    MINIMAP2_HOST(ch_reads, ch_host_ref)
    ch_versions = ch_versions.mix(MINIMAP2_HOST.out.versions.ifEmpty(null))

    // Relabel the host-removed reads as their own readset so ASSEMBLY can treat them exactly
    // like any other input. This is what makes "did removing host actually help?" a question
    // the pipeline answers rather than one it invites. Inert unless --host_filtered_assembly
    // is set, since MINIMAP2_HOST only writes the split FASTQs when it is.
    ch_host_filtered = MINIMAP2_HOST.out.host_removed
        .map { meta, fq ->
            tuple(
                [ id: "${meta.sample}_host_filtered".toString(),
                  sample: meta.sample,
                  readset: 'host_filtered' ],
                [ fq ]
            )
        }

    emit:
    kraken_report = KRAKEN2.out.report
    kraken_summary = KRAKEN2.out.summary
    host_summary  = MINIMAP2_HOST.out.summary
    host_removed  = MINIMAP2_HOST.out.host_removed
    host_filtered_reads = ch_host_filtered   // [ meta(readset: host_filtered), [ fastq ] ]
    // Kraken2's own report is parsed natively by MultiQC, so it needs no conversion step.
    // Both: the raw report for MultiQC's native kraken module, and our own bucket table,
    // which is the one that actually names what was found.
    // Host mapping joins them. It was previously emitted as `host_summary` and consumed by
    // nobody, so the numbers reached contamination/host_mapping/ and stopped there.
    for_multiqc   = KRAKEN2.out.report.map { meta, f -> f }
                        .mix(KRAKEN2.out.mqc)
                        .mix(MINIMAP2_HOST.out.mqc)
    versions      = ch_versions
}
