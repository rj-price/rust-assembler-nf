//
// ASSEMBLY_QC — the same evidence gathered for every candidate.
//
// This subworkflow is deliberately blind to where an assembly came from: it consumes the
// ASSEMBLIES channel and nothing else. That is what makes adding Hi-C haplotypes, Verkko, or
// the ALL/run1/run2 subsets free — they inherit the whole QC suite without touching this file.
//

include { GFASTATS          } from '../../modules/local/gfastats'
include { BUSCO             } from '../../modules/local/busco'
include { MERQURY           } from '../../modules/local/merqury'
include { MINIMAP2_ASSEMBLY } from '../../modules/local/minimap2_assembly'
include { COVERAGE_SUMMARY  } from '../../modules/local/coverage_summary'
include { HAPLOTYPE_COMPETITION } from '../../modules/local/haplotype_competition'
include { FCS_GX            } from '../../modules/local/fcs_gx'
include { FCS_GX_CLEAN      } from '../../modules/local/fcs_gx_clean'
include { TELOMERES         } from '../../modules/local/telomeres'
include { MITO_ALIGN        } from '../../modules/local/mito_screen'
include { MITO_CALL         } from '../../modules/local/mito_screen'

workflow ASSEMBLY_QC {

    take:
    ch_assemblies   // [ meta(id,sample,assembler,assembly_type,readset), fasta, gfa ]
    ch_reads        // [ meta, [ fastqs ] ]
    ch_meryl_db     // [ meta, meryl_db ]
    ch_yield_map    // value: [ sample: total bases ] — empty map when READ_QC did not run

    main:
    ch_versions = Channel.empty()

    // Assembly-level contamination (the second pass; read-level was the first). This runs
    // FIRST, on the raw assemblies, because everything below measures what it leaves.
    ch_fcs_gx_db = params.fcs_gx_db
        ? file(params.fcs_gx_db, checkIfExists: false)
        : file("${projectDir}/assets/NO_FILE_fcs_gx_db")

    FCS_GX(ch_assemblies, ch_fcs_gx_db)
    ch_versions = ch_versions.mix(FCS_GX.out.versions.first().ifEmpty(null))

    // ...and the one place that ACTS on a contamination call. FCS-GX's report is optional
    // output (no findings, or the tool disabled), so join on remainder and substitute the
    // per-slot placeholder — the cleaning script treats "no report" as "nothing to remove".
    ch_to_clean = ch_assemblies
        .map { meta, fasta, gfa -> tuple(meta, fasta) }
        .join(FCS_GX.out.report, remainder: true)
        .map { meta, fasta, report ->
            tuple(meta, fasta, report ?: file("${projectDir}/assets/NO_FILE_fcs_report"))
        }
        .filter { meta, fasta, report -> fasta != null }

    FCS_GX_CLEAN(ch_to_clean)
    ch_versions = ch_versions.mix(FCS_GX_CLEAN.out.versions.first().ifEmpty(null))

    // Mitochondrial removal, downstream of contamination cleaning because it is a different
    // question: FCS-GX asks whether a sequence is from another organism, and the
    // mitochondrion is not. See modules/local/mito_screen.nf.
    ch_cleaned  = FCS_GX_CLEAN.out.fasta
    ch_mito_tsv = Channel.empty()
    ch_mito_json = Channel.empty()

    if (params.run_mito_screen) {
        ch_mito_reference = file(params.mito_reference, checkIfExists: true)

        MITO_ALIGN(FCS_GX_CLEAN.out.fasta, ch_mito_reference)
        ch_versions = ch_versions.mix(MITO_ALIGN.out.versions.first().ifEmpty(null))

        MITO_CALL(FCS_GX_CLEAN.out.fasta.join(MITO_ALIGN.out.hits))
        ch_versions = ch_versions.mix(MITO_CALL.out.versions.first().ifEmpty(null))

        ch_cleaned   = MITO_CALL.out.fasta
        ch_mito_tsv  = MITO_CALL.out.manifest
        ch_mito_json = MITO_CALL.out.json
    }

    // QC measures the CLEANED assembly, because that is the one carried forward -- into
    // NuclearPhaser, the scaffolder, and whatever gets submitted. Measuring the raw one
    // reported M. larici-populina hifiasm hap1 as 291 contigs / 136.03 Mb and flagged it
    // `oversized`, when 232 of those contigs were bacterial and the assembly actually handed
    // on was 59 contigs / 104.40 Mb. How much cleaning removed is still in the record, from
    // the cleaning JSONs below, so nothing is lost by the switch.
    //
    // join(remainder: true) falls back to the raw FASTA per assembly when cleaning is off or
    // failed for it, so a cleaning fault degrades a row instead of dropping it -- and the
    // record's qc_input column says which one was measured. With cleaning off, ch_cleaned is
    // empty and every tuple is exactly what it was before, so existing task hashes hold.
    ch_qc = ch_assemblies
        .join(ch_cleaned, remainder: true)
        .filter { meta, raw, gfa, clean -> raw != null }
        .map { meta, raw, gfa, clean -> tuple(meta, clean ?: raw, gfa) }

    // Contiguity and composition.
    GFASTATS(ch_qc)
    ch_versions = ch_versions.mix(GFASTATS.out.versions.first().ifEmpty(null))

    // Gene-space completeness. The full breakdown is kept — duplication may be real biology.
    // Nullable site path; see the note in subworkflows/local/contamination.nf.
    ch_busco_db = params.busco_db
        ? file(params.busco_db, checkIfExists: false)
        : file("${projectDir}/assets/NO_FILE_busco_db")

    BUSCO(ch_qc, ch_busco_db)
    ch_versions = ch_versions.mix(BUSCO.out.versions.first().ifEmpty(null))

    // Reference-free QV and k-mer completeness, reusing the meryl DB from KMER_ANALYSIS.
    //
    // Matched to the assembly's OWN sample. This used to be ch_meryl_db.first(), which takes
    // a single emission and broadcasts it: correct for a one-sample project, but on the
    // seven-isolate Pst run it scored all 35 assemblies against whichever sample's k-mers
    // happened to finish counting first. Being a race, it was not even stable between runs
    // -- gkPucStri1's hicanu QV read 70.34 in one run and 48.04 in the next, same assembly.
    //
    // combine(by: 0) rather than join(): many assemblies share one sample, and join()
    // consumes a key after its first match, so six of seven samples would silently vanish.
    // The '*' key carries a single project-wide DB supplied by --meryl_db, which is how
    // --qc_only runs: there is no per-sample entry to match, so every assembly falls back
    // to it, reproducing the old broadcast where that is genuinely what is wanted.
    ch_meryl_map = ch_meryl_db
        .map { meta_k, db -> [ (meta_k.sample ?: '*'): db ] }
        .reduce([:]) { acc, entry -> acc + entry }

    ch_asm_meryl = ch_qc
        .combine(ch_meryl_map)
        .map { meta, fasta, gfa, dbmap -> tuple(meta, fasta, gfa, dbmap[meta.sample] ?: dbmap['*']) }
        .filter { meta, fasta, gfa, db -> db != null }

    MERQURY(ch_asm_meryl)
    ch_versions = ch_versions.mix(MERQURY.out.versions.first().ifEmpty(null))

    // Read-to-assembly mapping — the basis for interpreting the dikaryotic structure.
    //
    // Same fix as above, and the same bug: ch_reads.first() mapped every assembly against
    // gkPucStri1's reads, so every depth, coverage mode and coverage interpretation in the
    // run described one sample's reads against another sample's contigs.
    ch_asm_reads = ch_qc
        .map { meta, fasta, gfa -> tuple(meta.sample, meta, fasta, gfa) }
        .combine(ch_reads.map { meta_r, fq -> tuple(meta_r.sample, fq) }, by: 0)
        .map { sample, meta, fasta, gfa, fq -> tuple(meta, fasta, gfa, fq) }

    MINIMAP2_ASSEMBLY(ch_asm_reads)
    ch_versions = ch_versions.mix(MINIMAP2_ASSEMBLY.out.versions.first().ifEmpty(null))

    // Attach each assembly's OWN sample yield, so the depth expectation is that sample's
    // rather than a project-wide average. combine() on a value channel broadcasts the map to
    // every assembly; the lookup falls back to the global param when the map is empty.
    ch_coverage_in = MINIMAP2_ASSEMBLY.out.depth
        .join(MINIMAP2_ASSEMBLY.out.coverage)
        .combine(ch_yield_map)
        .map { meta, depth, coverage, ymap ->
            tuple(meta, depth, coverage, ymap[meta.sample] ?: params.read_yield_bases)
        }

    COVERAGE_SUMMARY(ch_coverage_in)
    ch_versions = ch_versions.mix(COVERAGE_SUMMARY.out.versions.first().ifEmpty(null))

    // The same expectation, as a table the histogram plot can draw one line per panel from.
    ch_expected_depths = ch_coverage_in
        .map { meta, depth, coverage, yield_bases ->
            def exp = GenomeSize.expectedHaplotypeDepth(yield_bases, params.genome_size)
            exp ? "${meta.id}\t${exp}" : null
        }
        .filter { it != null }

    // Competitive haplotype mapping: hap1 and hap2 in ONE index, so reads must choose.
    //
    // Pairs are formed per (sample, assembler) from the assembly_type suffix, so this picks
    // up any future phased pair -- Verkko's, or Hi-C haplotypes -- without editing. Samples
    // with only one haplotype assembled produce no pair and are simply skipped.
    ch_hap_pairs = ch_qc
        .map { meta, fasta, gfa ->
            def m = (meta.assembly_type =~ /^(.*)hap([12])\b/)
            m ? tuple(tuple(meta.sample, meta.assembler, meta.readset, m[0][1]),
                      m[0][2] as Integer, meta, fasta)
              : null
        }
        .filter { it != null }
        .groupTuple()
        .filter { key, haps, metas, fastas -> haps.size() == 2 }
        .map { key, haps, metas, fastas ->
            // groupTuple preserves arrival order, not hap order, so sort explicitly rather
            // than trusting hap1 to land first.
            def byHap = [haps, fastas].transpose().sort { it[0] }.collect { it[1] }
            def (sample, assembler, readset, prefix) = key
            def meta = [ id       : "${sample}_${readset}_${assembler}_${prefix}hap_pair",
                         sample   : sample,
                         assembler: assembler,
                         readset  : readset ]
            tuple(sample, meta, byHap[0], byHap[1])
        }
        .combine(ch_reads.map { meta_r, fq -> tuple(meta_r.sample, fq) }, by: 0)
        .map { sample, meta, h1, h2, fq -> tuple(meta, h1, h2, fq) }

    HAPLOTYPE_COMPETITION(ch_hap_pairs)
    ch_versions = ch_versions.mix(HAPLOTYPE_COMPETITION.out.versions.first().ifEmpty(null))

    // Telomere placement: capped contig ends, and any telomere sitting in a contig INTERIOR
    // (a mis-join). Cheap enough to run on every candidate, and it is the only contiguity
    // measure here that can be scored against an expectation — a genome has a known number
    // of chromosome ends, whereas N50 can only be compared between assemblies.
    TELOMERES(ch_qc)
    ch_versions = ch_versions.mix(TELOMERES.out.versions.first().ifEmpty(null))

    //
    // Gather per-assembly evidence into one record for the summary table. Joins are on meta,
    // and the optional tools use ifEmpty so a disabled or failed step degrades the record
    // rather than dropping the assembly from the comparison entirely.
    //
    // The two cleaning JSONs ride along so the record can say what was measured (raw,
    // FCS-GX-cleaned, or also mito-free) and how much each step removed.
    ch_for_summary = ch_qc
        .map { meta, fasta, gfa -> tuple(meta, fasta) }
        .join(GFASTATS.out.stats,           remainder: true)
        .join(BUSCO.out.json,               remainder: true)
        .join(MERQURY.out.qv,               remainder: true)
        .join(COVERAGE_SUMMARY.out.json,    remainder: true)
        .join(TELOMERES.out.json,           remainder: true)
        .join(FCS_GX_CLEAN.out.json,        remainder: true)
        .join(ch_mito_json,                 remainder: true)

    emit:
    for_summary  = ch_for_summary
    cleaned      = ch_cleaned                    // [ meta, cleaned fasta ] — NuclearPhaser input
    mito_manifest = ch_mito_tsv                 // empty unless --run_mito_screen
    clean_json   = FCS_GX_CLEAN.out.json
    busco_tables = BUSCO.out.full_table          // per-gene table — NuclearPhaser input
    gfastats     = GFASTATS.out.stats            // [ meta, gfastats.txt ] — NuclearPhaser size gate
    for_multiqc  = BUSCO.out.summary.map { meta, f -> f }
    depth_hists  = COVERAGE_SUMMARY.out.histogram.map { meta, f -> f }
    expected_depths = ch_expected_depths         // "assembly_id\texpected depth" rows
    hap_competition = HAPLOTYPE_COMPETITION.out.json
    versions     = ch_versions
}
