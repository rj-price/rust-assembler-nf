//
// SCAFFOLDING — Hi-C scaffolding of PHASED haplotypes.
//
// Ordering and orienting contigs into chromosomes using the Hi-C contact map. No new sequence
// is created: a scaffolder joins what is already assembled and pads each join with Ns. That
// makes this the last step that can change contiguity without changing the evidence, and the
// first step whose output must never be confused with the contigs it came from — so scaffolds
// are published under scaffolding/ and are NOT fed back into the candidate comparison.
//
// WHY HAPLOTYPES ONLY. Scaffolding a collapsed dikaryotic primary into 18 chromosomes is not a
// well-posed question: the contigs of both nuclei are present, most Hi-C contacts are within a
// nucleus, and the scaffolder has no way to know which of two homologous contigs a contact
// belongs to. Phased haplotypes are one nucleus each, which is the case every scaffolder is
// built for. Both sources of haplotypes qualify: hifiasm's --h1/--h2 output and NuclearPhaser's
// Haplotype_0/Haplotype_1.
//
// WHAT THE PAPER DID, since this branch exists because of it. Duplessis et al. (2026) reached
// 18 scaffolds for M. larici-populina 98AG31 by REFERENCE-GUIDED scaffolding — RagTag against
// the v2 assembly, itself anchored to a genetic map from a selfing progeny. That is a
// different operation from de novo Hi-C scaffolding and is not reproduced here; see
// modules/local/yahs.nf. De novo Hi-C scaffolding is what they used for the other species in
// the same paper.
//
// Entirely inert without --hic_r1/--hic_r2 and --run_scaffolding.
//

include { HIC_ALIGN } from '../../modules/local/hic_align'
include { HIC_BAM   } from '../../modules/local/hic_bam'
include { YAHS      } from '../../modules/local/yahs'
include { HAPHIC    } from '../../modules/local/haphic'
include { GFASTATS as GFASTATS_SCAFFOLDS } from '../../modules/local/gfastats'
include { SCAFFOLD_RECORD  } from '../../modules/local/scaffold_summary'
include { SCAFFOLD_SUMMARY } from '../../modules/local/scaffold_summary'

workflow SCAFFOLDING {

    take:
    ch_candidates     // [ meta, fasta ] — cleaned assemblies from ASSEMBLY_QC
    ch_np_haplotypes  // [ meta, [ Haplotype_*.fasta ] ] — from PHASING, may be empty

    main:
    ch_versions = Channel.empty()

    def target_re = java.util.regex.Pattern.compile(params.scaffold_targets)

    // Source one: haplotypes an assembler produced directly, which reach here as ordinary
    // candidates and are selected by id.
    ch_from_assembler = ch_candidates.filter { meta, fasta -> meta.id ==~ target_re }

    // Source two: NuclearPhaser's haplotypes, which arrive as a LIST of FASTAs under the id of
    // the assembly they were phased from. Split into one entry per haplotype and give each an
    // id of its own, so downstream publishing and the scaffolds table can tell them apart.
    //
    // Empty files are dropped rather than passed on. NuclearPhaser exits 0 and writes empty
    // haplotype FASTAs when it fails — it did exactly that for the Verkko assembly in job
    // 34173414, crashing on a KeyError after every other candidate had phased — and an empty
    // FASTA would fail the scaffolder in a way that looks like a scaffolding problem.
    ch_from_nuclearphaser = ch_np_haplotypes
        .flatMap { meta, files ->
            (files instanceof List ? files : [files])
                .findAll { f -> f.name ==~ /Haplotype_\d+_final\.fasta/ && f.size() > 0 }
                .collect { f ->
                    def hap = (f.name =~ /Haplotype_(\d+)_final\.fasta/)[0][1]
                    tuple(meta + [id: "${meta.id}_np_hap${hap}"], f)
                }
        }

    ch_targets = ch_from_assembler.mix(ch_from_nuclearphaser)

    ch_hic = Channel.value(
        tuple([id: 'hic'], file(params.hic_r1 ?: "${projectDir}/assets/NO_FILE_hic_r1"),
                           file(params.hic_r2 ?: "${projectDir}/assets/NO_FILE_hic_r2"))
    )

    // Hi-C aligned against each haplotype in turn. PHASING aligns against primaries, and those
    // alignments are useless here: contact coordinates only mean something against the exact
    // contig set they were produced from, and a haplotype is a different contig set.
    HIC_ALIGN(ch_targets, ch_hic)
    ch_versions = ch_versions.mix(HIC_ALIGN.out.versions.first().ifEmpty(null))

    HIC_BAM(ch_targets.join(HIC_ALIGN.out.sam))
    ch_versions = ch_versions.mix(HIC_BAM.out.versions.first().ifEmpty(null))

    ch_scaffold_in = HIC_BAM.out.contigs.join(HIC_BAM.out.bam)

    // One scaffolder per run, chosen by --scaffolder. Both are wired so the two can be
    // compared on the same input by running twice with -resume; everything up to and including
    // HIC_BAM is then cached and only the scaffolding itself repeats.
    if (params.scaffolder == 'haphic') {
        HAPHIC(ch_scaffold_in)
        ch_scaffolds = HAPHIC.out.scaffolds
        ch_agp       = HAPHIC.out.agp
        ch_versions  = ch_versions.mix(HAPHIC.out.versions.first().ifEmpty(null))
    }
    else {
        YAHS(ch_scaffold_in)
        ch_scaffolds = YAHS.out.scaffolds
        ch_agp       = YAHS.out.agp
        ch_versions  = ch_versions.mix(YAHS.out.versions.first().ifEmpty(null))
    }

    // Contiguity of the scaffolds, which is the number this branch exists to produce — a
    // scaffold N50 and a scaffold count are what published assemblies report, and until now
    // this pipeline could only report contigs. GFASTATS takes [meta, fasta, gfa]; there is no
    // graph for a scaffolded assembly, hence the placeholder.
    //
    // The id is suffixed so these rows can never be mistaken for the contig-level stats of the
    // same assembly, which GFASTATS writes to the same directory -- they land side by side as
    // <id>.gfastats.txt and <id>_scaffolds.gfastats.txt, which is the comparison you want.
    ch_stats_in = ch_scaffolds.map { meta, fasta ->
        tuple(meta + [id: "${meta.id}_scaffolds"], fasta, file("${projectDir}/assets/NO_FILE_gfastats"))
    }

    GFASTATS_SCAFFOLDS(ch_stats_in)
    ch_versions = ch_versions.mix(GFASTATS_SCAFFOLDS.out.versions.first().ifEmpty(null))

    // The scaffolding result as a table. gfastats alone cannot state it: it reports every
    // scaffold, including the small unplaced ones, so '# scaffolds' is 55 where the answer is
    // 18. SCAFFOLD_RECORD measures the length distribution from the FASTA to separate the two.
    //
    // Joined on the SUFFIXED meta, not the original: ch_stats_in is where the '_scaffolds' id
    // is minted, and GFASTATS_SCAFFOLDS emits that same meta, so the two line up. Joining
    // ch_scaffolds instead would silently never match, because its meta still carries the
    // unsuffixed id.
    SCAFFOLD_RECORD(
        ch_stats_in.map { meta, fasta, gfa -> tuple(meta, fasta) }
                   .join(GFASTATS_SCAFFOLDS.out.stats)
    )
    ch_versions = ch_versions.mix(SCAFFOLD_RECORD.out.versions.first().ifEmpty(null))

    SCAFFOLD_SUMMARY(SCAFFOLD_RECORD.out.record.collect())
    ch_versions = ch_versions.mix(SCAFFOLD_SUMMARY.out.versions.first().ifEmpty(null))

    emit:
    scaffolds = ch_scaffolds
    agp       = ch_agp
    stats     = GFASTATS_SCAFFOLDS.out.stats
    summary   = SCAFFOLD_SUMMARY.out.tsv
    mqc       = SCAFFOLD_SUMMARY.out.mqc
    versions  = ch_versions
}
