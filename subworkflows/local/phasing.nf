//
// PHASING — Hi-C phasing of assemblies that were not phased by their assembler.
//
// The pipeline already supports Hi-C in one place: hifiasm's --h1/--h2, which turns its
// partially phased bp.hap1/bp.hap2 into fully phased hic.hap1/hic.hap2. That only helps
// hifiasm. Flye, HiCanu and Verkko cannot use Hi-C at all, so a Hi-C library bought for this
// project would improve exactly one of four candidates.
//
// NuclearPhaser closes that gap. It phases an assembly AFTER the fact, from Hi-C contacts
// plus gene and BUSCO synteny, so the same library applies to every candidate — and phasing a
// hifiasm primary with it gives an independent check on hifiasm's own phasing, which on a
// dikaryon is worth having.
//
// Entirely inert without --hic_r1/--hic_r2, --nuclearphaser_genes and --run_nuclearphaser.
//

include { PBLAT_GENES     } from '../../modules/local/pblat_genes'
include { HIC_ALIGN       } from '../../modules/local/hic_align'
include { HIC_PAIRS       } from '../../modules/local/hic_pairs'
include { HIC_MATRIX      } from '../../modules/local/hic_matrix'
include { NUCLEARPHASER   } from '../../modules/local/nuclearphaser'

workflow PHASING {

    take:
    ch_candidates    // [ meta, fasta ] — cleaned assemblies from ASSEMBLY_QC
    ch_busco_tables  // [ meta, busco full_table.tsv ]
    ch_gfastats      // [ meta, gfastats.txt ] of the same cleaned assemblies

    main:
    ch_versions = Channel.empty()

    // Phasing every candidate is wasteful: haplotype outputs are already phased, and unitig
    // graphs are not assemblies. The default targets primaries only, which is where the
    // question "can Hi-C separate the two nuclei here?" actually lives.
    def target_re = java.util.regex.Pattern.compile(params.nuclearphaser_targets)
    def by_name = ch_candidates.filter { meta, fasta -> meta.id ==~ target_re }

    // NuclearPhaser splits an assembly into two nuclei, so it needs one that HOLDS both. A
    // collapsed primary has one copy of each chromosome: nearly everything lands in one bin,
    // the other has no trans contacts, and NuclearPhaser dies on a ZeroDivisionError. On
    // M. larici-populina that was hifiasm primary (103.3 / 2.7 Mb), IPA primary and Flye --
    // three wasted Hi-C alignments and three red tasks. The threshold is size_flag's own
    // `collapsed` rule (below 0.75 x the dikaryon minimum), so the gate and the report agree.
    // An assembly with no gfastats is kept: a missing measurement is not evidence of collapse.
    def min_size = (GenomeSize.dikaryonMin(params.genome_size, params.dikaryon_size_frac_min) * 0.75) as long
    ch_targets = by_name
        .join(ch_gfastats, remainder: true)
        .filter { meta, fasta, stats -> fasta != null }
        .filter { meta, fasta, stats ->
            if (!params.nuclearphaser_skip_collapsed || stats == null) return true
            def line = stats.readLines().find { it.toLowerCase().startsWith('total scaffold length:') }
            if (line == null) return true
            def size = line.tokenize(':')[1].trim() as long
            if (size >= min_size) return true
            log.warn "NuclearPhaser: skipping ${meta.id} -- ${size} bp is below ${min_size} bp, " +
                     "so it holds one nucleus and there is nothing to phase " +
                     "(--nuclearphaser_skip_collapsed false to force)"
            return false
        }
        .map { meta, fasta, stats -> tuple(meta, fasta) }

    ch_hic = Channel.value(
        tuple([id: 'hic'], file(params.hic_r1 ?: "${projectDir}/assets/NO_FILE_hic_r1"),
                           file(params.hic_r2 ?: "${projectDir}/assets/NO_FILE_hic_r2"))
    )

    ch_genes = params.nuclearphaser_genes
        ? file(params.nuclearphaser_genes, checkIfExists: true)
        : file("${projectDir}/assets/NO_FILE_genes")

    // Gene synteny signal. PSL, from pblat rather than the unmaintained BioKanga — see the
    // module for why that substitution is safe.
    PBLAT_GENES(ch_targets, ch_genes)
    ch_versions = ch_versions.mix(PBLAT_GENES.out.versions.first().ifEmpty(null))

    // Hi-C contacts, built against each target assembly. Per-assembly, not once: contact
    // coordinates are meaningless against a different set of contigs.
    //
    // Three processes rather than one, because they have very different costs and very
    // different reasons to change. Re-binning the matrix at a new resolution, or moving the
    // MAPQ threshold, should not re-align the library.
    HIC_ALIGN(ch_targets, ch_hic)
    ch_versions = ch_versions.mix(HIC_ALIGN.out.versions.first().ifEmpty(null))

    HIC_PAIRS(HIC_ALIGN.out.sam.join(HIC_ALIGN.out.chrom_sizes))
    ch_versions = ch_versions.mix(HIC_PAIRS.out.versions.first().ifEmpty(null))

    HIC_MATRIX(HIC_PAIRS.out.pairs.join(HIC_ALIGN.out.chrom_sizes))
    ch_versions = ch_versions.mix(HIC_MATRIX.out.versions.first().ifEmpty(null))

    // NuclearPhaser needs the per-gene BUSCO table, which ASSEMBLY_QC already computes for
    // every candidate — no second BUSCO run. Joining on meta also means a candidate whose
    // BUSCO failed simply drops out of phasing rather than failing the task on a missing file.
    ch_np_in = ch_targets
        .join(PBLAT_GENES.out.psl)
        .join(ch_busco_tables)
        .join(HIC_MATRIX.out.contacts)

    NUCLEARPHASER(ch_np_in)
    ch_versions = ch_versions.mix(NUCLEARPHASER.out.versions.first().ifEmpty(null))

    emit:
    haplotypes     = NUCLEARPHASER.out.haplotypes
    phase_switches = NUCLEARPHASER.out.phase_switches
    contacts       = HIC_MATRIX.out.contacts
    hic_stats      = HIC_PAIRS.out.stats
    versions       = ch_versions
}
