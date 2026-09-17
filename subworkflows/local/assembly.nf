//
// ASSEMBLY — produce several candidates, deliberately.
//
// The point is not to pick a winner by contiguity. It is to ask whether independent
// assemblers agree on the structure of the genome, and to keep enough evidence (graphs,
// haplotype labels, provenance) to judge which result is biologically credible.
//
// Each assembler is its own process so that tweaking a QC plot never re-triggers a
// multi-day assembly.
//

include { HIFIASM   } from '../../modules/local/hifiasm'
include { FLYE      } from '../../modules/local/flye'
include { HICANU    } from '../../modules/local/hicanu'
include { VERKKO    } from '../../modules/local/verkko'
include { IPA       } from '../../modules/local/ipa'
include { GFA2FASTA } from '../../modules/local/gfa2fasta'

//
// Classify a hifiasm GFA by filename into the assembly_type vocabulary.
//
// This distinction is the whole point: hifiasm's HiFi-only bp.hap1/bp.hap2 are PARTIALLY
// phased, while Hi-C mode's hic.hap1/hic.hap2 are FULLY phased. Calling both "hap1" would
// throw away exactly the information that matters.
//
def classifyHifiasmGfa(String name) {
    if (name.contains('.hic.hap1')) return 'fully_phased_hap1'
    if (name.contains('.hic.hap2')) return 'fully_phased_hap2'
    if (name.contains('.hap1'))     return 'partially_phased_hap1'
    if (name.contains('.hap2'))     return 'partially_phased_hap2'
    if (name.contains('.a_ctg'))    return 'alternate'
    if (name.contains('.p_ctg'))    return 'primary'
    if (name.contains('.r_utg'))    return 'raw_unitig'
    if (name.contains('.p_utg'))    return 'processed_unitig'
    return 'unknown'
}

// Unitig graphs are retained as evidence but are not assembly candidates for QC.
def CONTIG_TYPES = [
    'primary', 'alternate',
    'partially_phased_hap1', 'partially_phased_hap2',
    'fully_phased_hap1', 'fully_phased_hap2'
]

workflow ASSEMBLY {

    take:
    ch_reads_combined   // [ meta(id,sample,readset), [ fastqs ] ]
    ch_reads_per_run    // [ meta(id,sample,run), fastq ]
    ch_extra_readsets   // [ meta(id,sample,readset), [ fastqs ] ] — derived datasets, may be empty

    main:
    ch_versions   = Channel.empty()
    ch_assemblies = Channel.empty()

    //
    // Build the set of read subsets to assemble.
    //
    // DECISION D9: the diagnostic ALL/run1/run2 comparison is wired but off by default,
    // because it triples the most expensive process. It answers "is each sequencing run
    // actually earning its place?".
    //
    ch_subsets = params.subset_assemblies
        ? ch_reads_combined.mix(
              ch_reads_per_run.map { meta, fq ->
                  tuple(
                      [ id: "${meta.sample}_${meta.run}".toString(),
                        sample: meta.sample,
                        readset: meta.run ],
                      [ fq ]
                  )
              }
          )
        : ch_reads_combined

    // Derived datasets — host-removed reads, and the quality/length-filtered dataset — join
    // the same channel. They are questions of the same shape as the run subsets ("does this
    // version of the data assemble better?"), so they get the same treatment: assembled by
    // hifiasm only, and carried through QC as ordinary candidates distinguished by readset.
    // The channel is empty unless the relevant flag is set, so this costs nothing by default.
    ch_subsets = ch_subsets.mix(ch_extra_readsets)

    //
    // Hi-C channel. Inert until --hic_r1/--hic_r2 are supplied; supplying them switches
    // hifiasm from partially-phased to fully-phased output with no other change.
    //
    ch_hic = params.hic_r1
        ? Channel.value(tuple([id: 'hic'], file(params.hic_r1), file(params.hic_r2)))
        : Channel.value(tuple([id: 'no_hic'], [], []))

    //
    // hifiasm — the central assembler.
    //
    HIFIASM(ch_subsets, ch_hic)
    ch_versions = ch_versions.mix(HIFIASM.out.versions.ifEmpty(null))

    // Fan the GFAs out into individually-typed assembly candidates.
    //
    // hifiasm writes a `.noseq.gfa` alongside each real graph: same topology, but with '*'
    // in place of every sequence. Those must be dropped BEFORE classification — their names
    // contain '.hap1'/'.p_ctg' too, so they would classify as duplicate candidates sharing a
    // meta.id with the real assembly, and a FASTA of '*' records is non-empty enough to slip
    // past a naive emptiness check.
    ch_hifiasm_gfa = HIFIASM.out.gfa
        .transpose()
        .filter { meta, gfa -> !gfa.name.contains('.noseq.') }
        .map { meta, gfa ->
            def type = classifyHifiasmGfa(gfa.name)
            tuple(
                [ id        : "${meta.sample}_${meta.readset}_hifiasm_${type}".toString(),
                  sample    : meta.sample,
                  readset   : meta.readset,
                  assembler : 'hifiasm',
                  assembly_type: type ],
                gfa
            )
        }
        .filter { meta, gfa -> meta.assembly_type in CONTIG_TYPES }

    GFA2FASTA(ch_hifiasm_gfa)
    ch_versions = ch_versions.mix(GFA2FASTA.out.versions.ifEmpty(null))

    ch_assemblies = ch_assemblies.mix(
        GFA2FASTA.out.fasta.join(ch_hifiasm_gfa).map { meta, fasta, gfa ->
            tuple(meta, fasta, gfa)
        }
    )

    //
    // Independent assemblers — a different approach to the same data.
    //
    // These take the COMBINED reads only, never the subsets or the derived readsets. The
    // ALL/run1/run2, filtered and host-filtered comparisons are questions about what a given
    // version of the DATA contributes, which hifiasm alone answers; putting them through Flye
    // and HiCanu as well would multiply the scarcest resource (900 GB big-memory node time)
    // for no extra information.
    //
    FLYE(ch_reads_combined)
    ch_versions = ch_versions.mix(FLYE.out.versions.ifEmpty(null))

    ch_assemblies = ch_assemblies.mix(
        FLYE.out.fasta.map { meta, fasta ->
            tuple(
                [ id        : "${meta.sample}_${meta.readset}_flye_primary".toString(),
                  sample    : meta.sample,
                  readset   : meta.readset,
                  assembler : 'flye',
                  assembly_type: 'primary' ],
                fasta, []
            )
        }
    )

    HICANU(ch_reads_combined)
    ch_versions = ch_versions.mix(HICANU.out.versions.ifEmpty(null))

    ch_assemblies = ch_assemblies.mix(
        HICANU.out.fasta.map { meta, fasta ->
            tuple(
                [ id        : "${meta.sample}_${meta.readset}_hicanu_primary".toString(),
                  sample    : meta.sample,
                  readset   : meta.readset,
                  assembler : 'hicanu',
                  assembly_type: 'primary' ],
                fasta, []
            )
        }
    )

    VERKKO(ch_reads_combined)
    ch_versions = ch_versions.mix(VERKKO.out.versions.ifEmpty(null))

    ch_assemblies = ch_assemblies.mix(
        VERKKO.out.fasta.map { meta, fasta ->
            tuple(
                [ id        : "${meta.sample}_${meta.readset}_verkko_primary".toString(),
                  sample    : meta.sample,
                  readset   : meta.readset,
                  assembler : 'verkko',
                  assembly_type: 'primary' ],
                fasta, []
            )
        }
    )

    // IPA emits a primary and an alternate contig set, so unlike Flye/HiCanu/Verkko it
    // contributes two candidates. Both go through QC: the alternate is where IPA puts the
    // sequence it judged to be the second haplotype, which is exactly the question a
    // dikaryotic rust assembly is asking.
    IPA(ch_reads_combined)
    ch_versions = ch_versions.mix(IPA.out.versions.ifEmpty(null))

    ch_assemblies = ch_assemblies.mix(
        IPA.out.fasta.map { meta, fasta ->
            tuple(
                [ id        : "${meta.sample}_${meta.readset}_ipa_primary".toString(),
                  sample    : meta.sample,
                  readset   : meta.readset,
                  assembler : 'ipa',
                  assembly_type: 'primary' ],
                fasta, []
            )
        }
    )

    ch_assemblies = ch_assemblies.mix(
        IPA.out.alternate.map { meta, fasta ->
            tuple(
                [ id        : "${meta.sample}_${meta.readset}_ipa_alternate".toString(),
                  sample    : meta.sample,
                  readset   : meta.readset,
                  assembler : 'ipa',
                  assembly_type: 'alternate' ],
                fasta, []
            )
        }
    )

    emit:
    assemblies = ch_assemblies   // [ meta(id,sample,assembler,assembly_type,readset), fasta, gfa ]
    versions   = ch_versions
}
