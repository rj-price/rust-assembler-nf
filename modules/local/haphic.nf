process HAPHIC {
    tag "${meta.id}"
    label 'process_medium'

    publishDir path: { "${params.outdir}/scaffolding/haphic/${meta.id}" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(contigs), path(fai), path(bam)

    output:
    tuple val(meta), path("${meta.id}.haphic.scaffolds.fa"), emit: scaffolds
    tuple val(meta), path("${meta.id}.haphic.agp")         , emit: agp
    path "${meta.id}.haphic.log"                           , emit: log
    path "versions.yml"                                    , emit: versions

    script:
    // The alternative to YaHS, and the one with the better provenance for this genome:
    // Duplessis et al. (2026) used HapHiC to scaffold M. allii-populina 12AY07, with the
    // expected chromosome count set to 36 for the combined diploid. It was built for
    // haplotype-phased and polyploid assemblies, where YaHS's assumption that every contig
    // belongs to one linear genome is wrong.
    //
    // THE CHROMOSOME COUNT IS NOT A HINT. HapHiC clusters contigs into exactly this many
    // groups, so passing the wrong number does not degrade the result gracefully -- it
    // produces confidently wrong chromosomes. Per HAPLOTYPE the number is 18 for
    // M. larici-populina, not 36; 36 is the dikaryon, and this pipeline scaffolds the
    // haplotypes separately. main.nf refuses to start without --scaffold_n_chromosomes rather
    // than guessing.
    //
    // No biocontainer exists (checked against the Galaxy depot). containers/haphic.def builds
    // one; --haphic_sif points at it.
    //
    // TESTED on mlp98AG31 hifiasm hap2 (43 contigs, 102.74 Mb, already near chromosome-level)
    // against the same Hi-C BAM YaHS used, at 18 chromosomes. It finishes, and the two agree
    // almost exactly -- 43 sequences and 102.74 Mb from both, N50 5.55 vs 5.54 Mb, identical
    // contigs in the top five bar one join. Independent corroboration of the 18 chromosomes,
    // and no reason to prefer HapHiC here.
    //
    // Expect this in the log on an input this contiguous:
    //     Parameter --nclusters (18) is greater than the number of clusters (11) after
    //     reassignment, try higher inflations
    // It is a WARNING, not an error, and on a near-complete assembly it is expected rather
    // than alarming: HapHiC's Markov clustering groups contigs into chromosomes, and when
    // most chromosomes are ALREADY single contigs there is nothing left to group. Raising
    // --max_inflation does not fix it -- tested at 8, which ran 70 clustering rounds and
    // still capped at the same number. HapHiC earns its keep on FRAGMENTED phased assemblies;
    // this pipeline's hifiasm Hi-C haplotypes are not that, which is why YaHS is the default.
    //
    // A related failure IS fatal and looks the same from a distance: if the alignments are
    // empty the clustering aborts with "Pipeline Aborted: Inflation recommendation failed".
    // Before blaming inflation, check how many records survived filter_bam below.
    def args = task.ext.args ?: ''
    """
    # HapHiC drives its own multi-stage pipeline and writes a cache under \$HOME. \$HOME is
    # read-only inside the container on the compute nodes -- the same failure VERKKO hit in
    # job 34035432 and IPA hit later. Point it at the task work directory.
    export HOME="\$PWD"
    export XDG_CACHE_HOME="\$PWD/.cache"
    export MPLCONFIGDIR="\$PWD/.mpl"
    mkdir -p "\$XDG_CACHE_HOME" "\$MPLCONFIGDIR"

    # HapHiC's own documentation asks for FILTERED alignments, and means it: its first
    # positional is described as "filtered Hi-C read alignments", and upstream's usage line
    # pipes every BAM through this utility before the pipeline sees it. HIC_BAM deliberately
    # applies no MAPQ filter -- it is shared with YAHS, which takes its own -q -- so the
    # filtering belongs here, on the copy this process uses.
    #
    #   MAPQ 1 keeps read pairs where BOTH ends are uniquely placed, which is the whole point
    #   of a Hi-C contact: a pair with one ambiguous end says nothing about which two contigs
    #   are adjacent.
    #
    #   --nm 3 drops pairs with a high edit distance. On a dikaryon that matters more than
    #   usual, because the commonest way to mis-map here is onto the OTHER haplotype's copy
    #   of the same locus, and those alignments carry the haplotype's differences as
    #   mismatches. Letting them through would create contacts between haplotypes.
    #
    # `haphic pipeline` accepts bam or pairs, not sam, so this goes back through samtools --
    # which is in containers/haphic.def for no other reason than this line.
    # --remove-singletons is NOT optional here, and the reason is upstream in this pipeline.
    # HIC_BAM drops unmapped, secondary and supplementary records with -F 0x904, which orphans
    # whichever mate survives; filter_bam then walks a name-sorted BAM expecting pairs, meets
    # an unpaired read and PANICS:
    #     thread 'main' panicked at src/main.rs:70:25:
    #     BAM may be coord-sorted or has singletons. Sort it by read name or try
    #     --remove-singletons
    # Measured on mlp98AG31 hap2: 124,727,986 records in, 102 out without this flag. The BAM
    # was name-sorted all along -- singletons, not sort order, were the problem.
    /opt/HapHiC/utils/filter_bam \\
        ${bam} \\
        1 \\
        --nm 3 \\
        --remove-singletons \\
        --threads ${task.cpus} \\
        | samtools view -b -@ ${task.cpus} -o ${meta.id}.filtered.bam -
    FILTER_STATUS=\${PIPESTATUS[0]}
    if [ "\$FILTER_STATUS" -ne 0 ]; then
        echo "filter_bam exited \$FILTER_STATUS for ${meta.id}" >&2
        exit "\$FILTER_STATUS"
    fi

    haphic pipeline \\
        ${contigs} \\
        ${meta.id}.filtered.bam \\
        ${params.scaffold_n_chromosomes} \\
        --threads ${task.cpus} \\
        --outdir haphic_out \\
        ${args} \\
        2>&1 | tee ${meta.id}.haphic.log

    # As in YAHS: `bash -ue` has no pipefail, so tee's exit status is the task's. Find the
    # build output by search rather than by a hard-coded stage directory -- HapHiC numbers its
    # stages, and the numbering shifts with the correction options.
    SCAF=\$(find haphic_out -name 'scaffolds.fa' | head -1)
    AGP=\$(find haphic_out -name 'scaffolds.agp' | head -1)
    if [ -z "\$SCAF" ] || [ ! -s "\$SCAF" ]; then
        echo "HapHiC produced no scaffolds for ${meta.id}; see ${meta.id}.haphic.log" >&2
        find haphic_out -type f >&2 || true
        exit 1
    fi

    cp "\$SCAF" ${meta.id}.haphic.scaffolds.fa
    cp "\$AGP"  ${meta.id}.haphic.agp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haphic: \$(haphic --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.haphic.scaffolds.fa ${meta.id}.haphic.agp ${meta.id}.haphic.log
    echo '"${task.process}": {haphic: stub}' > versions.yml
    """
}
