process YAHS {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/scaffolding/yahs/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(contigs), path(fai), path(bam)

    output:
    tuple val(meta), path("${meta.id}.yahs.scaffolds.fa"), emit: scaffolds
    tuple val(meta), path("${meta.id}.yahs.agp")         , emit: agp
    path "${meta.id}.yahs.log"                           , emit: log
    path "versions.yml"                                  , emit: versions

    script:
    // Hi-C scaffolding: order and orient contigs into chromosomes using the contact map. The
    // sequence does not change -- YaHS joins what is already assembled, padding each join with
    // a run of Ns -- so this cannot rescue a bad assembly and cannot inflate a good one.
    //
    // WHAT THIS DOES AND DOES NOT REPRODUCE. Duplessis et al. (2026) reached 18 scaffolds for
    // M. larici-populina by REFERENCE-GUIDED scaffolding: RagTag v2.1.0 against the v2
    // assembly, which was itself anchored to a genetic map. That is a different operation from
    // this one and carries the reference's assumptions with it. De novo Hi-C scaffolding is
    // what they used for the other species in that paper (HapHiC on M. allii-populina), and it
    // is what this pipeline offers, because a genetic-map-anchored reference is not something
    // most projects have.
    //
    // -q is the one lever that matters. A Hi-C read pair mapping ambiguously across a repeat
    // is worse than no pair at all: it invents a contact between contigs that are not
    // neighbours, and a scaffolder cannot tell that from a real long-range contact. Shares
    // params.hic_min_mapq with the phasing branch so one setting governs what counts as a
    // contact everywhere in the pipeline.
    //
    // No -e: the Hi-C libraries here are not from a single known restriction enzyme, and
    // guessing one wrong biases the assembly-error correction step. Without it YaHS uses its
    // own break-point detection, which is the documented default.
    def args = task.ext.args ?: ''
    """
    yahs \\
        -o ${meta.id}.yahs \\
        -q ${params.hic_min_mapq} \\
        ${args} \\
        ${contigs} \\
        ${bam} \\
        2>&1 | tee ${meta.id}.yahs.log

    # tee returns 0 whatever yahs did, and Nextflow runs this under `bash -ue` without
    # pipefail, so the pipeline's status is tee's. Check for the output instead -- the same
    # trap that let a crashed NuclearPhaser report success in job 34173414.
    if [ ! -s ${meta.id}.yahs_scaffolds_final.fa ]; then
        echo "YaHS produced no scaffolds for ${meta.id}; see ${meta.id}.yahs.log" >&2
        exit 1
    fi

    mv ${meta.id}.yahs_scaffolds_final.fa  ${meta.id}.yahs.scaffolds.fa
    mv ${meta.id}.yahs_scaffolds_final.agp ${meta.id}.yahs.agp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yahs: \$(yahs --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.yahs.scaffolds.fa ${meta.id}.yahs.agp ${meta.id}.yahs.log
    echo '"${task.process}": {yahs: stub}' > versions.yml
    """
}
