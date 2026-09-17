process HIFIASM {
    tag "${meta.id}"
    label 'process_assembly'

    // The assembly GRAPHS are first-class outputs, not throwaway intermediates — they carry
    // the haplotype structure that a flat FASTA discards.
    publishDir path: { "${params.outdir}/assembly/hifiasm/${meta.id}" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)
    tuple val(meta_hic), path(hic_r1), path(hic_r2)

    output:
    tuple val(meta), path("${meta.id}*.gfa")     , emit: gfa
    tuple val(meta), path("${meta.id}*.bin")     , emit: bin,  optional: true
    tuple val(meta), path("${meta.id}.hifiasm.log"), emit: log
    path "versions.yml"                          , emit: versions

    when:
    params.run_hifiasm

    script:
    // Hi-C mode changes the output vocabulary: HiFi-only gives PARTIALLY phased bp.hap1/hap2,
    // Hi-C gives FULLY phased hic.hap1/hap2. The subworkflow reads that distinction off the
    // filenames, so the two are never conflated in the report.
    def hic_args = hic_r1 ? "--h1 ${hic_r1} --h2 ${hic_r2}" : ''
    """
    hifiasm \\
        -o ${meta.id} \\
        -t ${task.cpus} \\
        -l ${params.hifiasm_purge_level} \\
        ${hic_args} \\
        ${params.hifiasm_extra} \\
        ${reads} \\
        2> ${meta.id}.hifiasm.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hifiasm: \$(hifiasm --version 2>&1 | head -1)
        hifiasm_purge_level: ${params.hifiasm_purge_level}
        hifiasm_mode: ${hic_r1 ? 'hi-c (fully phased)' : 'hifi-only (partially phased)'}
    END_VERSIONS
    """

    stub:
    // Mirrors hifiasm's real output set INCLUDING the .noseq.gfa companions, so that a stub
    // run catches any regression in the noseq filtering (they classify identically to the
    // real graphs and would otherwise create duplicate candidates).
    //
    // The INFIX is load-bearing. hifiasm names its graphs <id>.bp.* in HiFi-only mode and
    // <id>.hic.* in Hi-C mode, and subworkflows/local/assembly.nf reads the fully-phased vs
    // partially-phased distinction straight off that infix. This block used to write .bp.*
    // unconditionally, so a stub run could never produce a fully_phased_hap* candidate --
    // which silently made the entire SCAFFOLDING branch untestable. Its default
    // --scaffold_targets matches only fully_phased/np haplotypes, so ch_targets came out
    // empty, no scaffolder ran, and the run still reported success: an empty .collect()
    // emits nothing, so a branch that never executed looks exactly like one that passed.
    //
    // Filenames verified against the real published runs, not from memory: mlp98AG31 (Hi-C)
    // wrote .hic.{p_ctg,hap1.p_ctg,hap2.p_ctg,r_utg,p_utg}.gfa and rust (HiFi-only) wrote the
    // same five under .bp., each with a .noseq.gfa companion and no a_ctg in either mode.
    //
    // hic_r1 is [] when --hic_r1 is unset (assembly.nf builds the no_hic tuple from empty
    // lists, not NO_FILE placeholders), so this truth test matches the script block's.
    def infix = hic_r1 ? 'hic' : 'bp'
    """
    for t in p_ctg hap1.p_ctg hap2.p_ctg r_utg p_utg; do
        printf 'S\\tstub1\\tACGT\\tLN:i:4\\n' > ${meta.id}.${infix}.\$t.gfa
        printf 'S\\tstub1\\t*\\tLN:i:4\\n'    > ${meta.id}.${infix}.\$t.noseq.gfa
    done
    touch ${meta.id}.hifiasm.log
    echo '"${task.process}": {hifiasm: stub}' > versions.yml
    """
}
