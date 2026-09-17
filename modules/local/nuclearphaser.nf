process NUCLEARPHASER {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/phasing/nuclearphaser/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(fasta), path(gene_mapping), path(busco_table), path(contacts)

    output:
    tuple val(meta), path("out/Haplotype_*.fasta")        , emit: haplotypes, optional: true
    tuple val(meta), path("out/*Unphased*.fasta")         , emit: unphased,   optional: true
    tuple val(meta), path("out/*phase*")                  , emit: phase_switches, optional: true
    tuple val(meta), path("${meta.id}.nuclearphaser.log") , emit: log
    path "versions.yml"                                   , emit: versions

    script:
    // Phases an EXISTING assembly into two nuclear haplotypes from Hi-C contacts plus gene
    // and BUSCO synteny. Unlike hifiasm's --h1/--h2, which phases inside hifiasm, this works
    // on any candidate — so one Hi-C library can be applied to the Flye, HiCanu and Verkko
    // assemblies as well, and to a hifiasm primary as an independent check on hifiasm's own
    // phasing. It was designed on dikaryotic rusts, which is exactly this case.
    //
    // TWO THINGS TO KNOW BEFORE READING THE OUTPUT:
    //
    // 1. This is NuclearPhaser's FIRST pass. The published method is two passes with a MANUAL
    //    phase-switch correction in between: pass one reports contigs it believes are phase
    //    switched, a human inspects and breaks them, and the inputs are regenerated. The
    //    pipeline deliberately stops after pass one rather than automating that judgement.
    //    The phase-switch files are published for exactly that inspection; re-enter with
    //    --qc_only and a corrected assembly to run the second pass.
    //
    // 2. It expects an assembly already cleaned of contaminant contigs — which is what
    //    FCS_GX_CLEAN produces, and why that is the input wired to it here.
    //
    // Non-zero exit is contained by the errorStrategy in conf/base.config: this is an
    // experimental branch and must not sink a run that has already spent days on assemblies.
    // Contained is not the same as hidden, though -- see the exit-status note below.
    """
    ASM=${fasta}
    if [[ "\$ASM" == *.gz ]]; then
        gunzip -c "\$ASM" > asm.fa
        ASM=asm.fa
    fi

    mkdir -p out

    NuclearPhaser.py \\
        -g ${gene_mapping} \\
        -b ${busco_table} \\
        -c ${contacts} \\
        -f "\$ASM" \\
        -o out/ \\
        2>&1 | tee ${meta.id}.nuclearphaser.log
    NP_STATUS=\${PIPESTATUS[0]}

    # PIPESTATUS, not \$?, and this is not a nicety. Nextflow runs script blocks under
    # `bash -ue` WITHOUT pipefail, so the exit status of `NuclearPhaser.py | tee` is tee's,
    # which is 0 whatever NuclearPhaser did. In job 34173414 NuclearPhaser died on a Python
    # traceback --
    #     File "/opt/NuclearPhaser/GeneBinning.py", line 51, in genes_shared
    #     KeyError: 'contig-0000617'
    # -- for the Verkko assembly, and the task recorded COMPLETED, exit 0. It had already
    # created out/, so the guard below passed too, and empty Haplotype_*.fasta files were
    # published as though they were a result. A failure that reports success is worse than a
    # failure: it is a wrong answer with a green tick next to it.
    #
    # So: test the real status, and test that the haplotypes have CONTENT rather than merely
    # existing. The errorStrategy still keeps the run alive; what changes is that the run says
    # so. Note that on a failed task Nextflow publishes nothing, so the log stays in the task
    # work directory, which .nextflow.log names.
    if [ "\$NP_STATUS" -ne 0 ]; then
        echo "NuclearPhaser exited \$NP_STATUS for ${meta.id}; see ${meta.id}.nuclearphaser.log" >&2
        tail -20 ${meta.id}.nuclearphaser.log >&2 || true
        exit "\$NP_STATUS"
    fi

    if ! ls out/Haplotype_*.fasta >/dev/null 2>&1; then
        echo "NuclearPhaser produced no haplotype FASTAs for ${meta.id}." >&2
        echo "  Usual causes: too few genes mapped, or too few trans Hi-C contacts." >&2
        echo "  Check ${meta.id}.nuclearphaser.log before concluding the assembly is unphaseable." >&2
        exit 1
    fi

    if ! find out -name 'Haplotype_*.fasta' -size +0c | grep -q .; then
        echo "NuclearPhaser wrote only EMPTY haplotype FASTAs for ${meta.id} -- it exited 0" >&2
        echo "  but produced nothing. Check ${meta.id}.nuclearphaser.log." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nuclearphaser: \$(NuclearPhaser.py 2>&1 | grep -io 'NuclearPhaser [0-9.]*' | head -1 || echo 'unknown')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p out
    touch out/Haplotype_0_genephasing.fasta out/Haplotype_1_genephasing.fasta
    touch ${meta.id}.nuclearphaser.log
    echo '"${task.process}": {nuclearphaser: stub}' > versions.yml
    """
}
