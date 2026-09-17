//
// Mitochondrial screening — the gap the M. larici-populina validation run exposed.
//
// Duplessis et al. (2026) removed the mitochondrion from their Melampsora haplotypes as a
// separate step; this pipeline did not, so the mito contig is sitting inside every haplotype
// it has produced. It inflates assembly size, it is the single contig running at many times
// nuclear depth, and it is one of the short contigs that put our contig counts above theirs.
//
// This runs AFTER FCS_GX_CLEAN, on its output, and for a reason: FCS-GX asks whether a
// sequence belongs to a different ORGANISM, and a mitochondrion does not. It is the right
// sequence in the wrong assembly, and no amount of contamination screening will find it.
//
// Two processes rather than one, and not by choice. The Galaxy depot's blast image has no
// python3 and the stock python image has no blastn, so a single process would need a
// container built for this one step. Splitting costs a hits table on disk and buys using
// only images this pipeline already proves elsewhere.
//
// BLASTn, not minimap2, and that was NOT the first design. Measured on the real mlp98AG31
// haplotypes against the RefSeq mitochondrion database, minimap2 -x asm20 was actively
// dangerous: it called 31 contigs / 3.41 Mb of hap1 mitochondrial. A fungal mitogenome is
// ~50-90 kb on one contig. The false calls were tandem-repeat arrays chaining onto unrelated
// mitogenomes -- h1tg000148l is seven copies of an ~11.8 kb unit, each matching 1,650 bases
// of an 11,800 bp block against a Penicillium mitochondrial scaffold, about 14% identity.
// minimap2 has no identity floor and its long chained blocks hide that.
//
// BLASTn with DUST masking and -perc_identity 90 -- the paper's own parameters -- separates
// them cleanly. Aligned fraction of each contig, intervals merged:
//
//     h1tg000096c    46,847 bp   0.349   the mitochondrion (circular, GC 30.4% vs 41%
//                                        nuclear, depth 522x vs 207x, hits A. psidii)
//     h1tg000270c    93,694 bp   0.349   EXACTLY 2x its length: a concatemer of the same
//                                        molecule, same GC -- correctly also called
//     h1tg000148l    74,764 bp   0.041   }
//     h1tg000149l   406,491 bp   0.037   } the repeat arrays minimap2 called mitochondrial
//     h1tg000030l   184,068 bp   0.028   }
//     h1tg000024l 1,088,795 bp   0.012   a chromosome carrying NUMTs
//
// True positives at 0.349, best false positive at 0.041. The default threshold sits in that
// gap, and it is a MEASURED gap rather than a guessed one.
//

process MITO_ALIGN {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(fasta)
    path mito_reference

    output:
    tuple val(meta), path("${meta.id}.mito_hits.tsv"), emit: hits
    path "versions.yml"                              , emit: versions

    when:
    params.run_mito_screen

    script:
    // -dust yes and -perc_identity 90 are the paper's parameters, not a choice made here:
    //     "Mitochondrial contigs were identified based on high coverage and via BLASTn
    //      against the NCBI RefSeq mitochondrial genome database (>=90% identity, with DUST
    //      masking parameters)"                        -- Duplessis et al., G3 2026
    // DUST is the load-bearing half. Mitogenomes and fungal repeat arrays are both AT-rich
    // and low-complexity, so without masking the database matches sequence composition
    // rather than homology, which is exactly how minimap2 came to call 3.41 Mb of tandem
    // repeat mitochondrial.
    //
    // The ASSEMBLY is the query. The calling rule is "what fraction of this CONTIG is
    // mitochondrial" -- a fraction of the QUERY length -- so swapping query and subject
    // would silently invert it into "how much of the mitogenome is present", which says
    // nothing about which contig it lives on.
    //
    // max_target_seqs is generous because the database holds ~18,000 mitogenomes and the
    // nearest relative may be several genera away; for M. larici-populina there is no
    // Melampsora mitogenome in RefSeq at all and the best hit is Austropuccinia psidii.
    """
    ASM=${fasta}
    if [[ "\$ASM" == *.gz ]]; then
        gunzip -c "\$ASM" > asm.fa
        ASM=asm.fa
    fi

    REF=${mito_reference}
    if [[ "\$REF" == *.gz ]]; then
        gunzip -c "\$REF" > mito_ref.fa
        REF=mito_ref.fa
    fi

    makeblastdb -in "\$REF" -dbtype nucl -out mito_db > makeblastdb.log 2>&1

    blastn \\
        -query "\$ASM" \\
        -db mito_db \\
        -dust yes \\
        -perc_identity ${params.mito_min_identity} \\
        -num_threads ${task.cpus} \\
        -max_target_seqs 50 \\
        -outfmt '6 qseqid qlen qstart qend pident length sseqid' \\
        > ${meta.id}.mito_hits.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version | head -1 | sed 's/blastn: //')
    END_VERSIONS
    """

    stub:
    """
    : > ${meta.id}.mito_hits.tsv
    echo '"${task.process}": {blast: stub}' > versions.yml
    """
}

process MITO_CALL {
    tag "${meta.id}"
    label 'process_single'

    publishDir "${params.outdir}/derived/mito_screen", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(fasta), path(hits)

    output:
    tuple val(meta), path("${meta.id}.mito_free.fa.gz")  , emit: fasta
    tuple val(meta), path("${meta.id}.mito.fa.gz")       , emit: mito
    tuple val(meta), path("${meta.id}.mito_screen.tsv")  , emit: manifest
    tuple val(meta), path("${meta.id}.mito_screen.json") , emit: json
    path "versions.yml"                                  , emit: versions

    when:
    params.run_mito_screen

    script:
    // The calling rule, the NUMT problem and the reasoning behind both thresholds are
    // documented in bin/mito_screen.py. The safeguards are the ones FCS_GX_CLEAN uses: the
    // input assembly is never modified, the removed sequence is written out rather than
    // deleted, and every contig with ANY alignment is recorded with its call so a surprising
    // one can be argued with.
    """
    mito_screen.py \\
        --fasta ${fasta} \\
        --hits ${hits} \\
        --assembly-id ${meta.id} \\
        --min-aligned-frac ${params.mito_min_aligned_frac} \\
        --max-length ${params.mito_max_length} \\
        --out-fasta ${meta.id}.mito_free.fa.gz \\
        --out-mito ${meta.id}.mito.fa.gz \\
        --out-tsv ${meta.id}.mito_screen.tsv \\
        --out-json ${meta.id}.mito_screen.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo ">stub" | gzip -c > ${meta.id}.mito_free.fa.gz
    echo ">stub" | gzip -c > ${meta.id}.mito.fa.gz
    echo -e "assembly_id\\tcontig\\tlength\\taligned_bp\\taligned_fraction\\tgc\\tcall\\treason" > ${meta.id}.mito_screen.tsv
    echo '{"assembly_id": "${meta.id}", "mito_contigs": 0}' > ${meta.id}.mito_screen.json
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
