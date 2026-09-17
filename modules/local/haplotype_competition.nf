process HAPLOTYPE_COMPETITION {
    tag "${meta.id}"
    label 'process_high'

    publishDir path: { "${params.outdir}/assembly_qc/haplotype_competition/${meta.id}" }, mode: params.publish_dir_mode,
        saveAs: { fn -> (fn.endsWith('.bam') || fn.endsWith('.bai')) ? null : fn }

    input:
    tuple val(meta), path(hap1), path(hap2), path(reads)

    output:
    tuple val(meta), path("${meta.id}.haplotype_competition.tsv") , emit: tsv
    tuple val(meta), path("${meta.id}.haplotype_competition.json"), emit: json
    path "versions.yml"                                           , emit: versions

    when:
    params.run_haplotype_competition

    script:
    // Does a "phased" haplotype pair actually represent two nuclei?
    //
    // MINIMAP2_ASSEMBLY maps reads to each assembly ALONE, which cannot answer that: with
    // hap1 as the only reference, reads from both nuclei have nowhere else to go, so a
    // perfectly phased haplotype still draws ~2x its own depth. Every "collapsed/shared"
    // call this pipeline made on a half-genome assembly is an artefact of that setup.
    //
    // Here the two haplotypes compete for the same reads in one index. Mean depth still
    // cannot separate the cases -- hap1+hap2 is a whole genome either way, so the average is
    // ~1n regardless -- so the discriminator is MAPPING QUALITY, not depth:
    //
    //   truly phased    reads have one good home; divergence makes it unique; MAPQ high
    //   collapsed pair  the same locus twice; every read ties; minimap2 emits MAPQ 0
    //
    // Some MAPQ 0 is expected even from a correct pair (conserved and repetitive sequence is
    // genuinely identical between nuclei), so the number is read comparatively, across
    // samples, not against an absolute threshold.
    """
    set -o pipefail

    # Namespace the contigs before concatenating: hifiasm's h1tg/h2tg names happen to differ
    # already, but flye and hicanu haplotypes would collide and silently lose half the reads.
    for H in 1 2; do
        SRC=\$([ "\$H" = 1 ] && echo "${hap1}" || echo "${hap2}")
        if [[ "\$SRC" == *.gz ]]; then gunzip -c "\$SRC"; else cat "\$SRC"; fi \\
            | awk -v h="\$H" '/^>/ { sub(/^>/, ">hap" h "|"); print; next } { print }'
    done > combined.fa

    SORT_THREADS=\$(( ${task.cpus} / 2 ))
    if [ "\$SORT_THREADS" -lt 1 ]; then SORT_THREADS=1; fi

    minimap2 -ax map-hifi -t ${task.cpus} combined.fa ${reads} \\
        | samtools sort -@ "\$SORT_THREADS" -m 1G -o ${meta.id}.bam -
    samtools index -@ ${task.cpus} ${meta.id}.bam

    # Primary alignments only, and mapped: -F 0x904 drops secondary (0x100), supplementary
    # (0x800) and unmapped (0x4), so each read is counted once at the single place minimap2
    # committed to. Unmapped reads carry '*' as their reference and would otherwise appear as
    # a third "haplotype" and skew every read_share.
    samtools view -F 0x904 ${meta.id}.bam \\
        | awk -F'\\t' '
            {
                split(\$3, a, "|"); hap = a[1]
                n[hap]++; total++
                if (\$5 == 0)  zero[hap]++
                if (\$5 < 10)  low[hap]++
                sum[hap] += \$5
            }
            END {
                for (h in n)
                    printf "%s\\t%d\\t%.4f\\t%.4f\\t%.2f\\t%.4f\\n",
                        h, n[h], zero[h]/n[h], low[h]/n[h], sum[h]/n[h], n[h]/total
            }' \\
        | sort > mapq.tsv

    # Mean depth per haplotype, from the same competitive alignment.
    samtools depth -a ${meta.id}.bam \\
        | awk -F'\\t' '{ split(\$1, a, "|"); s[a[1]] += \$3; n[a[1]]++ }
                       END { for (h in s) printf "%s\\t%.2f\\n", h, s[h]/n[h] }' \\
        | sort > depth.tsv

    printf 'haplotype\\tprimary_reads\\tmapq0_frac\\tmapq_lt10_frac\\tmean_mapq\\tread_share\\tmean_depth\\n' \\
        > ${meta.id}.haplotype_competition.tsv
    # awk, not join(1): this is the minimap2/samtools container, which ships neither coreutils'
    # join nor a python interpreter. Same reason the JSON below is assembled by hand.
    awk -F'\\t' 'NR == FNR { d[\$1] = \$2; next }
                 { printf "%s\\t%s\\n", \$0, (\$1 in d ? d[\$1] : "NA") }' \\
        depth.tsv mapq.tsv >> ${meta.id}.haplotype_competition.tsv

    # JSON built with awk, not python: this container has no interpreter (and no join(1),
    # no gawk, no mawk -- its awk is busybox). Busybox awk parses NAME( as a function call,
    # so a ternary written inline after a variable, `rec (i > 1 ? ...)`, dies with "Call to
    # undefined function". Hence the plain if/else below: it is not stylistic, it is the only
    # form this awk accepts. Values arrive via -v rather than string interpolation so an id
    # can never break the program.
    awk -F'\\t' \\
        -v ID="${meta.id}" -v SAMPLE="${meta.sample}" -v ASM="${meta.assembler}" '
        NR == 1 { for (i = 1; i <= NF; i++) h[i] = \$i; next }
        {
            rec = "    {"
            for (i = 1; i <= NF; i++) {
                sep = ""
                if (i > 1) sep = ", "
                val = \$i
                if (i == 1) val = "\\"" \$i "\\""
                rec = rec sep "\\"" h[i] "\\": " val
            }
            n = n + 1
            recs[n] = rec "}"
        }
        END {
            printf "{\\n  \\"id\\": \\"%s\\",\\n  \\"sample\\": \\"%s\\",\\n", ID, SAMPLE
            printf "  \\"assembler\\": \\"%s\\",\\n  \\"haplotypes\\": [\\n", ASM
            for (i = 1; i <= n; i++) {
                comma = ","
                if (i == n) comma = ""
                printf "%s%s\\n", recs[i], comma
            }
            printf "  ]\\n}\\n"
        }' ${meta.id}.haplotype_competition.tsv > ${meta.id}.haplotype_competition.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    printf 'haplotype\\tprimary_reads\\tmapq0_frac\\tmapq_lt10_frac\\tmean_mapq\\tread_share\\tmean_depth\\n' \\
        > ${meta.id}.haplotype_competition.tsv
    printf 'hap1\\t1000\\t0.1000\\t0.1200\\t45.00\\t0.5000\\t17.50\\n' >> ${meta.id}.haplotype_competition.tsv
    printf 'hap2\\t1000\\t0.1000\\t0.1200\\t45.00\\t0.5000\\t17.50\\n' >> ${meta.id}.haplotype_competition.tsv
    echo '{"id": "${meta.id}", "haplotypes": []}' > ${meta.id}.haplotype_competition.json
    echo '"${task.process}": {minimap2: stub}' > versions.yml
    """
}
