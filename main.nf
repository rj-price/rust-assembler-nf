#!/usr/bin/env nextflow
/*
 * rust-assembler-nf
 *
 * Phased assembly of a dikaryotic rust (Pucciniales) genome from PacBio HiFi reads.
 * A research/decision pipeline: it produces several assembly candidates plus the evidence
 * to judge them, rather than one "final" FASTA.
 *
 * See README.md for usage and for how to interpret the output.
 */

nextflow.enable.dsl = 2

include { READ_QC       } from './subworkflows/local/read_qc'
include { KMER_ANALYSIS } from './subworkflows/local/kmer_analysis'
include { CONTAMINATION } from './subworkflows/local/contamination'
include { ASSEMBLY      } from './subworkflows/local/assembly'
include { ASSEMBLY_QC   } from './subworkflows/local/assembly_qc'
include { PHASING       } from './subworkflows/local/phasing'
include { SCAFFOLDING   } from './subworkflows/local/scaffolding'
include { REPORTING     } from './subworkflows/local/reporting'

/* ------------------------------------------------------------------------------------
 * Help
 * ---------------------------------------------------------------------------------- */

def helpMessage() {
    log.info """
    rust-assembler-nf v${workflow.manifest.version}

    Usage:
      nextflow run . --input samplesheet.csv --genome_size 525m \\
          --outdir /scratch/me/results -profile <site>

    Required:
      --input                  Samplesheet CSV with columns: sample,run,fastq
      --genome_size            Haploid size of ONE nucleus, e.g. 525m [${params.genome_size}]
      --fcs_gx_taxid           NCBI taxid of your species (FCS-GX is on by default;
                               or set --run_fcs_gx false)

    Key options (see nextflow.config for the full set):
      --outdir                 Output directory [${params.outdir}]
      --assembler_genome_size  Size given to Flye/HiCanu; defaults to --genome_size
      --read_yield_bases       Total sequenced bases; enables expected-depth scoring
      --run_hifiasm/_flye/_hicanu/_ipa/_verkko
                               Which assemblers to run [verkko off, rest on]
      --ipa_sif                Locally built IPA image (containers/ipa.def), needed with --run_ipa
      --subset_assemblies      Diagnostic ALL/run1/run2 hifiasm comparison [${params.subset_assemblies}]
      --qc_only                Re-run QC + reporting over existing assemblies, no rebuild
      --assemblies_from        Directory of assembly FASTAs to QC (with --qc_only)
      --meryl_db               Existing meryl DB, or a directory of per-sample .meryl DBs
      --yields_from            Directory of a previous run's *.seqkit_stats.tsv (with --qc_only)
      --hic_r1 / --hic_r2      Enables fully-phased hifiasm Hi-C assembly
      --host_reference         Optional host FASTA for contamination mapping
      --host_filtered_assembly Assemble host-removed reads too (needs --host_reference)
      --filtered_assembly      Assemble the quality/length-filtered dataset too
      --run_fcs_gx             Contamination screen with FCS-GX [${params.run_fcs_gx}]
      --fcs_gx_clean           Write contamination-cleaned assemblies [${params.fcs_gx_clean}]
      --run_mito_screen        Remove mitochondrial contigs from the cleaned assemblies [${params.run_mito_screen}]
      --mito_reference         Mitogenome FASTA to screen against, e.g. the NCBI RefSeq
                               mitochondrion release (required with --run_mito_screen)
      --run_nuclearphaser      Hi-C phasing of any assembly (needs --hic_r1/--hic_r2 and
                               --nuclearphaser_genes, --nuclearphaser_sif)
      --nuclearphaser_targets  Regex of assembly ids to phase [${params.nuclearphaser_targets}]
      --nuclearphaser_skip_collapsed
                               Skip targets too small to hold both nuclei [${params.nuclearphaser_skip_collapsed}]
      --run_scaffolding        Hi-C scaffolding of phased haplotypes (needs --hic_r1/--hic_r2)
      --scaffold_targets       Regex of assembly ids to scaffold [${params.scaffold_targets}]
      --scaffolder             yahs | haphic [${params.scaffolder}]
      --scaffold_n_chromosomes Expected chromosomes PER HAPLOTYPE (required by haphic)
      --scaffold_chromosome_min_length
                               Length at or above which a scaffold counts as chromosome-scale
                               in scaffold_summary.tsv. Reporting only — changes no sequence
                               and no decision. [${params.scaffold_chromosome_min_length}]
      --min_read_length        [${params.min_read_length}]  (0 = no filtering)
      --min_read_q             [${params.min_read_q}]  (0 = no filtering)
      --busco_lineage          [${params.busco_lineage}]
      --kraken2_confidence     Min fraction of a read's k-mers supporting a taxon
                               [${params.kraken2_confidence}]  (Kraken2's own default is 0.0)

    Site profiles supply cluster paths, databases and queue routing:
      gruffalo               Crop Diversity HPC (slurm + apptainer + its databases)
      slurm / apptainer      Building blocks for your own site profile

    Dataset profiles supply species/dataset values:
      mlp                    Melampsora larici-populina 98AG31, the published validation set
      map                    Melampsora allii-populina 12AY07, second validation set

    Development: test (smoke test), stub (wiring only), robust (unattended), debug
    """.stripIndent()
}

if (params.help) { helpMessage(); exit 0 }

/* ------------------------------------------------------------------------------------
 * Parameter validation — fail fast, before anything expensive is submitted
 * ---------------------------------------------------------------------------------- */

def validateParams() {
    def errors = []

    if (!params.input) {
        errors << "--input is required (samplesheet CSV with columns: sample,run,fastq)"
    }
    else if (!file(params.input).exists()) {
        errors << "--input samplesheet not found: ${params.input}"
    }

    // Genome size has no defensible default across rusts (~100 Mb to ~2 Gb), so it is
    // required rather than guessed. Parse it here so a typo fails now rather than four days
    // into an assembly.
    if (!params.genome_size) {
        errors << "--genome_size is required (haploid size of ONE nucleus, e.g. 525m)"
    }
    else {
        try { GenomeSize.parse(params.genome_size) }
        catch (Exception e) { errors << e.message }
    }
    if (params.assembler_genome_size) {
        try { GenomeSize.parse(params.assembler_genome_size) }
        catch (Exception e) { errors << e.message }
    }

    // Database paths are site-specific and default to null, so distinguish "not configured"
    // from "configured but missing" — they need different fixes.
    if (params.run_kraken2) {
        if (!params.kraken_db) {
            errors << "--kraken_db is not set (use a site profile, pass it, or --run_kraken2 false)"
        }
        else if (!file(params.kraken_db).exists()) {
            errors << "Kraken2 database not found: ${params.kraken_db}"
        }

        // Kraken2 accepts a confidence outside [0,1] on the command line and then classifies
        // nothing, which looks like a clean run with an empty report. Catch it here instead.
        def conf = params.kraken2_confidence
        if (!(conf instanceof Number) || conf < 0 || conf > 1) {
            errors << "--kraken2_confidence must be a number between 0 and 1 (got: ${conf})"
        }
    }

    if (params.run_busco) {
        if (!params.busco_db) {
            errors << "--busco_db is not set (use a site profile, pass it, or --run_busco false)"
        }
        else if (!file("${params.busco_db}/${params.busco_lineage}").exists()) {
            errors << "BUSCO lineage not found: ${params.busco_db}/${params.busco_lineage}"
        }
    }

    if (params.run_fcs_gx) {
        if (!params.fcs_gx_sif) {
            errors << "--fcs_gx_sif is not set (use a site profile, pass it, or --run_fcs_gx false)"
        }
        else if (!file(params.fcs_gx_sif).exists()) {
            errors << "FCS-GX image not found: ${params.fcs_gx_sif}"
        }
        if (!params.fcs_gx_taxid) {
            errors << "--fcs_gx_taxid is required with --run_fcs_gx (NCBI taxid of your species)"
        }
    }

    // Some clusters back up and file-count-cap \$HOME, where publishing an assembly tree
    // throttles job concurrency. Site profiles opt into this guard.
    def outdir_abs = file(params.outdir).toAbsolutePath().toString()
    if (params.forbid_home_outdir && outdir_abs.startsWith(System.getProperty('user.home'))) {
        errors << "--outdir must not be inside \$HOME (got: ${outdir_abs}). Use scratch."
    }

    if (params.host_reference && !file(params.host_reference).exists()) {
        errors << "--host_reference not found: ${params.host_reference}"
    }

    // Hi-C is all-or-nothing.
    if ((params.hic_r1 && !params.hic_r2) || (params.hic_r2 && !params.hic_r1)) {
        errors << "Hi-C requires BOTH --hic_r1 and --hic_r2"
    }

    // Derived-dataset assemblies each depend on the step that produces the dataset. Catch the
    // combination that silently produces nothing, rather than letting the run finish with a
    // readset the user expected and never got.
    if (params.host_filtered_assembly && !params.host_reference) {
        errors << "--host_filtered_assembly requires --host_reference (there is nothing to remove without one)"
    }
    if (params.filtered_assembly && !(params.min_read_length > 0 || params.min_read_q > 0)) {
        errors << "--filtered_assembly requires --min_read_length and/or --min_read_q (no filter set = no filtered dataset)"
    }

    // IPA needs a locally built image: every stock pbipa biocontainer is broken (BusyBox
    // sort, and snakemake 8 on the newest build). Caught here rather than at task time --
    // IPA sits behind nothing, so a missing image would otherwise surface only after the
    // scheduler had found a 64-core node.
    if (params.run_ipa && !params.ipa_sif) {
        errors << "--run_ipa requires --ipa_sif (build it once: apptainer build ipa.sif containers/ipa.def, or set --run_ipa false)"
    }

    // Mitochondrial screening is reference-based; without a mitogenome to align to there is
    // nothing to do. It also sits downstream of FCS_GX_CLEAN, so the cleaning it filters must
    // actually be running -- the same dependency NuclearPhaser has, for the same reason.
    if (params.run_mito_screen) {
        if (!params.mito_reference) {
            errors << "--run_mito_screen requires --mito_reference (a mitogenome FASTA, ideally same genus)"
        }
        else if (!file(params.mito_reference).exists()) {
            errors << "--mito_reference not found: ${params.mito_reference}"
        }
        if (!params.run_fcs_gx || !params.fcs_gx_clean) {
            errors << "--run_mito_screen filters the CLEANED assemblies, so it needs --run_fcs_gx and --fcs_gx_clean"
        }
        if (params.mito_min_identity <= 0 || params.mito_min_identity > 100) {
            errors << "--mito_min_identity is a BLASTn percentage, so it must be >0 and <=100, got ${params.mito_min_identity}"
        }
        if (params.mito_min_aligned_frac <= 0 || params.mito_min_aligned_frac > 1) {
            errors << "--mito_min_aligned_frac must be >0 and <=1 (it is a fraction of contig length), got ${params.mito_min_aligned_frac}"
        }
    }

    // NuclearPhaser: three separate prerequisites, and missing any one of them is a different
    // fix, so report them separately.
    if (params.run_nuclearphaser) {
        if (!params.hic_r1 || !params.hic_r2) {
            errors << "--run_nuclearphaser requires --hic_r1 and --hic_r2"
        }
        if (!params.nuclearphaser_genes) {
            errors << "--run_nuclearphaser requires --nuclearphaser_genes (FASTA of gene sequences for the synteny signal)"
        }
        else if (!file(params.nuclearphaser_genes).exists()) {
            errors << "--nuclearphaser_genes not found: ${params.nuclearphaser_genes}"
        }
        if (!params.nuclearphaser_sif) {
            errors << "--nuclearphaser_sif is not set. No biocontainer exists for NuclearPhaser; " +
                     "build the bundled definition once with: apptainer build nuclearphaser.sif containers/nuclearphaser.def"
        }
        else if (!file(params.nuclearphaser_sif).exists()) {
            errors << "NuclearPhaser image not found: ${params.nuclearphaser_sif}"
        }
        if (!params.run_busco) {
            errors << "--run_nuclearphaser needs BUSCO (it reads the per-gene full_table.tsv); do not set --run_busco false"
        }
        if (!params.run_fcs_gx || !params.fcs_gx_clean) {
            errors << "--run_nuclearphaser phases the CLEANED assemblies, so it needs --run_fcs_gx and --fcs_gx_clean"
        }
    }

    // Scaffolding. Fewer prerequisites than phasing, but the same principle: refuse at
    // validation rather than at the end of a run that has already spent days on assemblies.
    if (params.run_scaffolding) {
        if (!params.hic_r1 || !params.hic_r2) {
            errors << "--run_scaffolding requires --hic_r1 and --hic_r2 (it scaffolds from the Hi-C contact map)"
        }
        if (!(params.scaffolder in ['yahs', 'haphic'])) {
            errors << "--scaffolder must be 'yahs' or 'haphic', not '${params.scaffolder}'"
        }
        if (params.scaffolder == 'haphic') {
            // HapHiC clusters into exactly this many groups, so a wrong or absent value does
            // not degrade gracefully -- it produces confidently wrong chromosomes.
            if (!params.scaffold_n_chromosomes) {
                errors << "--scaffolder haphic requires --scaffold_n_chromosomes (expected chromosomes PER HAPLOTYPE; " +
                         "18 for M. larici-populina, not the 36 of the dikaryon)"
            }
            if (!params.haphic_sif) {
                errors << "--haphic_sif is not set. No biocontainer exists for HapHiC; " +
                         "build the bundled definition once with: apptainer build haphic.sif containers/haphic.def"
            }
            else if (!file(params.haphic_sif).exists()) {
                errors << "HapHiC image not found: ${params.haphic_sif}"
            }
        }
        if (!params.run_fcs_gx || !params.fcs_gx_clean) {
            errors << "--run_scaffolding scaffolds the CLEANED assemblies, so it needs --run_fcs_gx and --fcs_gx_clean"
        }
    }

    if (params.qc_only) {
        // In QC-only mode the assemblers are irrelevant; what matters is that there are
        // assemblies to read and a meryl DB to score them against.
        if (!params.assemblies_from) {
            errors << "--qc_only requires --assemblies_from <dir containing assembly FASTAs>"
        }
        else if (!file(params.assemblies_from).exists()) {
            errors << "--assemblies_from not found: ${params.assemblies_from}"
        }
        if (params.run_merqury && !params.meryl_db) {
            errors << "--qc_only needs --meryl_db <existing meryl DB> for Merqury (or --run_merqury false)"
        }
        if (params.meryl_db && !file(params.meryl_db).exists()) {
            errors << "--meryl_db not found: ${params.meryl_db}"
        }
        if (params.yields_from && !file(params.yields_from).exists()) {
            errors << "--yields_from not found: ${params.yields_from}"
        }
    }
    else if (!(params.run_hifiasm || params.run_flye || params.run_hicanu || params.run_verkko)) {
        errors << "No assembler enabled — nothing to do"
    }

    if (errors) {
        log.error "Parameter validation failed:\n  - " + errors.join("\n  - ")
        exit 1
    }
}

/* ------------------------------------------------------------------------------------
 * Samplesheet parsing
 * ---------------------------------------------------------------------------------- */

def parseSamplesheet(path) {
    Channel
        .fromPath(path)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.run || !row.fastq) {
                error "Samplesheet row missing sample/run/fastq: ${row}"
            }
            def fq = file(row.fastq)
            if (!fq.exists()) {
                error "FASTQ not found (samplesheet row '${row.run}'): ${row.fastq}"
            }
            def meta = [
                id     : "${row.sample}_${row.run}".toString(),
                sample : row.sample,
                run    : row.run
            ]
            tuple(meta, fq)
        }
}

/* ------------------------------------------------------------------------------------
 * QC-only re-entry: rebuild the ASSEMBLIES channel from FASTAs already on disk
 * ---------------------------------------------------------------------------------- */

//
// Recover (sample, readset, assembler, assembly_type) from a published filename.
//
// Two naming conventions are in play, because GFA2FASTA and the direct assemblers publish
// differently:
//   rust_all_hifiasm_partially_phased_hap1.fa.gz   <- <sample>_<readset>_<assembler>_<type>
//   rust_all.flye.fasta.gz                         <- <sample>_<readset>.<assembler>
// Rather than guess by position (readset and assembly_type are both multi-token), locate the
// assembler token, which is from a known closed set, and split around it.
//
def classifyAssemblyFile(String filename) {
    // Declared inside the function on purpose: a script-level `def` in Nextflow is a local of
    // the script's run method, so it is not in scope inside a function body.
    def KNOWN_ASSEMBLERS = ['hifiasm', 'flye', 'hicanu', 'verkko']

    def base = filename.replaceAll(/\.(fa|fasta)\.gz$/, '').replaceAll(/\.(fa|fasta)$/, '')

    // Deliberately NOT stripping a .cleaned suffix here. Doing so makes meta.id equal the
    // ORIGINAL assembly's id, so FCS_GX_CLEAN's "${meta.id}.cleaned.fa.gz" output collides
    // with the name of its own staged input -- and because Nextflow stages inputs as
    // symlinks, the redirect writes straight through the link and destroys the source file.
    // That is exactly what happened on 2026-09-08: a --stub-run over derived/cleaned_assemblies
    // truncated all 35 cleaned FASTAs to empty gzips. Keeping the suffix in the id keeps the
    // two filenames distinct; strip it downstream if a summary needs to be diffed by id.

    def assembler = null
    def type      = 'primary'
    def stem      = base

    // Dotted form: <stem>.<assembler>
    def dotted = base.tokenize('.')
    if (dotted.size() > 1 && dotted[-1] in KNOWN_ASSEMBLERS) {
        assembler = dotted[-1]
        stem      = dotted[0..-2].join('.')
    }
    else {
        def toks = base.tokenize('_')
        def idx  = toks.findIndexOf { it in KNOWN_ASSEMBLERS }
        if (idx < 0) return null
        assembler = toks[idx]
        stem      = toks[0..<idx].join('_')
        if (idx + 1 < toks.size()) type = toks[(idx + 1)..-1].join('_')
    }

    // stem is <sample>_<readset>. Most readsets are one token (all, filtered, run1, run2...)
    // but host_filtered is two, and taking the last token alone would split it into a sample
    // called '<sample>_host' with readset 'filtered' — a candidate that then fails to join
    // against anything. So match the multi-token readsets explicitly, longest first, and only
    // fall back to the last-token rule for the single-token ones.
    def MULTI_TOKEN_READSETS = ['host_filtered']

    def st = stem.tokenize('_')
    if (st.size() < 2) return null

    def readset = null
    def sample  = null
    for (rs in MULTI_TOKEN_READSETS) {
        def n = rs.tokenize('_').size()
        if (st.size() > n && st[-n..-1].join('_') == rs) {
            readset = rs
            sample  = st[0..(st.size() - n - 1)].join('_')
            break
        }
    }
    if (readset == null) {
        readset = st[-1]
        sample  = st[0..-2].join('_')
    }

    // Rebuild the id from its parts rather than reusing the filename, so that a dotted name
    // (rust_all.flye) yields the same id as the pipeline's own convention
    // (rust_all_flye_primary) and results stay comparable across runs.
    return [ id           : [sample, readset, assembler, type].join('_'),
             sample       : sample,
             readset      : readset,
             assembler    : assembler,
             assembly_type: type ]
}

def assembliesFromDir(String dir) {
    // Two patterns, not one: Nextflow's `**` matches one or more directories, so a lone
    // `**/*.fa.gz` misses assemblies sitting directly in the given directory -- which is
    // exactly how derived/cleaned_assemblies is laid out. checkIfExists is dropped because
    // one of the two patterns is legitimately empty in either layout; the emptiness check
    // moves to the channel instead, where it can still fail loudly.
    Channel
        .fromPath(["${dir}/*.{fa,fasta}.gz", "${dir}/**/*.{fa,fasta}.gz"])
        // FCS_GX_CLEAN writes the removed sequence next to the cleaned assembly as
        // <id>.contaminants.fa.gz. Pointing --assemblies_from at that directory would
        // otherwise QC the contaminant bins as though they were assemblies.
        .filter { fa -> !fa.name.contains('.contaminants.') }
        .ifEmpty { error "No assembly FASTAs (*.fa.gz / *.fasta.gz) found under: ${dir}" }
        .map { fa ->
            def meta = classifyAssemblyFile(fa.name)
            if (!meta) {
                error "Cannot classify assembly file (expected <sample>_<readset>_<assembler>_<type>): ${fa}"
            }
            tuple(meta, fa, [])
        }
}

/* ------------------------------------------------------------------------------------
 * Main workflow
 * ---------------------------------------------------------------------------------- */

workflow {

    validateParams()

    ch_versions = Channel.empty()

    // Per-run reads. Run identity is preserved through QC so the report can compare the
    // two runs' very different character (11 kb/Q40 vs 27.5 kb/Q29).
    ch_reads_per_run = parseSamplesheet(params.input)

    // Combined dataset used for assembly: one entry per sample, all runs pooled.
    ch_reads_combined = ch_reads_per_run
        .map { meta, fq -> tuple(meta.sample, fq) }
        .groupTuple()
        .map { sample, fqs ->
            tuple([ id: "${sample}_all".toString(), sample: sample, readset: 'all' ], fqs)
        }

    //
    // QC-ONLY BRANCH: assemblies already exist; gather evidence over them and stop.
    //
    // Deliberately skips READ_QC, KMER_ANALYSIS, CONTAMINATION and ASSEMBLY. The reads are
    // still parsed above because MINIMAP2_ASSEMBLY needs them, but nothing re-derives what
    // is already on disk.
    //
    if (params.qc_only) {
        log.info """
        QC-only mode: reading assemblies from ${params.assemblies_from}
          meryl DB : ${params.meryl_db ?: '(none — Merqury disabled)'}
          Assemblers, read QC, k-mer analysis and read-level contamination are SKIPPED.
        """.stripIndent()

        ch_assemblies = assembliesFromDir(params.assemblies_from)

        // --meryl_db may be one DB or a directory of them. A directory is how a multi-sample
        // project re-QCs: each <sample>_<readset>.k*.meryl is tagged with its own sample so
        // MERQURY matches assembly to k-mers. A single DB keeps meta.sample null, which
        // ASSEMBLY_QC reads as the project-wide '*' fallback and applies to everything --
        // right for a one-sample project, wrong for seven, hence the distinction.
        ch_meryl = Channel.empty()
        if (params.meryl_db) {
            def meryl_path = file(params.meryl_db, checkIfExists: true)
            if (meryl_path.isDirectory() && !meryl_path.name.endsWith('.meryl')) {
                ch_meryl = Channel
                    .fromPath("${params.meryl_db}/*.meryl", type: 'dir', checkIfExists: true)
                    .map { db -> tuple([ id: db.name, sample: db.name.tokenize('_')[0] ], db) }
            }
            else {
                ch_meryl = Channel.value(tuple([id: 'meryl'], meryl_path))
            }
        }

        // READ_QC does not run here, so yields are read back from a previous run's published
        // seqkit tables when --yields_from points at one. Without it the map is empty and
        // every assembly falls back to params.read_yield_bases, which is null on a project
        // whose samples differ enough in yield that one figure would misjudge six of them.
        ch_yield_map = params.yields_from
            ? Channel.fromPath("${params.yields_from}/*.seqkit_stats.tsv", checkIfExists: true)
                // The sample comes from the TSV's own name (<sample>_<readset>.seqkit_stats.tsv):
                // the table's `file` column holds the raw read filename, which on this project
                // is a PacBio well id and carries no sample identity at all.
                .map { tsv -> tuple(tsv.name.tokenize('_')[0], tsv) }
                .splitCsv(sep: '\t', header: true, elem: 1)
                .map { sample, row -> tuple(sample, (row['sum_len'].replaceAll(',', '') as BigInteger)) }
                .groupTuple()
                .map { sample, sums -> [ (sample): sums.sum() ] }
                .reduce([:]) { acc, entry -> acc + entry }
                .ifEmpty([:])
            : Channel.value([:])

        ASSEMBLY_QC(ch_assemblies, ch_reads_combined, ch_meryl, ch_yield_map)
        ch_versions = ch_versions.mix(ASSEMBLY_QC.out.versions)

        // Phasing is available here too: it consumes cleaned assemblies and BUSCO tables,
        // both of which ASSEMBLY_QC has just produced. This is also the re-entry point for
        // NuclearPhaser's second pass, after a human has broken the phase-switched contigs
        // its first pass reported.
        ch_np_haplotypes = Channel.empty()
        if (params.run_nuclearphaser) {
            PHASING(ASSEMBLY_QC.out.cleaned, ASSEMBLY_QC.out.busco_tables, ASSEMBLY_QC.out.gfastats)
            ch_versions = ch_versions.mix(PHASING.out.versions)
            ch_np_haplotypes = PHASING.out.haplotypes
        }

        // Scaffolding is independent of phasing: hifiasm's own hap1/hap2 are already phased
        // haplotypes and can be scaffolded without NuclearPhaser ever running. When phasing
        // HAS run, its haplotypes are scaffolded too.
        ch_scaffold_mqc = Channel.empty()
        if (params.run_scaffolding) {
            SCAFFOLDING(ASSEMBLY_QC.out.cleaned, ch_np_haplotypes)
            ch_versions = ch_versions.mix(SCAFFOLDING.out.versions)
            ch_scaffold_mqc = SCAFFOLDING.out.mqc
        }

        REPORTING(
            ASSEMBLY_QC.out.for_summary,
            ASSEMBLY_QC.out.depth_hists,
            ASSEMBLY_QC.out.expected_depths,
            ASSEMBLY_QC.out.for_multiqc.mix(ch_scaffold_mqc).collect(),
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

        return
    }

    //
    // STAGE: read QC. No filtering by default (DECISION D2).
    //
    READ_QC(ch_reads_per_run)
    ch_versions = ch_versions.mix(READ_QC.out.versions)

    //
    // STAGE: k-mer analysis. Runs before assembly because hifiasm's coverage-based
    // homozygous/heterozygous decisions are exactly what goes wrong on an odd genome.
    //
    KMER_ANALYSIS(ch_reads_combined)
    ch_versions = ch_versions.mix(KMER_ANALYSIS.out.versions)

    //
    // STAGE: read-level contamination. Assessed, never silently acted on.
    //
    CONTAMINATION(ch_reads_combined)
    ch_versions = ch_versions.mix(CONTAMINATION.out.versions)

    //
    // Derived datasets that are worth assembling in their own right. Both are empty unless
    // explicitly asked for, and both exist to make a question answerable rather than assumed:
    //
    //   host_filtered  did removing apparent host contamination actually improve anything?
    //   filtered       did the quality/length filter actually improve anything?
    //
    // Neither replaces the raw `all` readset — they are assembled ALONGSIDE it, which is the
    // only arrangement in which the comparison means anything (DECISION D6).
    //
    ch_extra_readsets = Channel.empty()

    if (params.host_filtered_assembly) {
        ch_extra_readsets = ch_extra_readsets.mix(CONTAMINATION.out.host_filtered_reads)
    }
    if (params.filtered_assembly) {
        ch_extra_readsets = ch_extra_readsets.mix(READ_QC.out.filtered_combined)
    }

    //
    // STAGE: assembly. Each assembler is an independent, separately-cacheable process.
    //
    ASSEMBLY(ch_reads_combined, ch_reads_per_run, ch_extra_readsets)
    ch_versions = ch_versions.mix(ASSEMBLY.out.versions)

    //
    // STAGE: assembly QC. Consumes the ASSEMBLIES channel blind to provenance, so any new
    // assembler (or Hi-C haplotype) inherits the full QC suite for free.
    //
    // Per-sample yields, folded into one map so ASSEMBLY_QC can look up each assembly's
    // own sample. reduce() gives a value channel, which combine() broadcasts rather than
    // consumes, and ifEmpty keeps the lookup total if SEQKIT_STATS produced nothing.
    ch_yield_map = READ_QC.out.sample_yield
        .map { sample, bases -> [ (sample): bases ] }
        .reduce([:]) { acc, entry -> acc + entry }
        .ifEmpty([:])

    ASSEMBLY_QC(
        ASSEMBLY.out.assemblies,
        ch_reads_combined,
        KMER_ANALYSIS.out.meryl_db,
        ch_yield_map
    )
    ch_versions = ch_versions.mix(ASSEMBLY_QC.out.versions)

    //
    // STAGE: Hi-C phasing of the finished assemblies. Unlike hifiasm's --h1/--h2, this
    // applies to every assembler's output. Inert unless --run_nuclearphaser.
    //
    ch_np_haplotypes = Channel.empty()
    if (params.run_nuclearphaser) {
        PHASING(ASSEMBLY_QC.out.cleaned, ASSEMBLY_QC.out.busco_tables, ASSEMBLY_QC.out.gfastats)
        ch_versions = ch_versions.mix(PHASING.out.versions)
        ch_np_haplotypes = PHASING.out.haplotypes
    }

    //
    // STAGE: Hi-C scaffolding of the phased haplotypes. Contigs -> chromosomes; no new
    // sequence, only order and orientation. Independent of --run_nuclearphaser, because
    // hifiasm's --h1/--h2 haplotypes are already phased. Inert unless --run_scaffolding.
    //
    ch_scaffold_mqc = Channel.empty()
    if (params.run_scaffolding) {
        SCAFFOLDING(ASSEMBLY_QC.out.cleaned, ch_np_haplotypes)
        ch_versions = ch_versions.mix(SCAFFOLDING.out.versions)
        ch_scaffold_mqc = SCAFFOLDING.out.mqc
    }

    //
    // STAGE: reporting. Machine-readable summaries for automated comparison.
    //
    // Read-level contamination and the k-mer genome model are included here deliberately:
    // both were previously computed, written to disk, and then absent from the report, which
    // made the two things every size judgement is scored against the two things a reader had
    // to go hunting for. Kraken2's report needs no conversion — MultiQC parses it natively.
    //
    ch_multiqc = READ_QC.out.for_multiqc
        .mix(KMER_ANALYSIS.out.for_multiqc)
        .mix(CONTAMINATION.out.for_multiqc)
        .mix(ASSEMBLY_QC.out.for_multiqc)
        .mix(ch_scaffold_mqc)

    REPORTING(
        ASSEMBLY_QC.out.for_summary,
        ASSEMBLY_QC.out.depth_hists,
        ASSEMBLY_QC.out.expected_depths,
        ch_multiqc.collect(),
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}

workflow.onComplete {
    log.info """
    Pipeline ${workflow.success ? 'completed' : 'FAILED'}
      Duration : ${workflow.duration}
      Results  : ${params.outdir}
      Work dir : ${workflow.workDir}

    ${workflow.success ? "Check assembly_summary.tsv.\n    Judge candidates on size_flag + BUSCO breakdown + QV + coverage modes — never on N50." : ""}
    """.stripIndent()
}
