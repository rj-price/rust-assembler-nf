process ASSEMBLY_RECORD {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(fasta), path(gfastats), path(busco), path(qv), path(coverage), path(telomeres), path(fcs_clean), path(mito)

    output:
    path "${meta.id}.record.json", emit: record
    path "versions.yml"          , emit: versions

    script:
    // Optional inputs arrive as null when a QC step was disabled or failed; pass only what
    // exists so a missing tool degrades the record rather than dropping the assembly.
    def meta_json = groovy.json.JsonOutput.toJson([
        id           : meta.id,
        sample       : meta.sample,
        assembler    : meta.assembler,
        assembly_type: meta.assembly_type,
        readset      : meta.readset
    ])
    // Placeholders are per-slot (assets/NO_FILE_<slot>), so match on the prefix.
    def real = { f -> f && !f.name.startsWith('NO_FILE') }
    // Bounds are derived from --genome_size at TASK time rather than stored as absolute
    // byte counts, so they can never go stale against a revised genome-size estimate.
    def dik_min = GenomeSize.dikaryonMin(params.genome_size, params.dikaryon_size_frac_min)
    def dik_max = GenomeSize.dikaryonMax(params.genome_size, params.dikaryon_size_frac_max)
    def gfastats_arg = real.call(gfastats) ? "--gfastats ${gfastats}" : ''
    def busco_arg    = real.call(busco)    ? "--busco ${busco}"       : ''
    def qv_arg       = real.call(qv)       ? "--qv ${qv}"             : ''
    def cov_arg      = real.call(coverage) ? "--coverage ${coverage}" : ''
    def telo_arg     = real.call(telomeres) ? "--telomeres ${telomeres}" : ''
    def clean_arg    = real.call(fcs_clean) ? "--fcs-clean ${fcs_clean}" : ''
    def mito_arg     = real.call(mito)      ? "--mito ${mito}"           : ''
    """
    assembly_summary.py record \\
        --meta '${meta_json}' \\
        ${gfastats_arg} \\
        ${busco_arg} \\
        ${qv_arg} \\
        ${cov_arg} \\
        ${telo_arg} \\
        ${clean_arg} \\
        ${mito_arg} \\
        --dikaryon-min ${dik_min} \\
        --dikaryon-max ${dik_max} \\
        --out ${meta.id}.record.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"assembly_id": "${meta.id}"}' > ${meta.id}.record.json
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
