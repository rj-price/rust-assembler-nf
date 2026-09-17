process SMUDGEPLOT_CUTOFF {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(hist)

    output:
    tuple val(meta), path("${meta.id}.cutoffs.txt"), emit: cutoffs
    path "versions.yml"                            , emit: versions

    when:
    params.run_smudgeplot

    script:
    // Step 1 of 3. Split out from SMUDGEPLOT only because of a container boundary: the
    // smudgeplot image has no meryl and the meryl image has no smudgeplot, so the middle step
    // physically cannot run in the same container as the two ends.
    //
    // Cheap (seconds, reads only the histogram), so the split costs nothing but a scheduler
    // round-trip, and it makes the cutoffs a durable artefact rather than a shell variable.
    """
    L=\$(smudgeplot.py cutoff ${hist} L)
    U=\$(smudgeplot.py cutoff ${hist} U)

    # Guard the handoff: an empty or non-numeric cutoff would otherwise reach `meryl print` as
    # an empty bracket expression, which fails obscurely a step later.
    for v in "\$L" "\$U"; do
        case "\$v" in
            ''|*[!0-9]*)
                echo "ERROR: smudgeplot cutoff returned a non-integer (L='\$L' U='\$U')." >&2
                exit 1
                ;;
        esac
    done

    printf '%s\\t%s\\n' "\$L" "\$U" > ${meta.id}.cutoffs.txt
    echo "smudgeplot coverage cutoffs: L=\$L U=\$U"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smudgeplot: \$(smudgeplot.py --version 2>&1 | sed 's/smudgeplot //' | head -1)
    END_VERSIONS
    """

    stub:
    """
    printf '10\\t500\\n' > ${meta.id}.cutoffs.txt
    echo '"${task.process}": {smudgeplot: stub}' > versions.yml
    """
}
