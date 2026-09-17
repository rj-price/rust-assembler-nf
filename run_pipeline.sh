#!/usr/bin/env bash
#SBATCH --job-name=rust_assembly
#SBATCH --partition=long
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=14-00:00:00
#SBATCH --output=logs/nextflow_%j.out
#SBATCH --error=logs/nextflow_%j.err
#
# Nextflow driver for rust-assembler-nf.
#
# The driver itself must run on a compute node — never a head node. It is small (2 CPU,
# 8 GB) because it only orchestrates; the real work is submitted as separate SLURM jobs. It
# needs a long walltime because it must outlive the assemblies it is waiting on.
#
# Everything site- or dataset-specific is an environment variable with a default, so this
# script does not need editing to be reused:
#
#   NF_PROFILE   comma-separated Nextflow profiles      [gruffalo]
#   NF_INPUT     samplesheet CSV (repo-relative or abs) [assets/samplesheet.csv]
#   NF_OUTDIR    results directory                      [$SCRATCH/rust-assembler-nf/results]
#   NF_WORK      Nextflow work directory                [from the site profile]
#   NF_PIPELINE  what to run: this clone, or a GitHub    [this repo]
#                project + revision, e.g. "rj-price/rust-assembler-nf -r v1.0.0"
#   NF_CONDA_ENV conda env holding nextflow, used only   [nextflow]
#                when nextflow is not already on PATH
#   NF_ARGS      extra params, e.g. "--genome_size 525m"
#
# Usage:
#   sbatch run_pipeline.sh --genome_size 525m
#   NF_PROFILE=mlp,gruffalo sbatch run_pipeline.sh
#   sbatch run_pipeline.sh -resume
#
# Smoke test (after bin/make_test_data.sh):
#   NF_PROFILE=test,gruffalo sbatch run_pipeline.sh

set -euo pipefail

# NOTE: sbatch executes a COPY of this script from the SLURM spool directory, so
# ${BASH_SOURCE} does NOT point at the repo when running under SLURM. SLURM_SUBMIT_DIR is
# the directory sbatch was invoked from, which is the repo; fall back to BASH_SOURCE only
# when running the script directly.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO_DIR="${SLURM_SUBMIT_DIR}"
else
    REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

if [[ ! -f "${REPO_DIR}/main.nf" ]]; then
    echo "ERROR: main.nf not found in ${REPO_DIR} — submit with 'sbatch' from the repo root" >&2
    exit 1
fi

cd "${REPO_DIR}"
mkdir -p logs

NF_PROFILE="${NF_PROFILE:-gruffalo}"
NF_INPUT="${NF_INPUT:-assets/samplesheet.csv}"

# NF_INPUT may be repo-relative (the default) or an absolute path to a samplesheet kept with
# the data rather than with the code. Resolve it here, once: the run line used to prefix
# REPO_DIR unconditionally, which turned an absolute path into ${REPO_DIR}//mnt/... and made
# the pipeline fail on a samplesheet that was perfectly valid.
if [[ "${NF_INPUT}" != /* ]]; then
    NF_INPUT="${REPO_DIR}/${NF_INPUT#${REPO_DIR}/}"
fi
if [[ ! -f "${NF_INPUT}" ]]; then
    echo "ERROR: samplesheet not found: ${NF_INPUT}" >&2
    exit 1
fi

# Results must not land in $HOME on a cluster that caps it: assemblies and BAMs are large
# and numerous, and $HOME file-count limits throttle job concurrency. Prefer $SCRATCH.
if [[ -z "${NF_OUTDIR:-}" ]]; then
    if [[ -n "${SCRATCH:-}" ]]; then
        NF_OUTDIR="${SCRATCH%/}/rust-assembler-nf/results"
    else
        echo "ERROR: set NF_OUTDIR (or \$SCRATCH) — refusing to guess a results location" >&2
        exit 1
    fi
fi

# Apptainer's image cache. Site profiles set params.apptainer_cache for the pipeline itself;
# this covers anything Apptainer does outside Nextflow's control.
if [[ -z "${NXF_APPTAINER_CACHEDIR:-}" ]]; then
    if [[ -d "/mnt/apps/users/${USER}" ]]; then
        NXF_APPTAINER_CACHEDIR="/mnt/apps/users/${USER}/apptainer_cache"
    else
        NXF_APPTAINER_CACHEDIR="${SCRATCH:-${HOME}}/apptainer_cache"
    fi
fi
export NXF_APPTAINER_CACHEDIR
export APPTAINER_CACHEDIR="${NXF_APPTAINER_CACHEDIR}"
mkdir -p "${NXF_APPTAINER_CACHEDIR}"

[[ -n "${NF_WORK:-}" ]] && export NXF_WORK="${NF_WORK}"

# Keep the driver's own JVM small; tasks get their memory from SLURM, not from here.
export NXF_OPTS='-Xms1g -Xmx4g'

echo "Starting pipeline on $(hostname) at $(date)"
echo "  Pipeline: ${NF_PIPELINE:-${REPO_DIR}}"
echo "  Profile : ${NF_PROFILE}"
echo "  Input   : ${NF_INPUT}"
echo "  Outdir  : ${NF_OUTDIR}"
echo "  Work    : ${NXF_WORK:-<from profile>}"

# Use nextflow from PATH when there is one; otherwise from a conda env. --no-capture-output
# matters: without it conda buffers stdout until exit, leaving the log empty for days.
if command -v nextflow >/dev/null 2>&1; then
    NXF=(nextflow)
else
    NXF=(conda run --no-capture-output -n "${NF_CONDA_ENV:-nextflow}" nextflow)
fi

# shellcheck disable=SC2086  # NF_PIPELINE may carry "-r <tag>"
"${NXF[@]}" run ${NF_PIPELINE:-"${REPO_DIR}"} \
    --input "${NF_INPUT}" \
    --outdir "${NF_OUTDIR}" \
    -profile "${NF_PROFILE}" \
    -with-report \
    -with-trace \
    ${NF_ARGS:-} \
    "$@"

echo "Pipeline finished at $(date)"
echo
echo "Next steps:"
echo "  1. seff on the assembly jobs, then correct conf/base.config's resource table"
echo "  2. Review assembly_summary.tsv — judge on size_flag, BUSCO breakdown, QV and"
echo "     coverage modes. Never on N50."
