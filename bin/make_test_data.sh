#!/usr/bin/env bash
#
# Build the smoke-test inputs: a small head-subsample of each real run.
#
# Reads are only ever READ from the raw location (raw FASTQs are immutable); subsamples are
# written to a separate test_data directory.
#
# Usage:
#   RAW_DIR=/path/to/fastqs bash bin/make_test_data.sh [n_reads]
#
# Environment:
#   RAW_DIR    directory holding the raw FASTQs                    (required)
#   FASTQS     space-separated FASTQ names or globs within RAW_DIR [*.fastq.gz]
#   TEST_DIR   where subsamples are written                        [$RAW_DIR/test_data]
#   SAMPLE     sample name written into the samplesheet            [rust]

set -euo pipefail

# 50k reads/run turned out to be ~2.2 Gb of sequence and made Flye the long pole of the
# smoke test (~1 h). 15k keeps the whole test to a sensible length while still exercising
# every stage on genuine data.
N_READS="${1:-15000}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLE="${SAMPLE:-rust}"

if [[ -z "${RAW_DIR:-}" ]]; then
    echo "ERROR: set RAW_DIR to the directory holding your raw HiFi FASTQs" >&2
    echo "  e.g. RAW_DIR=\$SCRATCH/rust bash bin/make_test_data.sh" >&2
    exit 1
fi
TEST_DIR="${TEST_DIR:-${RAW_DIR}/test_data}"

# shellcheck disable=SC2206  # word splitting is the point: FASTQS may be a glob list
if [[ -n "${FASTQS:-}" ]]; then
    read -r -a NAMES <<< "${FASTQS}"
    SRCS=()
    for n in "${NAMES[@]}"; do SRCS+=( ${RAW_DIR}/${n} ); done
else
    SRCS=( "${RAW_DIR}"/*.fastq.gz )
fi

if [[ ! -e "${SRCS[0]}" ]]; then
    echo "ERROR: no FASTQs matched in ${RAW_DIR}" >&2
    exit 1
fi

mkdir -p "${TEST_DIR}"

# seqkit is not installed cluster-wide; use the container the pipeline itself uses.
SEQKIT_IMG='https://depot.galaxyproject.org/singularity/seqkit:2.8.2--h9ee0642_0'

# Cache images off $HOME where the cluster caps it — $HOME must stay small in both size and
# file count or SLURM throttles job concurrency.
if [[ -z "${APPTAINER_CACHEDIR:-}" ]]; then
    if [[ -d "/mnt/apps/users/${USER}" ]]; then
        APPTAINER_CACHEDIR="/mnt/apps/users/${USER}/apptainer_cache"
    else
        APPTAINER_CACHEDIR="${SCRATCH:-${HOME}}/apptainer_cache"
    fi
fi
export APPTAINER_CACHEDIR
mkdir -p "${APPTAINER_CACHEDIR}"

# Bind the filesystems the reads and outputs actually live on.
BINDS=()
for d in "${RAW_DIR}" "${TEST_DIR}"; do
    BINDS+=( -B "$(cd "${d}" && pwd)" )
done

echo "Building test subsamples (${N_READS} reads each) from ${RAW_DIR}..."

SHEET="${REPO_DIR}/assets/samplesheet_test.csv"
echo "sample,run,fastq" > "${SHEET}"

i=0
for src in "${SRCS[@]}"; do
    i=$((i + 1))
    run="run${i}"
    dest="${TEST_DIR}/test_${run}.fastq.gz"
    if [[ -s "${dest}" ]]; then
        echo "  exists, skipping: ${dest}"
    else
        echo "  $(basename "${src}") -> ${dest} (${N_READS} reads)"
        apptainer exec "${BINDS[@]}" "${SEQKIT_IMG}" \
            seqkit head -n "${N_READS}" "${src}" -o "${dest}"
    fi
    echo "${SAMPLE},${run},${dest}" >> "${SHEET}"
done

echo "Wrote ${SHEET}"
echo "Now run: nextflow run . -profile test,<site>"
echo
echo "NOTE: head-subsampling is not quality-representative. The first n reads of a PacBio"
echo "run can look far better than the run as a whole, so QC observations from the smoke"
echo "test say nothing about the real data. Use 'seqkit sample' if that matters."
