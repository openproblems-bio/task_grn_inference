#!/usr/bin/env bash
# Run the GeneRNBI agentic assistant.
#
# Usage:
#   bash run.sh                          # interactive mode (reads .env or env vars)
#   bash run.sh "your question here"     # single-query mode
#
# Keys can be provided in three ways (first found wins):
#   1. Environment variable:  OPENAI_API_KEY=sk-... bash run.sh
#   2. Shell export:          export OPENAI_API_KEY=sk-...  then  bash run.sh
#   3. .env file:             cp .env.template .env  and fill in your key
#
# To use the Singularity image instead of the conda env, set:
#   USE_SINGULARITY=1 bash run.sh [question]
# and optionally override the image path:
#   SINGULARITY_IMAGE=/path/to/agentic.sif USE_SINGULARITY=1 bash run.sh [question]

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ─── Singularity mode ─────────────────────────────────────────────────────────
if [[ "${USE_SINGULARITY:-0}" == "1" ]]; then
    SIF="${SINGULARITY_IMAGE:-${SCRIPT_DIR}/singularity/agentic.sif}"
    if [[ ! -f "$SIF" ]]; then
        echo "ERROR: Singularity image not found at ${SIF}" >&2
        echo "  Build it with:  bash ${SCRIPT_DIR}/singularity/build_and_push.sh --build-only" >&2
        echo "  Or download it: aws s3 cp s3://openproblems-data/resources/grn/images/agentic.sif ${SIF} --no-sign-request" >&2
        exit 1
    fi
    exec singularity run \
        --bind "${SCRIPT_DIR}:/agentic" \
        --bind "${SCRIPT_DIR}/../docs:/docs" \
        --env SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt \
        "$SIF" "$@"
fi

# ─── Conda mode (default) ─────────────────────────────────────────────────────
CONDA_ENV="genernbi_agentic"

if command -v conda &>/dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate "$CONDA_ENV"
fi

cd "$SCRIPT_DIR/src"
python run.py "$@"
