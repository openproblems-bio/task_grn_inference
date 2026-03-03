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

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_ENV="biomni_e1"

# Activate conda env if available
if command -v conda &>/dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate "$CONDA_ENV"
fi

cd "$SCRIPT_DIR/src"
python run.py "$@"
