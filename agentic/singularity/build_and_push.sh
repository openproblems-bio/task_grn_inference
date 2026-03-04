#!/usr/bin/env bash
# Build the agentic Singularity image and push it to S3.
#
# Usage:
#   bash build_and_push.sh              # build + push
#   bash build_and_push.sh --build-only # build only, skip S3 upload
#
# Prerequisites:
#   - apptainer (or singularity) installed
#   - AWS CLI configured with write access to s3://openproblems-data/

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE_NAME="agentic.sif"
DEF_FILE="${SCRIPT_DIR}/agentic.def"
OUTPUT_IMAGE="${SCRIPT_DIR}/${IMAGE_NAME}"
S3_DEST="s3://openproblems-data/resources/grn/images/${IMAGE_NAME}"

BUILD_ONLY=false
for arg in "$@"; do
    [[ "$arg" == "--build-only" ]] && BUILD_ONLY=true
done

# ─── Build ────────────────────────────────────────────────────────────────────
echo "Building ${IMAGE_NAME} from ${DEF_FILE} ..."
apptainer build --fakeroot "${OUTPUT_IMAGE}" "${DEF_FILE}"
echo "Build complete: ${OUTPUT_IMAGE}"

# ─── Push ─────────────────────────────────────────────────────────────────────
if [[ "${BUILD_ONLY}" == "false" ]]; then
    echo "Uploading ${IMAGE_NAME} to ${S3_DEST} ..."
    aws s3 cp "${OUTPUT_IMAGE}" "${S3_DEST}"
    echo "Upload complete."
    echo
    echo "Public download command (no credentials required):"
    echo "  aws s3 cp ${S3_DEST} . --no-sign-request"
fi
