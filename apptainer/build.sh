#!/usr/bin/env bash
# Build the eggnog-mapper Apptainer image.
#
# Usage:
#   ./build.sh [--local] [EMAPPER_VERSION] [DIAMOND_VERSION]
#
#   --local   Bundle the working directory as-is (includes uncommitted changes).
#             Useful for testing dev builds before committing.
#             Output is named eggnog-mapper-VERSION-dev.sif.
#
#   Default (no --local): bundle from git HEAD (committed state only).
#             Use this for release images.
#             Output is named eggnog-mapper-VERSION.sif.
#
# EMAPPER_VERSION defaults to the version in eggnogmapper/version.py.
# DIAMOND_VERSION defaults to 2.1.24.
#
# Requires: apptainer >= 1.1, git

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

# ── parse flags ──────────────────────────────────────────────────────────────
LOCAL=0
POSITIONAL=()
for arg in "$@"; do
    case "$arg" in
        --local) LOCAL=1 ;;
        *)       POSITIONAL+=("$arg") ;;
    esac
done
set -- "${POSITIONAL[@]+"${POSITIONAL[@]}"}"

EMAPPER_VERSION="${1:-$(grep "__VERSION__" "$REPO_ROOT/eggnogmapper/version.py" | cut -d"'" -f2)}"
DIAMOND_VERSION="${2:-2.1.24}"
SRC_ARCHIVE="/tmp/emapper-src.tar.gz"

if [[ "$LOCAL" == "1" ]]; then
    OUT="${SCRIPT_DIR}/eggnog-mapper-${EMAPPER_VERSION}-dev.sif"
    SRC_LABEL="working directory (uncommitted changes included)"
else
    OUT="${SCRIPT_DIR}/eggnog-mapper-${EMAPPER_VERSION}.sif"
    SRC_LABEL="git HEAD: $(git -C "$REPO_ROOT" log -1 --oneline)"
fi

echo "=== eggnog-mapper Apptainer build ==="
echo "  emapper : ${EMAPPER_VERSION}"
echo "  diamond : ${DIAMOND_VERSION}"
echo "  source  : ${SRC_LABEL}"
echo "  output  : ${OUT}"

# ── create source archive ────────────────────────────────────────────────────
if [[ "$LOCAL" == "1" ]]; then
    tar -czf "$SRC_ARCHIVE" \
        --exclude='.git' \
        --exclude='__pycache__' \
        --exclude='*.pyc' \
        --exclude='*.pyo' \
        --exclude='build' \
        --exclude='dist' \
        --exclude='*.egg-info' \
        --exclude='apptainer' \
        --exclude='*.sif' \
        --exclude='.v3-venv' \
        --exclude='emappertmp_*' \
        -C "$REPO_ROOT" .
else
    git -C "$REPO_ROOT" archive --format=tar.gz HEAD > "$SRC_ARCHIVE"
fi

echo "  archive : ${SRC_ARCHIVE} ($(du -sh "$SRC_ARCHIVE" | cut -f1))"
echo ""

# ── build ────────────────────────────────────────────────────────────────────
apptainer build \
    --build-arg EMAPPER_VERSION="${EMAPPER_VERSION}" \
    --build-arg DIAMOND_VERSION="${DIAMOND_VERSION}" \
    "${OUT}" \
    "${SCRIPT_DIR}/eggnog-mapper.def"

rm -f "$SRC_ARCHIVE"

echo ""
echo "Done: ${OUT}"
echo ""
echo "Test:  apptainer test ${OUT}"
echo "Run:   apptainer run ${OUT} --help"
echo ""
echo "Example with bind mounts:"
echo "  apptainer run \\"
echo "    --bind /path/to/eggnog-data:/data \\"
echo "    --bind \$(pwd):/work \\"
echo "    ${OUT} \\"
echo "    -i /work/query.fa --itype proteins -m diamond \\"
echo "    --data_dir /data --output my_run --output_dir /work --cpu 20"
