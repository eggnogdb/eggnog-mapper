#!/usr/bin/env bash
#
# release-db.sh — publish a versioned eggnog-mapper DB set with manifest+checksum.
#
# Steps: validate source-dir → compress → sha256 + manifest.json + SHA256SUMS
#        → upload (rsync or local) → re-fetch manifest from public URL and verify.
#
# Usage:
#   scripts/release-db.sh \
#     --db-version 5.0.2 \
#     --source-dir /path/to/data/e7/full/final/mapper/ \
#     --target rsync://eggdb@downloads.example.org:/srv/emapper/emapperdb-5.0.2/ \
#     --public-base-url https://downloads.eggnogdb.org/emapper/emapperdb-5.0.2/
#
# Targets supported:
#   rsync://user@host:/abs/path/         — rsync-over-ssh upload
#   local:/abs/path/                     — copy to a local dir (staging)
#
set -euo pipefail

DB_VERSION=""
SOURCE_DIR=""
TARGET=""
PUBLIC_BASE_URL=""
SCHEMA_VERSION="v7"
MIN_MAPPER_VERSION=""
DRY_RUN=0

err()  { echo "ERROR: $*" >&2; exit 1; }
note() { echo "==> $*"; }
maybe(){ if [[ $DRY_RUN -eq 1 ]]; then echo "DRY-RUN: $*"; else "$@"; fi; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --db-version)          DB_VERSION="$2"; shift 2 ;;
    --source-dir)          SOURCE_DIR="$2"; shift 2 ;;
    --target)              TARGET="$2"; shift 2 ;;
    --public-base-url)     PUBLIC_BASE_URL="$2"; shift 2 ;;
    --schema-version)      SCHEMA_VERSION="$2"; shift 2 ;;
    --min-mapper-version)  MIN_MAPPER_VERSION="$2"; shift 2 ;;
    --dry-run)             DRY_RUN=1; shift ;;
    -h|--help)
      sed -n '2,17p' "$0"; exit 0 ;;
    *) err "unknown flag: $1" ;;
  esac
done

[[ -z "$DB_VERSION" ]]       && err "--db-version required"
[[ -z "$SOURCE_DIR" ]]       && err "--source-dir required"
[[ -z "$TARGET" ]]           && err "--target required"
[[ -z "$PUBLIC_BASE_URL" ]]  && err "--public-base-url required (used for post-upload verify)"
[[ -d "$SOURCE_DIR" ]]       || err "--source-dir not a directory: $SOURCE_DIR"

# default min_mapper_version from this checkout's version.py
if [[ -z "$MIN_MAPPER_VERSION" ]]; then
  REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
  MIN_MAPPER_VERSION="$(grep -oP "__VERSION__\s*=\s*'\K[^']+" "$REPO_ROOT/eggnogmapper/version.py" 2>/dev/null || echo "3.0.0")"
fi

note "db-version=$DB_VERSION  schema=$SCHEMA_VERSION  min-mapper=$MIN_MAPPER_VERSION"
note "source-dir=$SOURCE_DIR"
note "target=$TARGET  dry-run=$DRY_RUN"

cd "$SOURCE_DIR"

# ---------- pre-flight: source files present ----------
REQUIRED=(eggnog.db eggnog.taxa.db eggnog_proteins.dmnd)
for f in "${REQUIRED[@]}"; do
  [[ -f "$f" ]] || err "missing required source file: $SOURCE_DIR/$f"
done

# Optional: traverse pkl ships alongside taxa.db
TRAVERSE_PKL=""
[[ -f "eggnog.taxa.db.traverse.pkl" ]] && TRAVERSE_PKL="eggnog.taxa.db.traverse.pkl"

# ---------- compress (idempotent — keep originals) ----------
note "compress: eggnog.db.gz"
[[ -f eggnog.db.gz && eggnog.db.gz -nt eggnog.db ]] || maybe gzip -k -f eggnog.db

note "compress: eggnog_proteins.dmnd.gz"
[[ -f eggnog_proteins.dmnd.gz && eggnog_proteins.dmnd.gz -nt eggnog_proteins.dmnd ]] || maybe gzip -k -f eggnog_proteins.dmnd

note "compress: eggnog.taxa.tar.gz"
TAXA_FILES=(eggnog.taxa.db)
[[ -n "$TRAVERSE_PKL" ]] && TAXA_FILES+=("$TRAVERSE_PKL")
maybe tar -czf eggnog.taxa.tar.gz "${TAXA_FILES[@]}"

# ---------- artifact list ----------
ARTIFACTS=(eggnog.db.gz eggnog.taxa.tar.gz eggnog_proteins.dmnd.gz)

# ---------- sha256 + sizes ----------
note "hash + size for ${#ARTIFACTS[@]} artifacts"
declare -A SHA SIZE
for a in "${ARTIFACTS[@]}"; do
  if [[ $DRY_RUN -eq 1 ]]; then
    SHA[$a]="<dry-run>"
    SIZE[$a]=0
  else
    [[ -f "$a" ]] || err "compressed artifact not produced: $a"
    SHA[$a]="$(sha256sum "$a" | awk '{print $1}')"
    SIZE[$a]="$(stat -c%s "$a")"
  fi
  echo "  $a  ${SIZE[$a]}  ${SHA[$a]}"
done

# ---------- write manifest.json ----------
note "write manifest.json"
BUILD_DATE="$(date -u +%Y-%m-%d)"
MANIFEST="manifest.json"
if [[ $DRY_RUN -eq 1 ]]; then
  echo "DRY-RUN: would write $MANIFEST"
else
  {
    echo '{'
    echo "  \"db_version\": \"$DB_VERSION\","
    echo "  \"schema_version\": \"$SCHEMA_VERSION\","
    echo "  \"min_mapper_version\": \"$MIN_MAPPER_VERSION\","
    echo "  \"build_date\": \"$BUILD_DATE\","
    echo '  "artifacts": {'
    for i in "${!ARTIFACTS[@]}"; do
      a="${ARTIFACTS[$i]}"
      printf '    "%s": {"size": %s, "sha256": "%s"}' "$a" "${SIZE[$a]}" "${SHA[$a]}"
      [[ $i -lt $((${#ARTIFACTS[@]}-1)) ]] && echo "," || echo
    done
    echo '  }'
    echo '}'
  } > "$MANIFEST"
fi

# ---------- write SHA256SUMS ----------
note "write SHA256SUMS"
if [[ $DRY_RUN -eq 1 ]]; then
  echo "DRY-RUN: would write SHA256SUMS"
else
  : > SHA256SUMS
  for a in "${ARTIFACTS[@]}"; do
    echo "${SHA[$a]}  $a" >> SHA256SUMS
  done
fi

# ---------- upload ----------
UPLOAD_FILES=("${ARTIFACTS[@]}" "$MANIFEST" "SHA256SUMS")

if [[ "$TARGET" == rsync://* ]]; then
  RSYNC_TARGET="${TARGET#rsync://}"
  note "upload: rsync → $RSYNC_TARGET"
  maybe rsync -av --progress "${UPLOAD_FILES[@]}" "$RSYNC_TARGET"
elif [[ "$TARGET" == local:* ]]; then
  LOCAL_TARGET="${TARGET#local:}"
  note "upload: local cp → $LOCAL_TARGET"
  maybe mkdir -p "$LOCAL_TARGET"
  maybe cp -v "${UPLOAD_FILES[@]}" "$LOCAL_TARGET/"
else
  err "unsupported --target backend: $TARGET (rsync:// or local: only)"
fi

# ---------- verify upload via public URL ----------
if [[ $DRY_RUN -eq 1 ]]; then
  note "skip verify: dry-run"
else
  note "verify: re-fetch manifest from $PUBLIC_BASE_URL"
  REMOTE_MANIFEST="$(mktemp)"
  if curl -fsSL "${PUBLIC_BASE_URL%/}/manifest.json" -o "$REMOTE_MANIFEST"; then
    if cmp -s "$MANIFEST" "$REMOTE_MANIFEST"; then
      note "verify OK: public manifest matches local"
    else
      err "public manifest differs from local — upload incomplete or stale CDN cache"
    fi
  else
    err "could not fetch ${PUBLIC_BASE_URL%/}/manifest.json — upload not visible yet"
  fi
  rm -f "$REMOTE_MANIFEST"
fi

note "done."
echo
[[ $DRY_RUN -eq 0 ]] && echo "Public manifest: ${PUBLIC_BASE_URL%/}/manifest.json"
