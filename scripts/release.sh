#!/usr/bin/env bash
#
# release.sh — cut a new eggnog-mapper release.
#
# Steps: pre-flight checks → build sdist → smoke install in fresh venv
# → tag → push → gh release create (draft by default).
#
# Usage:
#   scripts/release.sh --version 3.0.0
#   scripts/release.sh --version 3.0.0 --dry-run
#   scripts/release.sh --version 3.0.0 --publish    # not draft
#   scripts/release.sh --version 3.0.0 --no-push    # local-only build+tag
#
# Prerequisites: clean working tree on main/release branch, gh CLI
# authenticated, version bumped in eggnogmapper/version.py, CHANGELOG.md
# has the "## [vX.Y.Z]" section, all tests passing.
#
set -euo pipefail

VERSION=""
DRY_RUN=0
NO_PUSH=0
DRAFT=1   # default ON per agreed plan; --publish flips

err()  { echo "ERROR: $*" >&2; exit 1; }
note() { echo "==> $*"; }
maybe(){ if [[ $DRY_RUN -eq 1 ]]; then echo "DRY-RUN: $*"; else "$@"; fi; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --version)   VERSION="$2"; shift 2 ;;
    --dry-run)   DRY_RUN=1; shift ;;
    --no-push)   NO_PUSH=1; shift ;;
    --publish)   DRAFT=0; shift ;;
    --draft)     DRAFT=1; shift ;;
    -h|--help)
      sed -n '2,15p' "$0"; exit 0 ;;
    *) err "unknown flag: $1" ;;
  esac
done

[[ -z "$VERSION" ]] && err "--version X.Y.Z is required"
[[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || err "version must be X.Y.Z (got '$VERSION')"

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

note "release: $VERSION  draft=$DRAFT  dry-run=$DRY_RUN  no-push=$NO_PUSH"

# ---------- pre-flight ----------
note "pre-flight: working tree clean"
[[ -z "$(git status --porcelain)" ]] || err "working tree has uncommitted/untracked files"

note "pre-flight: branch sanity"
BRANCH="$(git rev-parse --abbrev-ref HEAD)"
case "$BRANCH" in
  main|master|v3-dev|release/*) ;;
  *) err "must be on main / master / v3-dev / release/* (on '$BRANCH')" ;;
esac

note "pre-flight: gh CLI authenticated"
gh auth status >/dev/null 2>&1 || err "gh CLI not authenticated. Run 'gh auth login'."

note "pre-flight: version.py says $VERSION"
PY_VERSION="$(grep -oP "__VERSION__\s*=\s*'\K[^']+" eggnogmapper/version.py)"
[[ "$PY_VERSION" == "$VERSION" ]] || err "eggnogmapper/version.py says '$PY_VERSION', expected '$VERSION'"

note "pre-flight: CHANGELOG.md section [v$VERSION] exists"
grep -qE "^## \[v?$VERSION" CHANGELOG.md || err "no '## [v$VERSION]' section in CHANGELOG.md"

note "pre-flight: tag v$VERSION does not exist"
if git rev-parse "v$VERSION" >/dev/null 2>&1; then
  err "tag v$VERSION already exists locally"
fi
if git ls-remote --tags origin "refs/tags/v$VERSION" 2>/dev/null | grep -q .; then
  err "tag v$VERSION already exists on origin"
fi

note "pre-flight: pytest"
maybe python -m pytest tests/unit/ eggnogmapper/annotator/tests/ -q

# ---------- build ----------
note "build: sdist"
maybe rm -rf dist/
maybe python -m build --sdist

if [[ $DRY_RUN -eq 0 ]]; then
  SDIST="$(ls dist/eggnog?mapper-${VERSION}.tar.gz 2>/dev/null || ls dist/eggnog_mapper-${VERSION}.tar.gz)"
  [[ -n "$SDIST" && -f "$SDIST" ]] || err "sdist not produced under dist/"
  note "sdist: $SDIST  ($(du -h "$SDIST" | cut -f1))"
fi

# ---------- smoke install ----------
note "smoke: fresh venv install + import"
SMOKE_VENV="$(mktemp -d)/smoke-venv"
maybe python -m venv "$SMOKE_VENV"
if [[ $DRY_RUN -eq 0 ]]; then
  "$SMOKE_VENV/bin/pip" install --quiet --upgrade pip wheel
  "$SMOKE_VENV/bin/pip" install --quiet "$SDIST"
  "$SMOKE_VENV/bin/python" -c "
import eggnogmapper
from eggnogmapper.annotator import encode_intlist, decode_intlist
from eggnogmapper.annotator.e7 import AnnotationEngine, EggnogDB, LineageFilter
from eggnogmapper.annotator import _codec, _collect_inner
from eggnogmapper.version import __VERSION__
assert __VERSION__ == '$VERSION', f'version mismatch: {__VERSION__}'
print('smoke install OK; version', __VERSION__)
"
  rm -rf "$(dirname "$SMOKE_VENV")"
fi

# ---------- tag ----------
note "tag: v$VERSION"
maybe git tag -a "v$VERSION" -m "Release v$VERSION"

# ---------- push ----------
if [[ $NO_PUSH -eq 1 ]]; then
  note "skip: --no-push set; not pushing branch or tag"
else
  note "push: $BRANCH and v$VERSION"
  maybe git push origin "$BRANCH"
  maybe git push origin "v$VERSION"
fi

# ---------- gh release ----------
if [[ $NO_PUSH -eq 1 ]]; then
  note "skip: --no-push set; not creating gh release"
else
  note "gh release: extracting CHANGELOG section"
  NOTES_FILE="$(mktemp)"
  awk -v ver="$VERSION" '
    BEGIN { capturing = 0 }
    /^## \[v?/ {
      if (capturing) { exit }
      if ($0 ~ "## \\[v?" ver) { capturing = 1 }
    }
    capturing { print }
  ' CHANGELOG.md > "$NOTES_FILE"
  [[ -s "$NOTES_FILE" ]] || err "extracted release notes are empty"

  DRAFT_FLAG=""
  [[ $DRAFT -eq 1 ]] && DRAFT_FLAG="--draft"
  note "gh release create v$VERSION (draft=$DRAFT)"
  maybe gh release create "v$VERSION" "$SDIST" \
    --title "v$VERSION" \
    --notes-file "$NOTES_FILE" \
    $DRAFT_FLAG
  rm -f "$NOTES_FILE"
fi

note "done."
echo
if [[ $DRY_RUN -eq 1 ]]; then
  echo "DRY RUN — no changes made."
elif [[ $DRAFT -eq 1 ]] && [[ $NO_PUSH -eq 0 ]]; then
  echo "Release created as DRAFT. Review at:"
  gh release view "v$VERSION" --json url --jq .url 2>/dev/null || \
    echo "  gh release view v$VERSION"
fi
