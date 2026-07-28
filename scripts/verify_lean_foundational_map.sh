#!/usr/bin/env bash
set -euo pipefail

# External formal certificate for the foundational Keller map. The source remains in Dean
# Cureton's repository and is fetched at an immutable audited revision. The audited source
# was released against a Lean release candidate, so this wrapper builds an exact worktree
# with the checked stable-toolchain overlay in verified/lean-foundational-v4.32.1.patch.
readonly REPOSITORY_URL="https://github.com/deancureton/jacobian.git"
readonly REVISION="0d4a9212d874226ad81ce5a926becddfa94e6a88"
readonly STABLE_RELEASE="v4.32.1"
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly REPOSITORY_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
readonly CHECKOUT_DIR="${LEAN_FOUNDATIONAL_DIR:-${REPOSITORY_ROOT}/.cache/lean-foundational}"
readonly BUILD_DIR="${LEAN_FOUNDATIONAL_BUILD_DIR:-${REPOSITORY_ROOT}/.cache/lean-foundational-stable}"
readonly OVERLAY_PATCH="${REPOSITORY_ROOT}/verified/lean-foundational-v4.32.1.patch"
readonly EXPECTED_CHANGED_PATHS=$'lake-manifest.json\nlakefile.toml\nlean-toolchain'
readonly EXPECTED_MANIFEST_BLOB="d8761427a061bd102867c6ab422e43b8a7442870"
readonly EXPECTED_LAKEFILE_BLOB="9c733c4355ee8655035f22f35a44715de6429c52"
readonly EXPECTED_TOOLCHAIN_BLOB="2b61957be5f4986fcc464cb4f62322f2d6ca3dfc"

prepare_only=false
if [[ "${1:-}" == "--prepare-only" ]]; then
  prepare_only=true
elif [[ $# -ne 0 ]]; then
  echo "usage: $0 [--prepare-only]" >&2
  exit 2
fi

if ! command -v git >/dev/null 2>&1; then
  echo "ERROR: git is required" >&2
  exit 1
fi

if [[ ! -f "${OVERLAY_PATCH}" ]]; then
  echo "ERROR: missing stable-toolchain overlay ${OVERLAY_PATCH}" >&2
  exit 1
fi

if [[ ! -d "${CHECKOUT_DIR}/.git" ]]; then
  if [[ -e "${CHECKOUT_DIR}" ]]; then
    echo "ERROR: ${CHECKOUT_DIR} exists but is not a git checkout" >&2
    exit 1
  fi
  mkdir -p "$(dirname "${CHECKOUT_DIR}")"
  git clone "${REPOSITORY_URL}" "${CHECKOUT_DIR}"
fi

actual_url="$(git -C "${CHECKOUT_DIR}" remote get-url origin)"
if [[ "${actual_url%.git}" != "${REPOSITORY_URL%.git}" ]]; then
  echo "ERROR: ${CHECKOUT_DIR} has unexpected origin ${actual_url}" >&2
  exit 1
fi

if [[ -n "$(git -C "${CHECKOUT_DIR}" status --porcelain)" ]]; then
  echo "ERROR: refusing to alter dirty external checkout ${CHECKOUT_DIR}" >&2
  exit 1
fi

if ! git -C "${CHECKOUT_DIR}" cat-file -e "${REVISION}^{commit}" 2>/dev/null; then
  git -C "${CHECKOUT_DIR}" fetch --depth 1 origin "${REVISION}"
fi

git -C "${CHECKOUT_DIR}" checkout --detach "${REVISION}"

if [[ -e "${BUILD_DIR}" ]]; then
  if [[ ! -d "${BUILD_DIR}" ]] ||
      ! git -C "${BUILD_DIR}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "ERROR: ${BUILD_DIR} exists but is not the expected generated worktree" >&2
    exit 1
  fi
  build_url="$(git -C "${BUILD_DIR}" remote get-url origin)"
  if [[ "${build_url%.git}" != "${REPOSITORY_URL%.git}" ]]; then
    echo "ERROR: ${BUILD_DIR} has unexpected origin ${build_url}" >&2
    exit 1
  fi
  build_revision="$(git -C "${BUILD_DIR}" rev-parse HEAD)"
  if [[ "${build_revision}" != "${REVISION}" ]]; then
    echo "ERROR: ${BUILD_DIR} is at unexpected revision ${build_revision}" >&2
    exit 1
  fi
else
  mkdir -p "$(dirname "${BUILD_DIR}")"
  git -C "${CHECKOUT_DIR}" worktree prune
  git -C "${CHECKOUT_DIR}" worktree add --detach "${BUILD_DIR}" "${REVISION}"
fi

changed_paths="$(git -C "${BUILD_DIR}" diff --name-only)"
if [[ -z "${changed_paths}" ]]; then
  git -C "${BUILD_DIR}" apply "${OVERLAY_PATCH}"
  changed_paths="$(git -C "${BUILD_DIR}" diff --name-only)"
fi

untracked_paths="$(git -C "${BUILD_DIR}" ls-files --others --exclude-standard)"
if [[ -n "${untracked_paths}" ]]; then
  echo "ERROR: ${BUILD_DIR} contains unexpected untracked files" >&2
  printf '%s\n' "${untracked_paths}" >&2
  exit 1
fi

if [[ "${changed_paths}" != "${EXPECTED_CHANGED_PATHS}" ]]; then
  echo "ERROR: ${BUILD_DIR} has changes outside the checked stable overlay" >&2
  git -C "${BUILD_DIR}" status --short >&2
  exit 1
fi

if [[ "$(git -C "${BUILD_DIR}" hash-object lake-manifest.json)" != "${EXPECTED_MANIFEST_BLOB}" ]] ||
    [[ "$(git -C "${BUILD_DIR}" hash-object lakefile.toml)" != "${EXPECTED_LAKEFILE_BLOB}" ]] ||
    [[ "$(git -C "${BUILD_DIR}" hash-object lean-toolchain)" != "${EXPECTED_TOOLCHAIN_BLOB}" ]]; then
  echo "ERROR: ${BUILD_DIR} does not match the checked ${STABLE_RELEASE} overlay" >&2
  exit 1
fi

git -C "${BUILD_DIR}" diff --check

if [[ "${prepare_only}" == true ]]; then
  echo "PREPARED: external certificate ${REVISION} with stable ${STABLE_RELEASE} overlay"
  exit 0
fi

if ! command -v lake >/dev/null 2>&1; then
  echo "ERROR: lake is required (install Lean via elan)" >&2
  exit 1
fi

(
  cd "${BUILD_DIR}"
  lake exe cache get
  lake build
)

actual_revision="$(git -C "${BUILD_DIR}" rev-parse HEAD)"
if [[ "${actual_revision}" != "${REVISION}" ]]; then
  echo "ERROR: expected ${REVISION}, built ${actual_revision}" >&2
  exit 1
fi

echo "PASS: Dean Cureton's Lean foundational-map certificate built at ${REVISION}"
echo "PASS: build used stable Lean/Mathlib ${STABLE_RELEASE} via the checked overlay"
