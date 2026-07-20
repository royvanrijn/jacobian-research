#!/usr/bin/env bash
set -euo pipefail

# External formal certificate for Claim C01. The source remains in Dean
# Cureton's repository and is fetched at an immutable audited revision.
readonly REPOSITORY_URL="https://github.com/deancureton/jacobian.git"
readonly REVISION="0d4a9212d874226ad81ce5a926becddfa94e6a88"
readonly CHECKOUT_DIR="${LEAN_C01_DIR:-.cache/lean-c01}"

if ! command -v git >/dev/null 2>&1; then
  echo "ERROR: git is required" >&2
  exit 1
fi

if ! command -v lake >/dev/null 2>&1; then
  echo "ERROR: lake is required (install Lean via elan)" >&2
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

(
  cd "${CHECKOUT_DIR}"
  lake exe cache get
  lake build
)

actual_revision="$(git -C "${CHECKOUT_DIR}" rev-parse HEAD)"
if [[ "${actual_revision}" != "${REVISION}" ]]; then
  echo "ERROR: expected ${REVISION}, built ${actual_revision}" >&2
  exit 1
fi

echo "PASS: Dean Cureton's Lean C01 certificate built at ${REVISION}"
