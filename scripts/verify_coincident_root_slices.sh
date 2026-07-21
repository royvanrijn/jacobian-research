#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
check_script="$repo_root/scripts/verify_coincident_root_slices.m2"

if command -v M2 >/dev/null 2>&1; then
    exec M2 --script "$check_script"
fi

if ! command -v docker >/dev/null 2>&1; then
    echo "Macaulay2 is unavailable: install M2 or Docker." >&2
    exit 1
fi

image_name="jacobian-macaulay2:1.24.11"
docker build --tag "$image_name" "$repo_root/scripts/macaulay2"
exec docker run --rm --volume "$repo_root:/work:ro" \
    "$image_name" /work/scripts/verify_coincident_root_slices.m2

