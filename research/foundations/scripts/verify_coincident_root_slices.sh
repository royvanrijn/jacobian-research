#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "$script_dir/.." && pwd)"
check_script="$script_dir/verify_coincident_root_slices.m2"
boundary_script="$script_dir/verify_cancellation_boundary_intersections.m2"

if command -v M2 >/dev/null 2>&1; then
    M2 --script "$check_script"
    exec M2 --script "$boundary_script"
fi

if ! command -v docker >/dev/null 2>&1; then
    echo "Macaulay2 is unavailable: install M2 or Docker." >&2
    exit 1
fi

image_name="jacobian-macaulay2:1.24.11"
docker build --tag "$image_name" "$repo_root/scripts/macaulay2"
docker run --rm --volume "$repo_root:/work:ro" \
    "$image_name" /work/scripts/verify_coincident_root_slices.m2
exec docker run --rm --volume "$repo_root:/work:ro" \
    "$image_name" /work/scripts/verify_cancellation_boundary_intersections.m2
