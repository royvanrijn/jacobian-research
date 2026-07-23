#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"

cd "${repo_dir}"
.venv/bin/python scripts/verify_degree_six_ritt_boundary_atlas.py
Singular -q scripts/verify_degree_six_ritt_boundary_atlas.sing

