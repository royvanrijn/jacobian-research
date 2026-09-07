#!/usr/bin/env bash
set -euo pipefail

REPO="${Q24_REPO:-/Users/royvanrijn/Documents/jacobian-research}"
PRIME="100003"
DO_RESOLVE=0
DO_EXACT=0
RR_ONLY=0

for arg in "$@"; do
  case "$arg" in
    --resolve) DO_RESOLVE=1 ;;
    --exact) DO_EXACT=1 ;;
    --rr-only) RR_ONLY=1 ;;
    '' ) ;;
    *[!0-9]*) echo "Unknown argument: $arg" >&2; exit 2 ;;
    *) PRIME="$arg" ;;
  esac
done

cd "$REPO"

MOD_ARGS=(--repo "$REPO" --prime "$PRIME")
if [[ "$RR_ONLY" == 1 ]]; then MOD_ARGS+=(--rr-only); fi

echo "[q24] reduced 16-variable modular probe, p=$PRIME"
sage -python elkies-k3/scripts/probe_h92_q24_d12_reduced16.sage "${MOD_ARGS[@]}"

if [[ "$DO_EXACT" == 1 ]]; then
  EXACT_ARGS=(--repo "$REPO" --exact)
  if [[ "$RR_ONLY" == 1 ]]; then EXACT_ARGS+=(--rr-only); fi
  echo "[q24] reduced exact QQ probe"
  sage -python elkies-k3/scripts/probe_h92_q24_d12_reduced16.sage "${EXACT_ARGS[@]}"
fi

RESOLUTION="artifacts/local/elkies-k3/q24-i9star-resolution-mod-$PRIME.json"
if [[ "$DO_RESOLVE" == 1 ]]; then
  echo "[q24] rebuilding explicit I9* chart tree"
  sage -python elkies-k3/scripts/derive_h92_q24_i9star_resolution_modp.sage \
    --repo "$REPO" --prime "$PRIME"
fi

if [[ -f "$RESOLUTION" ]]; then
  echo "[q24] pulling (u^8,xP*y+yP*x) through the chart atlas"
  sage -python elkies-k3/scripts/probe_h92_q24_i9star_ideal_atlas_modp.sage \
    --repo "$REPO" --prime "$PRIME"
else
  echo "[q24] no resolution artifact; ideal-atlas probe skipped"
  echo "       rerun with --resolve"
fi
