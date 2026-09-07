#!/usr/bin/env bash
# Reproduce one complete alternate-Q80 singleton-twist Frobenius certificate.
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
label=${1:?usage: run_r17_singleton_toric_frobenius.sh LABEL [PRIME]}
prime=${2:-131}
toric_commit=74cda9e8148cd8e9a3928fc15a558c9a70b67cc1
toric_root="$root/artifacts/local/tools/ToricControlledReduction-$toric_commit"
sage_command=${SAGE_COMMAND:-sage}
sage_python=$(
  "$sage_command" -c 'import sys; print(sys.executable)'
)
sage_prefix=$(cd "$(dirname "$sage_python")/.." && pwd)

if [[ ! -d "$toric_root/.git" ]]; then
  git clone https://github.com/edgarcosta/ToricControlledReduction.git "$toric_root"
  git -C "$toric_root" checkout --detach "$toric_commit"
fi
if [[ $(git -C "$toric_root" rev-parse HEAD) != "$toric_commit" ]]; then
  echo "unexpected ToricControlledReduction checkout at $toric_root" >&2
  exit 1
fi
if [[ ! -x "$toric_root/build/examples/readfile.exe" ]]; then
  (
    cd "$toric_root"
    ./configure --disable-gdb --with-ntl="$sage_prefix" --with-gmp="$sage_prefix"
    make
  )
fi

local_dir="$root/artifacts/local/elkies-k3/r17-singleton-toric-frobenius/$label/p$prime"
mkdir -p "$local_dir"

"$sage_python" "$root/elkies-k3/scripts/export_r17_singleton_toric_frobenius_input.sage" \
  --label "$label" \
  --prime "$prime"

export LD_LIBRARY_PATH="$sage_prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
"$toric_root/build/examples/readfile.exe" \
  "$local_dir/toric-controlled-reduction.input" \
  "$local_dir/toric-controlled-reduction.output" \
  >"$local_dir/toric-controlled-reduction.log" 2>&1

"$sage_python" "$root/elkies-k3/scripts/certify_r17_product_toric_frobenius.sage" \
  --pair-key "$label" \
  --prime "$prime" \
  --toric-output "$local_dir/toric-controlled-reduction.output" \
  --toric-commit "$toric_commit" \
  --toric-executable "$toric_root/build/examples/readfile.exe" \
  --sage-prefix "$sage_prefix"

