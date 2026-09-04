#!/usr/bin/env bash
# Reproduce one complete product-twist Frobenius certificate at an extra prime.
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
pair_key=${1:-alternate-orbit-19bad:alternate-orbit-083ad}
prime=${2:-151}
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

tag=${pair_key//:/--}
local_dir="$root/artifacts/local/elkies-k3/r17-product-toric-frobenius/$tag/p$prime"
audit="$root/artifacts/generated-results/elkies-k3-r17-product-extra-prime-audits/$tag-p$prime-v1.json"
mkdir -p "$local_dir"

(
  cd "$root"
  "$sage_python" elkies-k3/scripts/audit_r17_product_twist_extra_prime.sage \
    --pair-key "$pair_key" \
    --prime "$prime"
)

"$sage_python" "$root/elkies-k3/scripts/export_r17_product_toric_frobenius_input.sage" \
  --pair-key "$pair_key" \
  --prime "$prime"

export LD_LIBRARY_PATH="$sage_prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
"$toric_root/build/examples/readfile.exe" \
  "$local_dir/toric-controlled-reduction.input" \
  "$local_dir/toric-controlled-reduction.output" \
  >"$local_dir/toric-controlled-reduction.log" 2>&1

"$sage_python" "$root/elkies-k3/scripts/certify_r17_product_toric_frobenius.sage" \
  --pair-key "$pair_key" \
  --prime "$prime" \
  --toric-output "$local_dir/toric-controlled-reduction.output" \
  --toric-commit "$toric_commit" \
  --toric-executable "$toric_root/build/examples/readfile.exe" \
  --sage-prefix "$sage_prefix" \
  --audit "$audit"
