#!/usr/bin/env bash
# Run the complete two-prime cohomological triage of the 17 selected products.
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
jobs=${JOBS:-3}
runner="$root/elkies-k3/scripts/run_r17_product_toric_frobenius.sh"
audit="$root/artifacts/generated-results/elkies-k3-r17-norm12-11952-product-twist-finite-field-bound-audit-v1.json"

python3 - "$audit" <<'PY' \
  | xargs -P "$jobs" -I{} "$runner" '{}' 131
import json
import sys
for record in json.load(open(sys.argv[1]))["targets"]:
    print(record["pair_key"])
PY

python3 - "$root" "$audit" <<'PY' \
  | xargs -P "$jobs" -r -I{} "$runner" '{}' 137
import json
from pathlib import Path
import sys
root = Path(sys.argv[1])
for record in json.load(open(sys.argv[2]))["targets"]:
    pair_key = record["pair_key"]
    tag = pair_key.replace(":", "--")
    path = root / "artifacts/generated-results" / (
        f"elkies-k3-r17-product-{tag}-p131-toric-frobenius-v1.json"
    )
    result = json.load(open(path))
    if result["bounds"]["geometric_twist_mw_rank_upper_bound"]:
        print(pair_key)
PY

sage -python \
  "$root/elkies-k3/scripts/certify_r17_all17_product_toric_frobenius_campaign.sage"

