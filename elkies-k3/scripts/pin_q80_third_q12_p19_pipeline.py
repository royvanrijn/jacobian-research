#!/usr/bin/env python3
"""Initialize or verify the fail-closed Q80 p=19 third-q12 hash manifest."""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "artifacts/generated-results/q80-third-q12-p19-pipeline-hash-manifest.json"

FILES = (
    ("worker", "canonical_horizontal", "elkies-k3/scripts/certify_q80_po0_rur_third_q12_modp.sage"),
    ("input", "canonical_horizontal", "artifacts/generated-results/q80-fixed-u-minus2-p19-height-shell-with-po1.json"),
    ("input", "canonical_horizontal", "artifacts/generated-results/q80-fixed-u-minus2-p19-po0-recursive-saturated-msolve/q80-third-q12-um2d1-p19-po0-polynomial-recursive-sign+1.solve"),
    ("output", "canonical_horizontal", "artifacts/generated-results/q80-fixed-u-minus2-p19-po0-rur-third-q12-modp.json"),
    ("worker", "independent_horizontal", "elkies-k3/scripts/certify_q80_third_q12_horizontal_mod19_quadratic.sage"),
    ("worker", "independent_horizontal", "elkies-k3/scripts/export_q80_third_q12_polynomial_closure_scheme.sage"),
    ("input", "independent_horizontal", "artifacts/generated-results/q80-third-q12-um2-p19-height-shell-complete.json"),
    ("input", "independent_horizontal", "artifacts/generated-results/q80-d7d5-mw5-height-lattice.json"),
    ("intermediate", "independent_horizontal", "artifacts/local/elkies-k3/q80-third-q12-um2-p19-polynomial-closure.ms"),
    ("intermediate", "independent_horizontal", "artifacts/generated-results/q80-third-q12-um2-p19-polynomial-closure-scheme.json"),
    ("output", "independent_horizontal", "artifacts/generated-results/q80-third-q12-um2-p19-quadratic-horizontal.json"),
    ("worker", "resolved_pencil", "elkies-k3/scripts/compile_q80_third_q12_um2_p19_resolved_pencil.sage"),
    ("output", "resolved_pencil", "artifacts/generated-results/q80-third-q12-um2-p19-resolved-pencil.json"),
    ("worker", "genus", "elkies-k3/scripts/verify_q80_third_q12_um2_p19_resolved_genus.sage"),
    ("input", "genus", "elkies-k3/scripts/analyze_q80_second_neighbor_chamber.sage"),
    ("input", "genus", "elkies-k3/data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt"),
    ("input", "genus", "elkies-k3/data/fibrations/kumar_q80_to_rootless_path.tsv"),
    ("output", "genus", "artifacts/generated-results/q80-third-q12-um2-p19-resolved-genus.json"),
)

EXPECTED_STATUSES = {
    "artifacts/generated-results/q80-fixed-u-minus2-p19-po0-rur-third-q12-modp.json":
        "PASS_EXACT_MODP2_THIRD_Q12_HORIZONTAL_FROBENIUS_ORBIT",
    "artifacts/generated-results/q80-third-q12-um2-p19-polynomial-closure-scheme.json":
        "PASS_EXACT_POLYNOMIAL_CLOSURE_SCHEME_EXPORTED",
    "artifacts/generated-results/q80-third-q12-um2-p19-quadratic-horizontal.json":
        "PASS_EXACT_QUADRATIC_THIRD_Q12_HORIZONTAL_MOD19",
    "artifacts/generated-results/q80-third-q12-um2-p19-resolved-pencil.json":
        "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_MOD19_QUADRATIC",
    "artifacts/generated-results/q80-third-q12-um2-p19-resolved-genus.json":
        "PASS_EXACT_THIRD_Q12_GENUS_ONE_BY_ADJUNCTION_MOD19_QUADRATIC",
}

REPRODUCE = (
    "sage -python elkies-k3/scripts/certify_q80_po0_rur_third_q12_modp.sage",
    "sage -python -c 'from pathlib import Path; from sage.all import *; "
    "from sage.repl.preparse import preparse_file; "
    "p=Path(\"elkies-k3/scripts/certify_q80_third_q12_horizontal_mod19_quadratic.sage\").resolve(); "
    "g=globals(); g.update({\"__file__\":str(p),\"__name__\":\"__main__\"}); "
    "exec(compile(preparse_file(p.read_text()),str(p),\"exec\"),g)'",
    "sage -python elkies-k3/scripts/compile_q80_third_q12_um2_p19_resolved_pencil.sage",
    "sage -python elkies-k3/scripts/verify_q80_third_q12_um2_p19_resolved_genus.sage",
    "python3 elkies-k3/scripts/pin_q80_third_q12_p19_pipeline.py",
)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def current_records():
    records = []
    for role, stage, relative in FILES:
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(relative)
        record = {
            "path": relative,
            "role": role,
            "stage": stage,
            "sha256": sha256(path),
            "size_bytes": path.stat().st_size,
        }
        expected_status = EXPECTED_STATUSES.get(relative)
        if expected_status is not None:
            payload = json.loads(path.read_text())
            if payload.get("status") != expected_status:
                raise ArithmeticError(
                    f"{relative}: expected status {expected_status}, got {payload.get('status')}"
                )
            record["status"] = expected_status
        records.append(record)
    return records


def verify_embedded_hashes(records):
    pinned = {record["path"]: record["sha256"] for record in records}
    checked = []
    for relative in EXPECTED_STATUSES:
        payload = json.loads((ROOT / relative).read_text())
        containers = []
        for key in ("inputs", "input"):
            value = payload.get(key)
            if isinstance(value, dict):
                containers.append(value)
        for container in containers:
            for item in container.values():
                if not isinstance(item, dict) or "path" not in item or "sha256" not in item:
                    continue
                path = Path(item["path"])
                if path.is_absolute():
                    try:
                        path = path.relative_to(ROOT)
                    except ValueError:
                        continue
                normalized = path.as_posix()
                if normalized in pinned and pinned[normalized] != item["sha256"]:
                    raise ArithmeticError(
                        f"{relative}: embedded hash mismatch for {normalized}"
                    )
                checked.append({"consumer": relative, "dependency": normalized})
    return checked


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--initialize",
    action="store_true",
    help="create the manifest once; refuses to overwrite an existing pin",
)
args = parser.parse_args()

records = current_records()
embedded = verify_embedded_hashes(records)
if args.initialize:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite pinned manifest {OUTPUT}")
    payload = {
        "schema": "elkies-k3.q80-third-q12-p19-pipeline-hash-manifest.v1",
        "status": "PASS_PINNED_Q80_THIRD_Q12_P19_PIPELINE",
        "specialization": {"u": "-2", "prime": 19, "extension": "r^2+12r+3"},
        "files": records,
        "embedded_dependency_hashes_checked": embedded,
        "reproduce": list(REPRODUCE),
        "refresh_policy": (
            "This manifest is immutable. If a pinned byte changes, investigate and create "
            "a versioned replacement rather than refreshing this file in place."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
else:
    if not OUTPUT.is_file():
        raise FileNotFoundError("manifest absent; run once with --initialize")
    payload = json.loads(OUTPUT.read_text())
    expected = {record["path"]: record for record in payload["files"]}
    actual = {record["path"]: record for record in records}
    if set(expected) != set(actual):
        raise ArithmeticError("pinned file set changed")
    for path, record in actual.items():
        if record["sha256"] != expected[path]["sha256"]:
            raise ArithmeticError(f"pinned hash changed: {path}")
        if record["size_bytes"] != expected[path]["size_bytes"]:
            raise ArithmeticError(f"pinned size changed: {path}")

print(
    f"Q80THIRDQ12P19PIN|files={len(records)}|embedded={len(embedded)}|"
    f"mode={'initialize' if args.initialize else 'verify'}|"
    "status=PASS_PINNED_Q80_THIRD_Q12_P19_PIPELINE",
    flush=True,
)
