#!/usr/bin/env sage-python
"""Aggregate the sharp Frobenius bounds for the two rank-55 singletons."""

import argparse
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUTS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-singleton-alternate-orbit-0fda0-p131-toric-frobenius-v1.json",
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-singleton-alternate-orbit-1037d-p157-toric-frobenius-v1.json",
)
AUDIT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-singleton-twist-finite-field-bound-audit-v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-rank55-singleton-toric-frobenius-bounds-v1.json"
)


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def build_payload() -> dict:
    targets = []
    expected = (
        ("alternate-orbit-0fda0", 131),
        ("alternate-orbit-1037d", 157),
    )
    for path, (label, prime) in zip(INPUTS, expected):
        record = json.loads(path.read_text())
        if (
            record.get("schema") != "elkies-k3.r17-singleton-toric-frobenius.v1"
            or record.get("character_kind") != "singleton"
            or record.get("labels") != [label]
            or int(record.get("prime")) != prime
            or record["H2"].get("degree") != 34
            or record["elliptic_L"].get("degree") != 24
            or record["elliptic_L"].get("independent_power_sum_check")
            != "PASS_AGAINST_STORED_FIBREWISE_N1_N2_AUDIT"
            or record["elliptic_L"]["weil_circle_check"].get("status")
            != "PASS_EXACT_REAL_ROOT_ISOLATION"
            or record["bounds"].get("trivial_lattice_rank") != 10
            or record["bounds"].get("geometric_twist_mw_rank_lower_bound") != 1
            or record["bounds"].get("geometric_twist_mw_rank_upper_bound") != 2
        ):
            raise ArithmeticError(f"singleton Frobenius gate failed for {label}")
        hits = record["elliptic_L"]["cyclotomic_hits_after_T_equals_pZ"]
        if len(hits) != 1 or (
            hits[0].get("order"), hits[0].get("multiplicity"), hits[0].get("total_degree")
        ) != (1, 2, 2):
            raise ArithmeticError(f"unexpected Tate factor for {label}")
        targets.append(
            {
                "label": label,
                "prime": prime,
                "geometric_mw_rank_interval": [1, 2],
                "normalized_tate_factor": "(Z-1)^2",
                "certificate": str(path.relative_to(ROOT)),
                "certificate_sha256": digest(path),
            }
        )
    return {
        "schema": "elkies-k3.r17-rank55-singleton-toric-frobenius-bounds.v1",
        "status": "PASS_GEOMETRIC_SINGLETON_RANK_INTERVALS_1_2",
        "targets": targets,
        "character_rank_interval": "17+[1,2]+[1,2]+0",
        "total_geometric_generic_rank_interval": [19, 21],
        "proof_boundary": (
            "The two singleton upper bounds are exact and unconditional, but the "
            "second possible direction is not excluded. Exact singleton rank one "
            "and the decomposition 17+1+1+0 remain UNKNOWN."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (Path(__file__).resolve(), AUDIT, *INPUTS)
        },
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--check", action="store_true")
parser.add_argument("--output", type=Path, default=OUTPUT)
args = parser.parse_args()
payload = build_payload()
if args.check:
    if json.loads(args.output.read_text()) != payload:
        raise ArithmeticError("stored singleton aggregate does not replay")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "R17RANK55SINGLETONS|bounds=1:2,1:2|characters=17+[1,2]+[1,2]+0|"
    f"status={payload['status']}",
    flush=True,
)

