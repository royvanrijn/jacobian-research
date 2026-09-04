#!/usr/bin/env sage-python
"""Aggregate the two-prime toric Frobenius campaign for all 17 products."""

import argparse
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AUDIT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-product-twist-finite-field-bound-audit-v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-all17-product-toric-frobenius-campaign-v1.json"
)
EXPORTER = ROOT / "elkies-k3/scripts/export_r17_product_toric_frobenius_input.sage"
PARSER = ROOT / "elkies-k3/scripts/parse_toric_controlled_reduction_output.py"
VERIFIER = ROOT / "elkies-k3/scripts/certify_r17_product_toric_frobenius.sage"
RUNNER = ROOT / "elkies-k3/scripts/run_r17_product_toric_frobenius.sh"
TORIC_COMMIT = "74cda9e8148cd8e9a3928fc15a558c9a70b67cc1"


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def artifact_path(pair_key: str, prime: int) -> Path:
    tag = pair_key.replace(":", "--")
    return (
        ROOT
        / "artifacts/generated-results"
        / f"elkies-k3-r17-product-{tag}-p{prime}-toric-frobenius-v1.json"
    )


def validate_reduction(path: Path, pair_key: str, prime: int) -> dict:
    record = json.loads(path.read_text())
    hits = record["elliptic_L"]["cyclotomic_hits_after_T_equals_pZ"]
    tate_degree = sum(int(hit["total_degree"]) for hit in hits)
    software = record["software"]
    if (
        record.get("schema") != "elkies-k3.r17-product-toric-frobenius.v1"
        or record.get("character_kind") != "product"
        or record.get("pair_key") != pair_key
        or int(record.get("prime")) != prime
        or record["H2"].get("degree") != 46
        or record["elliptic_L"].get("degree") != 28
        or record["elliptic_L"].get("independent_power_sum_check")
        != "PASS_AGAINST_STORED_FIBREWISE_N1_N2_AUDIT"
        or record["elliptic_L"]["weil_circle_check"].get("status")
        != "PASS_EXACT_REAL_ROOT_ISOLATION"
        or record["bounds"].get("trivial_lattice_rank") != 18
        or record["bounds"].get("geometric_twist_mw_rank_upper_bound")
        != tate_degree
        or software.get("ToricControlledReduction_commit") != TORIC_COMMIT
        or software.get("exporter_sha256") != digest(EXPORTER)
        or software.get("raw_output_parser_sha256") != digest(PARSER)
        or software.get("independent_verifier_sha256") != digest(VERIFIER)
        or software.get("runner_sha256") != digest(RUNNER)
    ):
        raise ArithmeticError(f"invalid product Frobenius artifact {path}")
    return {
        "prime": prime,
        "tate_degree_with_multiplicity": tate_degree,
        "geometric_mw_rank_upper_bound": tate_degree,
        "certificate": str(path.relative_to(ROOT)),
        "certificate_sha256": digest(path),
    }


def build_payload() -> dict:
    audit = json.loads(AUDIT.read_text())
    pair_keys = [record["pair_key"] for record in audit["targets"]]
    if len(pair_keys) != 17 or len(set(pair_keys)) != 17:
        raise ArithmeticError("campaign target set is not the pinned 17")
    targets = []
    zero_pairs = []
    survivor_pairs = []
    for pair_key in pair_keys:
        reductions = [validate_reduction(artifact_path(pair_key, 131), pair_key, 131)]
        if reductions[0]["tate_degree_with_multiplicity"]:
            reductions.append(
                validate_reduction(artifact_path(pair_key, 137), pair_key, 137)
            )
        best_upper = min(row["geometric_mw_rank_upper_bound"] for row in reductions)
        closed_rank_zero = best_upper == 0
        (zero_pairs if closed_rank_zero else survivor_pairs).append(pair_key)
        targets.append(
            {
                "pair_key": pair_key,
                "reductions": reductions,
                "best_geometric_mw_rank_upper_bound": best_upper,
                "geometric_rank_zero": closed_rank_zero,
            }
        )
    return {
        "schema": "elkies-k3.r17-all17-product-toric-frobenius-campaign.v1",
        "status": "PASS_COMPLETE_TWO_PRIME_TRIAGE_OF_ALL_17_PRODUCTS",
        "target_count": 17,
        "rank_zero_count": len(zero_pairs),
        "persistent_tate_survivor_count": len(survivor_pairs),
        "rank_zero_pair_keys": zero_pairs,
        "persistent_tate_survivor_pair_keys": survivor_pairs,
        "targets": targets,
        "proof_boundary": (
            "Every zero-Tate target has geometric product-twist rank zero. A "
            "persistent Tate factor at both reductions is only a cohomological "
            "survivor and does not prove a characteristic-zero section."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (Path(__file__).resolve(), AUDIT)
        },
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--check", action="store_true")
parser.add_argument("--output", type=Path, default=OUTPUT)
args = parser.parse_args()
payload = build_payload()
if args.check:
    if json.loads(args.output.read_text()) != payload:
        raise ArithmeticError("stored all-17 Frobenius campaign does not replay")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"R17ALL17TORIC|rank_zero={payload['rank_zero_count']}|"
    f"survivors={payload['persistent_tate_survivor_count']}|"
    f"status={payload['status']}",
    flush=True,
)
