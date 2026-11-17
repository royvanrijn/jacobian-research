#!/usr/bin/env python3
"""Check retained F2 checkpoint integrity and metadata without rebuilding a DAG.

This does not replay algebra, permutation enumeration or formal continuation.
Hashes certify retained bytes; the canonical proof notes determine their scope.
Missing files fail the audit and never trigger automatic reconstruction.
"""

from copy import deepcopy
import hashlib
import json
from pathlib import Path

from f2_formal_homotopy_metadata import boundary_text, comparison_record

ROOT = Path(__file__).resolve().parents[2]
RETAINED = {
    "jc2_f2_75_125_carrier_specializations.json": "00602d2c614ded350710c6c6c72f418352d153e17a7163ee1118f015fad7242a",
    "jc2_f2_75_125_global_attachment.json": "419c970e322b16e1bfb6403dc36b1a38b95eb9a52403def6b7ee067c42fe8ddc",
    "jc2_f2_a6_simple_spectator_gluing.json": "f91e4e3c187c089fa560bf4cade1604434c404fa92bf9a3a79c8a618164613be",
    "jc2_f2_modified_chart_bridge.json": "ac7dbc170cafbcf028079b9ccdb41afd78c333e3850abc0761f64cc056e7d7b8",
    "jc2_f2_modified_laurent_family.json": "bca206498c153e41a2f31344015df2ce63890f8b12228b6d0ba2c0970eb87c85",
    "jc2_f2_75_125_nonlinear_forcing.json": "381854d84a0377a14745280c3b74ac5376f78fa7a07906eac24a826fd566daed",
    "jc2_f2_75_125_modular_probe.json": "4a64dc8393d43047563c45f2760b49224bb724dcf80d20b0ed6f402580432b55",
    "jc2_f2_75_125_tangent_obstruction.json": "ce5b7fe491ff58b884a986aa6041b69fd3938c74c2593b6ab1e1cbf5cf1d7027",
    "jc2_f2_75_125_formal_homotopy.json": "ef65694180eef00fe24652e8ecc6ecfb4f120eccd093989059906a804e9a41ef",
    "jc2_f2_75_125_formal_homotopy_regular_gauge.json": "7d0a9818deba415a84f83be8265b1e05903d16b34d75af0d9e27dd22ed8d412f",
    "jc2_f2_75_125_formal_homotopy_mod61.json": "10b76a5f042f84b88e18003eb6152f56954aea3b0d7454f6f7abed29b6c8b68f",
}


def digest(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def main() -> None:
    records = {}
    for name, expected in RETAINED.items():
        raw = (ROOT / "artifacts/generated-results" / name).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == expected, f"retained checkpoint changed: {name}"
        records[name] = json.loads(raw)

    forcing = records["jc2_f2_75_125_nonlinear_forcing.json"]
    assert forcing["equation_ledger"]["total"] == 366
    coupled = forcing["coupled_Laurent_forcing"]
    assert coupled["divisibility_coordinate_count"] == 294
    assert coupled["pinned_quotient_cokernel_coordinate_count"] == 53
    assert coupled["full_Laurent_cokernel_coordinate_count"] == 347
    assert not records["jc2_f2_75_125_modular_probe.json"]["converged"]

    for suffix, prime, order, residuals in (
        ("", 31, 16, 62),
        ("_regular_gauge", 31, 16, 50),
        ("_mod61", 61, 8, 62),
    ):
        record = records[f"jc2_f2_75_125_formal_homotopy{suffix}.json"]
        assert record["field"]["prime"] == prime
        rho, y = record["field"]["rho"], record["field"]["y"]
        assert (rho * rho - 3 * rho + 1) % prime == 0
        assert (27 * y * y - 9 * y + 1) % prime == 0
        assert record["achieved_order"] == record["requested_order"] == order
        assert record["Jacobian_rank"] == 214
        assert record["Jacobian_cokernel_dimension"] == 153
        assert record["obstruction"] is None
        assert [step["order"] for step in record["steps"]] == list(range(1, order + 1))
        assert all(step["cokernel_projection_nonzero_count"] == 0 for step in record["steps"])
        series = record["variable_series"]
        assert digest(series) == record["variable_series_digest_sha256"]
        assert all(len(row) == order + 1 for row in series.values())
        point = {name: sum(row) % prime for name, row in series.items() if sum(row) % prime}
        target = record["truncated_lambda_one_evaluation"]
        assert digest(point) == target["point_digest_sha256"]
        assert target["nonzero_total"] == residuals == sum(target["nonzero_by_group"].values())
        assert not target["is_exact_modular_point"]
        assert (ROOT / record["source_circuit_artifact"]).is_file()

        # Exercise the two exact metadata aliases against real retained records.
        corrected = deepcopy(record)
        corrected["claim_boundary"] = boundary_text(prime)
        if suffix == "_regular_gauge":
            corrected["higher_order_gauge"]["name"] = "selected-coordinate-zero"
        assert comparison_record(record) == comparison_record(corrected)

    print(f"PASS: {len(records)} retained F2 hashes; three formal-jet metadata audits")
    print("The regular order-16 jet still has 50 nonzero equations at lambda=1.")
    print("No algebra, permutation enumeration or formal continuation was replayed.")


if __name__ == "__main__":
    main()
