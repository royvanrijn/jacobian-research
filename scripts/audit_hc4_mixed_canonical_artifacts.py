#!/usr/bin/env python3
"""Audit committed HC4MCP1--10 ledgers without rerunning their searches."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "artifacts" / "generated-results"
EXPECTED = {
    "HC4MCP1": (
        "hc4_mixed_canonical_pivot_search.json",
        "10574297d6ba72240e7e1c7acb217be0a545badc1eae022927fea64063bbbacd",
    ),
    "HC4MCP2": (
        "hc4_mixed_quadratic_words.json",
        "17402479fe5b60068459f1936665d970e2204e9cc6f78f7b7ad3b88e06a5d6f7",
    ),
    "HC4MCP3": (
        "hc4_mixed_quadratic_cubic_words.json",
        "5bdc44741ad499f9c6f63be90d2e8d3a4e02e9675160f62490be61dd5ab6c10c",
    ),
    "HC4MCP4": (
        "hc4_canonical_signed_quadratic_cubic_words.json",
        "0baf8823d8ed090356d756b0ba9689b91402ce5fc6c6ea9044de40077d2b1f3d",
    ),
    "HC4MCP5": (
        "hc4_symbolic_quadratic_cubic_words.json",
        "c7ebbc6345d3037a9fdbecb080ca152ed7390cd6db7fcc5dd36dd7b8ffde082d",
    ),
    "HC4MCP6": (
        "hc4_symbolic_cubic_quadratic_words.json",
        "1c33fdbcb2296efe04df6e6e86d79bd407793a1adc5142f89e11912b419a36d9",
    ),
    "HC4MCP7": (
        "hc4_mixed_quadratic_cubic_commutators.json",
        "95cd7757483cc71e97c8ed8925a0bce9e2d351794b2366a8ee83dab41f1ab359",
    ),
    "HC4MCP8": (
        "hc4_noncoordinate_coisotropic_scalar_gate.json",
        "a6a806d649c1af0c1ea2c26e01937817832de502daacda97a1a015faca472eb2",
    ),
    "HC4MCP9": (
        "hc4_nonlinear_unit_schur_blocks.json",
        "7f235d427e0cf63e3aeddf198d6ade72c5478ae90774d526fb3a5610dae9286e",
    ),
    "HC4MCP10": (
        "hc4_symbolic_unit_schur_classification.json",
        "e92465e4991e7635f07fcc70895995f5d0465a1c3b816a6c9a88643500865e30",
    ),
}


def load(entry_id: str) -> dict[str, object]:
    filename, expected_hash = EXPECTED[entry_id]
    path = GENERATED / filename
    actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    assert actual_hash == expected_hash, (entry_id, actual_hash, expected_hash)
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    rows = {entry_id: load(entry_id) for entry_id in EXPECTED}

    for entry_id in ("HC4MCP1", "HC4MCP2", "HC4MCP3", "HC4MCP4"):
        assert rows[entry_id]["status"] == "bounded_search", entry_id
        assert rows[entry_id]["scope"]["prime"] == 1_000_003, entry_id
    assert rows["HC4MCP4"]["scope"]["affine_subspaces"].startswith(
        "coordinate subspaces only"
    )

    for entry_id in ("HC4MCP5", "HC4MCP6"):
        assert rows[entry_id]["status"] == "exact_symbolic_finite_support_theorem"
        assert rows[entry_id]["summary"]["patterns"] == 54
        assert "longer words" in rows[entry_id]["scope"]["limitations"]
    assert rows["HC4MCP5"]["summary"]["parent_survivors"] == 0
    assert rows["HC4MCP6"]["summary"]["parent_survivors"] == 54
    assert rows["HC4MCP6"]["summary"]["rank_four_two_pivot_pairs"] == 102

    assert rows["HC4MCP7"]["claim_type"] == "bounded_exact_computation"
    assert rows["HC4MCP7"]["summary"]["parent_survivors_requiring_expansion"] == 0
    assert rows["HC4MCP8"]["status"] == "bounded_exact_search"
    assert rows["HC4MCP8"]["scope"]["full_descended_determinants"] == "not formed"
    assert rows["HC4MCP9"]["status"] == "exact_finite_box_search"
    assert rows["HC4MCP9"]["summary"]["unit_schur_blocks"] == 0
    assert "directions outside the declared coefficient box" in rows[
        "HC4MCP9"
    ]["scope"]["limitations"]

    symbolic = rows["HC4MCP10"]
    assert symbolic["status"] == "exact_symbolic_family_classification"
    assert symbolic["scope"]["b"] == "arbitrary nonzero complex coefficient"
    assert symbolic["summary"]["patterns"] == 54
    assert symbolic["summary"]["grassmann_charts"] == 810
    assert symbolic["summary"]["constant_two_plane_unit_schur_blocks"] == 0
    assert symbolic["summary"]["certificate_kinds"] == {
        "joint_cubic_plus_origin_determinant": 654,
        "sampled_determinant_unit_ideal": 156,
    }

    print(
        "PASS: 10 committed HC4MCP1--10 ledgers match exact hashes and retain "
        "their finite-box versus coefficient-uniform scopes; no SymPy import, "
        "Singular run, search, or artifact rewrite"
    )
    print(
        "SCOPE: HC4MCP10 removes constant coefficient/direction bounds only on "
        "the 54 HC4MCP6 resonance families; moving pivots, other supports, and "
        "longer canonical words remain open"
    )


if __name__ == "__main__":
    main()
