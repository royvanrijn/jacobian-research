#!/usr/bin/env python3
"""Verify the pinned dessin-first degree-(72,108) deformation certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
CAS = ROOT / "plane-jc" / "cas"
CORE = CAS / "jc2_degree108_belyi_deformations.py"
EXPECTED_CORE_SHA256 = "2930cfba18577b9986ca224d2dce8b40f65effae16cedd0b945932d6db2b36d6"
EXPECTED_CERTIFICATE_SHA256 = (
    "c319b662474c856644393cef2cd8c81f4f54cf93aff535bfe0d6e6a6fc84ee33"
)
assert hashlib.sha256(CORE.read_bytes()).hexdigest() == EXPECTED_CORE_SHA256
sys.path.insert(0, str(CAS))

from jc2_degree108_belyi_deformations import OUT, compile_report  # noqa: E402


def audit_existing_only(certificate: Path) -> None:
    """Validate the committed result graph without exact reconstruction."""

    actual_hash = hashlib.sha256(certificate.read_bytes()).hexdigest()
    if actual_hash != EXPECTED_CERTIFICATE_SHA256:
        raise RuntimeError(
            "pinned Belyi-deformation certificate bytes changed: "
            f"expected {EXPECTED_CERTIFICATE_SHA256}, got {actual_hash}"
        )
    pinned = json.loads(certificate.read_text(encoding="utf-8"))
    assert pinned["scope"] == "no-vertical-edge (72,108) / (8,28) Laurent residue only"
    assert "not a stand-alone proof" in pinned["claim_boundary"]

    field = pinned["coefficient_field"]
    assert field["field_of_moduli_degree"] == 5
    assert field["irreducible_over_Q"] is True
    assert field["galois_group"] == "S5"
    assert field["galois_group_order"] == 120

    dessins = pinned["dessins"]
    assert dessins["labelled_center_sets"] == 85
    assert dessins["rotation_orbits"] == 5
    representatives = dessins["representatives"]
    assert len(representatives) == 5
    center_sets = [tuple(row["center_set"]) for row in representatives]
    assert len(center_sets) == len(set(center_sets)) == 5
    for row in representatives:
        assert row["transitive"] is True
        assert row["automorphism_order"] == 1
        for key in ("sigma_0", "sigma_1", "sigma_infinity"):
            assert sorted(row[key]) == list(range(21))

    top = pinned["top_reconstruction"]
    assert top["passport"] == [
        [2] * 10 + [1],
        [3] * 7,
        [17, 1, 1, 1, 1],
    ]
    assert top["infinity_multiplicity"] == 17
    assert top["third_fiber_degree"] == 4
    assert top["squarefree_finite_fibers"] is True
    quotient_graph = ROOT / top["quotient_graph_source"]
    assert hashlib.sha256(quotient_graph.read_bytes()).hexdigest() == (
        top["quotient_graph_sha256"]
    )

    deformations = pinned["deformations"]
    assert deformations["parameter_names"] == ["p0", "p1", "p2", "p3", "p4"]
    expected_maps = {
        "BE": ([20, 19], 17, 2, 3),
        "CF": ([20, 21], 18, 3, 2),
        "G_derivative": ([20, 12], 12, 0, 8),
    }
    for label, (shape, rank, kernel, cokernel) in expected_maps.items():
        row = deformations["linear_maps"][label]
        assert row["shape"] == shape
        assert row["rank"] == rank
        assert row["kernel_dimension"] == kernel == shape[1] - rank
        assert row["cokernel_dimension_in_declared_target"] == cokernel == shape[0] - rank
    assert deformations["G_constant_kernel_dimension"] == 1
    assert deformations["G_cokernel_equations"] == 7
    assert deformations["J0_equations"] == 18
    assert deformations["final_equation_count"] == 25
    assert len(deformations["final_equation_term_counts"]) == 25

    terminal = pinned["terminal_open_check"]
    assert terminal == {
        "saturation": "<final equations, z*B_8-1>",
        "unit_ideal": True,
        "basis_size": 1,
        "singular_input_hash": (
            "496015c8103357d8010837ed85da050fc330fb8206d977c24cf7f157fcce51bb"
        ),
    }
    print(
        "JC2_DEGREE108_BELYI_COMMITTED_AUDIT_PASS "
        "(5 exact labelled representatives; quotient graph; deformation ranks; "
        "pinned saturated unit record; no reconstruction or Singular run)"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quick",
        action="store_true",
        help="skip replaying the terminal five-parameter Singular unit-ideal check",
    )
    parser.add_argument("--certificate", type=Path, default=OUT)
    parser.add_argument("--singular-timeout", type=int, default=300)
    parser.add_argument(
        "--audit-existing-only",
        action="store_true",
        help="validate the committed certificate without reconstruction or Singular",
    )
    arguments = parser.parse_args()

    certificate = arguments.certificate.resolve()
    if arguments.audit_existing_only:
        audit_existing_only(certificate)
        return

    pinned = json.loads(certificate.read_text(encoding="utf-8"))
    run_full = not arguments.quick
    replay = compile_report(run_full, arguments.singular_timeout)
    if not run_full:
        replay["terminal_open_check"] = pinned["terminal_open_check"]
    if replay != pinned:
        raise RuntimeError("pinned Belyi-deformation certificate differs from exact replay")

    assert pinned["coefficient_field"]["galois_group"] == "S5"
    assert pinned["dessins"]["rotation_orbits"] == 5
    assert pinned["top_reconstruction"]["passport"] == [
        [2] * 10 + [1],
        [3] * 7,
        [17, 1, 1, 1, 1],
    ]
    assert pinned["deformations"]["linear_maps"]["BE"]["kernel_dimension"] == 2
    assert pinned["deformations"]["linear_maps"]["CF"]["kernel_dimension"] == 3
    assert pinned["deformations"]["G_constant_kernel_dimension"] == 1
    assert pinned["terminal_open_check"]["unit_ideal"] is True
    print("JC2_DEGREE108_BELYI_DEFORMATIONS_PASS")
    print(f"MODE={'full' if run_full else 'quick'}")
    print(f"CERTIFICATE={certificate.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
