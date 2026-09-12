#!/usr/bin/env python3
"""Fast exact regressions for the reusable support-saturation compiler."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from jcsearch.support_saturation import (
    CertificateAssurance,
    CompilerOptions,
    ModulePresentation,
    NormalFiltration,
    PolynomialRing,
    SupportSaturationCompiler,
    SupportSaturationProblem,
)


GENERATED = ROOT / "artifacts" / "generated-results"


def torsion_calibration() -> None:
    # R/(x) direct-sum R/(x,y,z).  The second basis vector is precisely the
    # closed-point local-cohomology class.
    presentation = ModulePresentation(
        ring=PolynomialRing(("x", "y", "z")),
        rank=2,
        generators=(
            ("x", "0"),
            ("0", "x"),
            ("0", "y"),
            ("0", "z"),
        ),
        label="closed-point-torsion-calibration",
    )
    result = SupportSaturationCompiler().compile(
        presentation,
        ("x", "y", "z"),
        distinguished_class=("0", "1"),
        filtration=NormalFiltration(("x", "y"), (1, 2, 3)),
    )
    assert not result["saturation"]["equal_to_presentation"]
    assert result["local_cohomology"]["generator_count"] == 1
    assert result["local_cohomology"]["annihilator"] == ["z", "y", "x"]
    assert result["local_cohomology"]["annihilator_radical"] == [
        "z",
        "y",
        "x",
    ]
    assert result["local_cohomology"]["least_annihilating_exponent"] == 1
    assert result["associated_primes"]["primes"] == [
        ["x"],
        ["z", "y", "x"],
    ]
    assert result["distinguished_class"] == {
        "zero_in_F_mod_N": False,
        "belongs_to_local_cohomology": True,
        "boundary_annihilation_exponent": 1,
        "remainder": ["gen(2)"],
        "annihilator": ["z", "y", "x"],
        "annihilator_radical": ["z", "y", "x"],
    }
    assert all(
        transition["surjective"]
        for transition in result["finite_jets"]["transitions"]
    )
    uniform = result["finite_jets"]["uniform_exponent_test"]
    assert uniform["status"] == "certified_on_requested_finite_tower"
    assert not uniform["all_order_uniform_bound_certified"]


def regularity_calibration() -> None:
    # On R/(x), y is regular.  Exercise the depth-certificate path without a
    # saturation computation.
    presentation = ModulePresentation(
        ring=PolynomialRing(("x", "y")),
        rank=1,
        generators=(("x",),),
        label="regular-element-calibration",
    )
    compiler = SupportSaturationCompiler(
        CompilerOptions(
            associated_primes="regularity",
            saturation_strategy="regularity",
        )
    )
    result = compiler.compile(presentation, ("y",))
    assert result["saturation"]["strategy"] == "regularity"
    assert result["saturation"]["equal_to_presentation"]
    assert result["local_cohomology"]["zero"]
    assert result["regular_elements"]["candidate"] == "y"
    assert result["associated_primes"]["status"] == (
        "boundary_noncontainment_certified_by_regular_element"
    )


def multi_generator_saturation_calibration() -> None:
    # Saturation is all of R^2 and genuinely has two generators.  This
    # guards against the invalid legacy idiom `sat(N,I)[1]`, which selects
    # only the first generator of Singular's returned module.
    presentation = ModulePresentation(
        ring=PolynomialRing(("x", "y")),
        rank=2,
        generators=(
            ("x", "0"),
            ("y", "0"),
            ("0", "x"),
            ("0", "y"),
        ),
        label="whole-saturated-module-calibration",
    )
    result = SupportSaturationCompiler().compile(
        presentation,
        ("x", "y"),
    )
    assert result["saturation"]["saturated_presentation"] == [
        "gen(1)",
        "gen(2)",
    ]
    assert result["local_cohomology"]["generator_count"] == 2
    assert result["local_cohomology"]["associated_primes"] == [["y", "x"]]


def distinguished_witness_calibration() -> None:
    presentation = ModulePresentation(
        ring=PolynomialRing(("x", "y")),
        rank=1,
        generators=(("x",), ("y",)),
        label="distinguished-support-witness-calibration",
    )
    result = SupportSaturationCompiler(
        CompilerOptions(torsion_exponent_bound=2)
    ).compile_distinguished_support_witness(
        presentation,
        ("x", "y"),
        ("1",),
    )
    assert result["local_cohomology"]["status"] == (
        "certified_nonzero_by_distinguished_class"
    )
    assert result["distinguished_class"][
        "boundary_annihilation_exponent"
    ] == 1
    assert result["associated_primes"]["primes"] is None
    assert result["regular_elements"][
        "certifies_no_regular_element_in_boundary"
    ]
    assert result["distinguished_class"]["annihilator_radical"] == [
        "y",
        "x",
    ]


def shared_input_schema_calibration() -> None:
    payload = json.loads(
        (ROOT / "schemas" / "support_saturation_example.json").read_text()
    )
    problem = SupportSaturationProblem.from_mapping(payload)
    assert problem.to_mapping()["schema"] == SupportSaturationProblem.schema
    compiler = SupportSaturationCompiler(
        CompilerOptions(**payload["compiler_options"])
    )
    result = compiler.compile_problem(problem)
    assert result["problem"] == problem.to_mapping()
    assert result["ideals"] == {
        "support_ideal": ["x"],
        "completion_ideal": ["t"],
    }
    assert result["certificate_state"] == {
        "backend_arithmetic": "exact",
        "backend_characteristic": 0,
        "claim_assurance": "exact",
        "target_characteristic": 0,
        "characteristic_zero_lift": "not_needed",
        "note": "Exact characteristic-zero calibration.",
    }
    assert result["local_cohomology"]["annihilator_radical"] == ["x"]
    assert result["associated_prime_candidates"]["primes"] == [
        ["x"],
        ["y", "x"],
    ]
    assert result["problem"]["parameter_base_variables"] == ["y"]
    assert result["problem"]["normal_variables"] == ["t"]
    assert len(result["problem_sha256"]) == 64

    modular = SupportSaturationProblem(
        presentation=ModulePresentation(
            ring=PolynomialRing(("x", "y"), characteristic=32003),
            rank=1,
            generators=(("x",),),
        ),
        support_ideal=("y",),
        parameter_base_variables=("y",),
        normal_variables=("x",),
        assurance=CertificateAssurance(
            claim="modular",
            target_characteristic=0,
            note="No characteristic-zero lift is claimed.",
        ),
    )
    modular_result = SupportSaturationCompiler(
        CompilerOptions(
            associated_primes="regularity",
            saturation_strategy="regularity",
        )
    ).compile_problem(modular)
    state = modular_result["certificate_state"]
    assert state["backend_arithmetic"] == "exact"
    assert state["backend_characteristic"] == 32003
    assert state["claim_assurance"] == "modular"
    assert state["characteristic_zero_lift"] == "not_claimed"


def generated_artifact_regression() -> None:
    cubic_path = (
        GENERATED / "support_saturation_universal_cubic_symbols.json"
    )
    if not cubic_path.exists():
        return
    cubic = json.loads(cubic_path.read_text())
    assert cubic["case"] == "universal-cubic-homogeneous-symbol-atlas"
    assert set(cubic["certificates"]) == {
        "smooth",
        "nodal",
        "cuspidal",
        "line-transverse-conic",
        "line-tangent-conic",
        "triangle",
        "concurrent-lines",
        "double-line",
        "triple-line",
        "zero",
    }
    assert all(
        item["saturation"]["equal_to_presentation"]
        and item["local_cohomology"]["zero"]
        for item in cubic["certificates"].values()
    )

    frontier_path = (
        GENERATED / "support_saturation_cubic_annihilator_frontier.json"
    )
    if frontier_path.exists():
        frontier = json.loads(frontier_path.read_text())
        assert frontier["case"] == (
            "cubic-formal-gauge-annihilator-frontier"
        )
        strata = frontier["stratification"]
        assert strata["smooth"][
            "universal_24_parameter_cotangent_saturation"
        ] == "proved"
        assert strata["smooth"]["quartic_nongauge_dimension"] == 0
        singular_queue = strata["singular_squarefree"]["queue"]
        assert [row["annihilator_type"] for row in singular_queue] == [
            "(x)",
            "(x^2)",
            "(yz)",
            "(y^3)",
            "(xyz)",
            "(x^3)",
        ]
        assert [
            row["quartic_nongauge_dimension"] for row in singular_queue
        ] == [2, 4, 4, 6, 6, 8]
        degenerate_queue = strata["double_triple_line_and_zero"]["queue"]
        assert [row["symbol"] for row in degenerate_queue] == [
            "double-line",
            "triple-line",
            "zero",
        ]
        assert all(
            row["gate_before_saturation"]
            == "generically etale and Keller compatibility"
            for row in degenerate_queue
        )

    degree42_path = (
        GENERATED / "support_saturation_degree42_ritt_fiber_mod32003.json"
    )
    if degree42_path.exists():
        degree42 = json.loads(degree42_path.read_text())
        certificate = degree42["certificate"]
        assert degree42["characteristic"] == 32003
        assert degree42["base_specialization"] == {
            "sync42p_e1": 1,
            "sync42p_e2": 2,
            "sync42p_t": 3,
        }
        assert [
            jet["order"] for jet in certificate["finite_jets"]["jets"]
        ] == [6, 7]
        assert certificate["finite_jets"]["strategy"] == (
            "full_saturation"
        )
        assert [
            jet["boundary_annihilation_exponent"]
            for jet in certificate["finite_jets"]["jets"]
        ] == [1, 2]
        transition = certificate["finite_jets"]["transitions"][0]
        assert transition["surjective"]
        assert transition["distinguished_class_lifts"]
        support = degree42["untruncated_base"]["support_witness"]
        assert support["local_cohomology"]["status"] == (
            "certified_nonzero_by_distinguished_class"
        )
        assert support["distinguished_class"][
            "boundary_annihilation_exponent"
        ] == 2
        assert support["regular_elements"][
            "certifies_no_regular_element_in_boundary"
        ]

    plane_path = (
        GENERATED / "support_saturation_plane_jc_boundary_layer.json"
    )
    if plane_path.exists():
        plane = json.loads(plane_path.read_text())
        certificate = plane["certificate"]
        assert not certificate["local_cohomology"]["zero"]
        assert certificate["distinguished_class"][
            "belongs_to_local_cohomology"
        ]
        assert certificate["distinguished_class"][
            "boundary_annihilation_exponent"
        ] >= 1


def main() -> None:
    torsion_calibration()
    regularity_calibration()
    multi_generator_saturation_calibration()
    distinguished_witness_calibration()
    shared_input_schema_calibration()
    generated_artifact_regression()
    print("PASS: support-saturation compiler exact calibrations")
    print("PASS: local cohomology, associated primes, classes, and jets agree")


if __name__ == "__main__":
    main()
