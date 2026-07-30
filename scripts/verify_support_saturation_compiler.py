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
    CompilerOptions,
    ModulePresentation,
    NormalFiltration,
    PolynomialRing,
    SupportSaturationCompiler,
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
    }
    assert all(
        transition["surjective"]
        for transition in result["finite_jets"]["transitions"]
    )


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
    generated_artifact_regression()
    print("PASS: support-saturation compiler exact calibrations")
    print("PASS: local cohomology, associated primes, classes, and jets agree")


if __name__ == "__main__":
    main()
