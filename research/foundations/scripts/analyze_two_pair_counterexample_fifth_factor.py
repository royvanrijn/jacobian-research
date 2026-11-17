#!/usr/bin/env python3
"""Factor reconstruction for the SIC2C4 fifth obstruction on a P2 slice.

This reads exact component-wide obstruction artifacts on h3=h4=0.  Eight
effective samples reconstruct a projective Q(sqrt(41))-factor within a
quadratic ansatz, and a ninth sample is reserved for exact validation.  The
finite interpolation does not by itself prove a universal degree bound.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_ROOT = ROOT / "artifacts" / "generated-results"
OUTPUT = (
    ARTIFACT_ROOT
    / "two_pair_counterexample_fifth_factor_plane_research.json"
)

FIT_SAMPLES = [
    ("axis_h0", (1, 0, 0)),
    ("plane01_mix01", (1, 1, 0)),
    ("plane01_mix02", (1, 0, 1)),
    ("plane01_check111", (1, 1, 1)),
    ("plane01_fit120", (1, 2, 0)),
    ("plane01_fit102", (1, 0, 2)),
    ("plane01_check123", (1, 2, 3)),
    ("plane01_fit121", (1, 2, 1)),
]
CHECK_SAMPLE = ("plane01_check112", (1, 1, 2))


def artifact_path(tag: str) -> Path:
    return (
        ARTIFACT_ROOT
        / f"two_pair_counterexample_fifth_component_research_{tag}.json"
    )


def obstruction_ratio(tag: str) -> sp.Rational:
    data = json.loads(artifact_path(tag).read_text())
    obstruction = sp.expand(
        sp.sympify(data["uniform_obstruction"]["constant"])
    )
    rational_part = obstruction.coeff(sp.sqrt(41), 0)
    radical_part = obstruction.coeff(sp.sqrt(41), 1)
    assert radical_part != 0
    return sp.factor(rational_part / radical_part)


def interpolation_row(
    point: tuple[int, int, int],
    value: sp.Rational,
) -> list[sp.Expr]:
    h0, h1, h2 = point
    return [
        h0**2,
        h0 * h1,
        h0 * h2,
        -value * h0**2,
        -value * h0 * h1,
        -value * h0 * h2,
        -value * h1**2,
        -value * h1 * h2,
        -value * h2**2,
    ]


def primitive_integer_vector(vector: sp.Matrix) -> list[int]:
    denominator = sp.ilcm(*[sp.denom(value) for value in vector])
    integers = [
        int(sp.factor(value * denominator))
        for value in vector
    ]
    common = abs(sp.gcd_list(integers))
    integers = [int(value // common) for value in integers]
    if integers[0] < 0:
        integers = [-value for value in integers]
    return integers


def main() -> None:
    fit_values = [
        obstruction_ratio(tag)
        for tag, _ in FIT_SAMPLES
    ]
    matrix = sp.Matrix(
        [
            interpolation_row(point, value)
            for (_, point), value in zip(FIT_SAMPLES, fit_values)
        ]
    )
    assert matrix.rank() == 8
    kernel = matrix.nullspace()
    assert len(kernel) == 1
    coefficients = primitive_integer_vector(kernel[0])

    h0, h1, h2 = sp.symbols("h0 h1 h2")
    a0, a1, a2, b00, b01, b02, b11, b12, b22 = coefficients
    linear = a0 * h0 + a1 * h1 + a2 * h2
    quadratic = (
        b00 * h0**2
        + b01 * h0 * h1
        + b02 * h0 * h2
        + b11 * h1**2
        + b12 * h1 * h2
        + b22 * h2**2
    )
    obstruction = sp.expand(h0 * linear + sp.sqrt(41) * quadratic)
    norm = sp.factor(
        sp.expand(h0**2 * linear**2 - 41 * quadratic**2)
    )
    conic_matrix = sp.hessian(obstruction, (h0, h1, h2)) / 2
    conic_determinant = sp.expand(conic_matrix.det())
    conic_determinant_norm = sp.factor(
        sp.expand(
            conic_determinant
            * conic_determinant.xreplace({sp.sqrt(41): -sp.sqrt(41)})
        )
    )
    assert conic_determinant_norm != 0

    check_tag, check_point = CHECK_SAMPLE
    check_value = obstruction_ratio(check_tag)
    substitution = dict(zip((h0, h1, h2), check_point))
    predicted = sp.factor(
        (h0 * linear / quadratic).subs(substitution)
    )
    assert predicted == check_value

    infinity_discriminant = sp.factor(
        sp.discriminant(
            sp.Poly(quadratic.subs(h0, 0), h1),
            h1,
        )
    )
    linear_h0 = sp.solve(linear, h0)[0]
    linear_restriction_numerator = sp.factor(
        sp.together(
            quadratic.subs(h0, linear_h0)
        ).as_numer_denom()[0]
    )
    linear_discriminant = sp.factor(
        sp.discriminant(
            sp.Poly(linear_restriction_numerator, h1),
            h1,
        )
    )
    infinity_discriminant_coefficient = int(
        infinity_discriminant.subs(h2, 1)
    )
    linear_discriminant_coefficient = int(
        linear_discriminant.subs(h2, 1)
    )
    assert infinity_discriminant_coefficient < 0
    assert linear_discriminant_coefficient > 0
    assert not sp.integer_nthroot(
        linear_discriminant_coefficient, 2
    )[1]
    result = {
        "format": "two-pair-counterexample-fifth-factor-plane-v1",
        "slice": "h3=h4=0",
        "fit_samples": [
            {
                "tag": tag,
                "direction": [*point, 0, 0],
                "ratio": str(value),
            }
            for (tag, point), value in zip(FIT_SAMPLES, fit_values)
        ],
        "check_sample": {
            "tag": check_tag,
            "direction": [*check_point, 0, 0],
            "actual_ratio": str(check_value),
            "predicted_ratio": str(predicted),
            "exact_match": True,
        },
        "coefficient_vector": coefficients,
        "linear_factor": str(sp.factor(linear)),
        "quadratic_factor": str(sp.factor(quadratic)),
        "quadratic_field_obstruction": str(obstruction),
        "norm": str(norm),
        "norm_factorization_over_Q_sqrt_41": str(
            sp.factor(norm, extension=sp.sqrt(41))
        ),
        "quadratic_matrix_rank": (
            sp.hessian(quadratic, (h0, h1, h2)) / 2
        ).rank(),
        "quadratic_matrix_determinant": str(
            sp.factor(
                (
                    sp.hessian(quadratic, (h0, h1, h2)) / 2
                ).det()
            )
        ),
        "conic_matrix_determinant": str(conic_determinant),
        "conic_matrix_determinant_norm": str(conic_determinant_norm),
        "conic_is_smooth": True,
        "infinity_discriminant": str(infinity_discriminant),
        "linear_factor_discriminant": str(linear_discriminant),
        "reconstructed_conic_has_rational_point": False,
        "rational_point_certificate": (
            "For rational h, the reconstructed Theta=0 forces "
            "h0*A1=B2=0.  On h0=0 the binary restriction of B2 has "
            "negative discriminant; on A1=0 its positive discriminant is "
            "not a rational square."
        ),
        "conclusion": (
            "Within the tested quadratic projective ansatz, the selected "
            "augmented minor reconstructs as h0*A1+sqrt(41)*B2 on "
            "h3=h4=0 and matches an exact holdout sample.  The reconstructed "
            "smooth conic has no rational projective points.  A symbolic "
            "degree bound or universal identity, followed by checks of the "
            "other minors on the conic over Q(sqrt(41)), is still required "
            "for geometric uniform obstruction on the whole slice."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS fifth factor: quadratic reconstruction and ninth-point check")
    print(f"WROTE {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
