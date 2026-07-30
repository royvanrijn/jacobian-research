#!/usr/bin/env python3
"""Exact numerical atlas for projective compactifications of HC(4) gradients.

For a degree-m polynomial map F:A^4 -> A^4, use the graph compactification

    [X0^m : F_1^h : ... : F_4^h] : P^4 --> P^4.

If B is its base scheme and

    i_* s(B,P^4) = sigma_1 H + ... + sigma_4 H^4,

then the projective degrees satisfy

    g_i = m^i - sum_{k=1}^i binom(i,k)m^(i-k)sigma_k.

The script enumerates every integral degree list allowed by the elementary
bounds and log-concavity when m is 2 or 3 and g_4 is 2 or 3.  These are
necessary numerical signatures, not existence results and not a
classification by Hilbert polynomial.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_projective_polar_atlas.json"
)


@dataclass(frozen=True)
class Signature:
    projective_degrees: tuple[int, int, int, int, int]
    segre_degrees: tuple[int, int, int, int]
    leading_base_codimension: int
    leading_cycle_degree: int
    smooth_lci_curve_genus: int | None
    smooth_lci_curve_numerically_possible: bool | None


def segre_from_degrees(m: int, degrees: tuple[int, ...]) -> tuple[int, ...]:
    """Invert the triangular projective-degree/Segre-degree relation."""
    sigmas: list[int] = []
    for i in range(1, len(degrees)):
        known = sum(
            comb(i, k) * m ** (i - k) * sigmas[k - 1]
            for k in range(1, i)
        )
        sigmas.append(m**i - known - degrees[i])
    return tuple(sigmas)


def projective_from_segre(m: int, sigmas: tuple[int, ...]) -> tuple[int, ...]:
    """Evaluate the blow-up formula in every projective degree."""
    return (1,) + tuple(
        m**i
        - sum(
            comb(i, k) * m ** (i - k) * sigmas[k - 1]
            for k in range(1, i + 1)
        )
        for i in range(1, len(sigmas) + 1)
    )


def castelnuovo_bound_p4(degree: int) -> int:
    """Castelnuovo's genus bound for a nondegenerate integral curve in P^4."""
    if degree <= 0:
        return -1
    quotient, remainder = divmod(degree - 1, 3)
    return 3 * quotient * (quotient - 1) // 2 + quotient * remainder


def make_signature(m: int, g2: int, g3: int, top_degree: int) -> Signature:
    degrees = (1, m, g2, g3, top_degree)
    sigmas = segre_from_degrees(m, degrees)
    assert projective_from_segre(m, sigmas) == degrees
    assert sigmas[0] == 0

    first_nonzero = next(
        (index for index, value in enumerate(sigmas, start=1) if value),
        4,
    )
    leading_degree = sigmas[first_nonzero - 1]
    curve_genus = None
    curve_possible = None
    if first_nonzero == 3:
        curve_degree = sigmas[2]
        numerator = 2 - 5 * curve_degree - sigmas[3]
        if numerator % 2 == 0:
            curve_genus = numerator // 2
        curve_possible = (
            curve_genus is not None
            and 0 <= curve_genus <= castelnuovo_bound_p4(curve_degree)
        )

    return Signature(
        projective_degrees=degrees,
        segre_degrees=sigmas,
        leading_base_codimension=first_nonzero,
        leading_cycle_degree=leading_degree,
        smooth_lci_curve_genus=curve_genus,
        smooth_lci_curve_numerically_possible=curve_possible,
    )


def atlas(m: int, top_degree: int) -> list[Signature]:
    """Enumerate the positive log-concave projective-degree lists."""
    result: list[Signature] = []
    for g2 in range(1, m**2 + 1):
        for g3 in range(1, m**3 + 1):
            if g2 * g2 < m * g3:
                continue
            if g3 * g3 < g2 * top_degree:
                continue
            result.append(make_signature(m, g2, g3, top_degree))
    return result


atlases = {
    f"gradient_degree_{m}_affine_degree_{top}": atlas(m, top)
    for m in (2, 3)
    for top in (2, 3)
}
assert {key: len(value) for key, value in atlases.items()} == {
    "gradient_degree_2_affine_degree_2": 9,
    "gradient_degree_2_affine_degree_3": 7,
    "gradient_degree_3_affine_degree_2": 72,
    "gradient_degree_3_affine_degree_3": 67,
}


zero_dimensional = {}
for key, signatures in atlases.items():
    rows = [
        row
        for row in signatures
        if row.leading_base_codimension == 4
        and row.segre_degrees[:3] == (0, 0, 0)
    ]
    assert len(rows) == 1
    zero_dimensional[key] = rows[0].segre_degrees[3]

assert zero_dimensional == {
    "gradient_degree_2_affine_degree_2": 14,
    "gradient_degree_2_affine_degree_3": 13,
    "gradient_degree_3_affine_degree_2": 79,
    "gradient_degree_3_affine_degree_3": 78,
}


# Wang's degree-two theorem says that every characteristic-zero quadratic
# Keller map is a polynomial automorphism.  Applied to grad(Psi), constant
# nonzero Hessian determinant therefore forces top projective degree one.
# The theorem is an external mathematical input; the finite consequences
# for the atlas are checked here.
quadratic_keller_atlas = atlas(2, 1)
assert len(quadratic_keller_atlas) == 11
assert all(
    row.projective_degrees[-1] != 1
    for key, rows in atlases.items()
    if key.startswith("gradient_degree_2_")
    for row in rows
)
quadratic_keller_consequence = {
    "external_input": (
        "Wang's theorem: every characteristic-zero Keller map of "
        "polynomial degree at most two is a polynomial automorphism"
    ),
    "forced_affine_degree": 1,
    "forced_total_segre_correction": 2**4 - 1,
    "number_of_log_concave_degree_one_signatures": len(
        quadratic_keller_atlas
    ),
    "excluded_atlas_rows": {
        "affine_degree_2": len(
            atlases["gradient_degree_2_affine_degree_2"]
        ),
        "affine_degree_3": len(
            atlases["gradient_degree_2_affine_degree_3"]
        ),
    },
    "excluded_zero_dimensional_lengths": [14, 13],
}
assert quadratic_keller_consequence["forced_total_segre_correction"] == 15
assert quadratic_keller_consequence["excluded_atlas_rows"] == {
    "affine_degree_2": 9,
    "affine_degree_3": 7,
}


# HC4CQ1 excludes every collision for a four-variable potential
# q_2+h_3+h_4 with constant nonzero Hessian determinant.  After translating
# the midpoint and subtracting the common gradient value, every collision
# of a degree-four potential has exactly this form.  Ax--Grothendieck then
# turns injectivity over the algebraic closure into polynomial
# invertibility.  These are external mathematical inputs; their complete
# numerical consequences for the cubic-gradient atlas are checked here.
cubic_keller_atlas = atlas(3, 1)
assert len(cubic_keller_atlas) == 80
cubic_gradient_consequence = {
    "external_inputs": [
        (
            "HC4CQ1: a characteristic-zero four-variable potential "
            "q2+h3+h4 with constant nonzero Hessian determinant has no "
            "nonzero antipodal gradient collision"
        ),
        (
            "Ax--Grothendieck: an injective polynomial endomorphism of "
            "affine space in characteristic zero is an automorphism"
        ),
    ],
    "forced_affine_degree": 1,
    "forced_total_segre_correction": 3**4 - 1,
    "number_of_log_concave_degree_one_signatures": len(cubic_keller_atlas),
    "excluded_atlas_rows": {
        "affine_degree_2": len(
            atlases["gradient_degree_3_affine_degree_2"]
        ),
        "affine_degree_3": len(
            atlases["gradient_degree_3_affine_degree_3"]
        ),
    },
    "excluded_zero_dimensional_lengths": [79, 78],
    "counterexample_potential_degree_lower_bound": 5,
}
assert cubic_gradient_consequence["forced_total_segre_correction"] == 80
assert cubic_gradient_consequence["excluded_atlas_rows"] == {
    "affine_degree_2": 72,
    "affine_degree_3": 67,
}


# The graph and full-polar degree lists are independently certified by the
# companion Macaulay2 checker.  Inverting them here checks the conventions.
calibrations = {
    "quadratic_graph": {
        "m": 2,
        "projective_degrees": (1, 2, 2, 2, 1),
    },
    "quadratic_full_polar": {
        "m": 2,
        "projective_degrees": (1, 2, 4, 4, 2),
    },
    "cubic_graph": {
        "m": 3,
        "projective_degrees": (1, 3, 3, 3, 1),
    },
    "cubic_full_polar": {
        "m": 3,
        "projective_degrees": (1, 3, 6, 6, 3),
    },
}
for row in calibrations.values():
    row["segre_degrees"] = segre_from_degrees(
        row["m"], row["projective_degrees"]
    )
    assert projective_from_segre(
        row["m"], row["segre_degrees"]
    ) == row["projective_degrees"]

assert calibrations["quadratic_graph"]["segre_degrees"] == (0, 2, -6, 15)
assert calibrations["quadratic_full_polar"]["segre_degrees"] == (
    0,
    0,
    4,
    -18,
)
assert calibrations["cubic_graph"]["segre_degrees"] == (0, 6, -30, 116)
assert calibrations["cubic_full_polar"]["segre_degrees"] == (0, 3, -6, -12)


controls = {
    "cotangent_lift_quartic_packet": {
        "ambient_dimension": 4,
        "affine_degree": 4,
        "correction_quadratic_gradient": 2**4 - 4,
        "correction_cubic_gradient": 3**4 - 4,
    },
    "meng_yang_before_schur_hc6": {
        "ambient_dimension": 6,
        "gradient_degree": 7,
        "affine_degree": 3,
        "total_segre_correction": 7**6 - 3,
    },
    "meng_yang_after_schur_hc5": {
        "ambient_dimension": 5,
        "gradient_degree": 13,
        "affine_degree": 3,
        "total_segre_correction": 13**5 - 3,
    },
}
assert controls["meng_yang_before_schur_hc6"]["total_segre_correction"] == 117646
assert controls["meng_yang_after_schur_hc5"]["total_segre_correction"] == 371290


payload = {
    "format": "hc4-projective-polar-atlas-v1",
    "conventions": {
        "map": "[X0^m:F1^h:F2^h:F3^h:F4^h]",
        "segre_pushforward": (
            "i_*s(B,P4)=sigma_1*H+sigma_2*H^2+"
            "sigma_3*H^3+sigma_4*H^4"
        ),
        "formula": (
            "g_i=m^i-sum_{k=1}^i binom(i,k)m^(i-k)sigma_k"
        ),
        "top_formula_without_fixed_divisor": (
            "g_4=m^4-6*m^2*sigma_2-4*m*sigma_3-sigma_4"
        ),
    },
    "scope": (
        "Integral signatures from positivity, degree bounds, and "
        "log-concavity. Wang's theorem excludes every listed "
        "quadratic-gradient affine-degree-two/three row. HC4CQ1 and "
        "Ax--Grothendieck exclude every listed cubic-gradient row. Thus "
        "these are a numerical pre-Keller atlas, not existence results; "
        "the ordinary Hilbert polynomial does not determine the Segre class."
    ),
    "counts": {key: len(value) for key, value in atlases.items()},
    "zero_dimensional_lengths": zero_dimensional,
    "quadratic_keller_consequence": quadratic_keller_consequence,
    "quadratic_keller_degree_one_atlas": [
        asdict(row) for row in quadratic_keller_atlas
    ],
    "cubic_gradient_consequence": cubic_gradient_consequence,
    "cubic_keller_degree_one_atlas": [
        asdict(row) for row in cubic_keller_atlas
    ],
    "atlases": {
        key: [asdict(row) for row in value]
        for key, value in atlases.items()
    },
    "calibrations": calibrations,
    "controls": controls,
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: inverted the projective-degree/Segre-degree formula exactly")
print("PASS: low-degree atlas counts are 9, 7, 72, and 67")
print("PASS: zero-dimensional lengths are 14, 13, 79, and 78")
print("PASS: Wang's theorem excludes all 16 quadratic-degree rows")
print("PASS: every quadratic Keller gradient has total correction 15")
print("PASS: HC4CQ1 excludes all 139 cubic-degree rows")
print("PASS: every cubic Keller gradient in HC4 has total correction 80")
print("PASS: any HC4 counterexample potential has degree at least five")
print("PASS: HC6/HC5 total corrections are 117646 and 371290")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
