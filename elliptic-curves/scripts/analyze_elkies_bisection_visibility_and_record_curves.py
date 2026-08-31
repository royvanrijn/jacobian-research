#!/usr/bin/env python3
"""Resolve bisection visibility and run exact published-R17 j-recognition.

The specialization census already expresses every split bisection modulo the
generic seventeen in a deterministic public-complement basis.  This replay
row-reduces those classes over F_2 and records a canonical complement to the
visible subspace.

It also solves the exact specialization-recognition equations

    j_R17(t) = j_target

for the 2024 rank-29 record and ICARM curves 273, 302, and 398--400.  For
each equation, irreducibility of its primitive degree-24 polynomial is
certified by an irreducible reduction of the same degree over a finite field.
The historical rank-28 fibre is a positive control: its polynomial has the
expected rational factor 5471*t+9529 and an irreducible degree-23 cofactor.
The exclusion is unchanged by quadratic twisting because j is twist-invariant.
It does not exclude an isogenous curve, another elliptic fibration, or another
K3 family.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any, Sequence

import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_irreducible_p


ROOT = Path(__file__).resolve().parents[2]


PROTOCOL = "ELKIES2026BISECTIONVISIBILITY"
STATUS = "PASS_EXACT_BISECTION_VISIBILITY_AND_RECORD_PROVENANCE"

MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
TARGETS = ROOT / "elliptic-curves/data/elkies_2026_r17_j_recognition_targets.json"
SPECIALIZATIONS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_bisection_specialization_controls_v1.json"
)
PRIORITY = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-2026-bisection-equation-priority-full.json"
)
COLLISIONS = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-2026-equation-bisection-collisions-full-compact.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_bisection_visibility_record_curves_v1.json"
)
REPRODUCING_COMMAND = (
    ".venv/bin/python "
    "elliptic-curves/scripts/analyze_elkies_bisection_visibility_and_record_curves.py"
)

EXPECTED_WITNESS_PRIMES = {
    "rank29": 461,
    "curve273": 367,
    "curve302": 397,
    "curve398": 1009,
    "curve399": 83,
    "curve400": 157,
}
EXPECTED_CONTROL_COFACTOR_PRIME = 197


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rational_text(value: Fraction | int) -> str:
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def rows_rref_mod2(
    rows: Sequence[Sequence[int]], width: int
) -> tuple[list[list[int]], list[int]]:
    matrix = [[int(value) & 1 for value in row] for row in rows]
    if any(len(row) != width for row in matrix):
        raise ValueError("binary row has the wrong width")
    pivot_row = 0
    pivots: list[int] = []
    for column in range(width):
        selected = next(
            (index for index in range(pivot_row, len(matrix)) if matrix[index][column]),
            None,
        )
        if selected is None:
            continue
        matrix[pivot_row], matrix[selected] = matrix[selected], matrix[pivot_row]
        for index in range(len(matrix)):
            if index != pivot_row and matrix[index][column]:
                matrix[index] = [
                    left ^ right
                    for left, right in zip(matrix[index], matrix[pivot_row])
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    nonzero = [row for row in matrix if any(row)]
    return nonzero, pivots


def binary_rank(rows: Sequence[Sequence[int]], width: int) -> int:
    return len(rows_rref_mod2(rows, width)[1])


def unit_vector(width: int, index: int) -> list[int]:
    answer = [0] * width
    answer[index] = 1
    return answer


def visibility_record(fibre: dict[str, Any]) -> dict[str, Any]:
    public = fibre["public_complement"]
    width = int(public["dimension"])
    labels = list(public["ordered_basis_labels"])
    source_indices = list(public["source_point_indices_one_based"])
    rows = [
        list(hit["finite_quotient_class_modulo_generic_17"]["coordinates_over_f2"])
        for hit in fibre["hits"]
    ]
    rref, pivots = rows_rref_mod2(rows, width)
    nonpivots = [column for column in range(width) if column not in pivots]
    complement = [unit_vector(width, column) for column in nonpivots]
    if binary_rank([*rref, *complement], width) != width:
        raise AssertionError("visible rows and canonical complement do not span")
    if len(pivots) != fibre["split_class_span"]["dimension_modulo_generic_17"]:
        raise AssertionError("recomputed visibility rank changed")
    return {
        "parameter": fibre["parameter"],
        "known_complement_dimension": width,
        "split_bisection_count": fibre["split_bisection_count"],
        "visible_span_dimension": len(pivots),
        "invisible_quotient_dimension": width - len(pivots),
        "visible_rref_rows": rref,
        "visible_pivot_columns_zero_based": pivots,
        "canonical_complement_columns_zero_based": nonpivots,
        "canonical_complement_basis": complement,
        "canonical_complement_labels": [labels[column] for column in nonpivots],
        "canonical_complement_source_point_indices_one_based": [
            source_indices[column] for column in nonpivots
        ],
        "basis_convention": (
            "The canonical complement uses standard basis vectors in the "
            "nonpivot columns of the visible row space in the stored ordered "
            "public-complement basis. It is a deterministic complement, not "
            "an intrinsic splitting of the Mordell-Weil group."
        ),
    }


def primitive_positive(poly: sp.Poly) -> tuple[int, sp.Poly]:
    content, primitive = sp.polys.polytools.primitive(poly)
    content = int(content)
    if int(primitive.LC()) < 0:
        content = -content
        primitive = -primitive
    return content, primitive


def first_irreducible_reduction(poly: sp.Poly) -> tuple[int, list[int]]:
    coefficients = [int(value) for value in poly.all_coeffs()]
    for prime in sp.primerange(2, 5000):
        reduced = [value % prime for value in coefficients]
        if reduced[0] and gf_irreducible_p(reduced, int(prime), ZZ):
            return int(prime), reduced
    raise AssertionError("no finite-field irreducibility witness below 5000")


def weierstrass_invariants(ainvs: Sequence[str | int]) -> dict[str, int]:
    if len(ainvs) != 5:
        raise ValueError("a generalized Weierstrass model needs five a-invariants")
    a1, a2, a3, a4, a6 = map(int, ainvs)
    b2 = a1**2 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    discriminant = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    if c4**3 - c6**2 != 1728 * discriminant:
        raise AssertionError("generalized Weierstrass invariant identity failed")
    return {"c4": c4, "c6": c6, "discriminant": discriminant}


def target_invariants(targets: dict[str, Any]) -> dict[str, dict[str, Any]]:
    answer: dict[str, dict[str, Any]] = {}
    for target in targets["targets"]:
        invariants = weierstrass_invariants(target["ainvs"])
        answer[target["label"]] = {**target, **invariants}
    return answer


def recognition_polynomial(
    family_c4: sp.Expr,
    family_delta: sp.Expr,
    target_c4: int | sp.Rational,
    target_delta: int | sp.Rational,
    parameter: sp.Symbol,
) -> tuple[int, sp.Poly, Fraction]:
    j_rational = sp.Rational(target_c4) ** 3 / sp.Rational(target_delta)
    j_value = Fraction(int(sp.numer(j_rational)), int(sp.denom(j_rational)))
    equation = sp.Poly(
        sp.expand(family_c4**3 * j_value.denominator - j_value.numerator * family_delta),
        parameter,
        domain=sp.ZZ,
    )
    content, primitive = primitive_positive(equation)
    return content, primitive, j_value


def provenance_records(
    model: dict[str, Any], targets: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    parameter = sp.symbols("t")
    coefficient_a = sum(
        sp.Integer(value) * parameter**index
        for index, value in enumerate(model["A_coefficients_low_to_high"])
    )
    coefficient_b = sum(
        sp.Integer(value) * parameter**index
        for index, value in enumerate(model["B_coefficients_low_to_high"])
    )
    family_c4 = -48 * coefficient_a
    family_delta = -16 * (4 * coefficient_a**3 + 27 * coefficient_b**2)
    family_gcd = sp.gcd(sp.Poly(family_c4**3, parameter), sp.Poly(family_delta, parameter))
    if family_gcd.degree() != 0:
        raise AssertionError("published R17 j-map was not reduced")

    records: dict[str, Any] = {}
    for label, invariants in target_invariants(targets).items():
        c4 = invariants["c4"]
        discriminant = invariants["discriminant"]
        content, primitive, j_value = recognition_polynomial(
            family_c4, family_delta, c4, discriminant, parameter
        )
        witness_prime, reduced = first_irreducible_reduction(primitive)
        if witness_prime != EXPECTED_WITNESS_PRIMES[label]:
            raise AssertionError("irreducibility witness prime changed")
        if primitive.degree() != 24 or not primitive.is_irreducible:
            raise AssertionError("specialization-recognition polynomial changed")
        if not gf_irreducible_p(reduced, witness_prime, ZZ):
            raise AssertionError("stored modular witness is reducible")

        coefficients_low_to_high = [
            int(value) for value in reversed(primitive.all_coeffs())
        ]
        records[label] = {
            "icarm_id": invariants["icarm_id"],
            "source": invariants["source"],
            "source_sha256": invariants["source_sha256"],
            "ainvs": invariants["ainvs"],
            "c4": str(c4),
            "c6": str(invariants["c6"]),
            "discriminant": str(discriminant),
            "certified_rank_lower_bound": invariants["certified_rank_lower_bound"],
            "j_invariant": rational_text(j_value),
            "recognition_equation": "c4_R17(t)^3*Delta_target-c4_target^3*Delta_R17(t)=0",
            "primitive_polynomial_degree": primitive.degree(),
            "primitive_polynomial_content_removed": str(content),
            "primitive_coefficients_low_to_high": [
                str(value) for value in coefficients_low_to_high
            ],
            "primitive_coefficients_sha256": sha256(
                ",".join(str(value) for value in coefficients_low_to_high).encode()
            ).hexdigest(),
            "irreducible_mod_prime_witness": {
                "prime": witness_prime,
                "coefficients_high_to_low": reduced,
                "same_degree": len(reduced) - 1 == primitive.degree(),
                "irreducible": True,
            },
            "rational_affine_parameters": [],
            "parameter_at_infinity": False,
            "conclusion": (
                "No rational parameter in the published R17 chart has the "
                "same j-invariant, so this curve is not a rational fibre of "
                "that fibration, including after quadratic twisting."
            ),
        }

    control = targets["positive_control"]
    control_parameter = sp.Rational(control["parameter"])
    control_c4 = sp.Rational(family_c4.subs(parameter, control_parameter))
    control_delta = sp.Rational(family_delta.subs(parameter, control_parameter))
    content, primitive, control_j = recognition_polynomial(
        family_c4, family_delta, control_c4, control_delta, parameter
    )
    factors = sorted(
        sp.factor_list(primitive.as_expr())[1],
        key=lambda item: sp.degree(item[0], parameter),
    )
    if [(sp.degree(factor, parameter), exponent) for factor, exponent in factors] != [
        (1, 1),
        (23, 1),
    ]:
        raise AssertionError("rank-28 control factorization changed")
    linear = sp.Poly(factors[0][0], parameter, domain=sp.ZZ)
    cofactor = sp.Poly(factors[1][0], parameter, domain=sp.ZZ)
    if linear.as_expr() not in (5471 * parameter + 9529, -5471 * parameter - 9529):
        raise AssertionError("rank-28 control linear factor changed")
    cofactor_prime, cofactor_reduction = first_irreducible_reduction(cofactor)
    if cofactor_prime != EXPECTED_CONTROL_COFACTOR_PRIME:
        raise AssertionError("rank-28 control cofactor witness changed")
    control_record = {
        **control,
        "j_invariant": rational_text(control_j),
        "primitive_polynomial_degree": primitive.degree(),
        "primitive_polynomial_content_removed": str(content),
        "factor_degrees": [1, 23],
        "linear_factor": "5471*t+9529",
        "recovered_parameter": "-9529/5471",
        "irreducible_degree23_cofactor_mod_prime_witness": {
            "prime": cofactor_prime,
            "coefficients_high_to_low": cofactor_reduction,
            "irreducible": True,
        },
        "conclusion": "The exact positive control recovers its published rational parameter.",
    }

    family = {
        "equation": "y^2=x^3+A(t)*x+B(t)",
        "c4": "-48*A(t)",
        "discriminant": "-16*(4*A(t)^3+27*B(t)^2)",
        "j_map_reduced": True,
        "degrees_c4_discriminant": [
            sp.Poly(family_c4, parameter).degree(),
            sp.Poly(family_delta, parameter).degree(),
        ],
    }
    return family, control_record, records


def build_artifact() -> dict[str, Any]:
    model = json.loads(MODEL.read_text())
    targets = json.loads(TARGETS.read_text())
    specializations = json.loads(SPECIALIZATIONS.read_text())
    priority = json.loads(PRIORITY.read_text())
    collisions = json.loads(COLLISIONS.read_text())

    if priority["complete_enumeration"]["surviving_translation_orbits"] != 39120:
        raise AssertionError("bisection translation-orbit count changed")
    if priority["complete_enumeration"]["surviving_unoriented_norm_ten_representatives"] != 806238:
        raise AssertionError("norm-ten representative count changed")
    if collisions["lattice_orbit_coverage"]["status"] != "COMPLETE_EXACT_ORBIT_COVERAGE":
        raise AssertionError("equation atlas no longer has complete orbit coverage")
    if collisions["distinct_quadratic_extensions"] != 39120:
        raise AssertionError("quadratic-extension count changed")

    visibility = [visibility_record(fibre) for fibre in specializations["fibres"]]
    rank28 = next(item for item in visibility if item["parameter"] == "-9529/5471")
    curve394 = next(item for item in visibility if item["parameter"] == "3/8")
    if rank28["invisible_quotient_dimension"] != 10:
        raise AssertionError("rank-28 visibility deficit changed")
    if curve394["invisible_quotient_dimension"] != 0:
        raise AssertionError("curve 394 is no longer fully visible")

    family, control, records = provenance_records(model, targets)
    return {
        "schema": "elliptic-curves.elkies-2026-bisection-visibility-record-curves.v1",
        "status": STATUS,
        "claim": {
            "rank28": (
                "The complete bisection atlas sees a one-dimensional subspace "
                "of the eleven-dimensional known quotient; the stored ordered "
                "basis gives a canonical ten-dimensional complement."
            ),
            "curve394": (
                "The complete bisection atlas sees all four known directions "
                "beyond the generic seventeen."
            ),
            "record_curves": (
                "The exact degree-24 j-recognition polynomials for the 2024 "
                "rank-29 curve and ICARM curves 273, 302, and 398--400 are "
                "irreducible over Q, so none is a rational fibre of the "
                "published R17 fibration, including after quadratic twisting."
            ),
        },
        "bisection_mechanism_boundary": {
            "translation_orbits": 39120,
            "unoriented_norm_ten_representatives": 806238,
            "distinct_quadratic_extensions": 39120,
            "next_trace_shell_decision": "DO_NOT_CONSTRUCT_AS_A_NEW_BISECTION_SHELL",
            "reason": (
                "The certified atlas already quotients every section-nonnegative "
                "degree-two (-2)-curve by generic-section translation and sign. "
                "A higher-height translated trace represents the same bisection "
                "orbit, quadratic extension, and class modulo the generic "
                "seventeen; it cannot reveal one of the missing quotient classes."
            ),
            "remaining_search_scope": (
                "To target the ten missing rank-28 directions, change the "
                "geometric mechanism: use higher-degree multisections, residual "
                "Selmer covers, or direct quotient-targeted point searches."
            ),
        },
        "visibility": visibility,
        "published_r17_family": family,
        "rank28_positive_control": control,
        "record_curve_provenance": records,
        "claim_boundary": (
            "The visibility spaces are finite-quotient classes inside the known "
            "public complements, not Mordell-Weil upper bounds. The record-curve "
            "result excludes rational fibres only in the published R17 fibration. "
            "It is twist-stable because twisting preserves j, but it does not "
            "exclude isogenies, other elliptic fibrations, or a different family."
        ),
        "generation": {
            "command": REPRODUCING_COMMAND,
            "inputs": {
                relative(path): file_sha256(path)
                for path in (MODEL, TARGETS, SPECIALIZATIONS, PRIORITY, COLLISIONS)
            },
            "software": {
                "python": sys.version.split()[0],
                "sympy": sp.__version__,
            },
        },
    }


def canonical_json(data: dict[str, Any]) -> str:
    return json.dumps(data, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--check",
        action="store_true",
        help="recompute and compare with the pinned artifact",
    )
    args = parser.parse_args()

    artifact = build_artifact()
    serialized = canonical_json(artifact)
    if args.check:
        if not args.output.exists():
            raise SystemExit(f"missing pinned artifact: {args.output}")
        if args.output.read_text() != serialized:
            raise SystemExit("pinned visibility/provenance artifact is stale")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)

    rank28 = next(item for item in artifact["visibility"] if item["parameter"] == "-9529/5471")
    witnesses = artifact["record_curve_provenance"]
    print(
        f"{PROTOCOL}|rank28_visible={rank28['visible_span_dimension']}|"
        f"rank28_invisible={rank28['invisible_quotient_dimension']}|"
        f"rank29_modp={witnesses['rank29']['irreducible_mod_prime_witness']['prime']}|"
        f"curve273_modp={witnesses['curve273']['irreducible_mod_prime_witness']['prime']}|"
        f"curve302_modp={witnesses['curve302']['irreducible_mod_prime_witness']['prime']}|"
        f"curve398_modp={witnesses['curve398']['irreducible_mod_prime_witness']['prime']}|"
        f"control_factor={artifact['rank28_positive_control']['linear_factor']}|"
        f"status=PASS"
    )


if __name__ == "__main__":
    main()
