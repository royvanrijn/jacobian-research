#!/usr/bin/env python3
"""Exactly recognize ICARM 302, 351, and 356 in the rootless R17 family.

For each target, factor the degree-24 projective j-matching equation over Q,
retain all rational points of P^1(Q), test nonsingularity and twist class, and
transport the seventeen generic sections whenever a Q-isomorphism exists.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import gcd, isqrt, lcm
from pathlib import Path
import sys
import time

import sympy


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
MODEL = LOCAL / "q12o5867-smooth-rr-qq.json"
SECTIONS = LOCAL / "q12o5867-rootless-selected-basis-qq.json"
EC_CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(EC_CAS))

import icarm_curve302  # noqa: E402
import icarm_curve356  # noqa: E402


TARGETS = {
    302: {
        "ainvs": tuple(icarm_curve302.GENERAL_WEIERSTRASS_COEFFICIENTS),
        "source": "https://elliptic-rank.icarm.cloud/curve/302.json",
        "source_sha256": "2bd422f34d0566f79c53f6bc94c6ec5f9ff048cd88fcc300f53e0b397af2e67b",
        "rank_lower_bound": 31,
    },
    351: {
        "ainvs": tuple(map(Q, (
            1,
            -1,
            1,
            -250918915934128421340307896120883808444030,
            45913507617476434834864347970127139550348536547695769941315597,
        ))),
        "source": "https://elliptic-rank.icarm.cloud/curve/351.json",
        "source_sha256": "02c0de1801d0c925dd6e42204f8461e99595926e95221e91da0c09466a6f67fd",
        "rank_lower_bound": 25,
    },
    356: {
        "ainvs": tuple(icarm_curve356.GENERAL_WEIERSTRASS_COEFFICIENTS),
        "source": "https://elliptic-rank.icarm.cloud/curve/356.json",
        "source_sha256": "58afbc62dbb6e01b47266c90edcf0e09bb003bb6a558333422b332e42546e89e",
        "rank_lower_bound": 29,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=(GENERATED / "elkies-k3-h3-q12o5867-icarm-302-351-356-recognition.json"),
    )
    parser.add_argument("--check", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational_text(value: Fraction | int) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def hash_strings(values: list[Fraction | int]) -> str:
    raw = json.dumps(list(map(str, values)), separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def trim(poly: list[Fraction | int]) -> list[Fraction]:
    answer = list(map(Q, poly)) or [Q(0)]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def poly_add(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    answer = [Q(0)] * max(len(left), len(right))
    for index, value in enumerate(left):
        answer[index] += value
    for index, value in enumerate(right):
        answer[index] += value
    return trim(answer)


def poly_scale(poly: list[Fraction], scalar: Fraction | int) -> list[Fraction]:
    return trim([Q(scalar) * value for value in poly])


def poly_mul(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    answer = [Q(0)] * (len(left) + len(right) - 1)
    for i, x_value in enumerate(left):
        for j, y_value in enumerate(right):
            answer[i+j] += x_value * y_value
    return trim(answer)


def poly_pow(poly: list[Fraction], exponent: int) -> list[Fraction]:
    answer = [Q(1)]
    base = poly
    while exponent:
        if exponent & 1:
            answer = poly_mul(answer, base)
        base = poly_mul(base, base)
        exponent >>= 1
    return answer


def poly_div_exact(numerator: list[Fraction], denominator: list[Fraction]) -> list[Fraction]:
    numerator = trim(numerator)
    denominator = trim(denominator)
    if denominator == [0]:
        raise ZeroDivisionError
    quotient = [Q(0)] * max(1, len(numerator)-len(denominator)+1)
    while numerator != [0] and len(numerator) >= len(denominator):
        shift = len(numerator)-len(denominator)
        coefficient = numerator[-1]/denominator[-1]
        quotient[shift] = coefficient
        subtract = [Q(0)]*shift + [coefficient*value for value in denominator]
        numerator = poly_add(numerator, poly_scale(subtract, -1))
    if numerator != [0]:
        raise ArithmeticError("inexact polynomial division")
    return trim(quotient)


def to_sympy(poly: list[Fraction], variable: sympy.Symbol, domain: str = "QQ") -> sympy.Poly:
    coefficients = [sympy.Rational(value.numerator, value.denominator) for value in reversed(trim(poly))]
    return sympy.Poly.from_list(coefficients, gens=variable, domain=domain)


def from_sympy(poly: sympy.Poly) -> list[Fraction]:
    return trim([Q(int(value.p), int(value.q)) for value in reversed(poly.all_coeffs())])


def primitive_integer_polynomial(poly: list[Fraction]) -> list[int]:
    denominator = 1
    for value in poly:
        denominator = lcm(denominator, value.denominator)
    coefficients = [value.numerator * (denominator // value.denominator) for value in poly]
    content = 0
    for value in coefficients:
        content = gcd(content, abs(value))
    coefficients = [value // content for value in coefficients]
    if coefficients[-1] < 0:
        coefficients = [-value for value in coefficients]
    return coefficients


def primitive_integer_pair(
    numerator: list[Fraction], denominator: list[Fraction]
) -> tuple[list[int], list[int]]:
    values = numerator + denominator
    common_denominator = 1
    for value in values:
        common_denominator = lcm(common_denominator, value.denominator)
    integers = [value.numerator*(common_denominator//value.denominator) for value in values]
    content = 0
    for value in integers:
        content = gcd(content, abs(value))
    integers = [value//content for value in integers]
    numerator_ZZ = integers[:len(numerator)]
    denominator_ZZ = integers[len(numerator):]
    if denominator_ZZ[-1] < 0:
        numerator_ZZ = [-value for value in numerator_ZZ]
        denominator_ZZ = [-value for value in denominator_ZZ]
    return numerator_ZZ, denominator_ZZ


def polynomial_record(poly: list[Fraction | int], include_coefficients: bool = False) -> dict[str, object]:
    poly = trim(poly)
    record: dict[str, object] = {
        "degree": len(poly)-1,
        "coefficient_sha256_low_to_high": hash_strings(poly),
        "maximum_coefficient_bits": max(
            max(abs(value.numerator).bit_length(), value.denominator.bit_length()) for value in poly
        ),
    }
    if include_coefficients:
        record["coefficients_low_to_high"] = list(map(str, poly))
    return record


def generalized_invariants(ainvs: tuple[Fraction, ...]) -> dict[str, Fraction]:
    a1, a2, a3, a4, a6 = map(Q, ainvs)
    b2 = a1*a1+4*a2
    b4 = 2*a4+a1*a3
    b6 = a3*a3+4*a6
    b8 = a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
    c4 = b2*b2-24*b4
    c6 = -b2**3+36*b2*b4-216*b6
    discriminant = -b2*b2*b8-8*b4**3-27*b6**2+9*b2*b4*b6
    return {
        "b2": b2,
        "c4": c4,
        "c6": c6,
        "discriminant": discriminant,
        "j": c4**3/discriminant,
        "short_A": -c4/Q(48),
        "short_B": -c6/Q(864),
    }


def projective_value(poly: list[Fraction], weight: int, a: int, b: int) -> Fraction:
    return sum((value*a**i*b**(weight-i) for i, value in enumerate(poly)), Q(0))


def rational_square_root(value: Fraction) -> Fraction | None:
    if value < 0:
        return None
    numerator_root = isqrt(value.numerator)
    denominator_root = isqrt(value.denominator)
    if numerator_root**2 != value.numerator or denominator_root**2 != value.denominator:
        return None
    return Q(numerator_root, denominator_root)


def target_generalized_point(
    point: tuple[Fraction, Fraction], ainvs: tuple[Fraction, ...]
) -> tuple[Fraction, Fraction]:
    X, Y = point
    a1, a2, a3, _a4, _a6 = map(Q, ainvs)
    x_value = X-(a1*a1+4*a2)/12
    y_value = Y-(a1*x_value+a3)/2
    return x_value, y_value


def on_generalized_curve(
    point: tuple[Fraction, Fraction], ainvs: tuple[Fraction, ...]
) -> bool:
    x_value, y_value = point
    a1, a2, a3, a4, a6 = map(Q, ainvs)
    return (
        y_value*y_value+a1*x_value*y_value+a3*y_value
        == x_value**3+a2*x_value*x_value+a4*x_value+a6
    )


def main() -> None:
    args = parse_args()
    started = time.monotonic()
    model_raw = MODEL.read_bytes()
    sections_raw = SECTIONS.read_bytes()
    model = json.loads(model_raw)
    lifted = json.loads(sections_raw)
    assert model["status"] == "PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN"
    assert model["child"]["root_rank"] == 0
    assert lifted["status"] == "PASS_EXACT_QQ_Q12O5867_ROOTLESS_17_SELECTED_SECTIONS"
    assert lifted["lifted_section_count"] == 17

    A = trim(list(map(Q, model["child"]["minimal_A_coefficients_low_to_high"])))
    B = trim(list(map(Q, model["child"]["minimal_B_coefficients_low_to_high"])))
    assert len(A) == 9 and len(B) == 13
    section_points = []
    for row in lifted["sections"]:
        section = row["section"]
        X = trim(list(map(Q, section["x_coefficients_low_to_high"])))
        Y = trim(list(map(Q, section["y_coefficients_low_to_high"])))
        assert len(X) <= 5 and len(Y) <= 7
        assert poly_pow(Y, 2) == poly_add(poly_add(poly_pow(X, 3), poly_mul(A, X)), B)
        section_points.append((X, Y))
    assert len(section_points) == 17

    c4 = poly_scale(A, -48)
    c6 = poly_scale(B, -864)
    discriminant = poly_scale(poly_add(poly_scale(poly_pow(A, 3), 4), poly_scale(poly_pow(B, 2), 27)), -16)
    variable = sympy.symbols("u")
    common = from_sympy(to_sympy(poly_pow(c4, 3), variable).gcd(to_sympy(discriminant, variable)))
    j_numerator = poly_div_exact(poly_pow(c4, 3), common)
    j_denominator = poly_div_exact(discriminant, common)
    j_num_ZZ, j_den_ZZ = primitive_integer_pair(j_numerator, j_denominator)

    target_records = []
    for curve_id in (302, 356, 351):
        target = TARGETS[curve_id]
        invariants = generalized_invariants(target["ainvs"])
        assert invariants["discriminant"] != 0
        j_target = invariants["j"]
        matching = poly_add(
            poly_scale(poly_pow(c4, 3), j_target.denominator),
            poly_scale(discriminant, -j_target.numerator),
        )
        primitive = primitive_integer_polynomial(matching)
        primitive_sympy = sympy.Poly.from_list(list(reversed(primitive)), gens=variable, domain="ZZ")
        factor_unit, factorization = sympy.factor_list(primitive_sympy)
        product_check = sympy.Poly(factor_unit, variable, domain="ZZ")
        factor_records = []
        finite_roots: list[tuple[Fraction, int]] = []
        for factor, multiplicity in factorization:
            factor_poly = from_sympy(factor)
            product_check *= factor**multiplicity
            factor_record = polynomial_record(factor_poly)
            factor_record["multiplicity"] = multiplicity
            if len(factor_poly) == 2:
                root = -factor_poly[0]/factor_poly[1]
                finite_roots.append((root, multiplicity))
                factor_record["rational_root"] = rational_text(root)
            factor_records.append(factor_record)
        assert product_check == primitive_sympy
        infinity_multiplicity = 24-(len(primitive)-1)
        assert infinity_multiplicity >= 0
        roots = [(root.numerator, root.denominator, multiplicity) for root, multiplicity in finite_roots]
        if infinity_multiplicity:
            roots.append((1, 0, infinity_multiplicity))

        solution_records = []
        for a, b, multiplicity in roots:
            specialized_A = projective_value(A, 8, a, b)
            specialized_B = projective_value(B, 12, a, b)
            specialized_sections = [
                (projective_value(X, 4, a, b), projective_value(Y, 6, a, b))
                for X, Y in section_points
            ]
            assert all(
                y_value**2 == x_value**3+specialized_A*x_value+specialized_B
                for x_value, y_value in specialized_sections
            )
            specialized_delta = -16*(4*specialized_A**3+27*specialized_B**2)
            record: dict[str, object] = {
                "parameter_projective": [str(a), str(b)],
                "parameter": "infinity" if b == 0 else rational_text(Q(a, b)),
                "root_multiplicity": multiplicity,
                "specialized_A": rational_text(specialized_A),
                "specialized_B": rational_text(specialized_B),
                "specialized_discriminant": rational_text(specialized_delta),
                "nonsingular": specialized_delta != 0,
                "generic_sections_specialized_and_verified": len(specialized_sections),
            }
            if specialized_delta:
                target_A = invariants["short_A"]
                target_B = invariants["short_B"]
                assert specialized_A and specialized_B and target_A and target_B
                twist_parameter = (target_B/specialized_B)/(target_A/specialized_A)
                assert target_A == twist_parameter**2*specialized_A
                assert target_B == twist_parameter**3*specialized_B
                scale = rational_square_root(twist_parameter)
                record["twist_parameter_in_Qmod_squares"] = rational_text(twist_parameter)
                record["q_isomorphic"] = scale is not None
                if scale is not None:
                    short_points = [(scale**2*x, scale**3*y) for x, y in specialized_sections]
                    assert all(y*y == x**3+target_A*x+target_B for x, y in short_points)
                    general_points = [target_generalized_point(point, target["ainvs"]) for point in short_points]
                    assert all(on_generalized_curve(point, target["ainvs"]) for point in general_points)
                    record["isomorphism_family_to_target_short"] = {
                        "x_target_over_x_family": rational_text(scale**2),
                        "y_target_over_y_family": rational_text(scale**3),
                    }
                    record["transported_section_count"] = len(general_points)
                    record["transported_generalized_points"] = [
                        [rational_text(x), rational_text(y)] for x, y in general_points
                    ]
            solution_records.append(record)

        target_records.append({
            "curve_id": curve_id,
            "public_source": target["source"],
            "public_source_sha256": target["source_sha256"],
            "public_rank_lower_bound": target["rank_lower_bound"],
            "ainvs": list(map(rational_text, target["ainvs"])),
            "c4": rational_text(invariants["c4"]),
            "c6": rational_text(invariants["c6"]),
            "discriminant": rational_text(invariants["discriminant"]),
            "j": rational_text(j_target),
            "projective_matching_equation": "den(j_target)*c4_R17(a,b)^3-num(j_target)*Delta_R17(a,b)=0",
            "primitive_affine_polynomial": polynomial_record(list(map(Q, primitive))),
            "factorization_over_Q": factor_records,
            "primitive_affine_polynomial_irreducible_over_Q": (
                len(factorization) == 1
                and factorization[0][0].degree() == 24
                and factorization[0][1] == 1
            ),
            "infinity_root_multiplicity": infinity_multiplicity,
            "rational_projective_solution_count": len(solution_records),
            "rational_projective_solutions": solution_records,
        })

    status = (
        "PASS_EXACT_QQ_ICARM_FAMILY_RECOGNITION_MATCHES_FOUND"
        if any(row["rational_projective_solution_count"] for row in target_records)
        else "PASS_EXACT_QQ_ICARM_FAMILY_RECOGNITION_NO_RATIONAL_MATCHES"
    )
    payload = {
        "schema": "elkies-k3.h92-q12o5867-icarm-specialization-recognition.v1",
        "status": status,
        "inputs": {
            str(MODEL.relative_to(ROOT)): hashlib.sha256(model_raw).hexdigest(),
            str(SECTIONS.relative_to(ROOT)): hashlib.sha256(sections_raw).hexdigest(),
            str((EC_CAS / "icarm_curve302.py").relative_to(ROOT)): sha256(EC_CAS / "icarm_curve302.py"),
            str((EC_CAS / "icarm_curve356.py").relative_to(ROOT)): sha256(EC_CAS / "icarm_curve356.py"),
        },
        "software": {"python": sys.version.split()[0], "sympy": sympy.__version__},
        "rootless_family": {
            "equation": "y^2=x^3+A(u)*x+B(u)",
            "A": polynomial_record(A),
            "B": polynomial_record(B),
            "section_count": len(section_points),
            "all_generic_section_identities_verified": True,
            "projective_infinity_fibre": {
                "A": rational_text(A[8]),
                "B": rational_text(B[12]),
                "discriminant": rational_text(-16*(4*A[8]**3+27*B[12]**2)),
                "nonsingular": -16*(4*A[8]**3+27*B[12]**2) != 0,
                "j": rational_text(
                    (-48*A[8])**3/(-16*(4*A[8]**3+27*B[12]**2))
                ),
            },
            "reduced_j": {
                "formula": "c4(u)^3/Delta(u), c4=-48*A, Delta=-16*(4*A^3+27*B^2)",
                "cancelled_gcd_degree": len(common)-1,
                "primitive_integer_numerator": polynomial_record(list(map(Q, j_num_ZZ)), True),
                "primitive_integer_denominator": polynomial_record(list(map(Q, j_den_ZZ)), True),
            },
        },
        "targets": target_records,
        "proof_boundary": (
            "Exact projective j-recognition, nonsingularity, twist/isomorphism, and section "
            "transport for every rational solution. Absence of a rational projective root "
            "proves only that the named target is not a Q-rational specialization of this "
            "fixed parameterized family; it makes no statement about other families."
        ),
        "reproducing_command": (
            ".venv/bin/python elkies-k3/scripts/"
            "recognize_h92_q12o5867_icarm_specializations_qq.py"
        ),
        "runtime_seconds": time.monotonic()-started,
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True)+"\n"
    if args.check:
        if not args.output.exists():
            raise SystemExit(f"missing recognition artifact: {args.output}")
        existing = json.loads(args.output.read_text())
        existing.pop("runtime_seconds", None)
        payload.pop("runtime_seconds", None)
        if existing != payload:
            raise SystemExit("stale ICARM specialization-recognition artifact")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(
        "ICARMR17|"
        + "|".join(
            f"curve{row['curve_id']}={row['rational_projective_solution_count']}"
            for row in target_records
        )
        + f"|status={status}"
    )


if __name__ == "__main__":
    main()
