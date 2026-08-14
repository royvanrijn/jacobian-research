#!/usr/bin/env python3
"""Bounded exact rational-point sieve on the Fermigier degree-32 component.

For the pinned ``P13 x R20E1`` bidegree-(2,1) pencil, write ``D(c,k)`` for
the irreducible nonlinear factor of the branch discriminant.  This script
searches the declared projective box

    (C:K:D),  D >= 0,  gcd(C,K,D)=1,
    max(|C|,|K|,D) <= H,

where affine parameters are ``c=C/D, k=K/D`` and ``D>0`` is the unique sign
normalization.  Five fixed modular filters determine the only possible ``K``
in the box for each ``(C,D)`` before the primitive integral homogenization is
evaluated exactly.  The projective boundary and every rational intersection
with the five already-known discriminant lines are classified separately.

This is honest bounded negative coverage, not a proof that the nonlinear
component has no other rational points.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from functools import reduce
import hashlib
import itertools
import json
from math import gcd, lcm, prod
from pathlib import Path
import platform
import re
import shutil
import subprocess
from typing import Any, Iterable

import sympy as sp

from analyze_fermigier_bidegree21_pilot import (
    AFFINE_CANCELLATION_SLOPE,
    E22_T,
    RANK20_T,
    pencil_polynomial,
)
from analyze_fermigier_exceptional_transport import (
    factor_signature,
    hyperelliptic_genus,
    rational_text,
    squareclass_kernel_degree,
)


DEFAULT_HEIGHT = 1024
SIEVE_PRIMES = (17, 19, 23, 29, 31)
TOTAL_DEGREE = 32
PARAMETER_DEGREE = 4
EXPECTED_TERM_COUNT = 561
EXPECTED_PRIMITIVE_COEFFICIENT_SHA256 = (
    "5c60ea4247ddc7eb99f1cc6726c592e569be3ee7c1de0b74892e0bad252d6eda"
)
UPSTREAM_FACTOR_SHA256 = (
    "8c8c84159171629f514bd52f00091b61c0e8a3b806273765c35fc35b032d8799"
)
OUTPUT_RELATIVE = Path(
    "artifacts/generated-results/"
    "elliptic_fermigier_bidegree21_p13_r20e1_nonlinear_points_h1024.json"
)


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def result_digest(record: dict[str, Any]) -> str:
    stable = dict(record)
    stable.pop("generated_at_utc", None)
    stable.pop("result_sha256", None)
    return sha256_text(json.dumps(stable, sort_keys=True, separators=(",", ":")))


def _singular_program() -> str:
    """Return the exact extraction program used by the upstream pilot."""

    return r'''
ring r=0,(c,k,T),dp;
number t1=39508/39;
number t2=28917/10;
number x1=256714/39;
number x2=-8545;
poly a=(x1*(c*t1+1)-x2*(c*t2+1))/(t1-t2);
poly b=x1*(c*t1+1)-a*t1;
poly D=c*T+1;
poly N=a*T+b+k*(T-t1)*(T-t2);
poly f0=18103855887324900+102302344648*T^2-879500*T^4+T^6;
poly f1=-257843832010380-650709150*T^2+1860*T^4;
poly f2=1195214262641+1718550*T^2-2*T^4;
poly f3=-2051321790-1860*T^2;
poly f4=1149050+T^2;
poly P=f0*D^4+f1*N*D^3+f2*N^2*D^2+f3*N^3*D+f4*N^4;
poly rr=resultant(P,diff(P,T),T);
poly dd=rr/((c-k)^2*(c+k)^2);
poly known=(k-c)^2*(k+c)^2*(28917c+10)^12*(39508c+39)^12
           *(5899690c+732683k)^12;
poly nonlinear=dd/known;
if (dd-known*nonlinear!=0) { "@@ERROR_DIVISION"; }
"@@TOTAL_DEGREE"; deg(nonlinear);
"@@TERM_COUNT"; size(nonlinear);
poly q=nonlinear;
while(q!=0) {
  "@@TERM"; leadcoef(q); leadexp(q);
  q=q-lead(q);
}
'''


def _integer_after(marker: str, output: str) -> int:
    match = re.search(rf"^{re.escape(marker)}\n(-?\d+)$", output, re.MULTILINE)
    if match is None:
        raise AssertionError(f"missing Singular marker {marker}")
    return int(match.group(1))


def normalize_coefficients(
    rational_coefficients: dict[tuple[int, int], Fraction],
) -> dict[tuple[int, int], int]:
    """Clear denominators, divide content and choose a positive leading sign."""

    denominator = reduce(lcm, (value.denominator for value in rational_coefficients.values()), 1)
    integers = {
        monomial: value.numerator * (denominator // value.denominator)
        for monomial, value in rational_coefficients.items()
    }
    content = reduce(gcd, (abs(value) for value in integers.values() if value), 0)
    integers = {monomial: value // content for monomial, value in integers.items()}
    leading_monomial = max(integers, key=lambda item: (sum(item), item[0], item[1]))
    if integers[leading_monomial] < 0:
        integers = {monomial: -value for monomial, value in integers.items()}
    return integers


def coefficient_digest(coefficients: dict[tuple[int, int], int]) -> str:
    rows = [
        [i, j, str(coefficients[i, j])]
        for i, j in sorted(coefficients, reverse=True)
    ]
    return sha256_text(json.dumps(rows, separators=(",", ":")))


def extract_nonlinear_factor(*, timeout: int = 300) -> tuple[dict[tuple[int, int], int], dict[str, Any]]:
    """Recompute and extract the primitive degree-32 factor with Singular."""

    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required to extract the nonlinear factor")
    completed = subprocess.run(
        [executable, "-q"],
        input=_singular_program(),
        text=True,
        capture_output=True,
        timeout=timeout,
        check=True,
    )
    if completed.stderr.strip():
        raise RuntimeError(f"Singular wrote to stderr: {completed.stderr.strip()}")
    output = completed.stdout
    if "@@ERROR_DIVISION" in output or "   ?" in output:
        raise AssertionError(f"Singular nonlinear-factor extraction failed:\n{output}")
    term_pattern = re.compile(
        r"^@@TERM\n([^\n]+)\n(\d+),(\d+),(\d+)$", re.MULTILINE
    )
    rational: dict[tuple[int, int], Fraction] = {}
    for coefficient, c_degree, k_degree, t_degree in term_pattern.findall(output):
        if int(t_degree) != 0:
            raise AssertionError("the extracted discriminant factor depends on T")
        monomial = (int(c_degree), int(k_degree))
        if monomial in rational:
            raise AssertionError(f"duplicate monomial {monomial}")
        rational[monomial] = Fraction(coefficient.strip())
    total_degree = _integer_after("@@TOTAL_DEGREE", output)
    term_count = _integer_after("@@TERM_COUNT", output)
    if total_degree != TOTAL_DEGREE or term_count != EXPECTED_TERM_COUNT:
        raise AssertionError("the pinned nonlinear factor shape changed")
    if len(rational) != term_count:
        raise AssertionError(f"parsed {len(rational)} terms, expected {term_count}")
    coefficients = normalize_coefficients(rational)
    if max(i for i, _ in coefficients) != TOTAL_DEGREE:
        raise AssertionError("the c-degree changed")
    if max(j for _, j in coefficients) != TOTAL_DEGREE:
        raise AssertionError("the k-degree changed")
    primitive_digest = coefficient_digest(coefficients)
    if primitive_digest != EXPECTED_PRIMITIVE_COEFFICIENT_SHA256:
        raise AssertionError("the primitive nonlinear factor changed")
    metadata = {
        "total_degree": total_degree,
        "term_count": term_count,
        "degree_in_c": max(i for i, _ in coefficients),
        "degree_in_k": max(j for _, j in coefficients),
        "primitive_coefficient_sha256": primitive_digest,
        "upstream_singular_factor_text_sha256": UPSTREAM_FACTOR_SHA256,
    }
    return coefficients, metadata


def evaluate_homogeneous(
    coefficients: dict[tuple[int, int], int], C: int, K: int, D: int
) -> int:
    """Evaluate ``D^32 f(C/D,K/D)`` exactly, including ``D=0``."""

    powers_c = [1]
    powers_k = [1]
    powers_d = [1]
    for _ in range(TOTAL_DEGREE):
        powers_c.append(powers_c[-1] * C)
        powers_k.append(powers_k[-1] * K)
        powers_d.append(powers_d[-1] * D)
    return sum(
        coefficient * powers_c[i] * powers_k[j] * powers_d[TOTAL_DEGREE - i - j]
        for (i, j), coefficient in coefficients.items()
    )


def affine_mod_value(
    coefficients: dict[tuple[int, int], int], c_value: int, k_value: int, prime: int
) -> int:
    """Evaluate the affine factor modulo ``prime`` with nested Horner loops."""

    rows: list[list[int]] = [[] for _ in range(TOTAL_DEGREE + 1)]
    for i in range(TOTAL_DEGREE + 1):
        rows[i] = [coefficients.get((i, j), 0) % prime for j in range(TOTAL_DEGREE - i + 1)]
    answer = 0
    for i in range(TOTAL_DEGREE, -1, -1):
        row_value = 0
        for coefficient in reversed(rows[i]):
            row_value = (row_value * k_value + coefficient) % prime
        answer = (answer * c_value + row_value) % prime
    return answer


def modular_root_table(
    coefficients: dict[tuple[int, int], int], prime: int
) -> tuple[tuple[tuple[int, ...], ...], dict[str, Any]]:
    """Precompute allowed projective K residues for every pair (C,D)."""

    affine_roots = tuple(
        tuple(
            k_value
            for k_value in range(prime)
            if affine_mod_value(coefficients, c_value, k_value, prime) == 0
        )
        for c_value in range(prime)
    )
    top = {
        (i, j): value
        for (i, j), value in coefficients.items()
        if i + j == TOTAL_DEGREE
    }
    infinity_roots = tuple(
        r
        for r in range(prime)
        if sum(value * pow(r, j, prime) for (i, j), value in top.items()) % prime == 0
    )
    k_leading = coefficients.get((0, TOTAL_DEGREE), 0) % prime

    table: list[list[tuple[int, ...]]] = [
        [tuple() for _ in range(prime)] for _ in range(prime)
    ]
    for C in range(prime):
        for D in range(prime):
            if D:
                inverse_d = pow(D, -1, prime)
                c_affine = C * inverse_d % prime
                table[C][D] = tuple(D * root % prime for root in affine_roots[c_affine])
            elif C:
                table[C][D] = tuple(C * root % prime for root in infinity_roots)
            else:
                # At (C,D)=(0,0), the homogeneous value is a_k32*K^32.
                # A prime dividing a_k32 is still a sound (if weaker) filter.
                table[C][D] = (0,) if k_leading else tuple(range(prime))
    metadata = {
        "prime": prime,
        "affine_Fp_point_count": sum(len(row) for row in affine_roots),
        "affine_vertical_fibers_without_points": sum(not row for row in affine_roots),
        "maximum_roots_in_one_affine_vertical_fiber": max(map(len, affine_roots)),
        "projective_points_at_infinity": len(infinity_roots),
        "k_power_32_coefficient_nonzero": bool(k_leading),
    }
    return tuple(tuple(row) for row in table), metadata


def _centered_residue(residue: int, modulus: int, height: int) -> int | None:
    if residue <= height:
        return residue
    if modulus - residue <= height:
        return residue - modulus
    return None


def search_affine_box(
    coefficients: dict[tuple[int, int], int], *, height: int, primes: tuple[int, ...]
) -> tuple[list[tuple[int, int, int]], dict[str, Any]]:
    """Search every normalized affine projective triple in the declared box."""

    modulus = prod(primes)
    if modulus <= 2 * height:
        raise ValueError("the CRT modulus must exceed twice the height bound")
    tables_and_metadata = [modular_root_table(coefficients, prime) for prime in primes]
    tables = [item[0] for item in tables_and_metadata]
    root_metadata = [item[1] for item in tables_and_metadata]
    crt_basis = [
        (modulus // prime) * pow(modulus // prime, -1, prime) % modulus
        for prime in primes
    ]
    crt_contribution = [
        tuple(root * basis % modulus for root in range(prime))
        for prime, basis in zip(primes, crt_basis, strict=True)
    ]

    first_empty_rejections = [0 for _ in primes]
    combinations_tested = 0
    centered_candidates = 0
    primitive_candidates = 0
    exact_evaluations = 0
    integer_exact_evaluations = 0
    hits: list[tuple[int, int, int]] = []
    seen: set[tuple[int, int, int]] = set()

    for D in range(1, height + 1):
        residues_d = [D % prime for prime in primes]
        for C in range(-height, height + 1):
            root_lists: list[tuple[int, ...]] = []
            rejected = False
            for index, (prime, table) in enumerate(zip(primes, tables, strict=True)):
                roots = table[C % prime][residues_d[index]]
                if not roots:
                    first_empty_rejections[index] += 1
                    rejected = True
                    break
                root_lists.append(roots)
            if rejected:
                continue
            for roots in itertools.product(*root_lists):
                combinations_tested += 1
                residue = sum(
                    crt_contribution[index][root]
                    for index, root in enumerate(roots)
                ) % modulus
                K = _centered_residue(residue, modulus, height)
                if K is None:
                    continue
                centered_candidates += 1
                triple = (C, K, D)
                if triple in seen:
                    raise AssertionError("duplicate CRT candidate")
                seen.add(triple)
                if gcd(gcd(abs(C), abs(K)), D) != 1:
                    continue
                primitive_candidates += 1
                exact_evaluations += 1
                if D == 1:
                    integer_exact_evaluations += 1
                if evaluate_homogeneous(coefficients, C, K, D) == 0:
                    hits.append(triple)

    total_cd_pairs = (2 * height + 1) * height
    metadata = {
        "height": height,
        "definition": (
            "D>0, gcd(C,K,D)=1, max(|C|,|K|,D)<=H; "
            "c=C/D and k=K/D"
        ),
        "sieve_primes": list(primes),
        "crt_modulus": modulus,
        "crt_modulus_exceeds_twice_height": True,
        "modular_root_tables": root_metadata,
        "C_D_pairs_scanned": total_cd_pairs,
        "first_empty_root_list_rejections": first_empty_rejections,
        "C_D_pairs_reaching_crt": total_cd_pairs - sum(first_empty_rejections),
        "crt_root_combinations_tested": combinations_tested,
        "centered_K_candidates": centered_candidates,
        "primitive_exact_candidates": primitive_candidates,
        "exact_homogeneous_evaluations": exact_evaluations,
        "integer_parameter_exact_evaluations": integer_exact_evaluations,
        "exact_hits": len(hits),
        "integer_parameter_hits": sum(D == 1 for _, _, D in hits),
    }
    return hits, metadata


def primitive_projective(c_value: Fraction, k_value: Fraction) -> tuple[int, int, int]:
    denominator = lcm(c_value.denominator, k_value.denominator)
    C = c_value.numerator * (denominator // c_value.denominator)
    K = k_value.numerator * (denominator // k_value.denominator)
    common = gcd(gcd(abs(C), abs(K)), denominator)
    return C // common, K // common, denominator // common


def rational_roots(polynomial: sp.Poly) -> list[tuple[Fraction, int]]:
    """Return every rational root and its exact multiplicity."""

    answer: list[tuple[Fraction, int]] = []
    for factor, exponent in sp.factor_list(polynomial)[1]:
        if factor.degree() == 1:
            a, b = factor.all_coeffs()
            answer.append((Fraction(-int(b.p), int(b.q)) / Fraction(int(a.p), int(a.q)), exponent))
    return sorted(answer)


def polynomial_from_coefficients(
    coefficients: dict[tuple[int, int], int], c_expr: sp.Expr, k_expr: sp.Expr, variable: sp.Symbol
) -> sp.Poly:
    expression = sum(
        value * c_expr**i * k_expr**j
        for (i, j), value in coefficients.items()
    )
    numerator = sp.together(expression).as_numer_denom()[0]
    polynomial = sp.Poly(numerator, variable, domain=sp.QQ).primitive()[1]
    return -polynomial if polynomial.LC() < 0 else polynomial


def known_line_intersections(
    coefficients: dict[tuple[int, int], int]
) -> tuple[list[tuple[int, int, int]], list[dict[str, Any]]]:
    """Factor the restriction to every known line and retain rational points."""

    z = sp.symbols("z")
    slope = sp.Rational(AFFINE_CANCELLATION_SLOPE.numerator, AFFINE_CANCELLATION_SLOPE.denominator)
    line_specs = (
        ("k-c=0", z, z, (1, 1, 0)),
        ("k+c=0", z, -z, (1, -1, 0)),
        ("5899690*c+732683*k=0", z, slope * z, (732683, -5899690, 0)),
        ("28917*c+10=0", -sp.Rational(10, 28917), z, (0, 1, 0)),
        ("39508*c+39=0", -sp.Rational(39, 39508), z, (0, 1, 0)),
    )
    points: set[tuple[int, int, int]] = set()
    records: list[dict[str, Any]] = []
    for label, c_expr, k_expr, infinity_point in line_specs:
        restriction = polynomial_from_coefficients(coefficients, c_expr, k_expr, z)
        roots = rational_roots(restriction)
        affine_points = []
        for root, multiplicity in roots:
            c_value = Fraction(c_expr) if c_expr != z else root
            if k_expr == z:
                k_value = root
            elif k_expr == -z:
                k_value = -root
            else:
                k_value = Fraction(slope) * root
            triple = primitive_projective(c_value, k_value)
            points.add(triple)
            affine_points.append(
                {
                    "c": rational_text(c_value),
                    "k": rational_text(k_value),
                    "projective": list(triple),
                    "intersection_multiplicity_on_line": multiplicity,
                }
            )
        C, K, D = infinity_point
        at_infinity = evaluate_homogeneous(coefficients, C, K, D) == 0
        if at_infinity:
            common = gcd(abs(C), abs(K))
            normalized = (C // common, K // common, 0)
            if normalized[0] < 0 or (normalized[0] == 0 and normalized[1] < 0):
                normalized = tuple(-value for value in normalized)
            points.add(normalized)
        records.append(
            {
                "line": label,
                "restriction_degree": restriction.degree(),
                "factor_signature_over_QQ": [list(item) for item in factor_signature(restriction)],
                "rational_affine_intersection_count": len(affine_points),
                "rational_affine_intersections": affine_points,
                "rational_intersection_at_infinity": at_infinity,
            }
        )
    return sorted(points), records


def _squareclass_record(polynomial: sp.Poly) -> dict[str, Any]:
    factors = sp.factor_list(polynomial)
    kernel_degree = squareclass_kernel_degree(polynomial)
    return {
        "T_degree": polynomial.degree(),
        "factor_signature_over_QQ": [list(item) for item in factor_signature(polynomial)],
        "squareclass_kernel_degree": kernel_degree,
        "squareclass_genus": hyperelliptic_genus(kernel_degree),
        "genus_at_most_one": kernel_degree <= 4,
        "square_kernel": kernel_degree == 0,
        "factorization_unit": rational_text(Fraction(factors[0])),
    }


def line_memberships(C: int, K: int, D: int) -> list[str]:
    tests = (
        ("k-c=0", K - C),
        ("k+c=0", K + C),
        ("5899690*c+732683*k=0", 5899690 * C + 732683 * K),
        ("28917*c+10=0", 28917 * C + 10 * D),
        ("39508*c+39=0", 39508 * C + 39 * D),
    )
    return [label for label, value in tests if value == 0]


def classify_point(
    projective: tuple[int, int, int], pencil: sp.Poly, symbols: dict[str, sp.Expr]
) -> dict[str, Any]:
    """Specialize the cleared branch polynomial and its squareclass kernel."""

    C, K, D = projective
    T, c, k = symbols["T"], symbols["c"], symbols["k"]
    if D:
        c_value = sp.Rational(C, D)
        k_value = sp.Rational(K, D)
        specialized = sp.Poly(
            sp.expand(pencil.as_expr().subs({c: c_value, k: k_value})), T, domain=sp.QQ
        ).primitive()[1]
        parameter_text = {"c": rational_text(Fraction(C, D)), "k": rational_text(Fraction(K, D))}
    else:
        expression = 0
        for monomial, coefficient in pencil.terms():
            t_degree, c_degree, k_degree = monomial
            parameter_degree = c_degree + k_degree
            if parameter_degree == PARAMETER_DEGREE:
                expression += coefficient * T**t_degree * C**c_degree * K**k_degree
        specialized = sp.Poly(sp.expand(expression), T, domain=sp.QQ).primitive()[1]
        parameter_text = {"c": "projective infinity", "k": "projective infinity"}
    if specialized.LC() < 0:
        specialized = -specialized
    memberships = line_memberships(C, K, D)
    invalid = [
        label for label in memberships if label in {"28917*c+10=0", "39508*c+39=0"}
    ]
    record = {
        "projective": [C, K, D],
        "projective_height": max(abs(C), abs(K), D),
        "parameters": parameter_text,
        "known_linear_component_memberships": memberships,
        "endpoint_valid": not invalid,
        "invalid_endpoint_lines": invalid,
        **_squareclass_record(specialized),
    }
    record["new_section"] = record["square_kernel"] and record["endpoint_valid"]
    record["new_low_genus_base_change"] = record["genus_at_most_one"] and record["endpoint_valid"]
    return record


def projective_boundary(
    coefficients: dict[tuple[int, int], int]
) -> tuple[list[tuple[int, int, int]], dict[str, Any]]:
    """Globally factor the top binary form and list all rational directions."""

    z = sp.symbols("z")
    top = sp.Poly(
        sum(
            value * z**j
            for (i, j), value in coefficients.items()
            if i + j == TOTAL_DEGREE
        ),
        z,
        domain=sp.QQ,
    ).primitive()[1]
    roots = rational_roots(top)
    points = []
    root_records = []
    for root, multiplicity in roots:
        C, K, _ = primitive_projective(Fraction(1), root)
        point = (C, K, 0)
        points.append(point)
        root_records.append(
            {
                "K_over_C": rational_text(root),
                "projective": list(point),
                "multiplicity": multiplicity,
            }
        )
    c_zero_point = (0, 1, 0)
    c_zero_is_point = evaluate_homogeneous(coefficients, *c_zero_point) == 0
    if c_zero_is_point:
        points.append(c_zero_point)
    return sorted(set(points)), {
        "top_binary_form_degree": top.degree(),
        "factor_signature_over_QQ_after_C_equals_1": [list(item) for item in factor_signature(top)],
        "rational_K_over_C_roots": root_records,
        "C_equals_0_point_present": c_zero_is_point,
        "rational_projective_point_count": len(set(points)),
    }


def run(root: Path, *, height: int = DEFAULT_HEIGHT, timeout: int = 300) -> dict[str, Any]:
    if height <= 0:
        raise ValueError("height must be positive")
    coefficients, factor_metadata = extract_nonlinear_factor(timeout=timeout)
    affine_hits, search_metadata = search_affine_box(
        coefficients, height=height, primes=SIEVE_PRIMES
    )
    infinity_hits, infinity_metadata = projective_boundary(coefficients)
    intersection_points, intersection_records = known_line_intersections(coefficients)

    pencil, symbols = pencil_polynomial(Fraction(256714, 39), Fraction(-8545))
    all_points = sorted(set(affine_hits) | set(infinity_hits) | set(intersection_points))
    classifications = [classify_point(point, pencil, symbols) for point in all_points]
    affine_classifications = [record for record in classifications if record["projective"][2] > 0]
    bounded_classifications = [
        record for record in classifications if record["projective_height"] <= height
    ]
    low_genus = [record for record in classifications if record["new_low_genus_base_change"]]
    sections = [record for record in classifications if record["new_section"]]

    artifact: dict[str, Any] = {
        "schema_version": "elliptic-curves.fermigier-bidegree21-nonlinear-points.v1",
        "status": "complete exact bounded rational-point sieve",
        "claim_level": "bounded exact computation; not a global rational-point theorem",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "family_id": "fermigier-mestre-v1",
        "representative_pair": {
            "left": "P13 at literal_shift_T=39508/39",
            "right": "R20E1 at literal_shift_T=28917/10",
        },
        "nonlinear_factor": factor_metadata,
        "search_region": search_metadata,
        "projective_boundary": infinity_metadata,
        "known_linear_component_intersections": {
            "method": "exact QQ factorization of all five line restrictions",
            "lines": intersection_records,
            "distinct_rational_intersection_points": len(intersection_points),
        },
        "classified_points": classifications,
        "scope": {
            "affine_box_complete": True,
            "projective_boundary_globally_classified": True,
            "known_line_rational_intersections_globally_classified": True,
            "bounded_projective_points_classified": len(bounded_classifications),
            "all_rational_points_on_degree32_component_classified": False,
            "other_exceptional_pairs_classified": False,
            "not_claimed": (
                "No claim is made outside projective height H, except for the exact "
                "boundary and known-line intersection factorizations; no claim is made "
                "for the other 79 exceptional pairs."
            ),
        },
        "outcome": {
            "affine_degree32_hits_in_box": len(affine_hits),
            "rational_points_at_projective_infinity": len(infinity_hits),
            "distinct_rational_known_line_intersections": len(intersection_points),
            "valid_genus_at_most_one_points_found": len(low_genus),
            "new_sections": len(sections),
            "new_specializations": 0,
            "target_met": bool(sections),
        },
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
            "singular": subprocess.run(
                [shutil.which("Singular") or "Singular", "--version"],
                input="",
                text=True,
                capture_output=True,
                check=True,
            ).stdout.splitlines()[0],
        },
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves:elliptic-curves/cas .venv/bin/python "
            "elliptic-curves/cas/search_fermigier_bidegree21_nonlinear_points.py "
            f"--height {height}"
        ),
    }
    artifact["result_sha256"] = result_digest(artifact)
    return artifact


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=DEFAULT_HEIGHT)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--output", type=Path, default=root / OUTPUT_RELATIVE)
    args = parser.parse_args()
    if args.height != DEFAULT_HEIGHT and args.output == root / OUTPUT_RELATIVE:
        raise SystemExit("use --output when changing the pinned height")
    artifact = run(root, height=args.height, timeout=args.timeout)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps(artifact["outcome"], sort_keys=True))
    print(f"result_sha256={artifact['result_sha256']}")


if __name__ == "__main__":
    main()
