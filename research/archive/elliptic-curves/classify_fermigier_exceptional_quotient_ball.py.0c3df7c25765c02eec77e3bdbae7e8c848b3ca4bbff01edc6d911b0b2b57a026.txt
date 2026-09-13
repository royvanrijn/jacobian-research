#!/usr/bin/env python3
"""Classify the signed weight-two Fermigier exceptional quotient balls.

This is a finite structural calculation, not a specialization score sweep.
At each of two certified fibers it starts with the exceptional quotient basis,
enumerates every nonzero vector in ``{-1,0,1}^n`` of support at most two,
and maps the resulting group point back to the pointed Fermigier quartic.
Every inverse is checked in two independent coordinates:

* the elementary pointed-quartic/Weierstrass birational group law; and
* the canonical Fermigier model used by the finite-reduction certificates.

The mod-5 certificate matrix also verifies separately that every enumerated
combination remains outside the generic rank-12 span.  Finally, all 200 by
128 cross-anchor affine interpolants are classified exactly.  An irreducible
reduction modulo one recorded prime is used as a Gauss-lemma certificate of
irreducibility over ``QQ``; the exact ``QQ`` factorization is retained for any
polynomial for which no such witness is found.

The family polynomial uses the legacy literal shift ``T=2u``.  Both that
alias and the canonical adapter parameter ``u`` are recorded at each anchor.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import gcd, lcm
import multiprocessing as mp
from pathlib import Path
import platform
import sys
from typing import Any, Iterable, Sequence

import sympy as sp

from ecsearch.fermigier import quartic_point_to_canonical_point
from ecsearch.fermigier_near_miss import canonical_ratpoints_output
from ecsearch.fermigier_rank import (
    parse_ratpoints_output,
    section_and_point_cloud_differences,
    specialize_fermigier_rank_sections,
)
from ecsearch.rank_certification import (
    IndependenceCertificate,
    add_rational_points,
    matrix_rank_mod_prime,
    subtract_rational_points,
    verify_independence_certificate,
)
from fermigier_mestre import FermigierMestreFamily
from search_fermigier_rank22_accidental_slices import (
    FERMIGIER_BIVARIATE_COEFFICIENTS,
    poly_add,
    poly_evaluate,
    poly_multiply,
    published_accidental_points,
    select_reconstruction_convention,
)
from search_nagao_section7_auxiliary_jacobians import (
    translate_polynomial,
    weierstrass_add,
    weierstrass_multiply,
)
from triage_nagao_rank13_finalists import point_on_short_curve


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

FAMILY_ID = "fermigier-mestre-v1"
E22_U = Q(19754, 39)
E22_T = 2 * E22_U
RANK20_U = Q(28917, 20)
RANK20_T = 2 * RANK20_U
RANK20_EXCEPTIONAL_X = (
    Q(-8545),
    Q(23004, 5),
    Q(-8817, 10),
    Q(8183, 10),
    Q(76563, 10),
    Q(-69561, 20),
    Q(-431673, 70),
    Q(-408943, 110),
)
RANK20_SELECTED_INDICES = tuple(range(12)) + (12, 14, 28, 30, 40, 44, 57, 59)
RANK20_ARTIFACT_RELATIVE = Path(
    "artifacts/generated-results/elliptic-curves/fermigier_rank20_near_miss_v1.json"
)
RANK20_ARTIFACT_SHA256 = (
    "8416e835887236e9e4eafcb01384a710ce4f1be0628701a97f4a7d7a07fe63b1"
)
RANK22_ARTIFACT_RELATIVE = Path(
    "artifacts/generated-results/elliptic-curves/elliptic_fermigier_rank22_accidental_slices.json"
)
RANK22_ARTIFACT_SHA256 = (
    "3794f23d37685edcc1ad5c8279d48fcd247b1d31b7c88aa4939ca5efd40f79a7"
)
MODULAR_WITNESS_PRIMES = (
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29,
    31,
    37,
    41,
    43,
    47,
    53,
    59,
    61,
    67,
    71,
    73,
    79,
    83,
    89,
    97,
    101,
    103,
    107,
    109,
    113,
    127,
    131,
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_lines(lines: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for line in lines:
        digest.update((line + "\n").encode())
    return digest.hexdigest()


def rational_text(value: Fraction | sp.Rational) -> str:
    if isinstance(value, sp.Rational):
        value = Q(int(value.p), int(value.q))
    else:
        value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def point_record(point: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {"x": rational_text(point[0]), "y": rational_text(point[1])}


def projective_height(value: Fraction) -> int:
    value = Q(value)
    return max(abs(value.numerator), value.denominator)


@dataclass(frozen=True)
class PointedQuartic:
    """Elementary birational group model for one pointed quartic fiber."""

    literal_shift: Fraction
    base_point: tuple[Fraction, Fraction]
    quartic_coefficients: tuple[Fraction, ...]
    shifted_coefficients: tuple[Fraction, ...]
    weierstrass_coefficients: tuple[Fraction, ...]

    @classmethod
    def construct(cls, literal_shift: Fraction) -> "PointedQuartic":
        literal_shift = Q(literal_shift)
        base = tuple(
            Q(value)
            for value in FermigierMestreFamily.known_quartic_points(literal_shift)[0]
        )
        # Ascending x-power order is used by translate_polynomial.
        coefficients = tuple(
            reversed(FermigierMestreFamily.quartic_coefficients(literal_shift))
        )
        shifted = translate_polynomial(coefficients, base[0])
        if len(shifted) != 5 or shifted[0] != base[1] ** 2:
            raise AssertionError("the selected pointed-quartic origin is invalid")
        _, d_value, c_value, b_value, a_value = shifted
        q_value = base[1]
        weierstrass = (
            d_value / q_value,
            c_value - d_value**2 / (4 * q_value**2),
            2 * q_value * b_value,
            -4 * q_value**2 * a_value,
            a_value * (d_value**2 - 4 * q_value**2 * c_value),
        )
        return cls(literal_shift, base, coefficients, shifted, weierstrass)

    def quartic_value(self, x_value: Fraction) -> Fraction:
        return poly_evaluate(self.quartic_coefficients, Q(x_value))

    def forward(self, point: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
        x_value, ordinate = (Q(value) for value in point)
        if ordinate**2 != self.quartic_value(x_value):
            raise ValueError("the point missed the pointed quartic")
        u_value = x_value - self.base_point[0]
        if u_value == 0:
            raise ValueError("the pointed origin maps to infinity")
        _, d_value, c_value, _, _ = self.shifted_coefficients
        q_value = self.base_point[1]
        image_x = (
            2 * q_value * (ordinate + q_value) + d_value * u_value
        ) / u_value**2
        image_y = (
            4 * q_value**2 * (ordinate + q_value)
            + 2 * q_value * (d_value * u_value + c_value * u_value**2)
            - d_value**2 * u_value**2 / (2 * q_value)
        ) / u_value**3
        image = image_x, image_y
        if not point_on_short_curve(self.weierstrass_coefficients, image):
            raise AssertionError("the pointed-quartic forward map failed")
        return image

    def inverse(
        self, point: tuple[Fraction, Fraction] | None
    ) -> tuple[Fraction, Fraction] | None:
        if point is None:
            return self.base_point
        if not point_on_short_curve(self.weierstrass_coefficients, point):
            raise ValueError("the point missed the auxiliary Weierstrass curve")
        image_x, image_y = point
        if image_y == 0:
            return None
        _, d_value, c_value, _, _ = self.shifted_coefficients
        q_value = self.base_point[1]
        u_value = (
            4 * q_value**2 * (image_x + c_value) - d_value**2
        ) / (2 * q_value * image_y)
        if u_value == 0:
            return None
        ordinate = (
            (image_x * u_value**2 - d_value * u_value) / (2 * q_value)
            - q_value
        )
        answer = self.base_point[0] + u_value, ordinate
        if ordinate**2 != self.quartic_value(answer[0]):
            raise AssertionError("the pointed-quartic inverse map failed")
        return answer


def signed_weight_two_vectors(dimension: int) -> tuple[tuple[int, ...], ...]:
    """Return every signed support-one or support-two vector, without +/- quotient."""

    vectors: list[tuple[int, ...]] = []
    for first in range(dimension):
        for first_sign in (-1, 1):
            vector = [0] * dimension
            vector[first] = first_sign
            vectors.append(tuple(vector))
    for first in range(dimension):
        for second in range(first + 1, dimension):
            for first_sign, second_sign in product((-1, 1), repeat=2):
                vector = [0] * dimension
                vector[first] = first_sign
                vector[second] = second_sign
                vectors.append(tuple(vector))
    expected = 2 * dimension + 4 * dimension * (dimension - 1) // 2
    if len(vectors) != expected or len(set(vectors)) != expected:
        raise AssertionError("the signed weight-two vector ball changed")
    return tuple(vectors)


def vector_id(labels: Sequence[str], vector: Sequence[int]) -> str:
    return "_".join(
        f"{'p' if coefficient > 0 else 'm'}{label}"
        for label, coefficient in zip(labels, vector, strict=True)
        if coefficient
    )


def exact_combination(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    vector: Sequence[int],
) -> tuple[Fraction, Fraction] | None:
    answer = None
    for point, scalar in zip(points, vector, strict=True):
        if scalar:
            answer = weierstrass_add(
                coefficients,
                answer,
                weierstrass_multiply(coefficients, point, scalar),
            )
    return answer


def canonical_combination(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    vector: Sequence[int],
) -> tuple[Fraction, Fraction] | None:
    answer = None
    for point, scalar in zip(points, vector, strict=True):
        if scalar:
            addend = point
            if scalar == -1:
                # Multiplication by -1 through the same exact generalized law.
                addend = weierstrass_multiply(coefficients, point, -1)
                if addend is None:
                    raise AssertionError("an exceptional basis point became 2-torsion")
            answer = add_rational_points(coefficients, answer, addend)
    return answer


def certificate_separates_generic_coset(
    certificate: IndependenceCertificate,
    exceptional_dimension: int,
    vector: Sequence[int],
) -> bool:
    ell = certificate.relation_prime
    if len(vector) != exceptional_dimension:
        raise ValueError("the exceptional vector has the wrong dimension")
    matrix = []
    for row in certificate.rows:
        if len(row.logs) != 12 + exceptional_dimension:
            raise AssertionError("the certificate point ordering changed")
        exceptional_log = sum(
            scalar * row.logs[12 + index]
            for index, scalar in enumerate(vector)
        ) % ell
        matrix.append([*(value % ell for value in row.logs[:12]), exceptional_log])
    return matrix_rank_mod_prime(matrix, ell) == 13


def _oriented_basis_point(
    specialization: Any,
    candidates: Sequence[tuple[Fraction, Fraction]],
    target: tuple[Fraction, Fraction],
    x_value: Fraction,
) -> tuple[Fraction, Fraction]:
    answers = []
    for point in candidates:
        if point[0] != x_value:
            continue
        canonical = quartic_point_to_canonical_point(
            specialization.quartic_model, point
        )
        difference = subtract_rational_points(
            specialization.canonical_model,
            canonical,
            specialization.canonical_points[0],
        )
        if difference == target:
            answers.append(point)
    if len(answers) != 1:
        raise AssertionError(
            f"expected one oriented quartic preimage for {x_value}, found {len(answers)}"
        )
    return answers[0]


def e22_basis() -> tuple[
    Any,
    tuple[str, ...],
    tuple[tuple[Fraction, Fraction], ...],
    tuple[tuple[Fraction, Fraction], ...],
    IndependenceCertificate,
]:
    specialization = specialize_fermigier_rank_sections(E22_U)
    _, _, reconstruction, _ = select_reconstruction_convention()
    accidentals = published_accidental_points(reconstruction)
    if [label for label, _ in accidentals] != [
        "P6",
        *[f"P{index}" for index in range(13, 23)],
    ]:
        raise AssertionError("the E22 accidental ordering changed")
    all_points = parse_ratpoints_output(
        specialization.quartic_model,
        canonical_ratpoints_output(tuple(point[0] for _, point in accidentals)),
    )
    positive_points = all_points[::2]
    cloud = section_and_point_cloud_differences(specialization, positive_points)
    selected_indices = tuple(range(12)) + tuple(range(13, 23))
    labels = tuple(label for label, _ in accidentals[1:])
    quartic_basis = tuple(positive_points[index] for index in range(1, 11))
    canonical_basis = tuple(cloud[index] for index in selected_indices[12:])
    for point, target in zip(quartic_basis, canonical_basis, strict=True):
        oriented = _oriented_basis_point(
            specialization, all_points, target, point[0]
        )
        if oriented != point:
            raise AssertionError("the E22 positive-sheet orientation changed")

    # The pinned reconstruction artifact carries the independent mod-5 matrix
    # in its exact rank certificate.  Recompute the same matrix through the
    # existing selector to avoid trusting only serialized point labels.
    from ecsearch.rank_certification import select_independent_subset

    selected, certificate = select_independent_subset(
        specialization.canonical_model,
        cloud,
        relation_prime=5,
        maximum_reduction_prime=2_000,
    )
    if selected != selected_indices:
        raise AssertionError("the E22 independent subset changed")
    selected_points = tuple(cloud[index] for index in selected)
    verify_independence_certificate(
        specialization.canonical_model, selected_points, certificate
    )
    return specialization, labels, quartic_basis, canonical_basis, certificate


def rank20_basis(root: Path) -> tuple[
    Any,
    tuple[str, ...],
    tuple[tuple[Fraction, Fraction], ...],
    tuple[tuple[Fraction, Fraction], ...],
    IndependenceCertificate,
]:
    source = root / RANK20_ARTIFACT_RELATIVE
    if sha256_file(source) != RANK20_ARTIFACT_SHA256:
        raise AssertionError("the pinned rank-20 source artifact changed")
    artifact = json.loads(source.read_text())
    specialization = specialize_fermigier_rank_sections(RANK20_U)
    abscissas = tuple(Q(value) for value in artifact["bounded_search"]["abscissas"])
    searched = parse_ratpoints_output(
        specialization.quartic_model, canonical_ratpoints_output(abscissas)
    )
    cloud = section_and_point_cloud_differences(specialization, searched)
    selected_indices = tuple(artifact["point_cloud"]["selected_indices"])
    if selected_indices != RANK20_SELECTED_INDICES:
        raise AssertionError("the rank-20 independent subset changed")
    selected_points = tuple(cloud[index] for index in selected_indices)
    certificate = IndependenceCertificate.from_json_object(
        artifact["point_cloud"]["certificate"]
    )
    verify_independence_certificate(
        specialization.canonical_model, selected_points, certificate
    )
    labels = tuple(f"R20E{index}" for index in range(1, 9))
    canonical_basis = tuple(cloud[index] for index in selected_indices[12:])
    quartic_basis = tuple(
        _oriented_basis_point(specialization, searched, target, x_value)
        for target, x_value in zip(
            canonical_basis, RANK20_EXCEPTIONAL_X, strict=True
        )
    )
    return specialization, labels, quartic_basis, canonical_basis, certificate


def direction_ball(
    *,
    anchor_name: str,
    canonical_u: Fraction,
    literal_t: Fraction,
    specialization: Any,
    labels: Sequence[str],
    quartic_basis: Sequence[tuple[Fraction, Fraction]],
    canonical_basis: Sequence[tuple[Fraction, Fraction]],
    certificate: IndependenceCertificate,
) -> tuple[dict[str, Any], tuple[dict[str, Any], ...]]:
    auxiliary = PointedQuartic.construct(literal_t)
    # ecsearch retains the unnormalized remainder ordinate, whereas the CAS
    # family divides the quartic by (50616*T)^2.  The abscissa is unchanged.
    ordinate_scale = 50616 * literal_t
    normalized_basis = tuple(
        (point[0], point[1] / ordinate_scale) for point in quartic_basis
    )
    auxiliary_basis = tuple(auxiliary.forward(point) for point in normalized_basis)
    for source, image in zip(normalized_basis, auxiliary_basis, strict=True):
        if auxiliary.inverse(image) != source:
            raise AssertionError("a basis point failed its pointed-quartic round trip")

    records = []
    seen_canonical: set[tuple[Fraction, Fraction]] = set()
    seen_quartic: set[tuple[Fraction, Fraction]] = set()
    for vector in signed_weight_two_vectors(len(labels)):
        auxiliary_sum = exact_combination(
            auxiliary.weierstrass_coefficients, auxiliary_basis, vector
        )
        quartic_point = auxiliary.inverse(auxiliary_sum)
        if quartic_point is None:
            raise AssertionError("a certified exceptional combination had no affine inverse")
        if quartic_point in seen_quartic:
            raise AssertionError("two signed vectors produced the same oriented quartic point")
        seen_quartic.add(quartic_point)
        canonical_sum = canonical_combination(
            specialization.canonical_model, canonical_basis, vector
        )
        if canonical_sum is None:
            raise AssertionError("a certified exceptional combination became zero")
        if canonical_sum in seen_canonical:
            raise AssertionError("two signed vectors produced the same canonical point")
        seen_canonical.add(canonical_sum)
        canonical_image = quartic_point_to_canonical_point(
            specialization.quartic_model,
            (quartic_point[0], quartic_point[1] * ordinate_scale),
        )
        canonical_difference = subtract_rational_points(
            specialization.canonical_model,
            canonical_image,
            specialization.canonical_points[0],
        )
        if canonical_difference != canonical_sum:
            raise AssertionError("the auxiliary and canonical group combinations disagree")
        if not certificate_separates_generic_coset(
            certificate, len(labels), vector
        ):
            raise AssertionError("an exceptional combination entered the generic mod-5 span")
        records.append(
            {
                "direction_id": vector_id(labels, vector),
                "coefficient_vector": list(vector),
                "support_weight": sum(value != 0 for value in vector),
                "quartic_point": {
                    "x": rational_text(quartic_point[0]),
                    "z": rational_text(quartic_point[1]),
                },
                "quartic_x_projective_height": projective_height(quartic_point[0]),
                "canonical_group_point": point_record(canonical_sum),
                "exact_pointed_quartic_round_trip": True,
                "exact_canonical_group_relation": True,
                "mod5_outside_generic_span": True,
            }
        )
    records.sort(key=lambda row: row["direction_id"])
    lines = (
        f"{row['direction_id']}|{','.join(map(str, row['coefficient_vector']))}|"
        f"{row['quartic_point']['x']}|{row['quartic_point']['z']}|"
        f"{row['canonical_group_point']['x']}|{row['canonical_group_point']['y']}"
        for row in records
    )
    summary = {
        "anchor": anchor_name,
        "family_id": FAMILY_ID,
        "canonical_parameter_u": rational_text(canonical_u),
        "aliases": {"literal_shift_T": rational_text(literal_t)},
        "generic_basis_dimension": 12,
        "exceptional_basis_labels": list(labels),
        "exceptional_basis_dimension": len(labels),
        "coefficient_alphabet": [-1, 0, 1],
        "maximum_support_weight": 2,
        "global_sign_quotient": False,
        "signed_direction_count": len(records),
        "support_weight_histogram": dict(
            sorted(Counter(row["support_weight"] for row in records).items())
        ),
        "unique_oriented_quartic_point_count": len(seen_quartic),
        "unique_canonical_group_point_count": len(seen_canonical),
        "all_exact_group_relations_verified": True,
        "all_mod5_separated_from_generic_span": True,
        "direction_manifest_sha256": sha256_lines(lines),
        "minimum_quartic_x_projective_height": min(
            row["quartic_x_projective_height"] for row in records
        ),
        "maximum_quartic_x_projective_height": max(
            row["quartic_x_projective_height"] for row in records
        ),
        "records": records,
    }
    return summary, tuple(records)


def primitive_integer_coefficients(
    coefficients: Sequence[Fraction],
) -> tuple[int, ...]:
    common_denominator = 1
    for value in coefficients:
        common_denominator = lcm(common_denominator, Q(value).denominator)
    integers = [int(Q(value) * common_denominator) for value in coefficients]
    content = 0
    for value in integers:
        content = gcd(content, abs(value))
    if content == 0:
        raise ValueError("zero polynomial")
    integers = [value // content for value in integers]
    while len(integers) > 1 and integers[-1] == 0:
        integers.pop()
    if integers[-1] < 0:
        integers = [-value for value in integers]
    return tuple(integers)


def affine_slice_coefficients(
    left_x: Fraction, right_x: Fraction
) -> tuple[int, ...]:
    slope = (Q(right_x) - Q(left_x)) / (RANK20_T - E22_T)
    intercept = Q(left_x) - slope * E22_T
    linear = (intercept, slope)
    power: tuple[Fraction, ...] = (Q(1),)
    result: tuple[Fraction, ...] = (Q(0),)
    for coefficient_polynomial in FERMIGIER_BIVARIATE_COEFFICIENTS:
        result = poly_add(result, poly_multiply(coefficient_polynomial, power))
        power = poly_multiply(power, linear)
    if poly_evaluate(result, E22_T) != FermigierMestreFamily.quartic_value(
        E22_T, Q(left_x)
    ):
        raise AssertionError("an affine transport missed the E22 endpoint")
    if poly_evaluate(result, RANK20_T) != FermigierMestreFamily.quartic_value(
        RANK20_T, Q(right_x)
    ):
        raise AssertionError("an affine transport missed the rank-20 endpoint")
    return primitive_integer_coefficients(result)


def _mod_trim(polynomial: Sequence[int], prime: int) -> list[int]:
    answer = [value % prime for value in polynomial]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def _mod_divmod(
    numerator: Sequence[int], denominator: Sequence[int], prime: int
) -> tuple[list[int], list[int]]:
    remainder = _mod_trim(numerator, prime)
    denominator = _mod_trim(denominator, prime)
    if denominator == [0]:
        raise ZeroDivisionError("polynomial division by zero")
    quotient = [0] * max(1, len(remainder) - len(denominator) + 1)
    inverse_lead = pow(denominator[-1], -1, prime)
    while remainder != [0] and len(remainder) >= len(denominator):
        degree = len(remainder) - len(denominator)
        coefficient = remainder[-1] * inverse_lead % prime
        quotient[degree] = coefficient
        for index, value in enumerate(denominator):
            remainder[index + degree] = (
                remainder[index + degree] - coefficient * value
            ) % prime
        remainder = _mod_trim(remainder, prime)
    return _mod_trim(quotient, prime), remainder


def _mod_gcd(left: Sequence[int], right: Sequence[int], prime: int) -> list[int]:
    left = _mod_trim(left, prime)
    right = _mod_trim(right, prime)
    while right != [0]:
        _, remainder = _mod_divmod(left, right, prime)
        left, right = right, remainder
    inverse = pow(left[-1], -1, prime)
    return [(value * inverse) % prime for value in left]


def _mod_multiply_reduce(
    left: Sequence[int],
    right: Sequence[int],
    modulus: Sequence[int],
    prime: int,
) -> list[int]:
    product_coefficients = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            product_coefficients[left_index + right_index] = (
                product_coefficients[left_index + right_index]
                + left_value * right_value
            ) % prime
    return _mod_divmod(product_coefficients, modulus, prime)[1]


def _mod_power(
    base: Sequence[int], exponent: int, modulus: Sequence[int], prime: int
) -> list[int]:
    answer = [1]
    power = _mod_trim(base, prime)
    while exponent:
        if exponent & 1:
            answer = _mod_multiply_reduce(answer, power, modulus, prime)
        power = _mod_multiply_reduce(power, power, modulus, prime)
        exponent >>= 1
    return answer


def irreducible_mod_prime(coefficients: Sequence[int], prime: int) -> bool:
    """Rabin irreducibility test for an ascending integral polynomial."""

    polynomial = _mod_trim(coefficients, prime)
    degree = len(polynomial) - 1
    if degree <= 0 or coefficients[-1] % prime == 0:
        return False
    inverse = pow(polynomial[-1], -1, prime)
    polynomial = [(value * inverse) % prime for value in polynomial]
    x_polynomial = [0, 1]
    frobenius = x_polynomial
    checks = {degree // divisor for divisor in (2, 3) if degree % divisor == 0}
    for iteration in range(1, degree + 1):
        frobenius = _mod_power(frobenius, prime, polynomial, prime)
        if iteration in checks:
            difference_length = max(len(frobenius), 2)
            difference = [0] * difference_length
            for index, value in enumerate(frobenius):
                difference[index] = value
            difference[1] = (difference[1] - 1) % prime
            if len(_mod_gcd(polynomial, difference, prime)) > 1:
                return False
    return _mod_trim(frobenius, prime) == x_polynomial


def hyperelliptic_genus(squareclass_degree: int) -> int:
    return max(0, (squareclass_degree - 1) // 2)


def polynomial_sha256(coefficients: Sequence[int]) -> str:
    return sha256_lines((",".join(map(str, coefficients)),))


def _audit_affine_pair(
    row: tuple[str, str, str, str]
) -> dict[str, Any]:
    left_id, left_x_text, right_id, right_x_text = row
    coefficients = affine_slice_coefficients(Q(left_x_text), Q(right_x_text))
    degree = len(coefficients) - 1
    witness = next(
        (
            prime
            for prime in MODULAR_WITNESS_PRIMES
            if irreducible_mod_prime(coefficients, prime)
        ),
        None,
    )
    exact_factorization_used = witness is None
    if witness is not None:
        signature = ((degree, 1),)
    else:
        symbol = sp.symbols("T")
        polynomial = sp.Poly(
            sum(value * symbol**index for index, value in enumerate(coefficients)),
            symbol,
            domain=sp.ZZ,
        )
        signature = tuple(
            (factor.degree(), exponent)
            for factor, exponent in sp.factor_list(polynomial)[1]
        )
    squareclass_degree = sum(
        factor_degree
        for factor_degree, exponent in signature
        if exponent % 2
    )
    return {
        "left_direction_id": left_id,
        "right_direction_id": right_id,
        "degree": degree,
        "factor_signature": [list(item) for item in signature],
        "irreducible_mod_prime_witness": witness,
        "exact_QQ_factorization_fallback": exact_factorization_used,
        "squareclass_kernel_degree": squareclass_degree,
        "squareclass_genus": hyperelliptic_genus(squareclass_degree),
        "primitive_polynomial_sha256": polynomial_sha256(coefficients),
    }


def affine_transport_audit(
    left_records: Sequence[dict[str, Any]],
    right_records: Sequence[dict[str, Any]],
    *,
    workers: int,
) -> dict[str, Any]:
    if workers <= 0:
        raise ValueError("workers must be positive")
    rows = [
        (
            left["direction_id"],
            left["quartic_point"]["x"],
            right["direction_id"],
            right["quartic_point"]["x"],
        )
        for left in left_records
        for right in right_records
    ]
    if workers == 1:
        records = [_audit_affine_pair(row) for row in rows]
    else:
        start_method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
        with mp.get_context(start_method).Pool(workers) as pool:
            records = pool.map(_audit_affine_pair, rows, chunksize=32)
    records.sort(
        key=lambda record: (
            record["left_direction_id"], record["right_direction_id"]
        )
    )
    histogram = Counter(
        (
            record["degree"],
            tuple(map(tuple, record["factor_signature"])),
            record["squareclass_kernel_degree"],
            record["squareclass_genus"],
        )
        for record in records
    )
    witness_histogram = Counter(
        record["irreducible_mod_prime_witness"] for record in records
    )
    manifest_lines = (
        f"{record['left_direction_id']}|{record['right_direction_id']}|"
        f"{record['degree']}|{record['factor_signature']}|"
        f"{record['irreducible_mod_prime_witness']}|"
        f"{record['squareclass_kernel_degree']}|"
        f"{record['primitive_polynomial_sha256']}"
        for record in records
    )
    return {
        "definition": (
            "the unique affine x(T)=a*T+b through each signed exceptional "
            "quotient-ball point at both anchors"
        ),
        "left_direction_count": len(left_records),
        "right_direction_count": len(right_records),
        "pair_count": len(records),
        "classification_method": (
            "primitive integral polynomial; irreducible reduction modulo the "
            "recorded prime proves irreducibility over QQ, otherwise exact QQ factorization"
        ),
        "modular_witness_prime_set": list(MODULAR_WITNESS_PRIMES),
        "modular_witness_histogram": {
            str(key): value
            for key, value in sorted(
                witness_histogram.items(), key=lambda item: (-1 if item[0] is None else item[0])
            )
        },
        "exact_QQ_factorization_fallback_count": sum(
            record["exact_QQ_factorization_fallback"] for record in records
        ),
        "factor_squareclass_histogram": [
            {
                "degree": key[0],
                "factor_signature": [list(item) for item in key[1]],
                "squareclass_kernel_degree": key[2],
                "squareclass_genus": key[3],
                "count": count,
            }
            for key, count in sorted(histogram.items())
        ],
        "low_genus_candidates": [
            record for record in records if record["squareclass_genus"] <= 1
        ],
        "manifest_sha256": sha256_lines(manifest_lines),
        "records": records,
    }


def result_digest(artifact: dict[str, Any]) -> str:
    stable = {
        "anchors": artifact["anchors"],
        "direction_balls": artifact["direction_balls"],
        "affine_transport": artifact["affine_transport"],
        "outcome": artifact["outcome"],
        "sources": {
            key: value
            for key, value in artifact["sources"].items()
            if key != "script_sha256"
        },
    }
    return hashlib.sha256(
        json.dumps(stable, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def run(root: Path, *, workers: int) -> dict[str, Any]:
    rank22_source = root / RANK22_ARTIFACT_RELATIVE
    if sha256_file(rank22_source) != RANK22_ARTIFACT_SHA256:
        raise AssertionError("the pinned rank-22 source artifact changed")

    e22 = e22_basis()
    rank20 = rank20_basis(root)
    e22_summary, e22_records = direction_ball(
        anchor_name="E22",
        canonical_u=E22_U,
        literal_t=E22_T,
        specialization=e22[0],
        labels=e22[1],
        quartic_basis=e22[2],
        canonical_basis=e22[3],
        certificate=e22[4],
    )
    rank20_summary, rank20_records = direction_ball(
        anchor_name="rank20",
        canonical_u=RANK20_U,
        literal_t=RANK20_T,
        specialization=rank20[0],
        labels=rank20[1],
        quartic_basis=rank20[2],
        canonical_basis=rank20[3],
        certificate=rank20[4],
    )
    affine = affine_transport_audit(e22_records, rank20_records, workers=workers)
    outcome = {
        "signed_E22_exceptional_directions": len(e22_records),
        "signed_rank20_exceptional_directions": len(rank20_records),
        "cross_anchor_affine_interpolants": affine["pair_count"],
        "affine_low_genus_candidates": len(affine["low_genus_candidates"]),
        "new_base_changes": 0,
        "new_sections": 0,
        "new_specializations": 0,
        "target_met": False,
    }
    artifact: dict[str, Any] = {
        "schema_version": "elliptic-curves.fermigier-exceptional-quotient-ball.v1",
        "status": "complete finite exact quotient-ball classification",
        "claim_level": "exact computation; no new low-genus transport found",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "anchors": {
            "E22": {
                "family_id": FAMILY_ID,
                "canonical_parameter": {"u": rational_text(E22_U)},
                "aliases": {"literal_shift_T": rational_text(E22_T)},
                "certified_rank_lower_bound": 22,
            },
            "rank20": {
                "family_id": FAMILY_ID,
                "canonical_parameter": {"u": rational_text(RANK20_U)},
                "aliases": {"literal_shift_T": rational_text(RANK20_T)},
                "certified_rank_lower_bound": 20,
            },
        },
        "direction_balls": {"E22": e22_summary, "rank20": rank20_summary},
        "affine_transport": affine,
        "outcome": outcome,
        "sources": {
            "rank22_artifact": str(RANK22_ARTIFACT_RELATIVE),
            "rank22_artifact_sha256": RANK22_ARTIFACT_SHA256,
            "rank20_artifact": str(RANK20_ARTIFACT_RELATIVE),
            "rank20_artifact_sha256": RANK20_ARTIFACT_SHA256,
            "script_sha256": sha256_file(Path(__file__)),
        },
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
        },
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves:elliptic-curves/cas .venv/bin/python "
            "elliptic-curves/cas/classify_fermigier_exceptional_quotient_ball.py"
        ),
    }
    artifact["result_sha256"] = result_digest(artifact)
    return artifact


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts/generated-results/elliptic-curves/elliptic_fermigier_exceptional_quotient_ball.json",
    )
    args = parser.parse_args()
    artifact = run(root, workers=args.workers)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps(artifact["outcome"], sort_keys=True))
    print(f"result_sha256={artifact['result_sha256']}")


if __name__ == "__main__":
    main()
