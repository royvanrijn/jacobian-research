#!/usr/bin/env python3
"""Prove the arithmetic generic rank of the Fermigier--Mestre K3 surface.

The canonical parameter in this verifier is the repository adapter ``u``.
The literal symmetric shift in Fermigier's product is ``s=2*u``; that bridge
is checked coefficientwise through the binary-quartic invariants and through
the two primitive degree-20 discriminants.

The canonical elliptic surface has twenty geometric ``I1`` fibers and one
split ``I4`` fiber at infinity.  Its trivial lattice therefore has rank five.
The frozen exact finite-reduction certificate for the twelve rational section
differences is replayed from first principles, giving seventeen independent
Q-defined divisor classes.

At the good prime 41, this script counts the resolved elliptic K3 surface over
both F_41 and F_(41^2).  Removing the seventeen known ``+41`` eigenvalues
from H^2 leaves a reciprocal degree-five factor.  Exact Weil reconstruction
gives

  (X+41) (X^4 + 94 X^3 + 4428 X^2 + 158014 X + 41^4).

The quartic has no root ``41*zeta`` for a root of unity zeta, and the whole
residual factor has no root ``+41``.  Smooth proper specialization and the
cycle-class eigenvalue condition, without the Tate conjecture, then prove
that no thirteenth section is defined over Q(u).  Thus the arithmetic generic
Mordell--Weil rank is exactly 12.  Geometrically, the possible ``-41`` class
leaves the unconditional interval 12 <= rank over Qbar(u) <= 13.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import gcd
import os
from pathlib import Path
import platform
import sys
from typing import Any, Iterable, Sequence

import sympy as sp


SCRIPT_PATH = Path(__file__).resolve()
REPOSITORY_ROOT = SCRIPT_PATH.parents[2]
PROGRAM_ROOT = SCRIPT_PATH.parents[1]
sys.path.insert(0, str(PROGRAM_ROOT))
sys.path.insert(0, str(PROGRAM_ROOT / "cas"))

from ecsearch.fermigier import (  # noqa: E402
    FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS,
    FERMIGIER_LITERAL_DISCRIMINANT_FACTOR_COEFFICIENTS,
    FERMIGIER_REPORTED_PARAMETER,
    fermigier_canonical_coefficients,
)
from ecsearch.fermigier_rank import specialize_fermigier_rank_sections  # noqa: E402
from ecsearch.rank_certification import (  # noqa: E402
    IndependenceCertificate,
    matrix_rank_mod_prime,
    verify_independence_certificate,
)
from fermigier_mestre import FermigierMestreFamily  # noqa: E402


Q = Fraction
GOOD_PRIME = 41
KNOWN_SECTION_RANK = 12
TRIVIAL_LATTICE_RANK = 5
KNOWN_RATIONAL_DIVISOR_RANK = KNOWN_SECTION_RANK + TRIVIAL_LATTICE_RANK
EXPECTED_CHARACTER_SUMS = {1: 358, 2: 21_834}
EXPECTED_SURFACE_POINT_COUNTS = {1: 2_244, 2: 2_856_000}
EXPECTED_H2_TRACES = {1: 562, 2: 30_238}
EXPECTED_RESIDUAL_TRACES = {1: -135, 2: 1_661}
ROOT_OF_UNITY_ORDERS_OF_DEGREE_AT_MOST_FOUR = (1, 2, 3, 4, 5, 6, 8, 10, 12)
SECTION_CERTIFICATE = (
    REPOSITORY_ROOT
    / "artifacts/generated-results/elliptic-curves/fermigier_rank_certificates_v1.json"
)
EXPECTED_SECTION_CERTIFICATE_SHA256 = (
    "94fc64d7f1744f6a20a0396d32914cd36330107db2538e03ee95cc3e32927051"
)
DEFAULT_OUTPUT = (
    REPOSITORY_ROOT
    / "artifacts/generated-results/elliptic-curves/elliptic_fermigier_generic_rank_exact.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas:elliptic-curves .venv/bin/python "
    "elliptic-curves/cas/verify_fermigier_generic_rank_exact.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def canonical_polynomials() -> dict[str, sp.Poly]:
    """Return the exact canonical Weierstrass invariants over Z[u]."""

    u = sp.symbols("u")
    a1 = sp.Integer(1)
    a2 = -8 * u**4 + 1_718_550 * u**2 + 298_803_565_660
    a3 = sp.Integer(1)
    a4 = (
        -64 * u**8
        - 18_151_200 * u**6
        + 1_028_009_011_008 * u**4
        + 317_946_481_466_562_025 * u**2
        + 27_856_983_036_916_830_925_012
    )
    a6 = (
        512 * u**12
        + 35_222_400 * u**10
        - 31_029_122_010_320 * u**8
        + 573_566_223_236_700_400 * u**6
        + 109_667_821_527_431_621_677_482 * u**4
        + 14_321_756_340_366_264_921_294_086_000 * u**2
        + 829_998_138_277_457_737_118_423_455_411_406
    )
    b2 = sp.expand(a1**2 + 4 * a2)
    b4 = sp.expand(2 * a4 + a1 * a3)
    b6 = sp.expand(a3**2 + 4 * a6)
    b8 = sp.expand(
        a1**2 * a6
        + 4 * a2 * a6
        - a1 * a3 * a4
        + a2 * a3**2
        - a4**2
    )
    c4 = sp.expand(b2**2 - 24 * b4)
    c6 = sp.expand(-b2**3 + 36 * b2 * b4 - 216 * b6)
    discriminant = sp.expand(
        -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    )
    answer = {
        "a1": sp.Poly(a1, u, domain=sp.ZZ),
        "a2": sp.Poly(a2, u, domain=sp.ZZ),
        "a3": sp.Poly(a3, u, domain=sp.ZZ),
        "a4": sp.Poly(a4, u, domain=sp.ZZ),
        "a6": sp.Poly(a6, u, domain=sp.ZZ),
        "c4": sp.Poly(c4, u, domain=sp.ZZ),
        "c6": sp.Poly(c6, u, domain=sp.ZZ),
        "discriminant": sp.Poly(discriminant, u, domain=sp.ZZ),
    }
    expected_degrees = {
        "a1": 0, "a2": 4, "a3": 0, "a4": 8, "a6": 12,
        "c4": 8, "c6": 12, "discriminant": 20,
    }
    if {name: polynomial.degree() for name, polynomial in answer.items()} != expected_degrees:
        raise AssertionError("the canonical Fermigier K3 degrees changed")
    observed_delta = tuple(
        int(answer["discriminant"].coeff_monomial(u**degree))
        for degree in range(21)
    )
    if observed_delta != FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS:
        raise AssertionError("the canonical primitive discriminant changed")
    return answer


def exact_adapter_bridge(polynomials: dict[str, sp.Poly]) -> dict[str, Any]:
    """Check coefficientwise that literal shift s equals 2u."""

    u = polynomials["a2"].gens[0]
    literal_shift = 2 * u
    # FermigierMestreFamily stores the quartic after dividing the raw
    # remainder by (50616*s)^2.  Its invariant I is c4, while J is 2*c6.
    a = literal_shift**2 + 1_149_050
    b = -30 * (62 * literal_shift**2 + 68_377_393)
    c = -(
        2 * literal_shift**4
        - 1_718_550 * literal_shift**2
        - 1_195_214_262_641
    )
    d = 30 * (
        62 * literal_shift**4
        - 21_690_305 * literal_shift**2
        - 8_594_794_400_346
    )
    e = (
        literal_shift**6
        - 879_500 * literal_shift**4
        + 102_302_344_648 * literal_shift**2
        + 18_103_855_887_324_900
    )
    invariant_i = sp.expand(12 * a * e - 3 * b * d + c**2)
    invariant_j = sp.expand(
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d**2
        - 27 * b**2 * e
        - 2 * c**3
    )
    if sp.expand(invariant_i - polynomials["c4"].as_expr()) != 0:
        raise AssertionError("the normalized quartic I/c4 bridge changed")
    if sp.expand(invariant_j - 2 * polynomials["c6"].as_expr()) != 0:
        raise AssertionError("the normalized quartic J/c6 bridge changed")

    literal = sum(
        coefficient * literal_shift**degree
        for degree, coefficient in enumerate(
            FERMIGIER_LITERAL_DISCRIMINANT_FACTOR_COEFFICIENTS
        )
    )
    canonical = polynomials["discriminant"].as_expr()
    if sp.expand(literal - 16 * canonical) != 0:
        raise AssertionError("Psi(2u)=16*Phi(u) changed")

    # Independent exact API checks at several rational parameters protect the
    # symbolic transcription of the canonical coefficient polynomials.
    for parameter in (Q(1), Q(7, 3), FERMIGIER_REPORTED_PARAMETER):
        api = fermigier_canonical_coefficients(parameter)
        observed = tuple(
            Q(polynomial.eval(Q(parameter)))
            for polynomial in (
                polynomials["a1"], polynomials["a2"], polynomials["a3"],
                polynomials["a4"], polynomials["a6"],
            )
        )
        if observed != api:
            raise AssertionError("the canonical adapter API bridge changed")
        normalized_i, normalized_j = FermigierMestreFamily.binary_invariants(
            2 * parameter
        )
        if normalized_i != Q(polynomials["c4"].eval(parameter)):
            raise AssertionError("the rational I bridge replay failed")
        if normalized_j != 2 * Q(polynomials["c6"].eval(parameter)):
            raise AssertionError("the rational J bridge replay failed")

    return {
        "canonical_parameter": "u",
        "literal_symmetric_shift": "s=2*u",
        "canonical_model": "[1,a2(u),1,a4(u),a6(u)]",
        "normalized_quartic_invariants": {"I": "c4(C_u)", "J": "2*c6(C_u)"},
        "raw_quartic_scale": "(50616*s)^2=(101232*u)^2",
        "raw_invariant_identities": {
            "I_raw": "101232^4*u^4*c4(C_u)",
            "J_raw": "2*101232^6*u^6*c6(C_u)",
        },
        "primitive_discriminant_identity": "Psi(2*u)=16*Phi(u)",
        "all_identities_checked_coefficientwise": True,
        "factor_two_source_discrepancy_resolved_here": False,
    }


def fiber_inventory(polynomials: dict[str, sp.Poly]) -> dict[str, Any]:
    """Prove the 20 I1 plus split I4 fiber configuration over Qbar."""

    discriminant = polynomials["discriminant"]
    c4 = polynomials["c4"]
    if sp.gcd(discriminant, discriminant.diff()).degree() != 0:
        raise AssertionError("the characteristic-zero discriminant is not squarefree")
    if sp.gcd(discriminant, c4).degree() != 0:
        raise AssertionError("a finite discriminant root has c4=0")

    x = sp.symbols("X")
    infinity_cubic = sp.Poly(x**3 - 8 * x**2 - 64 * x + 512, x, domain=sp.ZZ)
    if sp.factor(infinity_cubic.as_expr()) != (x - 8) ** 2 * (x + 8):
        raise AssertionError("the infinity nodal cubic changed")
    infinity_c4 = 4_096
    leading_delta = int(discriminant.LC())
    if leading_delta != 1_322_860_664_586_240_000:
        raise AssertionError("the order-four infinity discriminant coefficient changed")
    # With v=1/u and x=v^-4 X, y=v^-6 Y, the scaled discriminant is
    # v^24*Phi(1/v), hence has order 24-20=4.  At the double root X=8,
    # Y^2=(X-8)^2(X+8) has tangent cone Y^2=16(X-8)^2.
    return {
        "surface_type": "elliptic K3 over P1_u",
        "fundamental_line_bundle_degree": 2,
        "weierstrass_coefficient_degrees": {
            "a1": 0, "a2": 4, "a3": 0, "a4": 8, "a6": 12,
        },
        "finite_discriminant": {
            "degree": 20,
            "squarefree_over_Q": True,
            "coprime_to_c4_over_Q": True,
            "geometric_fibers": "20 fibers of Kodaira type I1",
        },
        "infinity_fiber": {
            "scaled_special_cubic": str(infinity_cubic.as_expr()),
            "factorization": "(X-8)^2*(X+8)",
            "scaled_c4_at_infinity": infinity_c4,
            "scaled_discriminant_order": 4,
            "tangent_slopes_at_node": [-4, 4],
            "split_over_Q": True,
            "Kodaira_type": "I4",
        },
        "Euler_number": 24,
        "trivial_lattice": {
            "zero_section_and_fiber_rank": 2,
            "I4_root_lattice": "A3",
            "I4_component_rank": 3,
            "rank": TRIVIAL_LATTICE_RANK,
            "all_generators_defined_over_Q": True,
        },
    }


def fraction_point_digest(points: Iterable[tuple[Fraction, Fraction]]) -> str:
    return hashlib.sha256(
        "".join(f"{x}\t{y}\n" for x, y in points).encode()
    ).hexdigest()


def verify_twelve_sections() -> dict[str, Any]:
    """Replay the exact specialized certificate proving 12 generic sections."""

    actual_sha = sha256_file(SECTION_CERTIFICATE)
    if actual_sha != EXPECTED_SECTION_CERTIFICATE_SHA256:
        raise AssertionError("the frozen Fermigier section certificate changed")
    source = json.loads(SECTION_CERTIFICATE.read_text())
    generic = source["generic_sections"]
    parameter = Q(generic["adapter_parameter"])
    if parameter != FERMIGIER_REPORTED_PARAMETER:
        raise AssertionError("the section-certificate specialization changed")
    specialization = specialize_fermigier_rank_sections(parameter)
    if specialization.quartic_model.shift != 2 * parameter:
        raise AssertionError("the section specialization lost the exact s=2u bridge")
    if len(specialization.quartic_points) != 13:
        raise AssertionError("the thirteen rational quartic sections changed")
    points = specialization.section_differences
    if len(points) != KNOWN_SECTION_RANK:
        raise AssertionError("the section-difference population changed")
    certificate = IndependenceCertificate.from_json_object(generic["certificate"])
    verify_independence_certificate(
        specialization.canonical_model, points, certificate
    )
    matrix_rank = matrix_rank_mod_prime(
        tuple(row.logs for row in certificate.rows), certificate.relation_prime
    )
    if matrix_rank != KNOWN_SECTION_RANK:
        raise AssertionError("the section finite-reduction matrix lost full rank")
    return {
        "status": "exact finite-reduction independence certificate replayed",
        "adapter_parameter_used_only_for_specialization": str(parameter),
        "literal_shift": str(specialization.quartic_model.shift),
        "quartic_section_count": len(specialization.quartic_points),
        "section_difference_count": len(points),
        "section_difference_sha256": fraction_point_digest(points),
        "relation_prime": certificate.relation_prime,
        "torsion_witness": {
            "prime": certificate.torsion_witness_prime,
            "group_order": certificate.torsion_witness_group_order,
        },
        "certificate_primes": [row.prime for row in certificate.rows],
        "combined_exact_matrix_rank": matrix_rank,
        "certificate": certificate.to_json_object(),
        "generic_inference": (
            "the twelve points are specializations of rational sections; any generic "
            "integral relation would specialize, so exact specialized independence "
            "proves rank at least twelve over Q(u)"
        ),
        "certified_generic_Mordell_Weil_rank_lower_bound": KNOWN_SECTION_RANK,
    }


def _ascending_coefficients_mod(polynomial: sp.Poly, prime: int) -> tuple[int, ...]:
    variable = polynomial.gens[0]
    return tuple(
        int(polynomial.coeff_monomial(variable**degree)) % prime
        for degree in range(polynomial.degree() + 1)
    )


def _evaluate_prime(coefficients: Sequence[int], value: int, prime: int) -> int:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = (answer * value + coefficient) % prime
    return answer


def verify_good_reduction(
    polynomials: dict[str, sp.Poly], prime: int = GOOD_PRIME
) -> dict[str, Any]:
    """Check semistable good reduction of the resolved K3 at p=41."""

    if prime != GOOD_PRIME:
        raise ValueError("this exact certificate is pinned at p=41")
    reduced_delta = sp.Poly(polynomials["discriminant"], modulus=prime)
    reduced_c4 = sp.Poly(polynomials["c4"], modulus=prime)
    if reduced_delta.degree() != 20 or reduced_c4.degree() != 8:
        raise AssertionError("a leading invariant coefficient vanished modulo 41")
    if sp.gcd(reduced_delta, reduced_delta.diff()).degree() != 0:
        raise AssertionError("the finite discriminant is not squarefree modulo 41")
    if sp.gcd(reduced_delta, reduced_c4).degree() != 0:
        raise AssertionError("a finite mod-41 singular fiber is not I1")
    infinity_c4 = 4_096 % prime
    infinity_delta_lead = int(polynomials["discriminant"].LC()) % prime
    if infinity_c4 == 0 or infinity_delta_lead == 0:
        raise AssertionError("the infinity I4 fiber degenerated modulo 41")
    if (16 % prime) == 0 or (4 % prime) == 0:
        raise AssertionError("the infinity split node degenerated modulo 41")
    return {
        "prime": prime,
        "finite_discriminant_degree": 20,
        "finite_discriminant_squarefree": True,
        "finite_discriminant_coprime_to_c4": True,
        "finite_fibers": "20 geometrically simple I1 fibers",
        "infinity_scaled_c4_mod_p": infinity_c4,
        "infinity_scaled_discriminant_leading_coefficient_mod_p": (
            infinity_delta_lead
        ),
        "infinity_fiber": "split I4 with distinct roots 8,-8 and slopes +/-4",
        "good_reduction_of_resolved_K3": True,
    }


def finite_affine_character_sum(
    polynomials: dict[str, sp.Poly], prime: int, extension_degree: int
) -> int:
    """Return the exact character sum for all finite u and affine x.

    For the general model ``[1,a2,1,a4,a6]``, completing the square gives

      (2y+x+1)^2 = (x+1)^2 + 4(x^3+a2*x^2+a4*x+a6).

    Over F_(p^2), the quadratic character of a nonzero element is the
    Legendre symbol of its norm to F_p.
    """

    coefficients_a2 = _ascending_coefficients_mod(polynomials["a2"], prime)
    coefficients_a4 = _ascending_coefficients_mod(polynomials["a4"], prime)
    coefficients_a6 = _ascending_coefficients_mod(polynomials["a6"], prime)
    nonzero_squares = {value * value % prime for value in range(1, prime)}

    if extension_degree == 1:
        total = 0
        for parameter in range(prime):
            a2 = _evaluate_prime(coefficients_a2, parameter, prime)
            a4 = _evaluate_prime(coefficients_a4, parameter, prime)
            a6 = _evaluate_prime(coefficients_a6, parameter, prime)
            for x in range(prime):
                square_discriminant = (
                    (x + 1) ** 2
                    + 4 * (x**3 + a2 * x**2 + a4 * x + a6)
                ) % prime
                if square_discriminant:
                    total += 1 if square_discriminant in nonzero_squares else -1
        return total
    if extension_degree != 2:
        raise ValueError("the verifier counts only F_p and F_(p^2)")

    nonsquare = next(
        value
        for value in range(2, prime)
        if value not in nonzero_squares
    )
    elements = tuple(
        (real, imaginary)
        for real in range(prime)
        for imaginary in range(prime)
    )

    def add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
        return ((left[0] + right[0]) % prime, (left[1] + right[1]) % prime)

    def multiply(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
        return (
            (left[0] * right[0] + nonsquare * left[1] * right[1]) % prime,
            (left[0] * right[1] + left[1] * right[0]) % prime,
        )

    def scale(integer: int, value: tuple[int, int]) -> tuple[int, int]:
        return (integer * value[0] % prime, integer * value[1] % prime)

    def evaluate(
        coefficients: Sequence[int], value: tuple[int, int]
    ) -> tuple[int, int]:
        result = (0, 0)
        for coefficient in reversed(coefficients):
            result = add(multiply(result, value), (coefficient, 0))
        return result

    x_data = []
    for x in elements:
        square = multiply(x, x)
        cube = multiply(square, x)
        linear = add(x, (1, 0))
        x_data.append((x, square, cube, multiply(linear, linear)))

    total = 0
    for parameter in elements:
        a2 = evaluate(coefficients_a2, parameter)
        a4 = evaluate(coefficients_a4, parameter)
        a6 = evaluate(coefficients_a6, parameter)
        for x, square, cube, linear_square in x_data:
            right = add(add(add(cube, multiply(a2, square)), multiply(a4, x)), a6)
            square_discriminant = add(linear_square, scale(4, right))
            if square_discriminant == (0, 0):
                continue
            norm = (
                square_discriminant[0] ** 2
                - nonsquare * square_discriminant[1] ** 2
            ) % prime
            total += 1 if norm in nonzero_squares else -1
    return total


def surface_point_count(
    polynomials: dict[str, sp.Poly], prime: int, extension_degree: int
) -> dict[str, int]:
    """Count the smooth resolved K3 surface over F_(p^n)."""

    q = prime**extension_degree
    character_sum = finite_affine_character_sum(
        polynomials, prime, extension_degree
    )
    # The q finite fibers contribute q(q+1)+sum chi.  The split I4 fiber at
    # infinity is a cycle of four P1s and contributes 4(q+1)-4=4q points.
    point_count = q * q + 5 * q + character_sum
    h2_trace = point_count - 1 - q * q
    return {
        "field_size": q,
        "finite_affine_character_sum": character_sum,
        "split_I4_infinity_points": 4 * q,
        "surface_point_count": point_count,
        "h2_frobenius_trace": h2_trace,
    }


def reconstruct_h2_frobenius(
    prime: int, trace_degree_one: int, trace_degree_two: int
) -> dict[str, Any]:
    """Reconstruct the full H2 polynomial after removing 17 known classes."""

    residual_trace_one = trace_degree_one - KNOWN_RATIONAL_DIVISOR_RANK * prime
    residual_trace_two = (
        trace_degree_two - KNOWN_RATIONAL_DIVISOR_RANK * prime**2
    )
    if (
        residual_trace_one != EXPECTED_RESIDUAL_TRACES[1]
        or residual_trace_two != EXPECTED_RESIDUAL_TRACES[2]
    ):
        raise AssertionError("the residual Frobenius traces changed")

    sign_records = []
    valid_signs = []
    pair_square_sum = residual_trace_two + 3 * prime**2
    for sign in (-1, 1):
        pair_sum = residual_trace_one - sign * prime
        product_numerator = pair_sum**2 - pair_square_sum
        pair_product = None
        discriminant = None
        valid = product_numerator % 2 == 0
        if valid:
            pair_product = product_numerator // 2
            discriminant = pair_sum**2 - 4 * pair_product
            lower = pair_sum + 4 * prime
            upper = 4 * prime - pair_sum
            valid = (
                discriminant >= 0
                and lower >= 0
                and upper >= 0
                and discriminant <= lower**2
                and discriminant <= upper**2
            )
        sign_records.append(
            {
                "real_eigenvalue_sign": sign,
                "pair_sum_U_plus_V": pair_sum,
                "pair_square_sum_U2_plus_V2": pair_square_sum,
                "pair_product_U_times_V": pair_product,
                "pair_sum_discriminant": discriminant,
                "satisfies_exact_Weil_interval": valid,
            }
        )
        if valid:
            valid_signs.append((sign, pair_sum, pair_product, discriminant))
    if valid_signs != [(-1, -94, 1_066, 4_572)]:
        raise AssertionError(f"the unique Weil reconstruction changed: {valid_signs}")

    sign, pair_sum, pair_product, discriminant = valid_signs[0]
    x = sp.symbols("X")
    quartic = sp.Poly(
        x**4
        - pair_sum * x**3
        + (pair_product + 2 * prime**2) * x**2
        - prime**2 * pair_sum * x
        + prime**4,
        x,
        domain=sp.ZZ,
    )
    residual = sp.Poly(
        sp.expand((x - sign * prime) * quartic.as_expr()), x, domain=sp.ZZ
    )
    full_h2 = sp.Poly(
        sp.expand((x - prime) ** KNOWN_RATIONAL_DIVISOR_RANK * residual.as_expr()),
        x,
        domain=sp.ZZ,
    )
    if quartic.as_expr() != (
        x**4 + 94 * x**3 + 4_428 * x**2 + 158_014 * x + 2_825_761
    ):
        raise AssertionError("the residual reciprocal quartic changed")
    if residual.eval(prime) != 2_136_275_316:
        raise AssertionError("the residual +p exclusion changed")
    if full_h2.degree() != 22:
        raise AssertionError("the reconstructed H2 polynomial lost degree 22")
    coefficients = [int(value) for value in full_h2.all_coeffs()]
    trace_one_replay = -coefficients[1]
    trace_two_replay = trace_one_replay**2 - 2 * coefficients[2]
    if (trace_one_replay, trace_two_replay) != (
        trace_degree_one, trace_degree_two
    ):
        raise AssertionError("Newton trace replay of the H2 polynomial failed")

    z = sp.symbols("Z")
    normalized_quartic = sp.Poly(
        sp.expand(quartic.as_expr().subs(x, prime * z) / prime**4),
        z,
        domain=sp.QQ,
    )
    cyclotomic_gcds = []
    for order in ROOT_OF_UNITY_ORDERS_OF_DEGREE_AT_MOST_FOUR:
        cyclotomic = sp.Poly(
            sp.cyclotomic_poly(order, z), z, domain=sp.QQ
        )
        common = sp.gcd(normalized_quartic, cyclotomic)
        cyclotomic_gcds.append({"order": order, "gcd_degree": common.degree()})
        if common.degree() != 0:
            raise AssertionError("the residual quartic gained p times a root of unity")

    return {
        "known_Q_defined_divisor_eigenvalues_removed": KNOWN_RATIONAL_DIVISOR_RANK,
        "known_factor": f"(X-{prime})^{KNOWN_RATIONAL_DIVISOR_RANK}",
        "residual_traces": {
            "degree_1": residual_trace_one,
            "degree_2": residual_trace_two,
        },
        "reciprocal_Weil_reconstruction": sign_records,
        "unique_real_residual_eigenvalue": str(sign * prime),
        "pair_sums_U_V": ["-47+3*sqrt(127)", "-47-3*sqrt(127)"],
        "residual_quartic": str(quartic.as_expr()),
        "residual_characteristic_polynomial": str(residual.as_expr()),
        "residual_characteristic_polynomial_at_positive_p": int(
            residual.eval(prime)
        ),
        "normalized_residual_quartic": str(normalized_quartic.as_expr()),
        "cyclotomic_gcd_checks": cyclotomic_gcds,
        "residual_p_times_root_of_unity_multiplicity_upper_bound": 1,
        "residual_positive_p_eigenvalue_multiplicity": 0,
        "full_H2_characteristic_polynomial": str(full_h2.as_expr()),
        "full_H2_coefficients_descending": [str(value) for value in coefficients],
        "full_H2_coefficients_sha256": stable_digest(coefficients),
        "full_positive_p_eigenvalue_multiplicity": (
            KNOWN_RATIONAL_DIVISOR_RANK
        ),
        "Newton_trace_replay": {
            "degree_1": trace_one_replay,
            "degree_2": trace_two_replay,
        },
        "Weil_absolute_values_and_reciprocity_checked_exactly": True,
    }


def verify_generic_rank() -> dict[str, Any]:
    """Run the full exact bridge, section, point-count, and rank proof."""

    polynomials = canonical_polynomials()
    bridge = exact_adapter_bridge(polynomials)
    fibers = fiber_inventory(polynomials)
    sections = verify_twelve_sections()
    good_reduction = verify_good_reduction(polynomials)
    counts = {
        degree: surface_point_count(polynomials, GOOD_PRIME, degree)
        for degree in (1, 2)
    }
    for degree in (1, 2):
        if counts[degree]["finite_affine_character_sum"] != EXPECTED_CHARACTER_SUMS[degree]:
            raise AssertionError("the finite affine character sum changed")
        if counts[degree]["surface_point_count"] != EXPECTED_SURFACE_POINT_COUNTS[degree]:
            raise AssertionError("the resolved K3 point count changed")
        if counts[degree]["h2_frobenius_trace"] != EXPECTED_H2_TRACES[degree]:
            raise AssertionError("the H2 Frobenius trace changed")
    frobenius = reconstruct_h2_frobenius(
        GOOD_PRIME,
        counts[1]["h2_frobenius_trace"],
        counts[2]["h2_frobenius_trace"],
    )
    arithmetic_rank = (
        frobenius["full_positive_p_eigenvalue_multiplicity"]
        - TRIVIAL_LATTICE_RANK
    )
    geometric_upper = (
        KNOWN_RATIONAL_DIVISOR_RANK
        + frobenius["residual_p_times_root_of_unity_multiplicity_upper_bound"]
        - TRIVIAL_LATTICE_RANK
    )
    if arithmetic_rank != 12 or geometric_upper != 13:
        raise AssertionError("the Shioda--Tate rank conclusion changed")
    return {
        "model_bridge": bridge,
        "fiber_inventory": fibers,
        "exact_generic_sections": sections,
        "good_reduction": good_reduction,
        "point_counts": {str(degree): counts[degree] for degree in (1, 2)},
        "H2_frobenius": frobenius,
        "theorem": {
            "arithmetic_generic_Mordell_Weil_rank_over_Q_of_u": arithmetic_rank,
            "arithmetic_rank_status": "proved unconditionally",
            "geometric_generic_Mordell_Weil_rank_interval_over_Qbar_of_u": [
                KNOWN_SECTION_RANK, geometric_upper
            ],
            "geometric_rank_status": (
                "proved unconditional interval; this computation does not decide 12 versus 13"
            ),
            "Q_defined_Neron_Severi_plus_p_eigenspace_dimension": (
                KNOWN_RATIONAL_DIVISOR_RANK
            ),
            "geometric_Picard_rank_upper_bound": (
                TRIVIAL_LATTICE_RANK + geometric_upper
            ),
            "Tate_conjecture_assumed": False,
            "argument": (
                "the zero section, fiber, three split-I4 components, and twelve exact "
                "sections give 17 Q-defined classes. Any further Q(u)-section would, by "
                "Shioda--Tate and smooth proper specialization, add a +41 eigenvector in "
                "H2. The reconstructed polynomial has +41 multiplicity exactly 17. A "
                "geometric divisor can only contribute 41 times a root of unity; the "
                "residual factor permits only the single eigenvalue -41."
            ),
        },
        "conditional_inferences": [],
        "scope_limits": [
            "does not decide whether the geometric generic rank is 12 or 13",
            "does not determine the rank of every specialization",
            "does not itself produce a rank-21 or rank-30 curve",
            "does not resolve Fermigier's printed factor-two parameter discrepancy",
        ],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def exclusive_write(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")


def main() -> None:
    args = build_parser().parse_args()
    if args.output.exists():
        raise SystemExit("refusing to overwrite the Fermigier generic-rank artifact")
    result = verify_generic_rank()
    adapter_source = PROGRAM_ROOT / "ecsearch/fermigier.py"
    section_source = PROGRAM_ROOT / "ecsearch/fermigier_rank.py"
    certificate_source = PROGRAM_ROOT / "ecsearch/rank_certification.py"
    artifact = {
        "schema_version": 1,
        "status": "exact arithmetic generic-rank theorem",
        **result,
        "provenance": {
            "script_path": str(SCRIPT_PATH.relative_to(REPOSITORY_ROOT)),
            "script_sha256": sha256_file(SCRIPT_PATH),
            "canonical_adapter_source_sha256": sha256_file(adapter_source),
            "section_source_sha256": sha256_file(section_source),
            "finite_certificate_verifier_sha256": sha256_file(certificate_source),
            "frozen_section_certificate": str(
                SECTION_CERTIFICATE.relative_to(REPOSITORY_ROOT)
            ),
            "frozen_section_certificate_sha256": (
                EXPECTED_SECTION_CERTIFICATE_SHA256
            ),
            "reproducing_command": REPRODUCING_COMMAND,
        },
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
        },
    }
    artifact["result_sha256"] = stable_digest(
        {
            "bridge": artifact["model_bridge"],
            "fibers": artifact["fiber_inventory"],
            "sections": artifact["exact_generic_sections"],
            "good_reduction": artifact["good_reduction"],
            "counts": artifact["point_counts"],
            "frobenius": artifact["H2_frobenius"],
            "theorem": artifact["theorem"],
        }
    )
    exclusive_write(args.output, artifact)
    print("#S(F_41)=2244; #S(F_41^2)=2856000")
    print("proved arithmetic generic Mordell--Weil rank over Q(u) equals 12")
    print("proved 12 <= geometric generic Mordell--Weil rank over Qbar(u) <= 13")
    print(args.output)


if __name__ == "__main__":
    main()
