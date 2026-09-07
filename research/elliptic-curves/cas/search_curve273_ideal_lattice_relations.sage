#!/usr/bin/env sage
"""Search full cubic-ideal lattices for a curve-273 Selmer relation.

The older CRT probes restrict principal elements to ``m-theta``.  Once the
remaining large-prime modulus exceeds the real-root scale, that one-dimensional
slice becomes inefficient.  This script instead forms the product of the
declared degree-one prime ideals in the cubic 2-division field and enumerates
short elements in its full three-dimensional Minkowski lattice.

This is a bounded relation search, not a complete 2-Selmer computation.
Every reported closure is nevertheless checked by exact principal-ideal
factorization.
"""

from __future__ import annotations

import argparse
import json
from heapq import heappush, heapreplace
from itertools import combinations, product
from math import gcd, lcm
from pathlib import Path
import sys
import time

from sage.all import (
    ComplexField,
    GF,
    Matrix,
    NumberField,
    PolynomialRing,
    QQ,
    RealField,
    ZZ,
    prime_range,
    prod,
    vector,
)


sys.path.insert(0, str(Path(__file__).resolve().parent))

from icarm_curve273 import short_coefficients


PROTOCOL = "R30IDEAL"

S_RATIONAL = {
    2,
    3,
    5,
    7,
    13,
    31,
    41,
    47,
    53,
    67,
    379,
    4349,
    25721454817,
    97018222656318846556561979214040553412450110580812087282349817173780902099339117104673990259247421230916714670243202937,
}


def sage_q(value):
    numerator = value.numerator
    denominator = value.denominator
    if callable(numerator):
        numerator = numerator()
    if callable(denominator):
        denominator = denominator()
    return QQ(ZZ(numerator)) / QQ(ZZ(denominator))


def parse_ideal(text):
    q, residue = text.split(":", 1)
    return ZZ(q), ZZ(residue)


def canonical_sign(values):
    values = tuple(ZZ(value) for value in values)
    for value in values:
        if value:
            return values if value > 0 else tuple(-entry for entry in values)
    return values


def canonical_field_element_key(basis_coefficients, coordinates):
    """Canonical power-basis key for an element, modulo multiplication by -1.

    A single algebraic integer can occur in several factor-base-twisted ideal
    bases.  Keying the shortlist by ``(twist, coordinates)`` therefore wastes
    most of the factorization budget on exact duplicates.  Power-basis
    coordinates give a basis-independent key without factoring the norm.
    """
    coefficients = []
    for degree in range(3):
        coefficient = sum(
            (
                QQ(coordinates[index])
                * basis_coefficients[index][degree]
                if degree < len(basis_coefficients[index])
                else QQ(0)
            )
            for index in range(3)
        )
        coefficients.append(coefficient)
    for coefficient in coefficients:
        if coefficient:
            if coefficient < 0:
                coefficients = [-entry for entry in coefficients]
            break
    return tuple(coefficients)


def integer_norm_form(basis):
    """Return ``(denominator, monomials)`` for the ternary norm form."""
    ring = PolynomialRing(QQ, names=("u", "v", "w"))
    variables = ring.gens()
    multiplication = sum(
        (
            variables[index] * Matrix(ring, element.matrix())
            for index, element in enumerate(basis)
        ),
        Matrix(ring, 3, 3),
    )
    norm_form = multiplication.det()
    denominator = lcm(
        *(ZZ(coefficient.denominator()) for coefficient in norm_form.coefficients())
    )
    monomials = tuple(
        (tuple(exponents), ZZ(denominator * coefficient))
        for exponents, coefficient in norm_form.dict().items()
    )
    return ZZ(denominator), monomials, norm_form


def evaluate_norm(denominator, monomials, values):
    powers = tuple(
        (ZZ(1), value, value * value, value * value * value)
        for value in values
    )
    numerator = ZZ(0)
    for exponents, coefficient in monomials:
        term = coefficient
        for index, exponent in enumerate(exponents):
            term *= powers[index][exponent]
        numerator += term
    if numerator % denominator:
        raise ArithmeticError("integral ideal element acquired a fractional norm")
    return numerator // denominator


def minkowski_lll_transform(field, basis, precision, shape_shift):
    """LLL transformation for a determinant-one diagonal Minkowski twist."""
    real_field = RealField(precision)
    complex_field = ComplexField(precision)
    places = field.places(prec=precision)
    real_places = [place for place in places if place.codomain() is real_field]
    complex_places = [place for place in places if place.codomain() is complex_field]
    # Sage may construct equal-but-not-identical codomains, so fall back to
    # detecting the unique real embedding by the image of the generator.
    if len(real_places) != 1 or len(complex_places) != 1:
        real_places = []
        complex_places = []
        for place in places:
            image = place(field.gen())
            if image.parent().is_real_field():
                real_places.append(place)
            else:
                complex_places.append(place)
    if len(real_places) != 1 or len(complex_places) != 1:
        raise RuntimeError("expected signature (1,1)")

    real_place = real_places[0]
    complex_place = complex_places[0]
    two = real_field(2)
    real_weight = two**shape_shift
    complex_weight = two ** (-QQ(shape_shift) / 2)
    sqrt_two = two.sqrt()

    rows = []
    for element in basis:
        real_image = real_field(real_place(element))
        complex_image = complex_place(element)
        rows.append(
            (
                real_weight * real_image,
                complex_weight * sqrt_two * real_field(complex_image.real()),
                complex_weight * sqrt_two * real_field(complex_image.imag()),
            )
        )

    largest = max(abs(entry) for row in rows for entry in row)
    scale = two ** (precision - 48) / largest
    rounded = Matrix(
        ZZ,
        [[ZZ((scale * entry).round()) for entry in row] for row in rows],
    )
    if rounded.det() == 0:
        raise ArithmeticError("rounded Minkowski basis lost rank; increase precision")
    _, transform = rounded.LLL(transformation=True)
    if abs(transform.det()) != 1:
        raise ArithmeticError("LLL transformation is not unimodular")
    return transform


def strip_support(value, smooth_support_product, target_primes):
    value = abs(ZZ(value))
    for prime in target_primes:
        if value % prime:
            return None
        value //= prime
    return strip_primorial(value, smooth_support_product)


def strip_primorial(value, primorial):
    """Remove all powers of primes represented once in ``primorial``."""
    value = int(value)
    while value != 1:
        common = gcd(value, primorial)
        if common == 1:
            break
        value //= common
    return ZZ(value)


def degree_one_labels(defining_polynomial, alpha, prime):
    """Return roots ``r`` for degree-one ideals ``(p,theta-r)`` containing alpha."""
    residue_field = GF(prime)
    alpha_coefficients = tuple(QQ(value) for value in alpha.polynomial().list())
    if any(
        coefficient.denominator() % prime == 0
        for coefficient in alpha_coefficients
    ):
        raise ArithmeticError(
            f"power-basis denominator is not invertible modulo {prime}"
        )
    labels = []
    for root, _ in defining_polynomial.change_ring(residue_field).roots():
        value = sum(
            residue_field(coefficient) * root**index
            for index, coefficient in enumerate(alpha_coefficients)
        )
        if value == 0:
            labels.append(ZZ(root))
    return tuple(labels)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--target",
        action="append",
        help="degree-one target ideal q:theta_residue (repeatable)",
    )
    parser.add_argument(
        "--target-plan",
        type=Path,
        help=(
            "JSON plan from analyze_curve273_relation_pool.py; select one "
            "ranked ideal with --target-plan-rank"
        ),
    )
    parser.add_argument(
        "--target-plan-rank",
        type=int,
        help="one-based target rank to select from --target-plan",
    )
    parser.add_argument("--factor-base-bound", type=int, default=1_000_000)
    parser.add_argument("--radius", type=int, default=18)
    parser.add_argument(
        "--twist-factor-base-bound",
        type=int,
        default=0,
        help=(
            "also search the target ideal multiplied by every degree-one "
            "prime ideal above p <= this bound (zero disables this mode)"
        ),
    )
    parser.add_argument(
        "--twist-radius",
        type=int,
        default=4,
        help="coefficient radius for each factor-base-twisted ideal lattice",
    )
    parser.add_argument(
        "--twist-depth",
        type=int,
        choices=(1, 2),
        default=1,
        help=(
            "search individual factor-base prime-ideal twists, or those "
            "twists together with every unordered pair"
        ),
    )
    parser.add_argument("--preselect", type=int, default=2000)
    parser.add_argument(
        "--presieve-bound",
        type=int,
        default=1000,
        help="rank candidates after stripping all prime powers below this bound",
    )
    parser.add_argument("--factor-top", type=int, default=400)
    parser.add_argument("--factor-max-bits", type=int, default=170)
    parser.add_argument("--print-best", type=int, default=20)
    parser.add_argument("--precision", type=int, default=320)
    parser.add_argument(
        "--shape-shifts",
        default="-48,-24,0,24,48",
        help="comma-separated real-embedding base-2 twists",
    )
    args = parser.parse_args()

    if args.target_plan:
        if args.target:
            parser.error("--target and --target-plan are mutually exclusive")
        if args.target_plan_rank is None or args.target_plan_rank < 1:
            parser.error("--target-plan requires a positive --target-plan-rank")
        try:
            plan = json.loads(args.target_plan.read_text())
        except (OSError, json.JSONDecodeError) as error:
            parser.error(f"cannot read --target-plan: {error}")
        if plan.get("schema") != "elliptic-curves.bnf-free-large-prime-target-plan.v1":
            parser.error("--target-plan has an unsupported schema")
        selected = [
            item
            for item in plan.get("targets", [])
            if item.get("rank") == args.target_plan_rank
        ]
        if len(selected) != 1:
            parser.error(
                f"--target-plan contains no unique rank {args.target_plan_rank}"
            )
        target = selected[0]
        raw_targets = [
            f"{target['rational_prime']}:{target['residue']}"
        ]
    else:
        if args.target_plan_rank is not None:
            parser.error("--target-plan-rank requires --target-plan")
        raw_targets = args.target or []

    targets = tuple(parse_ideal(text) for text in raw_targets)
    if len(targets) < 1 or len(set(targets)) != len(targets):
        raise SystemExit("target prime-ideal labels must be nonempty and distinct")

    coefficients = short_coefficients()
    A = ZZ(sage_q(coefficients[3]))
    B = ZZ(sage_q(coefficients[4]))
    polynomial_ring = PolynomialRing(QQ, "x")
    x = polynomial_ring.gen()
    defining_polynomial = x**3 + A * x + B
    field = NumberField(defining_polynomial, "theta")
    theta = field.gen()

    target_ideals = []
    for q, residue in targets:
        if not q.is_prime(proof=True):
            raise SystemExit(f"target q={q} is not prime")
        if defining_polynomial(residue) % q:
            raise SystemExit(f"target {q}:{residue} is not a root")
        prime_ideal = field.ideal(q, theta - residue)
        if prime_ideal.norm() != q or not prime_ideal.is_prime():
            raise SystemExit(f"target {q}:{residue} is not a degree-one prime ideal")
        target_ideals.append(prime_ideal)

    target_norm = prod(q for q, _ in targets)
    ideal = prod(target_ideals, field.ideal(1))
    if ideal.norm() != target_norm:
        raise ArithmeticError("target ideal product has the wrong norm")
    basis = tuple(ideal.basis())
    if len(basis) != 3:
        raise ArithmeticError("expected a ternary ideal basis")

    denominator, monomials, norm_form = integer_norm_form(basis)
    shifts = tuple(
        ZZ(value.strip())
        for value in args.shape_shifts.split(",")
        if value.strip()
    )
    small_primes = tuple(
        ZZ(prime) for prime in prime_range(2, args.factor_base_bound + 1)
    )
    smooth_support_product = prod(
        small_primes
        + tuple(
            ZZ(prime)
            for prime in S_RATIONAL
            if prime > args.factor_base_bound
        ),
        ZZ(1),
    )
    presieve_primes = tuple(
        int(prime) for prime in prime_range(2, args.presieve_bound + 1)
    )
    presieve_primorial = prod(presieve_primes, 1)

    print(
        f"{PROTOCOL}|stage=input|targets="
        + ",".join(f"{q}:{residue}" for q, residue in targets)
        + f"|target_norm_bits={target_norm.nbits()}"
        + f"|field_disc_bits={abs(field.discriminant()).nbits()}"
        + f"|ideal_norm={ideal.norm()}|radius={args.radius}|shifts={','.join(map(str, shifts))}"
        + f"|twist_factor_base_bound={args.twist_factor_base_bound}"
        + f"|twist_radius={args.twist_radius}|twist_depth={args.twist_depth}",
        flush=True,
    )
    print(f"{PROTOCOL}|stage=norm_form|form={norm_form}", flush=True)

    ideal_entries = [("none", ideal, basis, denominator, monomials, args.radius)]
    if args.twist_factor_base_bound:
        if args.twist_factor_base_bound < 2 or args.twist_radius < 1:
            raise SystemExit("twist bound and radius must be positive")
        factor_prime_ideals = []
        for rational_prime in prime_range(2, args.twist_factor_base_bound + 1):
            for root, _ in defining_polynomial.change_ring(GF(rational_prime)).roots():
                root = ZZ(root)
                prime_ideal = field.ideal(rational_prime, theta - root)
                if prime_ideal.norm() != rational_prime or not prime_ideal.is_prime():
                    continue
                factor_prime_ideals.append((f"{rational_prime}:{root}", prime_ideal))

        twist_products = list(factor_prime_ideals)
        if args.twist_depth == 2:
            twist_products.extend(
                (f"{left_label}+{right_label}", left_ideal * right_ideal)
                for (left_label, left_ideal), (right_label, right_ideal) in combinations(
                    factor_prime_ideals, 2
                )
            )

        for twist_label, twist_ideal in twist_products:
            twisted_ideal = ideal * twist_ideal
            twisted_basis = tuple(twisted_ideal.basis())
            twisted_denominator, twisted_monomials, _ = integer_norm_form(
                twisted_basis
            )
            ideal_entries.append(
                (
                    twist_label,
                    twisted_ideal,
                    twisted_basis,
                    twisted_denominator,
                    twisted_monomials,
                    args.twist_radius,
                )
            )
        print(
            f"{PROTOCOL}|stage=twists|count={len(ideal_entries)-1}"
            + f"|bound={args.twist_factor_base_bound}|radius={args.twist_radius}"
            + f"|depth={args.twist_depth}",
            flush=True,
        )

    # A bounded max-heap keyed by the exact residual after one declared copy
    # of each target rational prime and removal of tiny prime powers.  This
    # avoids selecting only the smallest raw norms: smoothness, not magnitude,
    # is the actual objective.  Original ideal-basis coordinates dedupe vectors
    # obtained from different Minkowski shapes.
    heap = []
    retained = {}
    enumerated = 0
    started = time.monotonic()

    basis_by_tag = {}
    for tag, entry_ideal, entry_basis, entry_denominator, entry_monomials, radius in ideal_entries:
        basis_by_tag[tag] = entry_basis
        entry_basis_coefficients = tuple(
            tuple(QQ(value) for value in element.polynomial().list())
            for element in entry_basis
        )
        for shift in shifts:
            transform = minkowski_lll_transform(
                field, entry_basis, args.precision, shift
            )
            if tag == "none" or args.twist_factor_base_bound <= 50:
                print(
                    f"{PROTOCOL}|stage=lll|twist={tag}|shift={shift}|transform="
                    + ";".join(",".join(map(str, row)) for row in transform.rows()),
                    flush=True,
                )
            for local in product(range(-radius, radius + 1), repeat=3):
                if local == (0, 0, 0):
                    continue
                original = canonical_sign(vector(ZZ, local) * transform)
                coordinate_gcd = gcd(
                    gcd(abs(int(original[0])), abs(int(original[1]))),
                    abs(int(original[2])),
                )
                if coordinate_gcd != 1:
                    continue
                identity = canonical_field_element_key(
                    entry_basis_coefficients, original
                )
                if identity in retained:
                    continue
                norm = abs(
                    evaluate_norm(entry_denominator, entry_monomials, original)
                )
                if norm == 0 or norm % target_norm:
                    raise ArithmeticError("enumerated element escaped the target ideal")
                quotient = norm // target_norm
                score = strip_primorial(quotient, presieve_primorial)
                key = (-score, -quotient, tag, original, identity)
                retained[identity] = (score, quotient)
                if len(heap) < args.preselect:
                    heappush(heap, key)
                elif (score, quotient) < (-heap[0][0], -heap[0][1]):
                    *_, evicted_identity = heapreplace(heap, key)
                    retained.pop(evicted_identity, None)
                else:
                    retained.pop(identity, None)
                enumerated += 1
        print(
            f"{PROTOCOL}|stage=enumerate|twist={tag}|enumerated={enumerated}"
            + f"|retained={len(heap)}|seconds={time.monotonic()-started:.3f}",
            flush=True,
        )

    candidates = sorted(
        (-negative_score, -negative_quotient, tag, coordinates)
        for negative_score, negative_quotient, tag, coordinates, _ in heap
    )
    print(
        f"{PROTOCOL}|stage=factor|candidates={min(args.factor_top, len(candidates))}"
        + f"|best_presieved_bits={candidates[0][0].nbits()}"
        + f"|best_quotient_bits={candidates[0][1].nbits()}",
        flush=True,
    )

    factored = 0
    norm_smooth = 0
    closures = 0
    records = []
    target_rational_primes = {q for q, _ in targets}
    for _, quotient, tag, coordinates in candidates[: args.factor_top]:
        basis = basis_by_tag[tag]
        alpha = sum(
            (coordinates[index] * basis[index] for index in range(3)),
            field(0),
        )
        norm = abs(ZZ(alpha.norm()))
        residual = strip_support(
            norm,
            smooth_support_product,
            tuple(q for q, _ in targets),
        )
        if residual is None or residual.nbits() > args.factor_max_bits:
            continue
        raw_factorization = list(residual.factor(proof=False)) if residual != 1 else []
        if prod(prime**exponent for prime, exponent in raw_factorization) != residual:
            continue
        if any(not ZZ(prime).is_prime(proof=True) for prime, _ in raw_factorization):
            continue
        factored += 1
        large = tuple(
            (ZZ(prime), ZZ(exponent))
            for prime, exponent in raw_factorization
            if prime > args.factor_base_bound and prime not in S_RATIONAL
        )
        external_large = tuple(
            item for item in large if item[0] not in target_rational_primes
        )
        if external_large:
            records.append(
                (
                    len(external_large),
                    max(prime.nbits() for prime, _ in external_large),
                    sum(prime.nbits() for prime, _ in external_large),
                    tag,
                    coordinates,
                    norm,
                    external_large,
                )
            )
            continue

        norm_smooth += 1
        principal = field.ideal(alpha)
        target_valuations = tuple(principal.valuation(prime_ideal) for prime_ideal in target_ideals)
        exact_factorization = tuple(principal.factor())
        extra_large = tuple(
            (prime_ideal, exponent)
            for prime_ideal, exponent in exact_factorization
            if exponent % 2
            and prime_ideal.smallest_integer() > args.factor_base_bound
            and prime_ideal.smallest_integer() not in S_RATIONAL
            and prime_ideal not in target_ideals
        )
        closed = all(value % 2 for value in target_valuations) and not extra_large
        closures += int(closed)
        if closed or not args.twist_factor_base_bound:
            print(
                f"{PROTOCOL}|relation|status={'EXACT_LP_CLOSURE' if closed else 'SMOOTH_NORM_ONLY'}"
                + f"|twist={tag}|coordinates={','.join(map(str, coordinates))}|norm={norm}"
                + f"|target_valuations={','.join(map(str, target_valuations))}"
                + f"|extra_large_ideals={extra_large}"
                + f"|ideal_factorization={exact_factorization}",
                flush=True,
            )

    # For recursive descent, one 163-bit LP is usually worse than two LPs of
    # 49 and 74 bits: either smaller ideal can be forced cheaply in the next
    # stage.  Rank nonclosures by the hardest residual prime first, then by
    # total residual size.  Keep support count only as the third discriminator.
    recursive_records = sorted(
        records,
        key=lambda record: (
            record[1],
            record[2],
            record[0],
            tuple(abs(value) for value in record[4]),
        ),
    )

    for large_count, max_bits, total_bits, tag, coordinates, norm, large in recursive_records[
        : args.print_best
    ]:
        basis = basis_by_tag[tag]
        alpha = sum(
            (coordinates[index] * basis[index] for index in range(3)),
            field(0),
        )
        labels = tuple(
            (prime, degree_one_labels(defining_polynomial, alpha, prime), exponent)
            for prime, exponent in large
        )
        print(
            f"{PROTOCOL}|best|ranking=recursive_descent"
            + f"|large_count={large_count}|max_bits={max_bits}"
            + f"|total_bits={total_bits}|twist={tag}"
            + f"|coordinates={','.join(map(str, coordinates))}"
            + f"|alpha={alpha}|norm={norm}|large={large}|labels={labels}",
            flush=True,
        )

    print(
        f"{PROTOCOL}|stage=summary|enumerated={enumerated}|preselected={len(candidates)}"
        + f"|factored={factored}|norm_smooth={norm_smooth}|closures={closures}"
        + f"|seconds={time.monotonic()-started:.3f}",
        flush=True,
    )
    print(
        f"{PROTOCOL}|status={'EXACT_LP_CLOSURE' if closures else 'BOUNDED_NO_CLOSURE'}",
        flush=True,
    )


if __name__ == "__main__":
    main()
