#!/usr/bin/env python3
"""Exhaust the small-field q=80 CM boundary for polynomial sections.

The q=80 ambient has fibers ``I1*``, ``I4``, and ``IV*`` at ``0``, ``1``,
and infinity.  Its four normalized parameters are ``(d,p,q,e)``.  This
script enumerates that ambient over a small prime field, retains the exact
one-``A1`` boundary (the residual discriminant quintic has gcd degree one
with its derivative), and exhausts every polynomial section

    x = x0 + ... + x4*T^4,    y = y0 + ... + y6*T^6.

It is a bounded finite-field discovery calculation.  A hit is not by itself
an identification with the characteristic-zero discriminant-24 CM anchor;
component labels, pair intersections, lifting, and rational reconstruction
remain separate gates.
"""

from __future__ import annotations

import argparse
from itertools import product
from math import comb
import time


PROTOCOL = "Q80CMSECT"


def trim(poly, modulus):
    values = [value % modulus for value in poly]
    while values and values[-1] == 0:
        values.pop()
    return values


def add(left, right, modulus):
    out = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        out[index] += value
    for index, value in enumerate(right):
        out[index] += value
    return trim(out, modulus)


def scale(poly, scalar, modulus):
    return trim([scalar * value for value in poly], modulus)


def multiply(left, right, modulus, cap=None):
    if not left or not right:
        return []
    size = len(left) + len(right) - 1
    if cap is not None:
        size = min(size, cap)
    out = [0] * size
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            index = left_index + right_index
            if index >= size:
                break
            out[index] += left_value * right_value
    return trim(out, modulus)


def power(poly, exponent, modulus, cap=None):
    out = [1]
    base = poly
    while exponent:
        if exponent & 1:
            out = multiply(out, base, modulus, cap)
        exponent >>= 1
        if exponent:
            base = multiply(base, base, modulus, cap)
    return out


def derivative(poly, modulus):
    return trim([index * poly[index] for index in range(1, len(poly))], modulus)


def evaluate(poly, value, modulus):
    result = 0
    for coefficient in reversed(poly):
        result = (result * value + coefficient) % modulus
    return result


def divmod_poly(numerator, denominator, modulus):
    numerator = trim(numerator, modulus)
    denominator = trim(denominator, modulus)
    if not denominator:
        raise ZeroDivisionError("polynomial division by zero")
    quotient = [0] * max(0, len(numerator) - len(denominator) + 1)
    inverse_lead = pow(denominator[-1], -1, modulus)
    while numerator and len(numerator) >= len(denominator):
        shift = len(numerator) - len(denominator)
        coefficient = numerator[-1] * inverse_lead % modulus
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            numerator[index + shift] -= coefficient * value
        numerator = trim(numerator, modulus)
    return trim(quotient, modulus), numerator


def gcd_poly(left, right, modulus):
    left = trim(left, modulus)
    right = trim(right, modulus)
    while right:
        _, remainder = divmod_poly(left, right, modulus)
        left, right = right, remainder
    if not left:
        return []
    return scale(left, pow(left[-1], -1, modulus), modulus)


def resolved_pair_intersection(left, right, modulus):
    """Compute P.Q=(P-Q).O for two polynomial sections on the K3 model."""
    x_left, y_left = map(list, left)
    x_right, y_right = map(list, right)
    difference = add(x_left, scale(x_right, -1, modulus), modulus)
    y_sum = add(y_left, y_right, modulus)
    numerator = add(
        power(y_sum, 2, modulus),
        scale(
            multiply(
                power(difference, 2, modulus),
                add(x_left, x_right, modulus),
                modulus,
            ),
            -1,
            modulus,
        ),
        modulus,
    )
    first = gcd_poly(difference, numerator, modulus)
    numerator_once, remainder = divmod_poly(numerator, first, modulus)
    assert not remainder
    cancellation = gcd_poly(difference, numerator_once, modulus)
    reduced_denominator, remainder = divmod_poly(difference, cancellation, modulus)
    assert not remainder
    cancellation_square = power(cancellation, 2, modulus)
    reduced_numerator, remainder = divmod_poly(
        numerator, cancellation_square, modulus
    )
    assert not remainder
    finite = max(0, len(reduced_denominator) - 1)
    excess = (len(reduced_numerator) - 1) - 2 * finite - 4
    if excess > 0 and excess % 2:
        raise AssertionError("odd infinity excess in pair-intersection gate")
    return finite + (excess // 2 if excess > 0 else 0)


def marked_cm24_short_triples(sections, d1, d2, d3, modulus):
    """Apply inverse-label and resolved pair gates to coarse CM24 triples."""
    hits = []
    for indices in product(d1, d2, d3):
        if len(set(indices)) != 3:
            continue
        representatives = [sections[index - 1] for index in indices]
        for signs in product((1, -1), repeat=3):
            oriented = [
                (x, tuple(value * sign % modulus for value in y))
                for (x, y), sign in zip(representatives, signs)
            ]
            y2 = [section[1][2] for section in oriented]
            y4 = [section[1][4] for section in oriented]
            # D1,D2 use the same D5 spinor class.  D2,D3 use the E6 class
            # inverse to D1.  These leading exceptional coordinates record
            # the orientations exactly.
            if not (y2[0] and y2[0] == y2[1]):
                continue
            if not (y4[0] and y4[1] == -y4[0] % modulus and y4[2] == y4[1]):
                continue
            pairs = tuple(
                resolved_pair_intersection(oriented[left], oriented[right], modulus)
                for left, right in ((0, 1), (0, 2), (1, 2))
            )
            if pairs == (0, 0, 0):
                hits.append((indices, signs, pairs))
    return tuple(hits)


def invert_matrix(matrix, modulus):
    size = len(matrix)
    augmented = [
        [value % modulus for value in row]
        + [1 if row_index == column else 0 for column in range(size)]
        for row_index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(
            row for row in range(column, size) if augmented[row][column] % modulus
        )
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        inverse = pow(augmented[column][column], -1, modulus)
        augmented[column] = [value * inverse % modulus for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column] % modulus
            if factor:
                augmented[row] = [
                    (left - factor * right) % modulus
                    for left, right in zip(augmented[row], augmented[column])
                ]
    return [row[size:] for row in augmented]


def matrix_vector(matrix, vector, modulus):
    return [
        sum(left * right for left, right in zip(row, vector)) % modulus
        for row in matrix
    ]


def square_roots_polynomial(poly, modulus):
    """Return the two sign roots of an exact square in F_p[T], or empty."""
    poly = trim(poly, modulus)
    if not poly:
        return ([0],)
    valuation = next(index for index, value in enumerate(poly) if value % modulus)
    if valuation % 2:
        return ()
    shift = valuation // 2
    lead = poly[valuation] % modulus
    lead_roots = [value for value in range(modulus) if value * value % modulus == lead]
    roots = []
    target_degree = (len(poly) - 1) // 2
    for lead_root in lead_roots:
        root = [0] * (target_degree + 1)
        root[shift] = lead_root
        inverse_twice_lead = pow(2 * lead_root, -1, modulus)
        for offset in range(1, target_degree - shift + 1):
            degree = valuation + offset
            known = 0
            for left in range(shift + 1, shift + offset):
                right = degree - left
                if 0 <= right < len(root):
                    known += root[left] * root[right]
            root[shift + offset] = (
                (poly[degree] if degree < len(poly) else 0) - known
            ) * inverse_twice_lead % modulus
        if trim(multiply(root, root, modulus), modulus) == poly:
            roots.append(trim(root, modulus))
    return tuple(roots)


def q80_surface(d, parameter_p, parameter_q, e, modulus, interpolation_inverse):
    d2 = d * d % modulus
    r = (-3 * d2 + 3 - parameter_p - parameter_q) % modulus
    coefficient_a = [0, 0, -3, parameter_p, parameter_q, r]

    # Taylor jets A(1+s), followed by the cubic binomial truncation of
    # 2*d^3*(A/(-3*d^2))^(3/2).
    jets_a = [
        sum(
            coefficient_a[degree] * comb(degree, jet)
            for degree in range(jet, len(coefficient_a))
        ) % modulus
        for jet in range(4)
    ]
    if jets_a[0] != -3 * d2 % modulus:
        raise AssertionError("the I4 normalization A(1)=-3*d^2 failed")
    denominator = (-3 * d2) % modulus
    inverse_denominator = pow(denominator, -1, modulus)
    u = [0] + [jets_a[index] * inverse_denominator % modulus for index in range(1, 4)]
    branch = [1]
    branch = add(branch, scale(u, 3 * pow(2, -1, modulus), modulus), modulus)
    branch = add(
        branch,
        scale(multiply(u, u, modulus, 4), 3 * pow(8, -1, modulus), modulus),
        modulus,
    )
    branch = add(
        branch,
        scale(power(u, 3, modulus, 4), -pow(16, -1, modulus), modulus),
        modulus,
    )
    branch += [0] * (4 - len(branch))
    branch = [2 * pow(d, 3, modulus) * value % modulus for value in branch[:4]]

    fixed = [
        (2 * comb(3, jet) + e * comb(8, jet)) % modulus for jet in range(4)
    ]
    right = [(left - known) % modulus for left, known in zip(branch, fixed)]
    b1, b2, b3, b4 = matrix_vector(interpolation_inverse, right, modulus)
    coefficient_b = [0, 0, 0, 2, b1, b2, b3, b4, e]
    return coefficient_a, coefficient_b


def section_x_classes(coefficient_a, coefficient_b, modulus, max_x_degree=4):
    sections = []
    for short_x_tuple in product(range(modulus), repeat=max_x_degree + 1):
        x_tuple = tuple(short_x_tuple) + (0,) * (4 - max_x_degree)
        x_poly = trim(x_tuple, modulus)
        right = add(
            add(power(x_poly, 3, modulus), multiply(coefficient_a, x_poly, modulus), modulus),
            coefficient_b,
            modulus,
        )
        roots = square_roots_polynomial(right, modulus)
        if roots:
            sections.append((tuple(x_tuple), tuple(roots[0] + [0] * (7 - len(roots[0])))))
    return sections


def coarse_cm24_short_profiles(
    sections, d, double_root, coefficient_a, coefficient_b, modulus
):
    """Classify the three polynomial CM24 basis profiles coarsely.

    In (A1,A3,D5,E6) order they are

        D1=(1,1,1,1), D2=(0,0,1,2), D3=(1,1,0,2).

    This records identity versus nonidentity only; inverse additive labels,
    tangent branches, and resolved pair intersections are later gates.
    """
    i2_node = repeated_cubic_root(
        coefficient_a, coefficient_b, double_root, modulus
    )
    buckets = {"D1": [], "D2": [], "D3": []}
    targets = {
        "D1": (True, True, True, True),
        "D2": (False, False, True, True),
        "D3": (True, True, False, True),
    }
    for index, (x_coefficients, y_coefficients) in enumerate(sections, 1):
        x_poly = trim(x_coefficients, modulus)
        y_poly = trim(y_coefficients, modulus)
        pattern = (
            evaluate(x_poly, double_root, modulus) == i2_node
            and evaluate(y_poly, double_root, modulus) == 0,
            evaluate(x_poly, 1, modulus) == d % modulus
            and evaluate(y_poly, 1, modulus) == 0,
            (x_coefficients[0], y_coefficients[0]) == (0, 0),
            (x_coefficients[4], y_coefficients[6]) == (0, 0),
        )
        for label, target in targets.items():
            if pattern == target:
                buckets[label].append(index)
    return tuple(buckets[label] for label in ("D1", "D2", "D3"))


def repeated_cubic_root(coefficient_a, coefficient_b, fiber, modulus):
    a_value = evaluate(coefficient_a, fiber, modulus)
    b_value = evaluate(coefficient_b, fiber, modulus)
    roots = [
        value
        for value in range(modulus)
        if (value**3 + a_value * value + b_value) % modulus == 0
        and (3 * value * value + a_value) % modulus == 0
    ]
    if len(roots) != 1:
        raise AssertionError("the residual I2 cubic has no unique node")
    return roots[0]


def coarse_generic_profiles(
    sections, d, double_root, coefficient_a, coefficient_b, modulus
):
    """Filter the two polynomial directions that continue generically.

    The transported generic basis specializes with profiles

        P1 = (1,1,1,1),  P2 = (0,0,0,0)

    at ``(A1,A3,D5,E6)``.  This routine checks only identity versus
    nonidentity.  It deliberately leaves the deeper A3/D5/E6 labels and all
    pair intersections to the exact blowup audit.
    """
    i2_node = repeated_cubic_root(
        coefficient_a, coefficient_b, double_root, modulus
    )
    p1 = []
    p2 = []
    for index, (x_coefficients, y_coefficients) in enumerate(sections, 1):
        x_poly = trim(x_coefficients, modulus)
        y_poly = trim(y_coefficients, modulus)
        hits_d5_cusp = (x_coefficients[0], y_coefficients[0]) == (0, 0)
        hits_e6_cusp = (x_coefficients[4], y_coefficients[6]) == (0, 0)
        hits_i4_node = (
            evaluate(x_poly, 1, modulus) == d % modulus
            and evaluate(y_poly, 1, modulus) == 0
        )
        hits_i2_node = (
            evaluate(x_poly, double_root, modulus) == i2_node
            and evaluate(y_poly, double_root, modulus) == 0
        )
        pattern = (hits_i2_node, hits_i4_node, hits_d5_cusp, hits_e6_cusp)
        if pattern == (True, True, True, True):
            p1.append(index)
        if pattern == (False, False, False, False):
            p2.append(index)
    return tuple(p1), tuple(p2)


def one_pole_sections(
    coefficient_a,
    coefficient_b,
    modulus,
    double_root,
    d,
):
    """Exhaust the one-pole x-numerator chart for the transported P3."""
    i2_node = repeated_cubic_root(
        coefficient_a, coefficient_b, double_root, modulus
    )
    records = []
    poles = list(range(modulus)) + [None]
    for pole in poles:
        h = [1] if pole is None else [(-pole) % modulus, 1]
        h4 = power(h, 4, modulus)
        h6 = power(h, 6, modulus)
        for numerator_tuple in product(range(modulus), repeat=7):
            numerator = trim(numerator_tuple, modulus)
            if pole is not None and evaluate(numerator, pole, modulus) == 0:
                continue
            right = add(
                add(
                    power(numerator, 3, modulus),
                    multiply(multiply(coefficient_a, numerator, modulus), h4, modulus),
                    modulus,
                ),
                multiply(coefficient_b, h6, modulus),
                modulus,
            )
            roots = square_roots_polynomial(right, modulus)
            if not roots:
                continue
            denominator_at = lambda value: evaluate(h, value, modulus)
            if pole in (1, double_root):
                continue

            def value_at(poly, value, exponent):
                denominator = pow(denominator_at(value), exponent, modulus)
                return evaluate(poly, value, modulus) * pow(denominator, -1, modulus) % modulus

            root = roots[0]
            x_i2 = value_at(numerator, double_root, 2)
            y_i2 = value_at(root, double_root, 3)
            x_i4 = value_at(numerator, 1, 2)
            y_i4 = value_at(root, 1, 3)
            hits_i2 = x_i2 == i2_node and y_i2 == 0
            hits_i4 = x_i4 == d % modulus and y_i4 == 0

            if pole == 0:
                d5_identity = True
            else:
                x_zero = value_at(numerator, 0, 2)
                y_zero = value_at(root, 0, 3)
                d5_identity = (x_zero, y_zero) != (0, 0)

            if pole is None:
                e6_identity = True
            else:
                numerator_full = numerator + [0] * (7 - len(numerator))
                root_full = root + [0] * (10 - len(root))
                e6_identity = (numerator_full[6], root_full[9]) != (0, 0)

            if (hits_i2, hits_i4, d5_identity, e6_identity) != (
                True,
                True,
                True,
                True,
            ):
                continue
            records.append(
                (
                    "infinity" if pole is None else str(pole),
                    tuple(numerator_tuple),
                    tuple(root + [0] * (10 - len(root))),
                )
            )
    return tuple(records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=5)
    parser.add_argument("--min-sections", type=int, default=2)
    parser.add_argument("--progress-every", type=int, default=50)
    parser.add_argument("--scan-one-pole", action="store_true")
    parser.add_argument("--max-x-degree", type=int, choices=(2, 3, 4), default=4)
    parser.add_argument(
        "--profile-mode", choices=("generic", "cm24-short"), default="generic"
    )
    args = parser.parse_args()

    prime = args.prime
    if prime in (2, 3) or any(prime % divisor == 0 for divisor in range(2, int(prime**0.5) + 1)):
        raise SystemExit("--prime must be an odd prime different from 3")

    interpolation = [
        [comb(degree, jet) % prime for degree in range(4, 8)]
        for jet in range(4)
    ]
    interpolation_inverse = invert_matrix(interpolation, prime)
    fixed_discriminant = multiply(
        [0] * 7 + [1], power([-1, 1], 4, prime), prime
    )

    t0 = time.monotonic()
    tested = 0
    cm_boundary = 0
    hits = 0
    print(
        f"{PROTOCOL}|stage=start|prime={prime}|ambient={(prime - 1) * prime**3}"
        f"|section_x_space={prime**(args.max_x_degree + 1)}|"
        f"max_x_degree={args.max_x_degree}|profile_mode={args.profile_mode}|"
        f"min_sections={args.min_sections}",
        flush=True,
    )

    for d in range(1, prime):
        for parameter_p, parameter_q, e in product(range(prime), repeat=3):
            tested += 1
            coefficient_a, coefficient_b = q80_surface(
                d, parameter_p, parameter_q, e, prime, interpolation_inverse
            )
            discriminant = add(
                scale(power(coefficient_a, 3, prime), 4, prime),
                scale(power(coefficient_b, 2, prime), 27, prime),
                prime,
            )
            residual, remainder = divmod_poly(discriminant, fixed_discriminant, prime)
            if remainder or len(residual) != 6:
                continue
            gcd = gcd_poly(residual, derivative(residual, prime), prime)
            if len(gcd) != 2:
                continue
            double_root = -gcd[0] * pow(gcd[1], -1, prime) % prime
            # The extra A1 must be a separate I2, not an enhancement of the
            # already pinned fibers at zero or one.
            if double_root in (0, 1):
                continue
            cm_boundary += 1
            sections = section_x_classes(
                coefficient_a, coefficient_b, prime, args.max_x_degree
            )
            if len(sections) < args.min_sections:
                continue
            if args.profile_mode == "generic":
                coarse_p1, coarse_p2 = coarse_generic_profiles(
                    sections, d, double_root, coefficient_a, coefficient_b, prime
                )
                profile_text = (
                    f"coarse_P1={','.join(map(str, coarse_p1))}"
                    f"|coarse_P2={','.join(map(str, coarse_p2))}"
                )
                if not coarse_p1 or not coarse_p2:
                    continue
            else:
                short_d1, short_d2, short_d3 = coarse_cm24_short_profiles(
                    sections, d, double_root, coefficient_a, coefficient_b, prime
                )
                profile_text = (
                    f"coarse_D1={','.join(map(str, short_d1))}"
                    f"|coarse_D2={','.join(map(str, short_d2))}"
                    f"|coarse_D3={','.join(map(str, short_d3))}"
                )
                if not short_d1 or not short_d2 or not short_d3:
                    continue
                marked_short = marked_cm24_short_triples(
                    sections, short_d1, short_d2, short_d3, prime
                )
                profile_text += f"|marked_triples={len(marked_short)}"
                if not marked_short:
                    continue
            hits += 1
            print(
                f"{PROTOCOL}|stage=hit|index={hits}|d={d}|p={parameter_p}"
                f"|q={parameter_q}|e={e}|double_root={double_root}"
                f"|section_x_classes={len(sections)}|{profile_text}",
                flush=True,
            )
            for index, (x_poly, y_poly) in enumerate(sections, 1):
                print(
                    f"{PROTOCOL}|section={index}|x={','.join(map(str, x_poly))}"
                    f"|y={','.join(map(str, y_poly))}",
                    flush=True,
                )
            print(
                f"{PROTOCOL}|surface_A={','.join(map(str, coefficient_a))}"
                f"|surface_B={','.join(map(str, coefficient_b))}",
                flush=True,
            )
            if args.profile_mode == "cm24-short":
                for triple_index, (indices, signs, pairs) in enumerate(marked_short, 1):
                    print(
                        f"{PROTOCOL}|marked_triple={triple_index}|indices={indices}|"
                        f"signs={signs}|pairs={pairs}",
                        flush=True,
                    )
            if args.scan_one_pole:
                one_pole = one_pole_sections(
                    coefficient_a,
                    coefficient_b,
                    prime,
                    double_root,
                    d,
                )
                print(
                    f"{PROTOCOL}|stage=one_pole|count={len(one_pole)}",
                    flush=True,
                )
                for index, (pole, numerator, root) in enumerate(one_pole, 1):
                    print(
                        f"{PROTOCOL}|one_pole={index}|pole={pole}"
                        f"|N={','.join(map(str, numerator))}"
                        f"|M={','.join(map(str, root))}",
                        flush=True,
                    )

            if tested % args.progress_every == 0:
                print(
                    f"{PROTOCOL}|stage=progress|tested={tested}|cm_boundary={cm_boundary}"
                    f"|hits={hits}|seconds={time.monotonic() - t0:.3f}",
                    flush=True,
                )

    print(
        f"{PROTOCOL}|stage=done|prime={prime}|tested={tested}|cm_boundary={cm_boundary}"
        f"|hits={hits}|seconds={time.monotonic() - t0:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
