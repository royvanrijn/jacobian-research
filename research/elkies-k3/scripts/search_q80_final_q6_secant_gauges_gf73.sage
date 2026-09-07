#!/usr/bin/env sage
"""Search exact two-generator q6 modules from fifth-child secant values.

The corrected q6 divisor does not contain the constant function: three local
affine corrections share only two global fiber copies.  Hence a gauge of the
form ``1, chord-h(s)`` cannot represent its Riemann--Roch space.  This script
uses the 27 known polynomial sections Q_i of the pair23 fifth and their chord
values

    mu_i=(y(Q_i)+y(S))/(x(Q_i)-x(S)).

For each pair it tests the genuine two-generator gauge

    R=(chord-mu_i)/(mu_j-mu_i).

Every test is exact over GF(73).  Genus-one hits are converted to their
unit-corrected binary-quartic Jacobians and classified by ADE root data.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, FunctionField, PolynomialRing, gcd


ROOT = Path(__file__).resolve().parents[2]
FIFTH = ROOT / "artifacts/generated-results/q80-deforming-fifth-pair23-gf73.json"
SECTIONS = ROOT / "artifacts/generated-results/q80-deforming-fifth-polynomial-sections-gf73.json"
Q6 = ROOT / "artifacts/generated-results/q80-deforming-fifth-q6-horizontal-candidates-gf73.json"
KNOWN_HASHES = {
    FIFTH: "23fc49bce2618a6d3c5f5e18ded34b4ffbee220be83523ae250bf7774a91db14",
    SECTIONS: "c53804dc81b8f6573e4d0818851cef2612135127205d01b5ee4dbd2ae226c48a",
    Q6: "fcd61f89daab0a68785a006e6b10dc3829b1c30c24243b67a4e1b80c7d6e6e09",
}
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--section-index", type=int, default=0)
parser.add_argument("--all-sections", action="store_true")
parser.add_argument("--mw-norm-bound", type=int, default=1)
parser.add_argument("--max-radical-degree", type=int)
parser.add_argument("--pair-start", type=int, default=0)
parser.add_argument("--pair-end", type=int)
parser.add_argument("--write-artifact", action="store_true")
arguments = parser.parse_args()

payloads = {}
for path, expected in KNOWN_HASHES.items():
    content = path.read_bytes()
    assert hashlib.sha256(content).hexdigest() == expected
    payloads[path] = json.loads(content)
fifth = payloads[FIFTH]
polynomial_sections = payloads[SECTIONS]["sections"]
q6_sections = payloads[Q6]["q6_candidates"]
assert len(polynomial_sections) == 27 and len(q6_sections) == 12
assert all(row["I2_labels"] == [0, 0, 0] for row in q6_sections)
assert all(row["I0star_correction"] == "1" for row in q6_sections)

finite = GF(73, impl="modn")
base = FunctionField(finite, "R")
R = base.gen()
old_ring = PolynomialRing(base, "s")
s = old_ring.gen()
old_field = old_ring.fraction_field()
A = old_ring(fifth["A_coefficients_low_to_high"])
B = old_ring(fifth["B_coefficients_low_to_high"])


def polynomial(coefficients):
    return old_ring([base(value) for value in coefficients])


def fraction(numerator_key, denominator_key, row):
    return old_field(polynomial(row[numerator_key]))/old_field(
        polynomial(row[denominator_key])
    )


def kodaira_data(ord_a, ord_b, ord_delta):
    if ord_a == 0 or ord_b == 0:
        n = int(ord_delta)
        return n-1, n*(n-1), n, n, f"I{n}"
    if ord_delta == 2:
        return 0, 0, 1, 2, "II"
    if ord_delta == 3:
        return 1, 2, 2, 3, "III"
    if ord_delta == 4:
        return 2, 6, 3, 4, "IV"
    if ord_delta == 6 and ord_a >= 2 and ord_b >= 3:
        return 4, 24, 4, 6, "I0*"
    if ord_delta >= 7 and ord_a == 2 and ord_b == 3:
        n = int(ord_delta-6)
        rank = n+4
        return rank, 2*rank*(rank-1), 4, n+6, f"I{n}*"
    if ord_delta == 8:
        return 6, 72, 3, 8, "IV*"
    if ord_delta == 9:
        return 7, 126, 2, 9, "III*"
    if ord_delta == 10:
        return 8, 240, 1, 10, "II*"
    raise ArithmeticError((ord_a, ord_b, ord_delta))


def classify_quartic(quartic, twist):
    coefficients = list(quartic.list())+[base(0)]*5
    e, d, c, b, a = coefficients[:5]
    invariant_I = 12*a*e-3*b*d+c**2
    invariant_J = 72*a*c*e+9*b*c*d-27*a*d**2-27*b**2*e-2*c**3
    jacobian_A = twist**2*(-27*invariant_I)
    jacobian_B = twist**3*(-27*invariant_J)
    delta_core = twist**6*(4*invariant_I**3-invariant_J**2)
    if (
        jacobian_A.denominator() != 1
        or jacobian_B.denominator() != 1
        or delta_core.denominator() != 1
    ):
        return {"status": "nonintegral_raw_jacobian"}
    new_A = jacobian_A.numerator()
    new_B = jacobian_B.numerator()
    new_Delta = delta_core.numerator()
    finite_scalings = []
    for factor, _ in gcd(new_A, new_B).factor():
        order = min(new_A.valuation(factor)//4, new_B.valuation(factor)//6)
        if order <= 0:
            continue
        new_A //= factor**(4*order)
        new_B //= factor**(6*order)
        new_Delta //= factor**(12*order)
        finite_scalings.append((str(factor), int(order)))
    if new_A.degree() > 8 or new_B.degree() > 12 or new_Delta.degree() > 24:
        return {
            "status": "nonminimal_at_infinity",
            "degrees": [new_A.degree(), new_B.degree(), new_Delta.degree()],
        }
    root_rank = root_count = 0
    root_determinant = 1
    euler = 0
    signature = []
    for factor, exponent in new_Delta.factor():
        data = kodaira_data(
            int(new_A.valuation(factor)),
            int(new_B.valuation(factor)),
            int(exponent),
        )
        rank, count, determinant, local_euler, kind = data
        degree = int(factor.degree())
        root_rank += degree*rank
        root_count += degree*count
        root_determinant *= determinant**degree
        euler += degree*local_euler
        signature.append((str(factor), degree, int(exponent), kind))
    infinity_orders = (
        8-new_A.degree(), 12-new_B.degree(), 24-new_Delta.degree()
    )
    infinity_kind = "smooth"
    if infinity_orders[2] > 0:
        rank, count, determinant, local_euler, infinity_kind = kodaira_data(
            *infinity_orders
        )
        root_rank += rank
        root_count += count
        root_determinant *= determinant
        euler += local_euler
    if euler != 24:
        return {"status": "wrong_euler", "euler": euler}
    return {
        "status": "classified",
        "degrees": [new_A.degree(), new_B.degree(), new_Delta.degree()],
        "finite_scalings": [list(row) for row in finite_scalings],
        "finite_signature": [list(row) for row in signature],
        "infinity_orders": list(infinity_orders),
        "infinity_fiber": infinity_kind,
        "root_data": [root_rank, root_count, root_determinant],
        "CM24_MW_rank": 18-root_rank,
        "A_coefficients_low_to_high": list(map(int, new_A.list())),
        "B_coefficients_low_to_high": list(map(int, new_B.list())),
        "Delta_coefficients_low_to_high": list(map(int, new_Delta.list())),
    }


selected = q6_sections if arguments.all_sections else [q6_sections[arguments.section_index]]


def cheap_specialized_radical_degree(polynomial_value):
    """Screen fixed square factors after a good specialization of R."""
    specialized_ring = PolynomialRing(finite, "z")
    z = specialized_ring.gen()
    for value in (2, 3, 5, 7, 11):
        coefficients = []
        for coefficient in polynomial_value.list():
            if coefficient.denominator()(finite(value)) == 0:
                break
            coefficients.append(
                coefficient.numerator()(finite(value))
                / coefficient.denominator()(finite(value))
            )
        else:
            specialized = specialized_ring(coefficients)
            radical = specialized//gcd(specialized, specialized.derivative())
            return int(radical.degree())
    return None

# Build a larger exact section pool when requested.  Norm is measured in the
# A2 form a^2-ab+b^2, so heights are norm/2.  Bound one reproduces the pinned
# 27 polynomial sections; bounds three and four add the next two MW shells.
if arguments.mw_norm_bound == 1:
    section_pool = []
    for known in polynomial_sections:
        section_pool.append((
            {"polynomial_section_indices": [int(known["index"])]},
            old_field(polynomial(known["X_coefficients_low_to_high"])),
            old_field(polynomial(known["Y_coefficients_low_to_high"])),
        ))
else:
    curve = EllipticCurve(
        old_field, [0, 0, 0, old_field(A), old_field(B)]
    )
    pinned_points = tuple(
        curve(
            old_field(polynomial(row["X_coefficients_low_to_high"])),
            old_field(polynomial(row["Y_coefficients_low_to_high"])),
        )
        for row in polynomial_sections
    )
    torsion = (curve(0),)+tuple(point for point in pinned_points if point[1] == 0)
    nontorsion = tuple(point for point in pinned_points if point[1] != 0)
    nontorsion_set = set(nontorsion)
    P = nontorsion[0]
    Q = next(
        point for point in nontorsion
        if P+point in nontorsion_set and P-point not in nontorsion_set
    )
    coefficient_bound = arguments.mw_norm_bound+1
    pool_by_point = {}
    for left in range(-coefficient_bound, coefficient_bound+1):
        for right in range(-coefficient_bound, coefficient_bound+1):
            norm = left**2-left*right+right**2
            if norm > arguments.mw_norm_bound:
                continue
            for torsion_index, torsion_point in enumerate(torsion):
                point = left*P+right*Q+torsion_point
                if point.is_zero():
                    continue
                pool_by_point.setdefault(point, {
                    "MW_coefficients": [left, right],
                    "MW_norm": norm,
                    "torsion_index": torsion_index,
                })
    section_pool = [
        (label, point[0], point[1])
        for point, label in pool_by_point.items()
    ]

all_results = []
for section in selected:
    Xs = fraction(
        "X_numerator_coefficients_low_to_high",
        "X_denominator_coefficients_low_to_high",
        section,
    )
    Ys = fraction(
        "Y_numerator_coefficients_low_to_high",
        "Y_denominator_coefficients_low_to_high",
        section,
    )
    assert Ys**2 == Xs**3+old_field(A)*Xs+old_field(B)
    mu_to_labels = {}
    for label, Xq, Yq in section_pool:
        assert Yq**2 == Xq**3+old_field(A)*Xq+old_field(B)
        if Xq == Xs:
            continue
        mu = (Yq+Ys)/(Xq-Xs)
        mu_to_labels.setdefault(mu, []).append(label)
    mu_rows = tuple(mu_to_labels.items())
    genus_one = []
    tests = 0
    pair_ordinal = 0
    for left in range(len(mu_rows)):
        mu_left, left_labels = mu_rows[left]
        for right in range(left+1, len(mu_rows)):
            mu_right, right_labels = mu_rows[right]
            current_ordinal = pair_ordinal
            pair_ordinal += 1
            if current_ordinal < arguments.pair_start:
                continue
            if (
                arguments.pair_end is not None
                and current_ordinal >= arguments.pair_end
            ):
                continue
            tests += 1
            chord = mu_left+old_field(R)*(mu_right-mu_left)
            branch = (
                chord**4-6*Xs*chord**2-8*Ys*chord
                -3*Xs**2-4*old_field(A)
            )
            square_class = branch.numerator()*branch.denominator()
            if arguments.max_radical_degree is not None:
                cheap_degree = cheap_specialized_radical_degree(square_class)
                if (
                    cheap_degree is not None
                    and cheap_degree > arguments.max_radical_degree
                ):
                    continue
            radical = square_class//gcd(square_class, square_class.derivative())
            radical_degree = int(radical.degree())
            if (
                arguments.max_radical_degree is not None
                and radical_degree > arguments.max_radical_degree
            ):
                continue
            factorization = square_class.factor()
            odd_part = old_ring.one()
            for factor, exponent in factorization:
                if int(exponent) % 2:
                    odd_part *= factor
            odd_degree = int(odd_part.degree())
            if odd_degree != 4:
                continue
            quartic = odd_part.monic()
            twist = base(factorization.unit())
            classification = classify_quartic(quartic, twist)
            row = {
                "section_index": int(section["index"]),
                "left_section_labels": left_labels,
                "right_section_labels": right_labels,
                "radical_degree": radical_degree,
                "cover_factor_degrees_exponents": [
                    [int(factor.degree()), int(exponent)]
                    for factor, exponent in factorization
                ],
                "odd_degree": odd_degree,
                "classification": classification,
            }
            genus_one.append(row)
            print(
                "Q80FINALQ6SECANT|"
                f"S={section['index']}|left={left_labels}|right={right_labels}|"
                f"factors={tuple(tuple(item) for item in row['cover_factor_degrees_exponents'])}|"
                f"classification={classification.get('root_data', classification['status'])}",
                flush=True,
            )
    targets = [
        row for row in genus_one
        if row["classification"].get("root_data") == [15, 66, 800]
    ]
    summary = {
        "section_index": int(section["index"]),
        "MW_norm_bound": arguments.mw_norm_bound,
        "max_radical_degree": arguments.max_radical_degree,
        "section_pool_size": len(section_pool),
        "distinct_secant_values": len(mu_rows),
        "pair_tests": tests,
        "total_pairs": pair_ordinal,
        "pair_start": arguments.pair_start,
        "pair_end": arguments.pair_end,
        "genus_one_hits": len(genus_one),
        "target_hits": len(targets),
        "results": genus_one,
    }
    all_results.append(summary)
    print(
        "Q80FINALQ6SECANT|"
        f"S={section['index']}|values={len(mu_rows)}|tests={tests}|"
        f"genus_one={len(genus_one)}|targets={len(targets)}|"
        "status=PASS_EXACT_BOUNDED_SECANT_SEARCH",
        flush=True,
    )

target_rows = [
    row
    for summary in all_results
    for row in summary["results"]
    if row["classification"].get("root_data") == [15, 66, 800]
]
if arguments.write_artifact:
    output = {
        "schema": "q80-final-q6-secant-gauge-search-gf73-v1",
        "status": "exact_bounded_finite_field_search",
        "prime": 73,
        "source_artifacts": [
            {"path": str(path.relative_to(ROOT)), "sha256": digest}
            for path, digest in KNOWN_HASHES.items()
        ],
        "results": all_results,
        "target_hit_count": len(target_rows),
        "rank_claim": None,
        "reproduce": (
            "sage elkies-k3/scripts/search_q80_final_q6_secant_gauges_gf73.sage "
            + ("--all-sections " if arguments.all_sections else f"--section-index {arguments.section_index} ")
            + "--write-artifact"
        ),
    }
    output_path = ROOT / "artifacts/local/q80-final-q6-secant-gauge-search-gf73.json"
    encoded = json.dumps(output, indent=2, sort_keys=True, default=int)+"\n"
    output_path.write_text(encoded)
    print(
        "Q80FINALQ6SECANT|"
        f"artifact={output_path}|sha256={hashlib.sha256(encoded.encode()).hexdigest()}|"
        f"targets={len(target_rows)}|status=PASS_ARTIFACT_WRITE",
        flush=True,
    )

assert all_results
