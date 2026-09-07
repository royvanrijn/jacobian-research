#!/usr/bin/env sage
"""Search polynomial sections on the compensated pair23 fifth over GF(73).

The exact fifth model has fibers ``I0*+3I4+3I2`` and Weierstrass degrees
``(4,9)``.  A polynomial section has ``deg(x)<=4`` and ``deg(y)<=6``.  If it
meets four finite nonidentity components, its quartic x-coordinate is fixed
up to one scalar by the four nodal x-values.  We enumerate all four-fiber
supports and all 73 remaining scalars, then test the Weierstrass right-hand
side for being a polynomial square.

This is an exact finite-field search.  It does not by itself prove that the
listed sections lift to characteristic zero.
"""

import hashlib
import json
import argparse
from itertools import combinations
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--match-q6", action="store_true")
parser.add_argument("--write-artifact", action="store_true")
arguments = parser.parse_args()
if arguments.write_artifact:
    arguments.match_q6 = True
SOURCE = (
    ROOT / "artifacts/generated-results/"
    "q80-deforming-fifth-pair23-gf73.json"
)
KNOWN_SOURCE_SHA256 = (
    "23fc49bce2618a6d3c5f5e18ded34b4ffbee220be83523ae250bf7774a91db14"
)
source_bytes = SOURCE.read_bytes()
assert hashlib.sha256(source_bytes).hexdigest() == KNOWN_SOURCE_SHA256
source = json.loads(source_bytes)
assert source["root_data"] == [16, 66, 2048]

finite = GF(73, impl="modn")
polynomial_ring = PolynomialRing(finite, "s")
s = polynomial_ring.gen()
x_ring = PolynomialRing(finite, "x")
x = x_ring.gen()

A = polynomial_ring(source["A_coefficients_low_to_high"])
B = polynomial_ring(source["B_coefficients_low_to_high"])
Delta = -finite(16)*(4*A**3+27*B**2)
assert (A.degree(), B.degree(), Delta.degree()) == (4, 9, 18)


def linear_root(factor):
    assert factor.degree() == 1
    return -factor[0]/factor[1]


i2_roots = tuple(sorted(
    (linear_root(factor) for factor, exponent in Delta.factor() if exponent == 2),
    key=int,
))
i4_roots = tuple(sorted(
    (linear_root(factor) for factor, exponent in Delta.factor() if exponent == 4),
    key=int,
))
assert tuple(map(int, i2_roots)) == (1, 56, 72)
assert tuple(map(int, i4_roots)) == (0, 64, 65)


def node_x(root):
    cubic = x**3+A(root)*x+B(root)
    common = cubic.gcd(cubic.derivative())
    assert common.degree() == 1
    return -common[0]/common[1]


all_roots = i2_roots+i4_roots
nodes = {root: node_x(root) for root in all_roots}


def polynomial_square_roots(polynomial):
    """Return all polynomial square roots of a polynomial of degree <=12."""
    assert polynomial.degree() <= 12
    if polynomial == 0:
        return (polynomial_ring.zero(),)
    shift = next(value for value in finite if polynomial(value) != 0)
    shifted = polynomial(s+shift)
    constant = shifted[0]
    if not constant.is_square():
        return ()
    roots = []
    for first in constant.sqrt(all=True):
        coefficients = [first]
        for degree in range(1, 7):
            known = sum(
                coefficients[left]*coefficients[degree-left]
                for left in range(1, degree)
            )
            coefficients.append((shifted[degree]-known)/(2*first))
        candidate_shifted = polynomial_ring(coefficients)
        if candidate_shifted**2 == shifted:
            roots.append(candidate_shifted(s-shift))
    return tuple(roots)


def constrained_x_family(selected_roots):
    interpolation = polynomial_ring.lagrange_polynomial(
        [(root, nodes[root]) for root in selected_roots]
    )
    vanishing = polynomial_ring.one()
    for root in selected_roots:
        vanishing *= s-root
    assert interpolation.degree() < len(selected_roots)
    return interpolation, vanishing


raw_candidates = []
for selected_roots in combinations(all_roots, 4):
    interpolation, vanishing = constrained_x_family(selected_roots)
    assert vanishing.degree() == 4
    for leading in finite:
        X = interpolation+leading*vanishing
        for Y in polynomial_square_roots(X**3+A*X+B):
            actual_support = tuple(root for root in all_roots if X(root) == nodes[root])
            raw_candidates.append((X, Y, actual_support))

# The same section is found once for every four-subset of its actual support.
unique = {}
for X, Y, support in raw_candidates:
    key = (tuple(X.list()), tuple(Y.list()))
    unique[key] = (X, Y, support)
candidates = tuple(unique.values())

print(
    "Q80FIFTHPOLYSEARCH|prime=73|"
    f"I2_nodes={tuple((int(r), int(nodes[r])) for r in i2_roots)}|"
    f"I4_nodes={tuple((int(r), int(nodes[r])) for r in i4_roots)}|"
    f"tests={15*73}|raw_hits={len(raw_candidates)}|"
    f"sections={len(candidates)}",
    flush=True,
)
for index, (X, Y, support) in enumerate(candidates):
    support_i2 = tuple(int(root) for root in support if root in i2_roots)
    support_i4 = tuple(int(root) for root in support if root in i4_roots)
    print(
        f"Q80FIFTHPOLYSEARCH|index={index}|"
        f"I2_support={support_i2}|I4_support={support_i4}|"
        f"X={tuple(map(int, X.list()))}|Y={tuple(map(int, Y.list()))}",
        flush=True,
    )


if arguments.match_q6:
    function_field = polynomial_ring.fraction_field()
    curve = EllipticCurve(
        function_field,
        [0, 0, 0, function_field(A), function_field(B)],
    )
    points = tuple(
        curve(function_field(X), function_field(Y))
        for X, Y, _ in candidates
    )
    torsion = (curve(0),)+tuple(
        point for point in points if point[1] == 0
    )
    assert len(torsion) == 4
    assert all(2*point == curve(0) for point in torsion)
    nontorsion = tuple(point for point in points if point[1] != 0)
    assert len(nontorsion) == 24
    nontorsion_set = set(nontorsion)

    # The 24 sections are the four 2-torsion translates of the six minimal
    # vectors in the CM Mordell-Weil lattice with Gram
    # [[1/2,-1/4],[-1/4,1/2]].  Pick adjacent minimal vectors P,Q, detected
    # by P+Q again being minimal; they form a free basis with that Gram.
    P = nontorsion[0]
    Q = next(
        point for point in nontorsion
        if P+point in nontorsion_set and P-point not in nontorsion_set
    )
    assert P-Q not in nontorsion_set

    def section_pole(point):
        assert not point.is_zero()
        x_coordinate = point[0]
        numerator_degree = x_coordinate.numerator().degree()
        denominator_degree = x_coordinate.denominator().degree()
        assert denominator_degree % 2 == 0
        infinity_excess = max(0, numerator_degree-denominator_degree-4)
        assert infinity_excess % 2 == 0
        return denominator_degree//2+infinity_excess//2

    def hits_node(point, root):
        if point.is_zero():
            return False
        x_coordinate, y_coordinate = point[0], point[1]
        if x_coordinate.denominator()(root) == 0:
            return False
        return (
            x_coordinate(root) == nodes[root]
            and y_coordinate(root) == 0
        )

    # At I4, an odd component label stays nonidentity after doubling, while
    # label 2 doubles to the identity.  This finds a label-1 reference up to
    # inversion; the correction k(4-k)/4 is orientation independent.
    i4_references = {}
    for root in i4_roots:
        i4_references[root] = next(
            point for point in nontorsion
            if hits_node(point, root) and hits_node(2*point, root)
        )

    def i4_label(point, root):
        reference = i4_references[root]
        answers = tuple(
            multiplier for multiplier in range(4)
            if not hits_node(point-multiplier*reference, root)
        )
        assert len(answers) == 1
        return answers[0]

    norm_seven_coefficients = tuple(
        (left, right)
        for left in range(-4, 5)
        for right in range(-4, 5)
        if left**2-left*right+right**2 == 7
    )
    assert len(norm_seven_coefficients) == 12
    raw_q6 = []
    for left, right in norm_seven_coefficients:
        for torsion_point in torsion:
            point = left*P+right*Q+torsion_point
            if point.is_zero() or section_pole(point) != 1:
                continue
            i2_labels = tuple(
                1 if hits_node(point, root) else 0 for root in i2_roots
            )
            i4_labels = tuple(i4_label(point, root) for root in i4_roots)
            i4_correction = sum(
                QQ(label*(4-label))/4 for label in i4_labels
            )
            # The transported q6 section is identity at all three I2 fibers,
            # meets endpoints of two I4 fibers, and is nonidentity at I0*.
            # Thus the finite multiplicative correction is 3/2 and the D4
            # correction is 1, for the required total 5/2.  The previous
            # two-I2/identity-I0* filter selected the wrong torsion lifts.
            if sum(i2_labels) != 0 or i4_correction != QQ(3)/2:
                continue
            raw_q6.append((point, left, right, torsion_point, i2_labels, i4_labels))

    q6_candidates = {}
    for row in raw_q6:
        point = row[0]
        key = (
            tuple(point[0].numerator().list()),
            tuple(point[0].denominator().list()),
            tuple(point[1].numerator().list()),
            tuple(point[1].denominator().list()),
        )
        q6_candidates[key] = row
    assert q6_candidates

    def rational_coefficients(value):
        return (
            tuple(map(int, value.numerator().list())),
            tuple(map(int, value.denominator().list())),
        )

    print(
        "Q80FIFTHPOLYSEARCH|"
        f"torsion_order={len(torsion)}|minimal_sections={len(nontorsion)}|"
        f"norm7_vectors={len(norm_seven_coefficients)}|"
        f"q6_candidates={len(q6_candidates)}",
        flush=True,
    )
    for index, row in enumerate(q6_candidates.values()):
        point, left, right, torsion_point, i2_labels, i4_labels = row
        print(
            f"Q80FIFTHPOLYSEARCH|q6_index={index}|basis_coefficients={(left, right)}|"
            f"torsion_zero={int(torsion_point.is_zero())}|P.O=1|height=7/2|"
            f"I2_labels={i2_labels}|I4_labels={i4_labels}|"
            f"X_num_den={rational_coefficients(point[0])}|"
            f"Y_num_den={rational_coefficients(point[1])}",
            flush=True,
        )
    if arguments.write_artifact:
        output_rows = []
        for index, row in enumerate(q6_candidates.values()):
            point, left, right, torsion_point, i2_labels, i4_labels = row
            output_rows.append({
                "index": index,
                "basis_coefficients": [left, right],
                "torsion_zero": bool(torsion_point.is_zero()),
                "section_P_dot_O": 1,
                "height": "7/2",
                "I2_labels": list(i2_labels),
                "I4_labels": list(i4_labels),
                "I0star_correction": "1",
                "X_numerator_coefficients_low_to_high": list(map(
                    int, point[0].numerator().list()
                )),
                "X_denominator_coefficients_low_to_high": list(map(
                    int, point[0].denominator().list()
                )),
                "Y_numerator_coefficients_low_to_high": list(map(
                    int, point[1].numerator().list()
                )),
                "Y_denominator_coefficients_low_to_high": list(map(
                    int, point[1].denominator().list()
                )),
            })
        output = {
            "schema": "q80-deforming-fifth-q6-horizontal-candidates-gf73-v2",
            "status": "exact_finite_field_lattice_profile_filtered_section_search",
            "prime": 73,
            "source_artifact": str(SOURCE.relative_to(ROOT)),
            "source_sha256": KNOWN_SOURCE_SHA256,
            "torsion_order": len(torsion),
            "minimal_polynomial_sections": len(nontorsion),
            "height_7_over_2_vectors": len(norm_seven_coefficients),
            "target_component_profile": {
                "I2_labels": [0, 0, 0],
                "I4": "two endpoint labels, total correction 3/2",
                "I0star_correction": "1",
            },
            "q6_candidates": output_rows,
            "rank_claim": None,
            "reproduce": (
                "sage elkies-k3/scripts/"
                "search_q80_deforming_fifth_polynomial_sections_gf73.sage "
                "--match-q6 --write-artifact"
            ),
        }
        output_path = (
            ROOT / "artifacts/generated-results/"
            "q80-deforming-fifth-q6-horizontal-candidates-gf73.json"
        )
        encoded = json.dumps(output, indent=2, sort_keys=True, default=int)+"\n"
        output_path.write_text(encoded)
        print(
            "Q80FIFTHPOLYSEARCH|"
            f"artifact={output_path}|"
            f"sha256={hashlib.sha256(encoded.encode()).hexdigest()}|"
            "status=PASS_ARTIFACT_WRITE",
            flush=True,
        )

        # Preserve the complete finite set separately.  The q6 candidate
        # artifact above remains stable; this companion artifact supplies the
        # rational secant values needed to impose the final vertical module.
        all_section_rows = []
        for index, (X, Y, support) in enumerate(candidates):
            point = points[index]
            all_section_rows.append({
                "index": index,
                "torsion": bool(point[1] == 0),
                "I2_support": [
                    int(root) for root in support if root in i2_roots
                ],
                "I4_support": [
                    int(root) for root in support if root in i4_roots
                ],
                "X_coefficients_low_to_high": list(map(int, X.list())),
                "Y_coefficients_low_to_high": list(map(int, Y.list())),
            })
        all_sections_output = {
            "schema": "q80-deforming-fifth-polynomial-sections-gf73-v1",
            "status": "complete_exact_node_constrained_polynomial_section_set",
            "prime": 73,
            "source_artifact": str(SOURCE.relative_to(ROOT)),
            "source_sha256": KNOWN_SOURCE_SHA256,
            "section_count": len(all_section_rows),
            "nonzero_two_torsion_count": sum(
                row["torsion"] for row in all_section_rows
            ),
            "sections": all_section_rows,
            "rank_claim": None,
            "reproduce": (
                "sage elkies-k3/scripts/"
                "search_q80_deforming_fifth_polynomial_sections_gf73.sage "
                "--match-q6 --write-artifact"
            ),
        }
        all_sections_path = (
            ROOT / "artifacts/generated-results/"
            "q80-deforming-fifth-polynomial-sections-gf73.json"
        )
        all_sections_encoded = json.dumps(
            all_sections_output, indent=2, sort_keys=True, default=int
        )+"\n"
        all_sections_path.write_text(all_sections_encoded)
        print(
            "Q80FIFTHPOLYSEARCH|"
            f"all_sections_artifact={all_sections_path}|"
            f"sha256={hashlib.sha256(all_sections_encoded.encode()).hexdigest()}|"
            "status=PASS_ALL_SECTIONS_ARTIFACT_WRITE",
            flush=True,
        )
    print(
        "Q80FIFTHPOLYSEARCH|status=PASS_Q6_HORIZONTAL_CANDIDATES",
        flush=True,
    )

assert candidates
print("Q80FIFTHPOLYSEARCH|status=PASS_EXACT_SEARCH", flush=True)
