#!/usr/bin/env sage -python
"""Extract finite-extension points from an msolve grevlex basis."""

import argparse
import json
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


parser = argparse.ArgumentParser()
parser.add_argument("--gb", type=Path, required=True)
parser.add_argument("--prime", type=int, required=True)
parser.add_argument("--modulus", default="z^2+1")
parser.add_argument("--seed", type=Path)
parser.add_argument(
    "--seed-index",
    type=int,
    default=0,
    help="zero-based record index when the seed file contains multiple MW3SEED lines",
)
parser.add_argument("--output", type=Path)
parser.add_argument("--y-mode", choices=("recursive", "sparse"), default="recursive")
parser.add_argument("--basis-marking", choices=("original", "resolved"), default="original")
parser.add_argument(
    "--variable-rotation",
    type=int,
    default=0,
    help="cyclically rotate the triangular-decomposition variable order",
)
args = parser.parse_args()

base_polynomials = PolynomialRing(GF(args.prime), "z")
z = base_polynomials.gen()
extension = GF(args.prime**2, "z", modulus=base_polynomials(args.modulus))
names = ["c", "c0_inverse", "c1_inverse", "x1", "x2", "x3", "x4"]
if args.y_mode == "recursive":
    names += ["y1", "y1_inverse"]
else:
    names += [f"y{index}" for index in range(1, 8)]
    names += ["y0_tangent_inverse", "y1_tangent_inverse", "yi_tangent_inverse"]
rotation = args.variable_rotation % len(names)
names = tuple(names[rotation:] + names[:rotation])
ring = PolynomialRing(extension, names, order="degrevlex")
text = args.gb.resolve().read_text()
body = text[text.index("[") + 1:text.rindex("]")]
polynomials = [ring(item.replace("^", "**")) for item in body.split(",\n") if item.strip()]
ideal = ring.ideal(polynomials)
quotient_dimension = ideal.vector_space_dimension()
print(
    f"NS0024P4GB|basis={len(polynomials)}|dimension={ideal.dimension()}"
    f"|quotient={quotient_dimension}",
    flush=True,
)
solutions = ideal.variety(ring=extension)
print(f"NS0024P4VARIETY|solutions={len(solutions)}", flush=True)
for solution in solutions[:5]:
    print("|".join(f"{name}={solution[ring(name)]}" for name in names), flush=True)

if args.seed is not None:
    seed_records = [line for line in args.seed.resolve().read_text().splitlines() if line.strip()]
    if args.seed_index < 0 or args.seed_index >= len(seed_records):
        raise SystemExit("--seed-index is outside the seed file")
    seed_fields = {}
    for item in seed_records[args.seed_index].strip().split("|")[1:]:
        key, value = item.split("=", 1)
        seed_fields[key] = value

    def seed_values(key):
        return [extension(int(value)) for value in seed_fields[key].split(",")]

    function_polynomials = PolynomialRing(extension, "t")
    t = function_polynomials.gen()
    function_field = function_polynomials.fraction_field()
    A = function_polynomials(seed_values("A"))
    B = function_polynomials(seed_values("B"))
    r1 = extension(int(seed_fields["r1"]))
    ri = extension(int(seed_fields["ri"]))
    curve = EllipticCurve(function_field, [0, 0, 0, A, B])

    def seed_point(index):
        X = function_polynomials(seed_values(f"P{index}X"))
        Y = function_polynomials(seed_values(f"P{index}Y"))
        return curve(function_field(X), function_field(Y))

    first = [seed_point(index) for index in (1, 2, 3)]

    def finite_value(value, support):
        numerator, denominator = value.numerator(), value.denominator()
        if denominator(support) == 0:
            return None
        return numerator(support) / denominator(support)

    def hits_node(point, fibre):
        if point.is_zero():
            return False
        if fibre < 2:
            x_value = finite_value(point[0], extension(fibre))
            y_value = finite_value(point[1], extension(fibre))
            return x_value == (r1 if fibre else 1) and y_value == 0
        x_num, x_den = point[0].numerator(), point[0].denominator()
        y_num, y_den = point[1].numerator(), point[1].denominator()
        x_excess = x_num.degree() - x_den.degree()
        y_excess = y_num.degree() - y_den.degree()
        x_value = 0 if x_excess < 4 else x_num.leading_coefficient() / x_den.leading_coefficient() if x_excess == 4 else None
        y_value = 0 if y_excess < 6 else y_num.leading_coefficient() / y_den.leading_coefficient() if y_excess == 6 else None
        return x_value == ri and y_value == 0

    def component_label(point, reference, order, fibre):
        labels = [multiplier for multiplier in range(order) if not hits_node(point - multiplier * reference, fibre)]
        return labels[0] if len(labels) == 1 else -1

    def intersection(left, right):
        difference = left - right
        if difference.is_zero():
            return -2
        numerator, denominator = difference[0].numerator(), difference[0].denominator()
        degree = denominator.degree() + max(0, numerator.degree() - denominator.degree() - 4)
        assert degree % 2 == 0
        return degree // 2

    if args.basis_marking == "original":
        declared_profiles = ((1,0,0), (2,1,3), (2,1,1), (1,1,1))
        expected_relative_profiles = declared_profiles
        expected_intersections = (
            (-2,1,2,1), (1,-2,0,1), (2,0,-2,1), (1,1,1,-2)
        )
    else:
        declared_profiles = ((6,0,0), (2,1,1), (4,2,0), (6,4,3))
        # component_label(Q,P4) returns Q's label in units of P4, not the
        # absolute root-chain label.  P4=(6,4,3) is invertible in all three
        # component groups, giving the following relative coordinates.
        expected_relative_profiles = ((1,0,0), (5,4,3), (3,3,0), (1,1,1))
        expected_intersections = (
            (-2,1,1,1), (1,-2,0,2), (1,0,-2,2), (1,2,2,-2)
        )
    accepted = []
    marking_counts = Counter()
    for solution in solutions:
        c = solution[ring("c")]
        h = t - c
        xs = [solution[ring(f"x{index}")] for index in range(1, 5)]
        x5 = r1 * (1-c)**2 - c**2 - ri - sum(xs)
        X = function_polynomials([c**2] + xs + [x5, ri])
        rhs = X**3 + A*X*h**4 + B*h**6
        if args.y_mode == "recursive":
            y = [extension.zero(), solution[ring("y1")]]
            y1_inverse = solution[ring("y1_inverse")]
            for degree in range(3, 10):
                known = sum(y[left]*y[degree-left] for left in range(2, degree-1))
                y.append((rhs[degree]-known)*y1_inverse/2)
        else:
            y = [extension.zero()] + [
                solution[ring(f"y{index}")] for index in range(1, 8)
            ]
            y.append(-sum(y[1:]))
        Y = function_polynomials(y + [extension.zero()])
        assert Y**2 == rhs and sum(y[1:]) == 0
        fourth = curve(function_field(X/h**2), function_field(Y/h**3))
        points = first + [fourth]
        profiles = tuple(
            tuple(component_label(point, fourth, order, fibre) for fibre, order in enumerate((7,5,4)))
            for point in points
        )
        intersections = tuple(tuple(intersection(left, right) for right in points) for left in points)
        marking_counts[(profiles, intersections)] += 1
        if profiles == expected_relative_profiles and intersections == expected_intersections:
            accepted.append((solution, X, Y, h))
    for rank, ((profiles, intersections), count) in enumerate(marking_counts.most_common(5), 1):
        print(
            "NS0024P4MARKINGCLASS|rank={}|count={}|relative_profiles={}|gram={}".format(
                rank, count, profiles, intersections
            ),
            flush=True,
        )
    print(f"NS0024P4MARKING|accepted={len(accepted)}", flush=True)
    if not accepted:
        raise SystemExit("no residue-algebra point has the exact MW4 marking")

    def encoded(value):
        coefficients = list(extension(value).polynomial())
        return [int(coefficients[index]) if index < len(coefficients) else 0 for index in range(2)]

    selected, X, Y, h = accepted[0]
    payload = {
        "schema": "elkies-k3.lattice-foundry-ns0024-mw4-point-modp.v1",
        "status": "PASS_EXACT_MW4_MARKED_POINT_OVER_QUADRATIC_EXTENSION",
        "prime": args.prime,
        "extension": {"generator": "z", "modulus": args.modulus},
        "basis_marking": (
            "resolved_component_depth_recommendation"
            if args.basis_marking == "resolved" else "original"
        ),
        "mw3_seed": str(args.seed),
        "mw3_seed_index": args.seed_index,
        "P4": {
            "X_coefficients_low_to_high": [encoded(value) for value in X],
            "Y_coefficients_low_to_high": [encoded(value) for value in Y],
            "H_coefficients_low_to_high": [encoded(value) for value in h],
        },
        "component_profiles_I7_I5_I4": [list(item) for item in declared_profiles],
        "section_intersection_gram": [list(item) for item in expected_intersections],
        "residue_algebra_dimension": int(quotient_dimension),
        "quadratic_extension_solutions": len(solutions),
        "exact_marking_solutions": len(accepted),
        "proof_boundary": "This is an exact finite-field marked point, not yet a one-dimensional component or a characteristic-zero reconstruction.",
    }
    if args.output is not None:
        args.output.resolve().write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
