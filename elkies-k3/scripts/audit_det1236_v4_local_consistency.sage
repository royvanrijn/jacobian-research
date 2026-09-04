#!/usr/bin/env sage
"""Exhaustive local audit of the determinant-1236 V4 tower.

The asserted tower has common elliptic quotient

    E = 618f1: Y^2 + X*Y = X^3 - 185*X + 1401.

The known genus-two quotient is the quadratic cover with squareclass

    h = (X-4)/54,

whose Prym is 618e1.  If the other two character spaces are respectively
618a1 x 618b1 and 618c1 x 618d1, a genus-three double cover z^2=b of E and
the complementary cover z^2=h*b must have those two Prym traces, in either
order.

Every genus-three quadratic extension of F_5(E) has a representative in
L(6*O): multiply its squareclass by a square so that

    div(b) = (four branch points) + 2*Q - 6*O.

By default this script enumerates every point of P(L(6*O))(F_5), retains
exactly the classes with branch degree four, and compares both degree-two
character sums.  At an even zero (notably the auxiliary double zero 2*Q),
the character contribution is computed from the leading local unit rather
than from the zero point value.  An optional good odd prime can be passed as
the first argument for diagnostics.  Degree two is insensitive to all
quadratic constant twists.

The result is a consistency audit, not an arithmetic exclusion of the
marked curve.  Local compatibility neither constructs the characteristic-
zero squareclass nor decides any rational lift.
"""

import json
import sys
from itertools import combinations, product
from pathlib import Path


arguments = list(sys.argv[1:])
write_artifact = "write" in arguments
arguments = [value for value in arguments if value != "write"]
use_alternate_elliptic_base = "alternate-e" in arguments
arguments = [value for value in arguments if value != "alternate-e"]
P = 5 if not arguments else ZZ(arguments[0])
assert len(arguments) <= 1
assert P.is_prime() and P not in (2, 3, 103)
FIELD = GF(P)
EXTENSION = GF(P**2, "z")
if use_alternate_elliptic_base:
    BASE_LABEL = "618e1"
    BASE_A_INVARIANTS = [0, -90, 0, 3969, 157464]
    KNOWN_COVER_PRYM = "618f1"
    KNOWN_COVER_SQUARECLASS = "X"
    KNOWN_H_ROOT = 0
    KNOWN_H_SCALE = 1
else:
    BASE_LABEL = "618f1"
    BASE_A_INVARIANTS = [1, 0, 0, -185, 1401]
    KNOWN_COVER_PRYM = "618e1"
    KNOWN_COVER_SQUARECLASS = "(X-4)/54"
    KNOWN_H_ROOT = 4
    KNOWN_H_SCALE = 54
E_RATIONAL = EllipticCurve(QQ, BASE_A_INVARIANTS)
E = EllipticCurve(EXTENSION, BASE_A_INVARIANTS)
POLYNOMIAL_RING = PolynomialRing(FIELD, "X")
X = POLYNOMIAL_RING.gen()
A1, A2, A3, A4, A6 = map(FIELD, BASE_A_INVARIANTS)
WEIERSTRASS_RIGHT_SIDE = X**3 + A2*X**2 + A4*X + A6
X_MAP_DISCRIMINANT = (A1*X+A3)**2 + 4*WEIERSTRASS_RIGHT_SIDE
KNOWN_H_NUMERATOR = X-FIELD(KNOWN_H_ROOT)
B_POLYNOMIAL = 1944*X**6 + 441*X**4 - 90*X**2 + 9
POLE_WEIGHTS = [0, 2, 3, 4, 5, 6]
LOCAL_PRECISION = 7

assert E_RATIONAL.has_good_reduction(P)
assert B_POLYNOMIAL.degree() == 6
assert B_POLYNOMIAL.gcd(B_POLYNOMIAL.derivative()) == 1


def projective_vectors(field, vector_dimension):
    for pivot in range(vector_dimension):
        for tail in product(field, repeat=vector_dimension-pivot-1):
            yield vector(field, [0]*pivot + [1] + list(tail))


def quadratic_character(value):
    if not value:
        return 0
    return 1 if value.is_square() else -1


def power_trace_two(curve, prime):
    return ZZ(curve.ap(prime))**2 - 2*prime


def pair_target(first, second):
    return -power_trace_two(EllipticCurve(first), P) - power_trace_two(
        EllipticCurve(second), P
    )


def pole_data(coefficients):
    nonzero_indices = [
        index for index, value in enumerate(coefficients) if value
    ]
    pole_degree = max(POLE_WEIGHTS[index] for index in nonzero_indices)
    leading_index = POLE_WEIGHTS.index(pole_degree)
    return pole_degree, coefficients[leading_index]


def coefficient_polynomials(coefficients):
    a_value = (
        coefficients[0]
        + coefficients[1]*X
        + coefficients[3]*X**2
        + coefficients[5]*X**3
    )
    b_value = coefficients[2] + coefficients[4]*X
    return a_value, b_value


def norm_from_polynomials(a_value, b_value):
    return (
        a_value**2
        - (A1*X+A3)*a_value*b_value
        - WEIERSTRASS_RIGHT_SIDE*b_value**2
    )


def norm_polynomial(coefficients):
    return norm_from_polynomials(*coefficient_polynomials(coefficients))


def branch_degree_from_polynomials(a_value, b_value, pole_degree):
    # A common factor g(X) of A and B vanishes at both points above a
    # non-ramified root of g.  Its square in Norm(b) would erase both odd
    # valuations if one inspected only the norm parity.  Separate it first.
    common = a_value.gcd(b_value).monic()
    primitive_a = a_value // common
    primitive_b = b_value // common
    primitive_norm = norm_from_polynomials(primitive_a, primitive_b)
    common_factors = {
        factor.monic(): int(multiplicity)
        for factor, multiplicity in common.factor()
    }
    primitive_factors = {
        factor.monic(): int(multiplicity)
        for factor, multiplicity in primitive_norm.factor()
    }
    affine_degree = 0
    for factor in set(common_factors).union(primitive_factors):
        common_multiplicity = common_factors.get(factor, 0)
        primitive_multiplicity = primitive_factors.get(factor, 0)
        if (X_MAP_DISCRIMINANT % factor).is_zero():
            # Over a ramified x-value, g(x) has even valuation 2*ord(g).
            fibre_contribution = primitive_multiplicity % 2
        else:
            # At the two points over an unramified x-value, a primitive
            # A+B*Y can vanish at at most one point.  Combine its parity with
            # the common g(X) parity before counting; adding the two branch
            # counts separately would double-count their overlap.
            fibre_contribution = (
                (common_multiplicity + primitive_multiplicity) % 2
                + common_multiplicity % 2
            )
        affine_degree += factor.degree()*fibre_contribution
    return affine_degree + pole_degree % 2


def branch_degree(coefficients):
    pole_degree, _ = pole_data(coefficients)
    return branch_degree_from_polynomials(
        *coefficient_polynomials(coefficients), pole_degree
    )


def product_branch_degree(coefficients):
    pole_degree, _ = pole_data(coefficients)
    a_value, b_value = coefficient_polynomials(coefficients)
    return branch_degree_from_polynomials(
        KNOWN_H_NUMERATOR*a_value,
        KNOWN_H_NUMERATOR*b_value,
        pole_degree+2,
    )


def local_basis_expansion(point):
    """Expand the L(6O) basis in a uniformizer through order six."""
    series_ring = PowerSeriesRing(
        EXTENSION, "t", default_prec=LOCAL_PRECISION
    )
    t_value = series_ring.gen()
    x_zero, y_zero = point[0], point[1]
    partial_y = 2*y_zero + EXTENSION(A1)*x_zero + EXTENSION(A3)
    if partial_y:
        x_series = series_ring(x_zero) + t_value
        y_series = series_ring(y_zero)
        derivative = partial_y
        solve_for_y = True
    else:
        y_series = series_ring(y_zero) + t_value
        x_series = series_ring(x_zero)
        derivative = (
            EXTENSION(A1)*y_zero
            - 3*x_zero**2
            - 2*EXTENSION(A2)*x_zero
            - EXTENSION(A4)
        )
        assert derivative
        solve_for_y = False

    def equation_value(x_argument, y_argument):
        return (
            y_argument**2
            + EXTENSION(A1)*x_argument*y_argument
            + EXTENSION(A3)*y_argument
            - x_argument**3
            - EXTENSION(A2)*x_argument**2
            - EXTENSION(A4)*x_argument
            - EXTENSION(A6)
        )

    for degree in range(1, LOCAL_PRECISION):
        error = equation_value(x_series, y_series)[degree]
        correction = -error/derivative
        if solve_for_y:
            y_series += correction*t_value**degree
        else:
            x_series += correction*t_value**degree
    final_error = equation_value(x_series, y_series)
    assert all(
        final_error[degree] == 0
        for degree in range(LOCAL_PRECISION)
    )
    basis = [
        series_ring(1),
        x_series,
        y_series,
        x_series**2,
        x_series*y_series,
        x_series**3,
    ]
    return [
        [basis[index][degree] for index in range(6)]
        for degree in range(LOCAL_PRECISION)
    ]


def local_order_and_leading(coefficients, expansion_rows):
    for order, row in enumerate(expansion_rows):
        leading = sum(
            coefficients[index]*row[index] for index in range(6)
        )
        if leading:
            return order, leading
    raise ArithmeticError(
        "nonzero element of L(6O) vanished to order greater than six"
    )


LOCAL_DATA = []
for local_point in E:
    if local_point.is_zero():
        continue
    local_rows = local_basis_expansion(local_point)
    h_coefficients = [
        -EXTENSION(KNOWN_H_ROOT), EXTENSION(1), 0, 0, 0, 0
    ]
    h_order, h_leading = local_order_and_leading(
        h_coefficients, local_rows
    )
    LOCAL_DATA.append(
        (
            local_point,
            local_rows,
            h_order,
            h_leading/EXTENSION(KNOWN_H_SCALE),
        )
    )


def character_sums(coefficients):
    lifted = [EXTENSION(value) for value in coefficients]
    pole_degree, leading = pole_data(lifted)
    if pole_degree == 0:
        sum_b = quadratic_character(leading)
    elif pole_degree % 2:
        sum_b = 0
    else:
        sum_b = quadratic_character(leading)

    product_pole_degree = pole_degree + 2
    sum_product = (
        0
        if product_pole_degree % 2
        else quadratic_character(leading / EXTENSION(KNOWN_H_SCALE))
    )
    for _, expansion_rows, h_order, h_leading in LOCAL_DATA:
        b_order, b_leading = local_order_and_leading(
            lifted, expansion_rows
        )
        if b_order % 2 == 0:
            sum_b += quadratic_character(b_leading)
        product_order = b_order + h_order
        if product_order % 2 == 0:
            sum_product += quadratic_character(b_leading*h_leading)
    return int(sum_b), int(sum_product)


target_ab = int(pair_target("618a1", "618b1"))
target_cd = int(pair_target("618c1", "618d1"))
known_cover_target = int(
    -power_trace_two(EllipticCurve(KNOWN_COVER_PRYM), P)
)
pair_partitions = {
    "ab|cd": (target_ab, target_cd),
    "ac|bd": (
        int(pair_target("618a1", "618c1")),
        int(pair_target("618b1", "618d1")),
    ),
    "ad|bc": (
        int(pair_target("618a1", "618d1")),
        int(pair_target("618b1", "618c1")),
    ),
}

# A bridge diagnostic: allow every dimension-two block in the complete
# level-618 newspace, not just pairs among 618a1,...,618d1.  This tests
# whether the local tower points to a different four-dimensional Prym.
dimension_two_blocks = []
for left, right in combinations("abcdg", 2):
    dimension_two_blocks.append(
        {
            "label": "%s%s" % (left, right),
            "support": [left, right],
            "target": int(pair_target("618%s1" % left, "618%s1" % right)),
            "kind": "two_rational_elliptic_factors",
        }
    )

newspace = ModularSymbols(618, 2, sign=1).cuspidal_subspace().new_subspace()
for orbit_index, factor_space in enumerate(newspace.decomposition()):
    if factor_space.dimension() != 2:
        continue
    eigenform = factor_space.q_eigenform(P+1, names="u%s" % orbit_index)
    ap_value = eigenform[P]
    coefficient_field_polynomial = eigenform.base_ring().defining_polynomial()
    signs = []
    for atkin_prime in (2, 3, 103):
        operator = factor_space.atkin_lehner_operator(atkin_prime).matrix()
        signs.append(int(ZZ(operator[0, 0])))
    dimension_two_blocks.append(
        {
            "label": "newform_orbit_%s" % orbit_index,
            "support": ["newform_orbit_%s" % orbit_index],
            "target": int(-(ap_value**2-2*P).trace()),
            "kind": "quadratic_coefficient_field_newform",
            "classical_signs_W2_W3_W103": signs,
            # Sage may choose different primitive generators on repeated
            # runs.  The squarefree polynomial discriminant identifies the
            # same quadratic field deterministically.
            "coefficient_field_squarefree_discriminant": int(
                ZZ(coefficient_field_polynomial.discriminant()).squarefree_part()
            ),
        }
    )

known_h_coefficients = vector(
    FIELD, [-KNOWN_H_ROOT, 1, 0, 0, 0, 0]
)
assert branch_degree(known_h_coefficients) == 2

# Unit-test the character convention on the already explicit B -> E cover.
known_cover_sum = quadratic_character(
    EXTENSION(1)/EXTENSION(KNOWN_H_SCALE)
)
for _, _, h_order, h_leading in LOCAL_DATA:
    if h_order % 2 == 0:
        known_cover_sum += quadratic_character(h_leading)
assert known_cover_sum == known_cover_target

# A literal square must contribute +1 at every place, including its double
# zeros.  Multiplication by that square must leave the known h character sum
# unchanged.  This guards against the tempting but incorrect use of chi(b(P))
# at an even zero of b.
test_root = FIELD(10)
test_square = vector(
    FIELD, [test_root**2, -2*test_root, 0, 1, 0, 0]
)
test_square_sum, test_square_times_h_sum = character_sums(test_square)
assert test_square_sum == E.cardinality()
assert test_square_times_h_sum == known_cover_sum

projective_count = 0
branch_four_count = 0
both_branch_four_count = 0
unfiltered_compatible_by_partition = {name: [] for name in pair_partitions}
signature_counts = {}
both_branch_signature_counts = {}
compatible_by_partition = {name: [] for name in pair_partitions}
for coefficients in projective_vectors(FIELD, 6):
    projective_count += 1
    signature = character_sums(coefficients)
    for name, targets in pair_partitions.items():
        if signature in (targets, tuple(reversed(targets))):
            unfiltered_compatible_by_partition[name].append(
                [int(value) for value in coefficients]
            )
    if branch_degree(coefficients) != 4:
        continue
    branch_four_count += 1
    signature_counts[signature] = signature_counts.get(signature, 0) + 1
    if product_branch_degree(coefficients) != 4:
        continue
    both_branch_four_count += 1
    both_branch_signature_counts[signature] = (
        both_branch_signature_counts.get(signature, 0) + 1
    )
    for name, targets in pair_partitions.items():
        if signature in (targets, tuple(reversed(targets))):
            compatible_by_partition[name].append(
                [int(value) for value in coefficients]
            )

assert projective_count == (P**6-1)//(P-1)
claimed_compatible = compatible_by_partition["ab|cd"]

if claimed_compatible:
    status = "LOCALLY_COMPATIBLE"
    conclusion = (
        "At p=%s the displayed B/E squareclass and asserted ab|cd V4 "
        "factor partition have local genus-three classes. This local "
        "diagnostic neither constructs the characteristic-zero cover nor "
        "decides a rational lift."
    ) % P
else:
    status = "EXACT_UPSTREAM_INCONSISTENCY"
    conclusion = (
        "The displayed B/E squareclass and the asserted ab|cd V4 factor "
        "partition admit no genus-three quadratic cover over F_%s in "
        "either order. The other pair partitions are diagnostic only. "
        "This is not an arithmetic exclusion of C_1236; the upstream "
        "quotient/factor identification must be reconciled before lifting."
    ) % P

compatible_dimension_four_partitions = []
for left, right in combinations(dimension_two_blocks, 2):
    if set(left["support"]).intersection(right["support"]):
        continue
    targets = (left["target"], right["target"])
    compatible_count = sum(
        both_branch_signature_counts.get(signature, 0)
        for signature in set((targets, tuple(reversed(targets))))
    )
    if compatible_count:
        compatible_dimension_four_partitions.append(
            {
                "blocks": [left["label"], right["label"]],
                "targets": list(targets),
                "compatible_classes": int(compatible_count),
            }
        )

payload = {
    "schema": "elkies-k3.det1236-v4-local-consistency.v1",
    "status": status,
    "prime": P,
    "field_extension_degree": 2,
    "elliptic_base": BASE_LABEL,
    "known_cover_squareclass": KNOWN_COVER_SQUARECLASS,
    "known_cover_prym": KNOWN_COVER_PRYM,
    "known_cover_character_sum": int(known_cover_sum),
    "known_cover_target": known_cover_target,
    "claimed_complementary_pryms": [
        ["618a1", "618b1"],
        ["618c1", "618d1"],
    ],
    "degree_two_targets": {"ab": target_ab, "cd": target_cd},
    "all_pair_partitions": {
        name: {
            "left_target": targets[0],
            "right_target": targets[1],
            "compatible_classes": compatible_by_partition[name],
        }
        for name, targets in pair_partitions.items()
    },
    "unfiltered_compatible_class_counts": {
        name: len(values)
        for name, values in unfiltered_compatible_by_partition.items()
    },
    "projective_classes_enumerated": int(projective_count),
    "branch_degree_four_classes": int(branch_four_count),
    "both_quadratic_characters_branch_degree_four_classes": int(
        both_branch_four_count
    ),
    "distinct_signature_pairs": int(len(signature_counts)),
    "distinct_both_branch_four_signature_pairs": int(
        len(both_branch_signature_counts)
    ),
    "both_branch_four_signature_histogram": [
        {
            "signature": list(signature),
            "classes": int(count),
        }
        for signature, count in sorted(both_branch_signature_counts.items())
    ],
    "dimension_two_block_targets": dimension_two_blocks,
    "compatible_dimension_four_partitions": compatible_dimension_four_partitions,
    "compatible_classes": claimed_compatible,
    "quadratic_twist_boundary": (
        "degree-two character sums are unchanged by quadratic constant twists"
    ),
    "conclusion": conclusion,
}
rendered = json.dumps(payload, indent=2, sort_keys=True, default=int)+"\n"
if write_artifact:
    repository_root = Path.cwd().resolve()
    if not (repository_root / "elkies-k3").is_dir():
        raise RuntimeError("run this certificate from the repository root")
    if use_alternate_elliptic_base:
        raise ValueError("alternate-base diagnostics do not have a canonical artifact")
    suffix = "" if P == 5 else "-p%s" % P
    output = repository_root / (
        "artifacts/generated-results/"
        "elkies-k3-det1236-v4-local-consistency%s-v1.json" % suffix
    )
    output.write_text(rendered)
print(rendered, end="")
