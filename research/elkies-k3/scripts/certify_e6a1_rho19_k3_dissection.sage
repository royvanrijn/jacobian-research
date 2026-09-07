#!/usr/bin/env sage-python
"""Certify the lattice and boundary anatomy of the new E6+A1 rho-19 K3.

status: ACTIVE_PROOF
claim: generic NS/T, saturation, singular boundaries, and first chamber obstruction
inputs: rational-surface quadratic-rank-search certificate
outputs: elkies-k3-e6a1-rho19-k3-dissection-v1.json

This companion to ``certify_rational_surface_quadratic_rank_search.sage``
pins the integral Neron--Severi lattice, its discriminant form and saturation,
the generic transcendental lattice, the positive frame used by elliptic-
neighbor searches, the descent from k to s=k^2, and four forced singular-K3
specializations.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import product
from pathlib import Path

from sage.all import (
    BinaryQF_reduced_representatives,
    CartanMatrix,
    EllipticCurve,
    IntegralLattice,
    NumberField,
    PolynomialRing,
    QQ,
    QuadraticForm,
    ZZ,
    ceil,
    gcd,
    identity_matrix,
    lcm,
    matrix,
    pari,
    prod,
    span,
    vector,
    xgcd,
    zero_matrix,
)


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
DEFAULT_OUTPUT = GEN / "elkies-k3-e6a1-rho19-k3-dissection-v1.json"
SOURCE_CERTIFICATE = GEN / "elkies-k3-rational-surface-quadratic-rank-search-v1.json"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def integer_rows(value):
    return [[int(entry) for entry in row] for row in matrix(ZZ, value).rows()]


def matrix_digest(value):
    encoded = ";".join(
        ",".join(str(entry) for entry in row)
        for row in matrix(ZZ, value).rows()
    )
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def rational_rows(value):
    return [[str(entry) for entry in row] for row in matrix(QQ, value).rows()]


def compact_height(curve, point, multiplier):
    """Compute the K3 height after clearing all fibre component groups."""

    coordinate = (multiplier * point)[0]
    numerator = coordinate.numerator()
    denominator = coordinate.denominator()
    if not denominator.is_square():
        raise ArithmeticError("x-coordinate denominator is not a square")
    intersection = max(
        denominator.degree() // 2,
        ceil((numerator.degree() - 4) / 2),
    )
    return {
        "multiplier": int(multiplier),
        "x_numerator_degree": int(numerator.degree()),
        "x_denominator_degree": int(denominator.degree()),
        "intersection_with_zero": int(intersection),
        "height": str((4 + 2 * intersection) / QQ(multiplier**2)),
    }


def discriminant_form(gram):
    """Return Smith orders and a Gram matrix for discriminant generators."""

    gram = matrix(ZZ, gram)
    smith, left, _ = gram.smith_form()
    indices = [
        index
        for index, entry in enumerate(smith.diagonal())
        if abs(entry) > 1
    ]
    orders = tuple(abs(ZZ(smith[index, index])) for index in indices)
    target = left.inverse()[:, indices]
    dual = gram.inverse() * target
    form = dual.transpose() * gram * dual
    return orders, form


def element_order(element, orders):
    return lcm(
        [
            order // gcd(order, ZZ(coordinate))
            for coordinate, order in zip(element, orders)
        ]
    )


def quadratic_value(element, form):
    item = vector(QQ, element)
    return item * form * item


def bilinear_value(left, right, form):
    return vector(QQ, left) * form * vector(QQ, right)


def congruent(left, right, modulus):
    return (left - right) / modulus in ZZ


def discriminant_forms_are_isometric(first, second):
    """Brute-force the small finite quadratic modules occurring here."""

    first_orders, first_form = discriminant_form(first)
    second_orders, second_form = discriminant_form(second)
    if first_orders != second_orders:
        return False, None
    elements = list(product(*[range(int(order)) for order in second_orders]))
    candidates = []
    for index, order in enumerate(first_orders):
        candidates.append(
            [
                element
                for element in elements
                if element_order(element, second_orders) == order
                and congruent(
                    quadratic_value(element, second_form),
                    first_form[index, index],
                    2,
                )
            ]
        )
    for images in product(*candidates):
        if any(
            not congruent(
                bilinear_value(images[i], images[j], second_form),
                first_form[i, j],
                1,
            )
            for i in range(len(first_orders))
            for j in range(i)
        ):
            continue
        generated = {
            tuple(
                sum(coefficients[j] * images[j][i] for j in range(len(images)))
                % second_orders[i]
                for i in range(len(second_orders))
            )
            for coefficients in product(
                *[range(int(order)) for order in first_orders]
            )
        }
        if len(generated) == prod(first_orders):
            return True, images
    return False, None


def generic_ns_gram():
    """Integral divisor Gram for O,F,2E6,A3,P0,P1."""

    gram = zero_matrix(ZZ, 19)
    gram[0, 0] = -2
    gram[0, 1] = gram[1, 0] = 1
    e6 = CartanMatrix(["E", 6])
    a3 = CartanMatrix(["A", 3])
    gram[2:8, 2:8] = -e6
    gram[8:14, 8:14] = -e6
    gram[14:17, 14:17] = -a3
    p0, p1 = 17, 18
    for section in (p0, p1):
        gram[section, section] = -2
        gram[1, section] = gram[section, 1] = 1
    gram[p0, p1] = gram[p1, p0] = 1
    gram[p0, 2] = gram[2, p0] = 1
    gram[p0, 8] = gram[8, p0] = 1
    gram[p0, 15] = gram[15, p0] = 1
    gram[p1, 15] = gram[15, p1] = 1
    return gram


def split_primitive_isotropic(ns, fibre):
    """Split a primitive U from an even NS lattice."""

    ns = matrix(ZZ, ns)
    fibre = vector(ZZ, fibre)
    if fibre * ns * fibre != 0:
        raise ValueError("fibre is not isotropic")
    current = ZZ(0)
    mate_entries = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(ns * fibre):
        if value == 0:
            continue
        divisor, old_scale, new_scale = xgcd(current, ZZ(value))
        mate_entries = [old_scale * entry for entry in mate_entries]
        mate_entries[index] += new_scale
        current = divisor
    if current == -1:
        mate_entries = [-entry for entry in mate_entries]
    if current not in (1, -1):
        raise ValueError("fibre is not primitive")
    mate = vector(ZZ, mate_entries)
    mate -= (mate * ns * mate // 2) * fibre
    if mate * ns * fibre != 1 or mate * ns * mate != 0:
        raise ArithmeticError("failed to normalize a hyperbolic mate")
    complement = matrix(
        ZZ, [fibre * ns, mate * ns]
    ).right_kernel().basis_matrix()
    frame = -(complement * ns * complement.transpose())
    return frame


def root_data(frame):
    """Return exact (rank,count,root-lattice determinant)."""

    frame = matrix(ZZ, frame)
    minimum = pari(frame).qfminim(2)
    count = int(minimum[0])
    if count == 0:
        return (0, 0, 1)
    half = [vector(ZZ, column) for column in matrix(ZZ, minimum[2]).columns()]
    basis = matrix(ZZ, half + [-item for item in half]).row_module().basis_matrix()
    gram = basis * frame * basis.transpose()
    return (int(basis.rank()), count, abs(int(gram.det())))


def branch_collision_ns(q_meets_new_a1):
    """NS Gram at a smooth-branch/I1 collision before saturation checks."""

    gram = zero_matrix(ZZ, 20)
    gram[0, 0] = -2
    gram[0, 1] = gram[1, 0] = 1
    e6 = CartanMatrix(["E", 6])
    a3 = CartanMatrix(["A", 3])
    gram[2:8, 2:8] = -e6
    gram[8:14, 8:14] = -e6
    gram[14:17, 14:17] = -a3
    gram[17, 17] = -2
    p0, p1 = 18, 19
    for section in (p0, p1):
        gram[section, section] = -2
        gram[1, section] = gram[section, 1] = 1
    gram[p0, p1] = gram[p1, p0] = 1
    gram[p0, 2] = gram[2, p0] = 1
    gram[p0, 8] = gram[8, p0] = 1
    gram[p0, 15] = gram[15, p0] = 1
    gram[p1, 15] = gram[15, p1] = 1
    if q_meets_new_a1:
        gram[p1, 17] = gram[17, p1] = 1
    return gram


def torsion_special_ns(order):
    """Build the saturated NS at s=0 (order 2) or s=-4/3 (order 3)."""

    if order == 2:
        root = CartanMatrix(["E", 7])
        tail = CartanMatrix(["A", 3])
        root_size = 7
        tail_size = 3
        torsion_root_node = 6
        torsion_tail_node = 1
        section_tail_node = 1
    elif order == 3:
        root = CartanMatrix(["E", 6])
        tail = CartanMatrix(["A", 5])
        root_size = 6
        tail_size = 5
        torsion_root_node = 0
        torsion_tail_node = 1
        section_tail_node = 2
    else:
        raise ValueError("only the order-two and order-three boundary points occur")

    dimension = 2 + 2 * root_size + tail_size + 1
    gram = zero_matrix(ZZ, dimension)
    gram[0, 0] = -2
    gram[0, 1] = gram[1, 0] = 1
    first = 2
    second = first + root_size
    tail_start = second + root_size
    section = dimension - 1
    gram[first:second, first:second] = -root
    gram[second:tail_start, second:tail_start] = -root
    gram[tail_start:section, tail_start:section] = -tail
    gram[section, section] = -2
    gram[1, section] = gram[section, 1] = 1
    gram[section, tail_start + section_tail_node] = 1
    gram[tail_start + section_tail_node, section] = 1

    glue = [QQ(0)] * dimension
    glue[first:second] = [
        -entry for entry in root.inverse().column(torsion_root_node)
    ]
    glue[second:tail_start] = [
        -entry for entry in root.inverse().column(torsion_root_node)
    ]
    glue[tail_start:section] = [
        -entry for entry in tail.inverse().column(torsion_tail_node)
    ]
    glue = vector(QQ, glue)
    if glue * gram * glue != -4:
        raise ArithmeticError("torsion glue norm changed")
    if not all(value in ZZ for value in glue * gram):
        raise ArithmeticError("torsion glue is not in the dual lattice")
    if not all(value in ZZ for value in order * glue):
        raise ArithmeticError("torsion glue has the wrong order")
    basis = span(
        list(identity_matrix(QQ, dimension).rows()) + [glue], ZZ
    ).basis_matrix()
    saturated = basis * gram * basis.transpose()
    if not all(value in ZZ for value in saturated.list()):
        raise ArithmeticError("torsion saturation is not integral")
    return saturated.change_ring(ZZ), glue


def specialized_family(field, k_value):
    """Return the K3 curve and its two displayed sections over field(t)."""

    k_value = field(k_value)
    D = 3 * k_value**2 - 4
    c = 2 * k_value / D
    lam = -(k_value**2 - 4) * D / 4
    p = 4 / D**2
    q = k_value**2 / D
    s_coefficient = 8 / D**3
    v = 2 * (k_value**2 + 2) / D**2
    polynomial_ring = PolynomialRing(field, "t")
    t = polynomial_ring.gen()
    function_field = polynomial_ring.fraction_field()
    H = 1 - t**2
    curve = EllipticCurve(
        function_field,
        [
            H**3 * (lam - 3 * H),
            H**4 * (c**2 * lam**2 + lam * H - 2 * H**2),
        ],
    )
    p0 = curve(-H**2, c * lam * H**2)
    p1 = curve(
        p * lam**2 + q * lam * H - H**2,
        lam * t * (s_coefficient * lam**2 + v * lam * H),
    )
    return {
        "curve": curve,
        "P0": p0,
        "P1": p1,
        "t": t,
        "c": c,
        "lambda": lam,
    }


def specialization_points_at_t3(record):
    curve = record["curve"]
    field = curve.base_field().base_ring()
    t3 = field(3)
    specialized = EllipticCurve(
        field,
        [coefficient(t3) for coefficient in curve.a_invariants()[3:]],
    )
    p0 = specialized(record["P0"][0](t3), record["P0"][1](t3))
    p1 = specialized(record["P1"][0](t3), record["P1"][1](t3))
    return specialized, p0, p1


def no_divisibility(record, primes):
    curve, p0, p1 = specialization_points_at_t3(record)
    combinations = {
        "P0": p0,
        "P1": p1,
        "P0+P1": p0 + p1,
        "P0-P1": p0 - p1,
    }
    result = {}
    for prime in primes:
        result[str(prime)] = {
            label: bool(point.is_divisible_by(prime))
            for label, point in combinations.items()
        }
        if any(result[str(prime)].values()):
            raise ArithmeticError(f"unexpected divisibility by {prime}")
    return {
        "curve_discriminant": str(curve.discriminant()),
        "points": {
            label: [str(point[0]), str(point[1])]
            for label, point in combinations.items()
        },
        "is_divisible": result,
    }


if not SOURCE_CERTIFICATE.exists():
    raise FileNotFoundError(SOURCE_CERTIFICATE)


# Generic Neron--Severi lattice and its saturation.
generic_ns = generic_ns_gram()
if generic_ns.det() != 36:
    raise ArithmeticError("generic NS determinant changed")
if QuadraticForm(ZZ, generic_ns).signature() != -17:
    raise ArithmeticError("generic NS signature changed")
if any(generic_ns[index, index] % 2 for index in range(19)):
    raise ArithmeticError("generic NS lattice is not even")
generic_orders, generic_form = discriminant_form(generic_ns)
if generic_orders != (ZZ(3), ZZ(12)):
    raise ArithmeticError("generic NS discriminant group changed")
expected_form = matrix(QQ, [[QQ(4) / 3, QQ(1) / 3], [QQ(1) / 3, -QQ(1) / 4]])
for i in range(2):
    for j in range(2):
        modulus = 2 if i == j else 1
        if not congruent(generic_form[i, j], expected_form[i, j], modulus):
            raise ArithmeticError("generic NS discriminant form changed")

isotropic_elements = []
for first in range(3):
    for second in range(12):
        if first == 0 and second == 0:
            continue
        if congruent(
            quadratic_value((first, second), expected_form), QQ(0), 2
        ):
            isotropic_elements.append((first, second))
if isotropic_elements != [(0, 4), (0, 8), (1, 4), (2, 8)]:
    raise ArithmeticError("generic NS overlattice census changed")

generic_k1 = specialized_family(QQ, QQ(1))
generic_saturation = no_divisibility(generic_k1, (3,))


# Generic transcendental lattice and uniqueness in its genus.
generic_T = matrix(ZZ, [[0, 3, 0], [3, 0, 0], [0, 0, 4]])
if generic_T.det() != -36 or QuadraticForm(ZZ, generic_T).signature() != 1:
    raise ArithmeticError("generic transcendental candidate changed")
isometric, witness = discriminant_forms_are_isometric(-generic_ns, generic_T)
if not isometric:
    raise ArithmeticError("generic NS and T discriminant forms do not match")
genus = QuadraticForm(ZZ, generic_T).global_genus_symbol()
genus_representatives = genus.representatives()
if len(genus_representatives) != 1:
    raise ArithmeticError("generic transcendental genus is no longer one-class")


# Split the old U and expose the positive frame for neighbor searches.
old_fibre = vector(ZZ, [0, 1] + [0] * 17)
old_zero = vector(ZZ, [1] + [0] * 18)
constraints = matrix(
    ZZ,
    [generic_ns * old_fibre, generic_ns * (old_zero + old_fibre)],
)
frame_basis = constraints.right_kernel().basis_matrix()
positive_frame = -(frame_basis * generic_ns * frame_basis.transpose())
if positive_frame.det() != 36 or not positive_frame.is_positive_definite():
    raise ArithmeticError("generic positive frame changed")
short_vectors = IntegralLattice(positive_frame).short_vectors(3)
roots = short_vectors[2]
root_basis = matrix(ZZ, roots).row_module().basis_matrix()
root_gram = root_basis * positive_frame * root_basis.transpose()
if len(roots) != 156 or root_basis.rank() != 15 or root_gram.det() != 36:
    raise ArithmeticError("generic frame root system changed")

# Complete nominal smallest degree-two isotropic layer.  In split coordinates
# NS=U+(-M), a class D=e+2f+w is isotropic exactly when w^2=4, but D.O=-1.
# Thus O is fixed and D-O has old-fibre degree one: these are section pencils
# written in a misleading degree-two presentation, not genuine 2-neighbours.
# The child-frame recurrence below is retained as a regression for that exact
# obstruction.
all_short = IntegralLattice(positive_frame).short_vectors(5)
positive_roots = [
    vector(ZZ, item)
    for item in roots
    if next(entry for entry in item if entry) > 0
]
positive_root_set = {tuple(item) for item in positive_roots}
simple_roots = [
    item
    for item in positive_roots
    if not any(tuple(item - left) in positive_root_set for left in positive_roots)
]
if len(simple_roots) != 15:
    raise ArithmeticError("failed to recover the current simple roots")


def dominant_representative(item):
    item = vector(ZZ, item)
    for _ in range(10000):
        for root in simple_roots:
            pairing = item * positive_frame * root
            if pairing < 0:
                item -= pairing * root
                break
        else:
            return tuple(int(entry) for entry in item)
    raise RuntimeError("Weyl reduction did not terminate")


degree_two_orbits = sorted(
    {dominant_representative(item) for item in all_short[4]}
)
if len(all_short[4]) != 17688 or len(degree_two_orbits) != 14:
    raise ArithmeticError("smallest degree-two Weyl census changed")
degree_two_records = []
for orbit_index, representative in enumerate(degree_two_orbits):
    # Use the anti-dominant sign so intersections with the current simple
    # fibre components are nonnegative.
    divisor = (
        old_fibre
        + 2 * (old_zero + old_fibre)
        - vector(ZZ, representative) * frame_basis
    )
    if divisor * generic_ns * divisor != 0 or gcd(list(divisor)) != 1:
        raise ArithmeticError("degree-two orbit did not give a primitive isotropic class")
    zero_pairing = divisor * generic_ns * old_zero
    reduced_divisor = divisor - old_zero
    if zero_pairing != -1:
        raise ArithmeticError("nominal degree-two zero-section obstruction changed")
    if (
        reduced_divisor * generic_ns * reduced_divisor != 0
        or reduced_divisor * generic_ns * old_fibre != 1
    ):
        raise ArithmeticError("zero-section removal did not give a degree-one pencil")
    child_frame = split_primitive_isotropic(generic_ns, divisor)
    child_roots = root_data(child_frame)
    if child_roots != (15, 156, 36):
        raise ArithmeticError("smallest degree-two layer changed root data")
    degree_two_records.append(
        {
            "orbit": orbit_index,
            "dominant_w": list(representative),
            "isotropic_divisor_in_ns_basis": [int(entry) for entry in divisor],
            "intersection_with_old_zero": int(zero_pairing),
            "fixed_component": "O with multiplicity one",
            "reduced_degree_one_divisor_in_ns_basis": [
                int(entry) for entry in reduced_divisor
            ],
            "child_frame_digest": matrix_digest(child_frame),
            "child_root_data": list(child_roots),
        }
    )


# Descend the unmarked equation to s=k^2 and compute the full unordered pair
# of squared residual-I1 locations z_i=t_i^2.
S_RING = PolynomialRing(QQ, "s")
s = S_RING.gen()
S_FIELD = S_RING.fraction_field()
s = S_FIELD(s)
D_s = 3 * s - 4
c_squared = 4 * s / D_s**2
lambda_s = -(s - 4) * D_s / 4
residual_a = 27 * c_squared**2
residual_b = 54 * c_squared + 4
residual_c = -(108 * c_squared + 9)
z_sum = 2 + lambda_s * residual_b / residual_c
z_product = (
    (2 * s + 1) * (9 * s - 20) ** 2 / (9 * (3 * s + 4) ** 2)
)
z_discriminant = z_sum**2 - 4 * z_product
expected_z_sum = (
    27 * s**4 - 54 * s**3 - 126 * s**2 + 656 * s + 544
) / (9 * (3 * s + 4) ** 2)
expected_z_discriminant = (
    (s - 4) ** 2 * (s**2 + QQ(4) * s / 3 + QQ(16) / 9) ** 3
    / (9 * (s + QQ(4) / 3) ** 4)
)
if z_sum != expected_z_sum or z_discriminant != expected_z_discriminant:
    raise ArithmeticError("residual-fibre configuration formulas changed")


# The two smooth-branch/I1 collisions.  They are singular K3 surfaces because
# the new I2 root raises the Shioda--Tate lower bound to 20.
Z_RING = PolynomialRing(QQ, "z")
z = Z_RING.gen()
sqrt5_field = NumberField(z**2 - 5, "sqrt5")
sqrt5 = sqrt5_field.gen()
special_20_9 = specialized_family(sqrt5_field, 2 * sqrt5 / 3)
height_20_9_p0 = compact_height(special_20_9["curve"], special_20_9["P0"], 12)
height_20_9_p1 = compact_height(special_20_9["curve"], special_20_9["P1"], 12)
if (height_20_9_p0["height"], height_20_9_p1["height"]) != ("1/3", "5/2"):
    raise ArithmeticError("s=20/9 heights changed")
ns_20_9 = branch_collision_ns(True)
if ns_20_9.det() != -60:
    raise ArithmeticError("s=20/9 NS determinant changed")
saturation_20_9 = no_divisibility(special_20_9, (2,))
T_20_9 = matrix(ZZ, [[4, 2], [2, 16]])
isometric_20_9, witness_20_9 = discriminant_forms_are_isometric(-ns_20_9, T_20_9)
if not isometric_20_9:
    raise ArithmeticError("s=20/9 transcendental lattice changed")

sqrt_minus2_field = NumberField(z**2 + 2, "sqrt_minus2")
sqrt_minus2 = sqrt_minus2_field.gen()
special_minus_half = specialized_family(sqrt_minus2_field, sqrt_minus2 / 2)
height_minus_half_p0 = compact_height(
    special_minus_half["curve"], special_minus_half["P0"], 12
)
height_minus_half_p1 = compact_height(
    special_minus_half["curve"], special_minus_half["P1"], 12
)
if (height_minus_half_p0["height"], height_minus_half_p1["height"]) != (
    "1/3",
    "3",
):
    raise ArithmeticError("s=-1/2 heights changed")
ns_minus_half = branch_collision_ns(False)
if ns_minus_half.det() != -72:
    raise ArithmeticError("s=-1/2 NS determinant changed")
saturation_minus_half = no_divisibility(special_minus_half, (2, 3))
T_minus_half = matrix(ZZ, [[6, 0], [0, 12]])
isometric_minus_half, witness_minus_half = discriminant_forms_are_isometric(
    -ns_minus_half, T_minus_half
)
if not isometric_minus_half:
    raise ArithmeticError("s=-1/2 transcendental lattice changed")


# At s=0 the invariant section becomes 2-torsion and IV* enhances to III*.
special_zero = specialized_family(QQ, QQ(0))
if 2 * special_zero["P0"] != special_zero["curve"](0):
    raise ArithmeticError("s=0 section is not two-torsion")
height_zero_p1 = compact_height(special_zero["curve"], special_zero["P1"], 4)
if height_zero_p1["height"] != "3":
    raise ArithmeticError("s=0 free-section height changed")
ns_zero, glue_zero = torsion_special_ns(2)
if ns_zero.det() != -12:
    raise ArithmeticError("s=0 saturated NS determinant changed")
_, zero_p0, zero_p1 = specialization_points_at_t3(special_zero)
if zero_p1.is_divisible_by(2):
    raise ArithmeticError("s=0 free section is not saturated")
T_zero = matrix(ZZ, [[4, 2], [2, 4]])
isometric_zero, witness_zero = discriminant_forms_are_isometric(-ns_zero, T_zero)
if not isometric_zero:
    raise ArithmeticError("s=0 transcendental lattice changed")


# At s=-4/3 the I2 fibre absorbs an I1, becoming I3 before base change and
# I6 after it; the invariant section becomes 3-torsion.
sqrt_minus3_field = NumberField(z**2 + 3, "sqrt_minus3")
sqrt_minus3 = sqrt_minus3_field.gen()
special_minus_four_thirds = specialized_family(
    sqrt_minus3_field, 2 * sqrt_minus3 / 3
)
if 3 * special_minus_four_thirds["P0"] != special_minus_four_thirds["curve"](0):
    raise ArithmeticError("s=-4/3 section is not three-torsion")
height_minus_four_thirds_p1 = compact_height(
    special_minus_four_thirds["curve"], special_minus_four_thirds["P1"], 6
)
if height_minus_four_thirds_p1["height"] != "5/2":
    raise ArithmeticError("s=-4/3 free-section height changed")
ns_minus_four_thirds, glue_minus_four_thirds = torsion_special_ns(3)
if ns_minus_four_thirds.det() != -15:
    raise ArithmeticError("s=-4/3 saturated NS determinant changed")
T_minus_four_thirds = matrix(ZZ, [[4, 1], [1, 4]])
isometric_minus_four_thirds, witness_minus_four_thirds = (
    discriminant_forms_are_isometric(
        -ns_minus_four_thirds, T_minus_four_thirds
    )
)
if not isometric_minus_four_thirds:
    raise ArithmeticError("s=-4/3 transcendental lattice changed")


payload = {
    "schema": "elkies-k3.e6a1-rho19-dissection.v1",
    "status": "PASS_EXACT_GENERIC_NS_T_AND_FOUR_SINGULAR_K3_BOUNDARIES",
    "inputs": {relative(SOURCE_CERTIFICATE): digest(SOURCE_CERTIFICATE)},
    "parameter_descent": {
        "coarse_parameter": "s=k^2",
        "equation_coefficients": {
            "c_squared": str(c_squared),
            "lambda": str(lambda_s),
        },
        "section_galois_action": {
            "k_to_minus_k": "P0 -> -P0, P1 -> P1",
            "generic_mw_rank_over_QQ(s)(t)": 1,
            "generic_mw_rank_over_QQ(k)(t)": 2,
        },
        "residual_squared_locations": {
            "polynomial": "z^2-S(s)*z+P(s)",
            "sum": str(z_sum),
            "product": str(z_product),
            "discriminant": str(z_discriminant),
        },
    },
    "generic_k3": {
        "fibre_profile": "2IV*+I4+4I1",
        "basis": [
            "O",
            "F",
            *[f"E6a_{i}" for i in range(1, 7)],
            *[f"E6b_{i}" for i in range(1, 7)],
            *[f"A3_{i}" for i in range(1, 4)],
            "P0",
            "P1",
        ],
        "section_component_profiles": {
            "P0": "nonzero 4/3 component in both IV*, middle component in I4",
            "P1": "identity component in both IV*, middle component in I4",
            "P0_dot_P1": 1,
        },
        "ns_gram": integer_rows(generic_ns),
        "ns_signature": [1, 18],
        "ns_determinant": 36,
        "ns_discriminant_group": [3, 12],
        "ns_discriminant_form": rational_rows(expected_form),
        "nonzero_isotropic_elements": [list(item) for item in isotropic_elements],
        "possible_proper_even_overlattice_indices": [3],
        "saturation_certificate": generic_saturation,
        "transcendental_lattice": {
            "isometry_type": "U(3)+<4>",
            "gram": integer_rows(generic_T),
            "signature": [2, 1],
            "determinant": -36,
            "discriminant_form_isometry_witness": [
                list(item) for item in witness
            ],
            "genus_symbol": str(genus),
            "genus_class_count": len(genus_representatives),
        },
        "positive_frame_for_neighbor_search": {
            "gram": integer_rows(positive_frame),
            "determinant": 36,
            "minimum": 2,
            "root_count": len(roots),
            "root_rank": int(root_basis.rank()),
            "root_lattice_determinant": int(root_gram.det()),
            "root_type": "2E6+A3",
            "warning": (
                "This is an exact frame input, not a proof that a rootless "
                "elliptic neighbor exists."
            ),
        },
        "deck_involution": {
            "map": "t -> -t",
            "holomorphic_two_form_character": -1,
            "mw_action": "P0 -> P0, P1 -> -P1",
            "quotient": "the E6+A1 rational elliptic surface",
            "fixed_locus": "one genus-one curve plus two rational curves",
            "two_elementary_invariants_partially_pinned": {"r": 12, "a": 8},
        },
    },
    "forced_singular_k3_specializations": [
        {
            "s": "0",
            "fibre_profile": "2III*+I4+2I1",
            "root_type": "2E7+A3",
            "torsion": "Z/2",
            "free_mw_rank": 1,
            "free_height": height_zero_p1,
            "ns_determinant": -12,
            "transcendental_gram": integer_rows(T_zero),
            "transcendental_discriminant": 12,
            "cm_discriminant": -12,
            "saturation": "Q is not 2-divisible at t=3",
            "discriminant_form_isometry_witness": [
                list(item) for item in witness_zero
            ],
        },
        {
            "s": "-4/3",
            "fibre_profile": "2IV*+I6+2I1",
            "root_type": "2E6+A5",
            "torsion": "Z/3",
            "free_mw_rank": 1,
            "free_height": height_minus_four_thirds_p1,
            "ns_determinant": -15,
            "transcendental_gram": integer_rows(T_minus_four_thirds),
            "transcendental_discriminant": 15,
            "cm_discriminant": -15,
            "saturation": "automatic because no nontrivial square divides 15",
            "discriminant_form_isometry_witness": [
                list(item) for item in witness_minus_four_thirds
            ],
        },
        {
            "s": "20/9",
            "fibre_profile": "2IV*+I4+I2+2I1",
            "root_type": "2E6+A3+A1",
            "torsion": "trivial",
            "mw_rank": 2,
            "height_gram": [
                [height_20_9_p0["height"], "0"],
                ["0", height_20_9_p1["height"]],
            ],
            "ns_determinant": -60,
            "transcendental_gram": integer_rows(T_20_9),
            "transcendental_discriminant": 60,
            "cm_discriminant": -60,
            "saturation_certificate": saturation_20_9,
            "discriminant_form_isometry_witness": [
                list(item) for item in witness_20_9
            ],
        },
        {
            "s": "-1/2",
            "fibre_profile": "2IV*+I4+I2+2I1",
            "root_type": "2E6+A3+A1",
            "torsion": "trivial",
            "mw_rank": 2,
            "height_gram": [
                [height_minus_half_p0["height"], "0"],
                ["0", height_minus_half_p1["height"]],
            ],
            "ns_determinant": -72,
            "transcendental_gram": integer_rows(T_minus_half),
            "transcendental_discriminant": 72,
            "cm_discriminant": -72,
            "saturation_certificate": saturation_minus_half,
            "discriminant_form_isometry_witness": [
                list(item) for item in witness_minus_half
            ],
        },
    ],
    "other_boundary_loci": {
        "degenerate_not_k3": [
            {"s": "4", "reason": "lambda=0 and the Weierstrass cubic is singular"},
            {"s": "4/3", "reason": "the c,lambda chart has a pole"},
        ],
        "residual_I1_pair_to_II": {
            "equation": "9*s^2+12*s+16=0",
            "fibre_profile": "2IV*+I4+2II",
            "root_rank_change": 0,
            "picard_status": "not forced above 19 by Shioda--Tate",
        },
        "compactified_points_not_resolved_here": ["s=4/3", "s=infinity"],
    },
    "neighbor_search_boundary": {
        "proved": (
            "Every Jacobian elliptic fibration with section splits U from the "
            "same determinant-36 NS lattice; its positive frame has rank 17, "
            "determinant 36, and the pinned discriminant form."
        ),
        "not_proved": (
            "No rootless or low-root neighbor is asserted. Nefness, complete "
            "root enumeration in each child frame, and equation compilation "
            "remain separate gates."
        ),
        "complete_smallest_degree_two_layer": {
            "shape": "D=e+2*f+w with w^2=4 in U+(-M)",
            "norm_four_vector_count": len(all_short[4]),
            "weyl_orbit_count": len(degree_two_orbits),
            "orbit_records": degree_two_records,
            "result": (
                "Every class has D.O=-1 and reduces by the fixed zero section "
                "to old-fibre degree one. The repeated child root data "
                "(15,156,36) therefore describe section presentations, not "
                "genuine degree-two neighbors."
            ),
            "scope": (
                "Complete for this isotropic shape modulo the current "
                "2E6+A3 Weyl group; the physical zero-section obstruction is exact."
            ),
        },
        "nominal_degree_three_obstruction": {
            "shape": "D=e+3*f+w with w^2=6 in U+(-M)",
            "intersection_with_old_zero": -2,
            "fixed_component": "O with multiplicity two",
            "degree_after_removal": 1,
            "result": "The norm-six layer is likewise not a genuine degree-three neighbor.",
        },
        "first_genuine_shapes": {
            "degree_two": "D=2*e+2*f+w with w^2=8 and D.O=0",
            "degree_three": "D=3*e+3*f+w with w^2=18 and D.O=0",
        },
    },
}


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
output_path = arguments.output.resolve()
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != encoded:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(encoded)

print(
    "E6A1RHO19|NS=36|disc=Z3+Z12|T=U3+4|"
    "CM=-12,-15,-60,-72|status=PASS",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
