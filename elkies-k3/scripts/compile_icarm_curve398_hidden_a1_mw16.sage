#!/usr/bin/env sage-python
"""Recover curve 398 from an exact norm-eight A1/MW16 neighbour.

The source is the committed rootless alternate-Q80 equation.  Priority trace
16875 defines ``D=(2,2,w)=O+P_w``.  This program rebuilds its binary-quartic
Jacobian over QQ, proves the fibre configuration ``I2+22I1``, factors the
curve-398 j-equation, and checks the unique rational parameter is a Q-isomorphic
copy of curve 398.

It then enumerates the complete degree-one old-section shell, chooses one old
section as zero and a saturated MW16 basis, specializes those sixteen sections,
and identifies them by exact integer coordinates in the displayed public
rank-30 subgroup.  A separate redacted fixture contains only the short curve,
the sixteen generic points, and their exact generic height Gram; it is the only
input permitted to the subsequent blind discovery experiment.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from math import ceil, sqrt
from pathlib import Path
import sys

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, block_diagonal_matrix, matrix, pari, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
TABLE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
TARGET = ROOT / "elliptic-curves/cas/icarm_curve398.py"
CHORD = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
SCREEN = ROOT / "elkies-k3/scripts/screen_icarm_curve398_norm8_a1_fibrations.sage"
TRUTH_OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_curve398_hidden_a1_mw16_v1.json"
BLIND_OUTPUT = ROOT / "elliptic-curves/data/icarm_curve398_mw16_blind_input_v1.json"
PRIORITY_RANK = 16875
EXPECTED_PARAMETER = QQ(
    -273478312517509127154149830485048828022673347107308547939067553994727903425458545978043182638015899676311550557441827100822466901248
) / QQ(
    243076210150914055804756105904064536659703543720469425499709810733677965174759784940636972086422417178984090368085211
)

# The first row is the new zero.  The remaining rows are a saturated MW16
# basis selected from the complete 166-element old-degree-one shell.
SECTION_VECTORS = (
    (0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
    (0, -1, 0, 0, 0, 0, 0, 0, -1, 1, 0, 1, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, -1, 0, 1, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, 0),
    (0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 1, 1, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 1, 0, 1, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, -1, 0, -1, 0, 1, 0, 1, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, 1),
    (0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, -1, 0, 0, 0, 1),
    (0, 0, 1, -1, 0, -1, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
    (0, 1, 0, -1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, -1, 0, 0, 0, 1, 0, 0, 1, -1, 0, 0, 0, -1, 1, 0),
    (0, 0, 0, 0, 0, 0, 1, 0, -1, 1, 0, 0, -1, -1, 0, 0, 1),
    (-1, 1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 0, -1, 0, 0, 0, 1),
)


def load(name: str, path: Path):
    return SourceFileLoader(name, str(path)).load_module()


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def qtext(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def poly_record(poly):
    return [qtext(poly[index]) for index in range(poly.degree() + 1)] if poly else ["0"]


def point_record(point):
    return {"x": qtext(point[0]), "y": qtext(point[1])}


def rational_function_record(value):
    return {
        "numerator_coefficients_low_to_high": poly_record(value.numerator()),
        "denominator_coefficients_low_to_high": poly_record(value.denominator()),
    }


def section_class(v, height_gram):
    v = vector(ZZ, v)
    return vector(QQ, [(v * height_gram * v - 2) / 2, 1] + list(v))


def ns_intersection(left, right, height_gram):
    return left[0] * right[1] + left[1] * right[0] - vector(QQ, left[2:]) * height_gram * vector(QQ, right[2:])


def enumerate_degree_one_vectors(height_gram, trace_vector):
    """Enumerate z=trace-2v of norm 12 in the fixed parity coset."""

    dimension = height_gram.nrows()
    change = height_gram.LLL_gram().transpose()
    reduced = change * height_gram * change.transpose()
    parity_qq = trace_vector * change.inverse()
    if any(value not in ZZ for value in parity_qq):
        raise ArithmeticError("LLL parity transport is not integral")
    parity = tuple(int(value) % 2 for value in parity_qq)
    lower = matrix(QQ, dimension, dimension)
    diagonal = []
    for index in range(dimension):
        lower[index, index] = 1
        value = QQ(reduced[index, index]) - sum(
            lower[index, prior] ** 2 * diagonal[prior] for prior in range(index)
        )
        diagonal.append(value)
        for row in range(index + 1, dimension):
            lower[row, index] = (
                QQ(reduced[row, index])
                - sum(lower[row, prior] * lower[index, prior] * diagonal[prior] for prior in range(index))
            ) / value
    lower_float = [[float(lower[row, column]) for column in range(dimension)] for row in range(dimension)]
    diagonal_float = [float(value) for value in diagonal]
    coordinates = [0] * dimension
    answer = []
    floating_bound = 12.25

    def visit(index, used, exact_partial):
        if index < 0:
            if exact_partial == 12:
                z = vector(ZZ, coordinates) * change
                v = (trace_vector - z) / 2
                if any(value not in ZZ for value in v):
                    raise ArithmeticError("fixed-parity enumeration produced a half-integral section")
                answer.append(vector(ZZ, v))
            return
        center = sum(lower_float[row][index] * coordinates[row] for row in range(index + 1, dimension))
        radius = sqrt(max(0.0, (floating_bound - used) / diagonal_float[index]))
        lower_bound = ceil(-center - radius - 1.0e-10)
        upper_bound = int((-center + radius + 1.0e-10) // 1)
        start = lower_bound + ((parity[index] - lower_bound) % 2)
        cross = sum(reduced[index, row] * coordinates[row] for row in range(index + 1, dimension))
        for entry in range(start, upper_bound + 1, 2):
            cost = diagonal_float[index] * (entry + center) ** 2
            if used + cost > floating_bound + 1.0e-10:
                continue
            coordinates[index] = entry
            visit(
                index - 1,
                used + cost,
                exact_partial + reduced[index, index] * entry**2 + 2 * entry * cross,
            )

    visit(dimension - 1, 0.0, ZZ(0))
    if any(
        v * height_gram * v - trace_vector * height_gram * v != 1
        for v in answer
    ):
        raise ArithmeticError("degree-one coset enumeration failed exact replay")
    return tuple(answer)


def shioda_gram(vectors, zero, trace_vector, height_gram):
    zero_class = section_class(zero, height_gram)
    old_zero = vector(QQ, [-1, 1] + [0] * 17)
    if zero * height_gram * zero != 4:
        raise ArithmeticError("the chosen new zero does not meet the identity component")

    def data(v):
        source_class = section_class(v, height_gram)
        return source_class, ns_intersection(source_class, zero_class, height_gram), int(ns_intersection(source_class, old_zero, height_gram))

    records = [data(vector(ZZ, v)) for v in vectors]
    result = matrix(QQ, len(vectors), len(vectors))
    for left in range(len(vectors)):
        for right in range(len(vectors)):
            left_class, left_zero, left_component = records[left]
            right_class, right_zero, right_component = records[right]
            if left == right:
                value = 4 + 2 * left_zero - QQ(left_component) / 2
            else:
                value = (
                    2
                    + left_zero
                    + right_zero
                    - ns_intersection(left_class, right_class, height_gram)
                    - QQ(left_component * right_component) / 2
                )
            result[left, right] = value
    return result


def evaluate_rational(function, value):
    return function.numerator()(value) / function.denominator()(value)


def invert_mobius(function, value, ring):
    numerator, denominator = map(ring, (function.numerator(), function.denominator()))
    if numerator.degree() > 1 or denominator.degree() > 1:
        raise ArithmeticError("an old degree-one section did not induce a Mobius base map")
    n0, n1 = numerator[0], numerator[1]
    d0, d1 = denominator[0], denominator[1]
    bottom = value * d1 - n1
    if not bottom:
        raise ArithmeticError("selected specialization meets this section at quartic infinity")
    answer = (n0 - value * d0) / bottom
    if evaluate_rational(function, answer) != value:
        raise ArithmeticError("Mobius inversion failed")
    return QQ(answer)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth-output", type=Path, default=TRUTH_OUTPUT)
    parser.add_argument("--blind-output", type=Path, default=BLIND_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    sys.path.insert(0, str(ROOT / "elliptic-curves"))
    from latent_lattice.elliptic import EllipticCurve as LatentEllipticCurve
    from latent_lattice.pari import recover_exact_embedding

    screen = load("curve398_hidden_screen", SCREEN)
    chord = load("curve398_hidden_chord", CHORD)
    target = load("curve398_hidden_target", TARGET)
    model = json.loads(MODEL.read_text())
    table = screen.load_rows(TABLE)
    trace_row = table[PRIORITY_RANK - 1]
    trace_vector = vector(ZZ, screen.parse_vector(trace_row["section_basis_w"]))
    height_gram = matrix(ZZ, model["sections"]["height_gram"])
    if trace_vector * height_gram * trace_vector != 8:
        raise ArithmeticError("selected trace lost norm eight")

    old_ring = PolynomialRing(QQ, "t")
    old_field = old_ring.fraction_field()
    old_a = old_ring(model["weierstrass_model"]["A_coefficients_low_to_high"])
    old_b = old_ring(model["weierstrass_model"]["B_coefficients_low_to_high"])
    old_curve = EllipticCurve(old_field, [old_a, old_b])
    old_basis = tuple(
        old_curve(
            screen.polynomial_from_record(record["X"], old_ring, QQ),
            screen.polynomial_from_record(record["Y"], old_ring, QQ),
        )
        for record in model["sections"]["records"]
    )
    trace = sum(
        (coefficient * point for coefficient, point in zip(trace_vector, old_basis) if coefficient),
        old_curve(0),
    )
    frame = chord.trace_chord_frame(trace[0], trace[1], old_ring)
    h, nx, ny, m0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
    if h.degree() != 2:
        raise ArithmeticError("selected trace is not in the finite-pole chart")

    parameter_ring = PolynomialRing(QQ, "lambda")
    parameter = parameter_ring.gen()
    bivariate_ring = PolynomialRing(parameter_ring, "t")
    hh, nnx, nny, mm0 = map(bivariate_ring, (h, nx, ny, m0))
    slope_numerator = mm0 + parameter * hh**2
    numerator = (
        slope_numerator**4
        - 6 * slope_numerator**2 * nnx
        - 8 * slope_numerator * nny
        - 3 * nnx**2
        - 4 * bivariate_ring(old_a) * hh**4
    )
    quartic, remainder = numerator.quo_rem(hh**6)
    if remainder or quartic.degree() != 4:
        raise ArithmeticError("residual chord did not produce a binary quartic")
    invariant_i, invariant_j = screen.binary_quartic_invariants(quartic, parameter_ring)
    child_a, child_b = -27 * invariant_i, -27 * invariant_j
    child_delta = parameter_ring(-16 * (4 * child_a**3 + 27 * child_b**2))
    if [child_a.degree(), child_b.degree(), child_delta.degree()] != [8, 12, 22]:
        raise ArithmeticError("A1 child lost its K3 degree profile")
    if child_delta.gcd(child_delta.derivative()).degree() != 0:
        raise ArithmeticError("finite child discriminant is not squarefree")

    a1, a2, a3, a4, a6 = tuple(QQ(str(value)) for value in target.GENERAL_WEIERSTRASS_COEFFICIENTS)
    b2 = a1**2 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3**2 + 4 * a6
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    target_a, target_b = -27 * c4, -54 * c6
    comparison = child_a**3 * target_b**2 - target_a**3 * child_b**2
    factorization = comparison.factor()
    linear = [factor for factor, multiplicity in factorization if factor.degree() == 1 for unused in range(multiplicity)]
    if len(linear) != 1:
        raise ArithmeticError("selected trace does not have one rational curve-398 parameter")
    recovered_parameter = -linear[0][0] / linear[0][1]
    if recovered_parameter != EXPECTED_PARAMETER:
        raise ArithmeticError("curve-398 parameter changed")
    specialized_a = QQ(child_a(recovered_parameter))
    specialized_b = QQ(child_b(recovered_parameter))
    child_curve = EllipticCurve(QQ, [specialized_a, specialized_b])
    target_short = EllipticCurve(QQ, [target_a, target_b])
    if not child_curve.is_isomorphic(target_short):
        raise ArithmeticError("j-match is a nontrivial quadratic twist")
    child_to_target = child_curve.isomorphism_to(target_short)

    # The lattice certificate: D=O+P_w, a complete degree-one shell, and a
    # determinant-474 MW16 basis (|disc NS|/disc(A1)=948/2).
    fibre = vector(ZZ, [2, 2] + list(trace_vector))
    old_zero_class = vector(ZZ, [-1, 1] + [0] * 17)
    trace_class = section_class(trace_vector, height_gram)
    ns = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -height_gram)
    if fibre != old_zero_class + trace_class or fibre * ns * fibre or fibre * ns * old_zero_class:
        raise ArithmeticError("D=O+P trace decomposition failed")
    degree_one = enumerate_degree_one_vectors(height_gram, trace_vector)
    if len(degree_one) != 166:
        raise ArithmeticError(f"degree-one old-section count changed: {len(degree_one)}")
    selected = tuple(vector(ZZ, row) for row in SECTION_VECTORS)
    if any(value not in degree_one for value in selected):
        raise ArithmeticError("a frozen generic section left the complete degree-one shell")
    generic_gram = shioda_gram(selected[1:], selected[0], trace_vector, height_gram)
    if generic_gram.det() != 474 or generic_gram.rank() != 16:
        raise ArithmeticError("selected MW16 basis is not saturated")

    new_zero_class = section_class(selected[0], height_gram)
    mate = fibre + new_zero_class
    complement = matrix(ZZ, [list(fibre * ns), list(mate * ns)]).right_kernel_matrix()
    transport = matrix(QQ, [list(fibre), list(mate)] + [list(row) for row in complement])
    child_frame = -(complement * ns * complement.transpose())
    root_count = int(pari(matrix(ZZ, child_frame)).qfminim(2)[0])
    if abs(transport.det()) != 1 or child_frame.det() != 948 or root_count != 2:
        raise ArithmeticError("new frame is not primitive with one A1 root")

    # Specialize the chosen old curves to points on the binary quartic.
    fixed_m = m0 + recovered_parameter * h**2
    fixed_quartic = old_ring([QQ(quartic[index](recovered_parameter)) for index in range(5)])
    sum_x = old_ring((fixed_m**2 - nx) // h**2)
    quartic_points = []
    base_maps = []
    for section_vector in selected:
        source_point = sum(
            (coefficient * point for coefficient, point in zip(section_vector, old_basis) if coefficient),
            old_curve(0),
        )
        source_x, source_y = source_point[0], source_point[1]
        base_map = old_field(
            (((source_y + trace[1]) / (source_x - trace[0])) * h - m0) / h**2
        )
        old_parameter = invert_mobius(base_map, recovered_parameter, old_ring)
        x_value = QQ(source_x(old_parameter))
        y_value = QQ(source_y(old_parameter))
        w_value = (2 * x_value - QQ(sum_x(old_parameter))) / QQ(h(old_parameter))
        if w_value**2 != fixed_quartic(old_parameter):
            raise ArithmeticError("old section missed the specialized quartic")
        quartic_points.append((old_parameter, w_value))
        base_maps.append(base_map)

    # Point the quartic at the selected zero section and use the exact
    # pointed-quartic birational map to the raw invariant model.
    t0, w0 = quartic_points[0]
    shift_ring = PolynomialRing(QQ, "z")
    z = shift_ring.gen()
    shifted = shift_ring(fixed_quartic(t0 + z))
    ee, dd, cc, bb, aa = [QQ(shifted[index]) for index in range(5)]
    if ee != w0**2:
        raise ArithmeticError("pointed quartic constant term changed")
    a1g = dd / w0
    a2g = cc - dd**2 / (4 * w0**2)
    a3g = 2 * w0 * bb
    a4g = -4 * w0**2 * aa
    a6g = a2g * a4g
    b2g = a1g**2 + 4 * a2g
    b4g = a1g * a3g + 2 * a4g
    b6g = a3g**2 + 4 * a6g
    c4g = b2g**2 - 24 * b4g
    c6g = -b2g**3 + 36 * b2g * b4g - 216 * b6g
    if 81 * (-c4g / 48) != specialized_a or 729 * (-c6g / 864) != specialized_b:
        raise ArithmeticError("pointed quartic normalization missed the invariant model")

    public_curve = EllipticCurve(QQ, [a1, a2, a3, a4, a6])
    target_to_public = target_short.isomorphism_to(public_curve)
    public_to_target = public_curve.isomorphism_to(target_short)
    generic_public_points = []
    generic_child_points = []
    for old_parameter, w_value in quartic_points[1:]:
        zz = old_parameter - t0
        if not zz:
            raise ArithmeticError("two selected old sections meet at the quartic origin")
        x_general = (2 * w0 * (w_value + w0) + dd * zz) / zz**2
        y_general = (
            4 * w0**2 * (w_value + w0)
            + 2 * w0 * dd * zz
            + (2 * w0 * cc - dd**2 / (2 * w0)) * zz**2
        ) / zz**3
        child_point = child_curve(
            9 * (x_general + b2g / 12),
            27 * (y_general + (a1g * x_general + a3g) / 2),
        )
        public_point = target_to_public(child_to_target(child_point))
        generic_child_points.append(child_point)
        generic_public_points.append(public_point)

    latent_curve = LatentEllipticCurve(tuple(Fraction(str(value)) for value in (a1, a2, a3, a4, a6)))
    public_basis = tuple((Fraction(str(x)), Fraction(str(y))) for x, y in target.POINTS)
    generic_public_affine = tuple(
        (Fraction(str(point[0])), Fraction(str(point[1]))) for point in generic_public_points
    )
    embedding = recover_exact_embedding(
        latent_curve,
        public_basis,
        generic_public_affine,
        digits=150,
        timeout=300.0,
    )
    if len(embedding.columns) != 16:
        raise ArithmeticError("public subgroup embedding has the wrong rank")

    generic_records = []
    for index, (section_vector, base_map, quartic_point, public_point, column) in enumerate(
        zip(selected[1:], base_maps[1:], quartic_points[1:], generic_public_points, embedding.columns),
        start=1,
    ):
        source_class = section_class(section_vector, height_gram)
        generic_records.append(
            {
                "basis_index_one_based": index,
                "source_section_basis_coordinates": list(map(int, section_vector)),
                "source_section_height": int(section_vector * height_gram * section_vector),
                "source_section_intersection_with_D": int(ns_intersection(source_class, fibre, height_gram)),
                "new_fibre_component": "nonidentity-O" if source_class * ns * old_zero_class else "identity-P_w",
                "base_map_lambda_of_t": rational_function_record(base_map),
                "specialized_quartic_point": {"t": qtext(quartic_point[0]), "W": qtext(quartic_point[1])},
                "specialized_public_point": point_record(public_point),
                "coordinates_in_public_rank30_points": list(column),
                "exact_public_group_law_replay": True,
            }
        )

    truth = {
        "schema": "elliptic-curves.icarm-curve398-hidden-a1-mw16.v1",
        "status": "PASS_EXACT_HIDDEN_A1_MW16_FIBRATION_PARAMETER_AND_PUBLIC_SUBGROUP",
        "curve_id": 398,
        "source_chart": model["divisor"]["label"],
        "fibration": {
            "priority_rank": PRIORITY_RANK,
            "orbit_mask": int(trace_row["orbit_mask"]),
            "orbit_hex": trace_row["orbit_hex"],
            "trace_section_basis_w": list(map(int, trace_vector)),
            "divisor_class_in_U_plus_M_minus": list(map(int, fibre)),
            "divisor_identity": "D=(2,2,w)=O+P_w",
            "old_fibre_degree": 2,
            "old_zero_degree": 0,
            "finite_pole_degree": int(h.degree()),
            "equation": "Y^2=X^3+A(lambda)*X+B(lambda)",
            "A_coefficients_low_to_high": poly_record(child_a),
            "B_coefficients_low_to_high": poly_record(child_b),
            "degrees_A_B_Delta": [8, 12, 22],
            "infinity_orders_c4_c6_Delta": [0, 0, 2],
            "fibre_configuration": "I2 at infinity + 22 I1",
            "finite_discriminant_squarefree": True,
            "child_frame_determinant": int(child_frame.det()),
            "child_frame_norm_two_vector_count_signed": root_count,
        },
        "parameter_recovery": {
            "comparison_polynomial_degree": int(comparison.degree()),
            "factor_degrees_with_multiplicity": [
                [int(factor.degree()), int(multiplicity)] for factor, multiplicity in factorization
            ],
            "lambda": qtext(recovered_parameter),
            "specialized_child_short_coefficients": [qtext(specialized_a), qtext(specialized_b)],
            "isomorphic_to_curve398_over_Q": True,
            "child_to_curve398_short_isomorphism_u_r_s_t": [qtext(value) for value in child_to_target.tuple()],
        },
        "generic_mw16": {
            "complete_old_degree_one_section_count": len(degree_one),
            "zero_source_section_basis_coordinates": list(map(int, selected[0])),
            "height_gram": [[qtext(value) for value in row] for row in generic_gram.rows()],
            "height_gram_determinant": qtext(generic_gram.det()),
            "rank": 16,
            "saturated": True,
            "records": generic_records,
        },
        "public_rank30_embedding": {
            "orientation": "columns give each generic MW16 point in the ordered public 30-point list",
            "matrix_30_by_16_columns": [list(column) for column in embedding.columns],
            "maximum_absolute_coordinate": embedding.max_abs_coordinate,
            "nonzero_coordinate_count": embedding.nonzero_coordinates,
            "height_dual_numerical_residual_max": embedding.numerical_residual_max,
            "exact_group_law_replay": True,
        },
        "inputs": {relative(path): digest(path) for path in (MODEL, TABLE, TARGET, CHORD, SCREEN)},
        "software": {"sage_version": SAGE_VERSION, "pari_version": ".".join(map(str, pari.version()))},
        "proof_boundary": (
            "This certifies one exact A1/MW16 fibration, its curve-398 rational parameter, "
            "a saturated generic rank-16 section basis, and its exact embedding in the public "
            "rank-30 subgroup. It does not claim rank exactly 30 or saturation of the public subgroup."
        ),
        "reproducing_command": "sage -python elkies-k3/scripts/compile_icarm_curve398_hidden_a1_mw16.sage --check",
    }

    # Redacted input: no public complement points, public-coordinate columns,
    # or truth-artifact path are included.
    blind = {
        "schema": "elliptic-curves.icarm-curve398-mw16-blind-input.v1",
        "status": "PASS_REDACTED_GENERIC_MW16_INPUT",
        "curve_label": "curve398-calibration",
        "short_model": ["0", "0", "0", qtext(target_a), qtext(target_b)],
        "generic_points": [point_record(public_to_target(point)) for point in generic_public_points],
        "generic_height_gram": [[qtext(value) for value in row] for row in generic_gram.rows()],
        "generic_rank": 16,
        "redaction": {
            "public_rank30_fixture_loaded_by_search": False,
            "held_out_point_count": 14,
            "contains_public_embedding_coordinates": False,
            "permitted_information": "curve equation, sixteen specialized generic sections, exact generic MW16 height Gram",
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "source_input_hashes": {relative(path): digest(path) for path in (MODEL, TABLE, TARGET)},
        },
    }

    truth_text = json.dumps(truth, indent=2, sort_keys=True) + "\n"
    blind_text = json.dumps(blind, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.truth_output.read_text() != truth_text or args.blind_output.read_text() != blind_text:
            raise ArithmeticError("stored curve-398 A1/MW16 calibration artifacts differ from replay")
    else:
        args.truth_output.parent.mkdir(parents=True, exist_ok=True)
        args.blind_output.parent.mkdir(parents=True, exist_ok=True)
        args.truth_output.write_text(truth_text)
        args.blind_output.write_text(blind_text)
    print(
        f"CURVE398A1MW16|lambda={qtext(recovered_parameter)}|degree1={len(degree_one)}|"
        f"mw_rank=16|mw_det={qtext(generic_gram.det())}|public_embedding=PASS|"
        f"truth={relative(args.truth_output)}|blind={relative(args.blind_output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
