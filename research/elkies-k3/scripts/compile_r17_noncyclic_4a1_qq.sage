#!/usr/bin/env sage-python
"""Compile the new non-cyclic 4A1/MW13 two-neighbour on published R17.

status: ACTIVE_COMPILER
claim: exact 4I2+16I1 equation, saturated arithmetic MW13, and reverse R17 hop
inputs: published R17 model/sections and pinned relative-U/bridge/planner artifacts
outputs: artifacts/generated-results/elkies-k3-r17-noncyclic-4a1-direct-fibration-v1.json

The marked fibre and zero come from the exact relative-U/local-bridge
certificate.  The universal old-degree-two chord construction is then run
directly on Elkies's published R17 equation.  The resulting pointed quartic
is based at the certified physical zero, not at the old zero.  Finally the
script transports a saturated thirteen-section basis and verifies the
birational reverse hop to the literal published R17 equation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    EllipticCurve,
    PolynomialRing,
    QQ,
    QuadraticForm,
    ZZ,
    block_diagonal_matrix,
    matrix,
    pari,
    vector,
)
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
PINNED_GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
SHORT_GRAM = ROOT / "elkies-k3/data/lattice/short_vector_basis_gram.txt"
SHORT_COORDS = ROOT / "elkies-k3/data/lattice/short_vector_basis_coords.txt"
PUBLISHED = ROOT / "artifacts/generated-results/elkies-2026-published-r17-target.json"
LOCAL_MUTATION = ROOT / "artifacts/generated-results/elkies-k3-r17-local-bridge-mutation-v1.json"
PRIME_LOCAL = ROOT / "artifacts/generated-results/elkies-k3-prime-local-bridge-mutation-v1.json"
PLANNER = ROOT / "artifacts/generated-results/elkies-k3-marked-u-realization-planner-controls-v1.json"
HISTORICAL = ROOT / "artifacts/generated-results/elkies-k3-h3-pinned-r17-current-suffix-marking.json"
PHYSICAL_Q8 = ROOT / "artifacts/generated-results/elkies-k3-h3-q4o164-c8-q8o376-4a1-old_zero-frame.txt"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-noncyclic-4a1-direct-fibration-v1.json"

HYPERBOLIC = matrix(ZZ, [[0, 1], [1, 0]])

# These thirteen old sections were selected from the exact D-degree-one shell.
# Their images, together with the four primitive A1 roots, form a unimodular
# basis of the new frame.  Keeping this short list literal makes the replay a
# narrow exact certificate rather than a new ambient section search.
MW13_SOURCE_SECTIONS_SHORT = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [-1, 0, 0, 0, -1, 1, 2, 0, 0, 0, 0, 1, 0, -1, 1, 1, 1],
    [-1, 1, 1, 2, 1, -1, -3, 1, -1, 1, 0, -3, -1, 1, -1, 1, -1],
    [0, -1, 0, 0, -1, 2, 3, -1, 0, -1, 1, 1, 0, -2, 1, 1, 1],
    [-2, 1, 0, 1, -1, 1, 1, 0, 1, 1, 0, 0, -1, 0, 0, 1, 1],
    [-1, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1],
    [0, -1, 0, 0, 0, 2, 2, -1, -1, -1, 1, 1, 0, -2, 0, 1, 1],
    [0, -1, 0, 0, 0, 1, 1, 0, -1, 0, 1, 0, 0, -1, 0, 1, 0],
    [0, 0, 1, 1, 2, -1, -3, 1, -2, 0, 0, -3, 0, 1, -1, 0, -1],
    [-1, 0, -1, -1, -1, 0, 2, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1],
    [-1, 1, 1, 2, 2, -3, -5, 2, -1, 2, -1, -4, -1, 3, -2, 0, -2],
    [0, 0, 0, 0, 0, 1, 1, -1, -1, -1, 1, 0, 0, -2, 0, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def load_matrix(path: Path):
    return matrix(
        ZZ,
        [
            [ZZ(entry) for entry in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def polynomial_text(poly) -> list[str]:
    if not poly:
        return ["0"]
    return [rational_text(poly[index]) for index in range(poly.degree() + 1)]


def polynomial_content(poly):
    """Return the signed rational content of a nonzero QQ polynomial."""

    denominator = ZZ(1)
    for coefficient in poly:
        denominator = denominator.lcm(QQ(coefficient).denominator())
    numerator = ZZ(0)
    for coefficient in poly:
        numerator = numerator.gcd(ZZ(QQ(coefficient) * denominator))
    return QQ(abs(numerator) / denominator)


def rational_function_record(value):
    value = value.parent()(value)
    return {
        "numerator_coefficients_low_to_high": polynomial_text(value.numerator()),
        "denominator_coefficients_low_to_high": polynomial_text(value.denominator()),
    }


def polynomial_over_function_field_record(poly):
    return [rational_function_record(poly[index]) for index in range(poly.degree() + 1)]


def reconstruct_basis(ring, A, B, section_data):
    points = []
    for expected_index, record in enumerate(section_data["sections"]):
        assert int(record["basis_index"]) == expected_index
        x_coordinate = ring([QQ(value) for value in record["x_coefficients_low_to_high"]])
        if expected_index == 0:
            y_coordinate = ring([QQ(value) for value in record["y_coefficients_low_to_high"]])
        else:
            reference_x, reference_y = points[int(record["chord"]["reference_basis_index"])]
            slope = ring([QQ(value) for value in record["chord"]["slope_coefficients_low_to_high"]])
            y_coordinate = reference_y + slope * (x_coordinate - reference_x)
        assert y_coordinate**2 == x_coordinate**3 + A * x_coordinate + B
        points.append((x_coordinate, y_coordinate))
    return points


def evaluate_polynomial(poly, value):
    result = value.parent()(0)
    for coefficient in reversed(list(poly)):
        result = result * value + value.parent()(coefficient)
    return result


def evaluate_rational(function, value):
    return evaluate_polynomial(function.numerator(), value) / evaluate_polynomial(
        function.denominator(), value
    )


def invert_mobius(function, new_variable):
    numerator = function.numerator()
    denominator = function.denominator()
    if numerator.degree() > 1 or denominator.degree() > 1:
        raise ArithmeticError(f"expected a degree-one base map, obtained {function}")
    n0, n1 = numerator[0], numerator[1]
    d0, d1 = denominator[0], denominator[1]
    inverse = (n0 - new_variable * d0) / (new_variable * d1 - n1)
    if evaluate_rational(function, inverse) != new_variable:
        raise ArithmeticError("failed to invert a section Mobius base map")
    return inverse


def row_isometry(source, target):
    witness = QuadraticForm(ZZ, source).is_globally_equivalent_to(
        QuadraticForm(ZZ, target), return_matrix=True
    )
    if witness is False:
        return None
    witness = matrix(ZZ, witness)
    candidates = [witness, witness.transpose()]
    inverse = witness.inverse()
    if inverse.change_ring(ZZ) == inverse:
        inverse = inverse.change_ring(ZZ)
        candidates.extend([inverse, inverse.transpose()])
    for candidate in candidates:
        if candidate * source * candidate.transpose() == target:
            assert abs(candidate.det()) == 1
            return candidate
    raise ArithmeticError("unrecognized integral-isometry orientation")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    model = json.loads(MODEL.read_text())
    section_data = json.loads(SECTIONS.read_text())
    published = json.loads(PUBLISHED.read_text())
    local = json.loads(LOCAL_MUTATION.read_text())
    prime_local = json.loads(PRIME_LOCAL.read_text())
    planner = json.loads(PLANNER.read_text())
    historical = json.loads(HISTORICAL.read_text())

    assert local["status"] == "PASS_EXACT_R17_LOCAL_BRIDGE_MUTATION"
    assert prime_local["status"] == "PASS_EXACT_PRIME_LOCAL_BRIDGE_MUTATION_CLASSIFICATION"
    assert planner["status"] == "PASS_MARKED_U_PLANNER_CONTROLS_AND_R17_END_TO_END"
    forward_plan = planner["ordered_controls"][1]["result"]
    reverse_plan = planner["ordered_controls"][2]["result"]
    assert forward_plan["selected_candidate_id"] == "published-R17-degree2-new-noncyclic-4A1"
    assert reverse_plan["selected_candidate_id"] == "new-4A1-cheap-return-to-published-rootless"
    assert reverse_plan["target_request"]["target_frame_supplied"] is False

    pinned = load_matrix(PINNED_GRAM)
    short = load_matrix(SHORT_GRAM)
    short_coordinates = load_matrix(SHORT_COORDS)
    assert short_coordinates * pinned * short_coordinates.transpose() == short
    basis_change = matrix(ZZ, published["pinned_identification"]["basis_change_matrix"])

    Rt = PolynomialRing(QQ, "t")
    Kt = Rt.fraction_field()
    t = Rt.gen()
    Aold = Rt([QQ(value) for value in model["A_coefficients_low_to_high"]])
    Bold = Rt([QQ(value) for value in model["B_coefficients_low_to_high"]])
    Delta_old = -16 * (4 * Aold**3 + 27 * Bold**2)
    basis_xy = reconstruct_basis(Rt, Aold, Bold, section_data)
    Eold = EllipticCurve(Kt, [Aold, Bold])
    published_basis = [Eold(Kt(x), Kt(y)) for x, y in basis_xy]

    new_u_short = matrix(
        ZZ, local["r17_example"]["relative_U"]["target_U_basis_in_U_plus_short_R17"]
    )
    trace_short = vector(ZZ, new_u_short.row(0)[2:])
    zero_short = vector(ZZ, (new_u_short.row(1) - new_u_short.row(0))[2:])
    trace_pinned = trace_short * short_coordinates
    zero_pinned = zero_short * short_coordinates
    assert trace_short * short * trace_short == 12
    assert zero_short * short * zero_short == 4

    trace_published = trace_pinned * basis_change.transpose()
    zero_published = zero_pinned * basis_change.transpose()
    trace_point = sum(
        (coefficient * point for coefficient, point in zip(trace_published, published_basis)),
        Eold(0),
    )
    zero_point = sum(
        (coefficient * point for coefficient, point in zip(zero_published, published_basis)),
        Eold(0),
    )
    xP, yP = Kt(trace_point[0]), Kt(trace_point[1])
    x_zero, y_zero = Kt(zero_point[0]), Kt(zero_point[1])
    assert yP**2 == xP**3 + Aold * xP + Bold
    assert y_zero**2 == x_zero**3 + Aold * x_zero + Bold

    x_denominator = Rt(xP.denominator())
    if not x_denominator.is_square():
        raise ArithmeticError("trace x denominator is not a square")
    h = x_denominator.sqrt()
    h /= h.leading_coefficient()
    if xP.denominator() != h**2 or yP.denominator() != h**3 or h.degree() != 4:
        raise ArithmeticError("trace section has the wrong normalized pole divisor")
    Nx = Rt(xP * h**2)
    Ny = Rt(yP * h**3)
    if Nx.gcd(h) != 1:
        raise ArithmeticError("trace x numerator is not invertible modulo h")
    M0 = Rt((-Ny * Nx.inverse_mod(h**2)) % h**2)
    assert (M0 * Nx + Ny) % h**2 == 0

    # Universal nonzero-trace degree-two coefficient block: c=4 and k=-1,
    # hence deg(a)<=7, deg(b)<=1 and 2c=8 congruence rows.
    columns = [(t**degree * Nx) % h**2 for degree in range(8)]
    columns.extend([(-t**degree * Ny) % h**2 for degree in range(2)])
    rr_matrix = matrix(QQ, 8, 10, lambda i, j: columns[j][i])
    rr_kernel = rr_matrix.right_kernel_matrix()
    if rr_matrix.rank() != 8 or rr_kernel.nrows() != 2:
        raise ArithmeticError("degree-two Riemann--Roch kernel is not two-dimensional")
    ab = []
    for kernel_row in rr_kernel.rows():
        kernel_row = vector(QQ, kernel_row)
        a = Rt(list(kernel_row[:8]))
        b = Rt(list(kernel_row[8:]))
        assert (a * Nx - b * Ny) % h**2 == 0
        ab.append((a, b))
    (a0, b0), (a1, b1) = ab

    Ru = PolynomialRing(QQ, "u")
    Ku = Ru.fraction_field()
    u = Ru.gen()
    Stu = PolynomialRing(Ku, "t")
    Ftu = Stu.fraction_field()
    tt = Stu.gen()
    lift_t = lambda value: Stu([Ku(coefficient) for coefficient in Rt(value)])
    hh, NNx, NNy = map(lift_t, (h, Nx, Ny))
    AAold, BBold = map(lift_t, (Aold, Bold))
    aa0, bb0, aa1, bb1 = map(lift_t, (a0, b0, a1, b1))
    xxP, yyP = Ftu(NNx / hh**2), Ftu(NNy / hh**3)
    numerator_m = aa1 - u * aa0
    denominator_m = u * bb0 - bb1
    slope_m = Ftu(numerator_m / (denominator_m * hh))
    radical = (
        slope_m**4
        - 6 * xxP * slope_m**2
        - 8 * yyP * slope_m
        - 3 * xxP**2
        - 4 * AAold
    )
    radical_numerator = Stu(radical.numerator())
    radical_denominator = Stu(radical.denominator())
    square_factor = radical_numerator.gcd(radical_numerator.derivative()).monic()
    quartic, remainder = radical_numerator.quo_rem(square_factor**2)
    if remainder or quartic.degree() != 4 or quartic.gcd(quartic.derivative()).degree():
        raise ArithmeticError("residual chord did not produce a smooth generic quartic")
    if not radical_denominator.is_square():
        raise ArithmeticError("chord radical denominator is not a square")
    denominator_sqrt = radical_denominator.sqrt()
    assert radical == Ftu(quartic * (square_factor / denominator_sqrt) ** 2)

    def pencil_map(point):
        Xold, Yold = Kt(point[0]), Kt(point[1])
        L0 = a0 * (Xold * h**2 - Nx) + b0 * (Yold * h**3 + Ny)
        L1 = a1 * (Xold * h**2 - Nx) + b1 * (Yold * h**3 + Ny)
        return Kt(L1 / L0)

    def quartic_point_from_affine_source(point):
        old_base_map = pencil_map(point)
        t_section = Ku(invert_mobius(old_base_map, u))
        x_section = Ku(evaluate_rational(Kt(point[0]), t_section))
        y_section = Ku(evaluate_rational(Kt(point[1]), t_section))
        slope_section = evaluate_rational(slope_m, t_section)
        xP_section = evaluate_rational(xxP, t_section)
        yP_section = evaluate_rational(yyP, t_section)
        if y_section + yP_section != slope_section * (x_section - xP_section):
            raise ArithmeticError("source section misses the compiled chord")
        radical_root = 2 * x_section - (slope_section**2 - xP_section)
        W_section = radical_root * evaluate_polynomial(
            denominator_sqrt, t_section
        ) / evaluate_polynomial(square_factor, t_section)
        if W_section**2 != evaluate_polynomial(quartic, t_section):
            raise ArithmeticError("source section misses the compiled quartic")
        return old_base_map, Ku(t_section), Ku(W_section)

    def old_zero_quartic_point():
        old_base_map = Kt(b1 / b0)
        t_section = Ku(invert_mobius(old_base_map, u))
        if denominator_m.degree() != 1:
            raise ArithmeticError("old zero does not cut a linear quartic point")
        assert t_section == Ku(-denominator_m[0] / denominator_m[1])
        normalization = Ftu(denominator_sqrt / denominator_m**2)
        if Stu(normalization.numerator()).degree() or Stu(normalization.denominator()).degree():
            raise ArithmeticError("old-zero radical normalization depends on t")
        normalization_u = Ku(
            Stu(normalization.numerator())[0] / Stu(normalization.denominator())[0]
        )
        W_section = Ku(
            normalization_u
            * numerator_m(t_section) ** 2
            / (hh(t_section) ** 2 * square_factor(t_section))
        )
        assert W_section**2 == quartic(t_section)
        return old_base_map, t_section, W_section

    # The marked target zero is P_{-b1}, not the old zero.  Pointing here is
    # what keeps the non-cyclic relative-U marking literal at equation level.
    zero_base_map, t0, v0 = quartic_point_from_affine_source(zero_point)
    assert v0 and v0**2 == quartic(t0)

    Sz = PolynomialRing(Ku, "z")
    z = Sz.gen()
    shifted = Sz(quartic(t0 + z))
    ee, dd, cc, bbb, aaa = [Ku(shifted[index]) for index in range(5)]
    assert ee == v0**2
    a1g = dd / v0
    a2g = cc - dd**2 / (4 * v0**2)
    a3g = 2 * v0 * bbb
    a4g = -4 * v0**2 * aaa
    a6g = a2g * a4g
    b2g = a1g**2 + 4 * a2g
    b4g = a1g * a3g + 2 * a4g
    b6g = a3g**2 + 4 * a6g
    c4g = b2g**2 - 24 * b4g
    c6g = -b2g**3 + 36 * b2g * b4g - 216 * b6g

    e, d, c, b, a = [Ku(quartic[index]) for index in range(5)]
    invariant_I = 12 * a * e - 3 * b * d + c**2
    invariant_J = (
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d**2
        - 27 * b**2 * e
        - 2 * c**3
    )
    Araw = Ku(-27 * invariant_I)
    Braw = Ku(-27 * invariant_J)
    assert Ku(81 * (-c4g / 48)) == Araw
    assert Ku(729 * (-c6g / 864)) == Braw

    factors_A = list(Ru(Araw.denominator()).factor())
    factors_B = list(Ru(Braw.denominator()).factor())
    if len(factors_A) != 1 or factors_A[0][1] != 8:
        raise ArithmeticError("unexpected A denominator profile")
    ell = factors_A[0][0].monic()
    if len(factors_B) != 1 or factors_B[0][1] != 12 or factors_B[0][0].monic() != ell:
        raise ArithmeticError("unexpected B denominator profile")
    gauge = ell**2
    Achild = Ru(Araw * gauge**4)
    Bchild = Ru(Braw * gauge**6)
    Delta_child = Ru(-16 * (4 * Achild**3 + 27 * Bchild**2))
    if (Achild.degree(), Bchild.degree(), Delta_child.degree()) != (8, 12, 24):
        raise ArithmeticError("child equation lost the K3 degree profile")
    if Achild.gcd(Delta_child).degree():
        raise ArithmeticError("finite child fibres are not multiplicative")
    delta_factors = list(Delta_child.factor())
    factor_profile = sorted((factor.degree(), int(exponent)) for factor, exponent in delta_factors)
    if factor_profile != [(1, 2)] * 4 + [(16, 1)]:
        raise ArithmeticError(f"unexpected 4A1 discriminant factorization: {factor_profile}")

    # Publication normalization.  Send the three largest rational I2 bases to
    # 0, 1, infinity.  The fourth becomes -146234/269481.  After removing the
    # common rational square from the weighted coefficients, the remaining
    # twist class is 3 and the equation has integral coefficients 27*A0,54*B0.
    i2_roots_u = sorted(
        [
            QQ(-factor[0] / factor[1])
            for factor, exponent in delta_factors
            if factor.degree() == 1 and exponent == 2
        ],
        reverse=True,
    )
    if len(i2_roots_u) != 4:
        raise ArithmeticError("failed to recover the four rational I2 bases")
    i2_a, i2_b, i2_c, i2_d = i2_roots_u
    Rs = PolynomialRing(QQ, "s")
    Ks = Rs.fraction_field()
    s = Rs.gen()
    base_numerator = s * i2_c * (i2_b - i2_a) - i2_a * (i2_b - i2_c)
    base_denominator = s * (i2_b - i2_a) - (i2_b - i2_c)
    u_of_s = Ks(base_numerator / base_denominator)
    s_of_u = Ku(
        (u - i2_a) * (i2_b - i2_c) / ((u - i2_c) * (i2_b - i2_a))
    )
    assert evaluate_rational(s_of_u, u_of_s) == s
    fourth_i2_s = QQ(
        (i2_d - i2_a)
        * (i2_b - i2_c)
        / ((i2_d - i2_c) * (i2_b - i2_a))
    )
    assert fourth_i2_s == -QQ(146234) / 269481

    def homogenized_base_change(poly, degree):
        return Rs(
            sum(
                poly[index]
                * base_numerator**index
                * base_denominator ** (degree - index)
                for index in range(poly.degree() + 1)
            )
        )

    A_homogeneous = homogenized_base_change(Achild, 8)
    B_homogeneous = homogenized_base_change(Bchild, 12)
    A_content = polynomial_content(A_homogeneous)
    B_content = polynomial_content(B_homogeneous)
    A_primitive = Rs(A_homogeneous / A_content)
    B_primitive = Rs(B_homogeneous / B_content)
    assert all(QQ(value).denominator() == 1 for value in A_primitive)
    assert all(QQ(value).denominator() == 1 for value in B_primitive)
    assert A_content**3 / B_content**2 == QQ(27) / 4
    twist_parameter = QQ(3 * B_content / (2 * A_content))
    assert A_content == 3 * twist_parameter**2
    assert B_content == 2 * twist_parameter**3
    twist_square = QQ(twist_parameter / 3)
    if not twist_square.numerator().is_square() or not twist_square.denominator().is_square():
        raise ArithmeticError("publication twist does not have square class 3")
    publication_scale = QQ(
        twist_square.numerator().sqrt() / twist_square.denominator().sqrt()
    )
    A_publication = Rs(A_homogeneous / publication_scale**4)
    B_publication = Rs(B_homogeneous / publication_scale**6)
    assert A_publication == 27 * A_primitive
    assert B_publication == 54 * B_primitive
    Delta_publication = Rs(-16 * (4 * A_publication**3 + 27 * B_publication**2))
    if (A_publication.degree(), B_publication.degree(), Delta_publication.degree()) != (8, 12, 22):
        raise ArithmeticError("publication model has the wrong finite/infinite degree profile")
    if A_publication.gcd(Delta_publication).degree():
        raise ArithmeticError("publication model has an additive finite fibre")
    publication_delta_factors = list(Delta_publication.factor())
    publication_factor_profile = sorted(
        (factor.degree(), int(exponent))
        for factor, exponent in publication_delta_factors
    )
    if publication_factor_profile != [(1, 2)] * 3 + [(16, 1)]:
        raise ArithmeticError(
            f"unexpected normalized discriminant factorization: {publication_factor_profile}"
        )

    def quartic_to_child(t_section, W_section):
        z_section = t_section - t0
        x_general = (2 * v0 * (W_section + v0) + dd * z_section) / z_section**2
        y_general = (
            4 * v0**2 * (W_section + v0)
            + 2 * v0 * dd * z_section
            + (2 * v0 * cc - dd**2 / (2 * v0)) * z_section**2
        ) / z_section**3
        X_section = Ku(gauge**2) * 9 * (x_general + b2g / 12)
        Y_section = Ku(gauge**3) * 27 * (
            y_general + (a1g * x_general + a3g) / 2
        )
        if Y_section**2 != X_section**3 + Achild * X_section + Bchild:
            raise ArithmeticError("quartic point misses the child Weierstrass equation")
        return Ku(X_section), Ku(Y_section)

    def raw_child_to_publication(X_section, Y_section):
        X_substituted = evaluate_rational(X_section, u_of_s)
        Y_substituted = evaluate_rational(Y_section, u_of_s)
        X_publication = Ks(
            base_denominator**4 * X_substituted / publication_scale**2
        )
        Y_publication = Ks(
            base_denominator**6 * Y_substituted / publication_scale**3
        )
        if Y_publication**2 != (
            X_publication**3
            + A_publication * X_publication
            + B_publication
        ):
            raise ArithmeticError("section misses the publication model")
        return X_publication, Y_publication

    # Marking and ADE certificate in the pinned R17 basis.
    ns_pinned = block_diagonal_matrix(HYPERBOLIC, -pinned)
    new_u_pinned_rows = []
    for source_row in new_u_short.rows():
        new_u_pinned_rows.append(
            list(source_row[:2]) + list(vector(ZZ, source_row[2:]) * short_coordinates)
        )
    new_u_pinned = matrix(ZZ, new_u_pinned_rows)
    fibre = new_u_pinned.row(0)
    new_zero = new_u_pinned.row(1) - fibre
    assert fibre == vector(ZZ, [3, 2] + list(trace_pinned))
    assert new_zero == vector(ZZ, [1, 1] + list(zero_pinned))
    assert new_u_pinned * ns_pinned * new_u_pinned.transpose() == HYPERBOLIC
    frame_basis = (new_u_pinned * ns_pinned).right_kernel_matrix()
    transport = new_u_pinned.stack(frame_basis)
    frame = -(frame_basis * ns_pinned * frame_basis.transpose())
    if abs(transport.det()) != 1 or frame.det() != 948:
        raise ArithmeticError("target U does not split primitively")

    simple_short = matrix(
        ZZ, local["r17_example"]["geometric_gate"]["simple_components_in_ambient_NS"]
    )
    simple_pinned = matrix(
        ZZ,
        [
            list(component[:2])
            + list(vector(ZZ, component[2:]) * short_coordinates)
            for component in simple_short.rows()
        ],
    )
    root_coordinates = frame_basis.solve_left(simple_pinned).change_ring(ZZ)
    root_gram = root_coordinates * frame * root_coordinates.transpose()
    if root_gram != 2 * matrix.identity(ZZ, 4):
        raise ArithmeticError("target roots are not 4A1")
    if abs(root_coordinates.row_module(ZZ).index_in(root_coordinates.row_module(ZZ).saturation())) != 1:
        raise ArithmeticError("target root span is not primitive")
    if int(pari(frame).qfminim(2)[0]) != 8:
        raise ArithmeticError("target frame has unexpected extra roots")

    local_frame = matrix(ZZ, local["r17_example"]["target_frame"]["gram"])
    local_isometry = row_isometry(frame, local_frame)
    if local_isometry is None:
        raise ArithmeticError("compiled frame misses the local-mutation target")
    historical_basis = matrix(
        ZZ,
        historical["current_suffix_stages"]["current_4A1"]["basis_in_pinned_R17"],
    )
    historical_split = historical_basis * ns_pinned * historical_basis.transpose()
    historical_frame = -historical_split[2:, 2:]
    physical_q8_frame = load_matrix(PHYSICAL_Q8)
    if row_isometry(frame, historical_frame) is not None:
        raise ArithmeticError("compiled frame equals the historical H3 4A1 frame")
    if row_isometry(frame, physical_q8_frame) is not None:
        raise ArithmeticError("compiled frame equals the physical-q8 H3 4A1 frame")

    transport_inverse = transport.inverse()
    section_frame_rows = []
    shioda_rows = []
    section_records = []
    for basis_index, short_mw_entries in enumerate(MW13_SOURCE_SECTIONS_SHORT):
        short_mw = vector(ZZ, short_mw_entries)
        old_height = int(short_mw * short * short_mw)
        old_class = vector(ZZ, [(old_height - 2) // 2, 1] + list(short_mw * short_coordinates))
        assert fibre * ns_pinned * old_class == 1
        new_class = old_class * transport_inverse
        if new_class[1] != 1 or any(entry not in ZZ for entry in new_class):
            raise ArithmeticError("selected old section is not an integral child section")
        frame_row = vector(ZZ, new_class[2:])
        shioda_row = frame_row - (
            frame_row * frame * root_coordinates.transpose()
        ) * root_gram.inverse() * root_coordinates

        pinned_mw = short_mw * short_coordinates
        published_mw = pinned_mw * basis_change.transpose()
        source_point = sum(
            (coefficient * point for coefficient, point in zip(published_mw, published_basis)),
            Eold(0),
        )
        if source_point.is_zero():
            old_base_map, t_section, W_section = old_zero_quartic_point()
            source_kind = "published R17 old zero"
        else:
            old_base_map, t_section, W_section = quartic_point_from_affine_source(source_point)
            source_kind = "published R17 section"
        X_section, Y_section = quartic_to_child(t_section, W_section)
        X_publication, Y_publication = raw_child_to_publication(X_section, Y_section)
        publication_base_map = Kt(evaluate_rational(s_of_u, old_base_map))
        publication_t_section = Ks(evaluate_rational(t_section, u_of_s))

        section_frame_rows.append(frame_row)
        shioda_rows.append(shioda_row)
        section_records.append(
            {
                "basis_index": basis_index,
                "source": source_kind,
                "source_short_R17_coordinates": list(map(int, short_mw)),
                "source_pinned_R17_coordinates": list(map(int, pinned_mw)),
                "source_published_R17_coordinates": list(map(int, published_mw)),
                "source_old_height": old_height,
                "new_frame_coordinates": list(map(int, frame_row)),
                "publication_base_map_s_of_t": rational_function_record(publication_base_map),
                "old_base_t_as_function_of_s": rational_function_record(publication_t_section),
                "X": rational_function_record(X_publication),
                "Y": rational_function_record(Y_publication),
                "equation_verified": True,
            }
        )

    full_frame_generators = root_coordinates.stack(matrix(ZZ, section_frame_rows))
    if abs(full_frame_generators.det()) != 1:
        raise ArithmeticError("four roots plus thirteen sections do not span the frame")
    height_gram = matrix(QQ, shioda_rows) * frame * matrix(QQ, shioda_rows).transpose()
    if height_gram.rank() != 13 or height_gram.det() != QQ(237) / 4:
        raise ArithmeticError("MW13 height lattice is not saturated")

    # Exact generic inverse.  Starting with W^2=q(t,u), the pointed-quartic
    # inverse recovers t and W, then the residual chord recovers the literal
    # published R17 coordinates.  This certifies the reverse target-free hop
    # on a dense open, not merely equality of j-invariants.
    PW = PolynomialRing(Ftu, "Wbar")
    Wbar = PW.gen()
    quartic_field = PW.quotient(Wbar**2 - Ftu(quartic), names=("Wbar",))
    Wq = quartic_field.gen()
    zq = quartic_field(Ftu(tt - t0))
    xgq = (2 * v0 * (Wq + v0) + dd * zq) / zq**2
    ygq = (
        4 * v0**2 * (Wq + v0)
        + 2 * v0 * dd * zq
        + (2 * v0 * cc - dd**2 / (2 * v0)) * zq**2
    ) / zq**3
    inverse_z = (4 * v0**2 * (xgq + cc) - dd**2) / (2 * v0 * ygq)
    inverse_W = (xgq * inverse_z**2 - dd * inverse_z) / (2 * v0) - v0
    assert inverse_z == zq
    assert inverse_W == Wq
    radical_root_q = Wq * quartic_field(Ftu(square_factor / denominator_sqrt))
    old_x_q = (radical_root_q + slope_m**2 - xxP) / 2
    old_y_q = slope_m * (old_x_q - xxP) - yyP
    assert old_y_q**2 == old_x_q**3 + AAold * old_x_q + BBold
    L0q = aa0 * (old_x_q * hh**2 - NNx) + bb0 * (old_y_q * hh**3 + NNy)
    L1q = aa1 * (old_x_q * hh**2 - NNx) + bb1 * (old_y_q * hh**3 + NNy)
    assert L1q / L0q == u

    publication_delta_factor_records = [
        {
            "coefficients_low_to_high": polynomial_text(factor),
            "degree": int(factor.degree()),
            "multiplicity": int(exponent),
        }
        for factor, exponent in publication_delta_factors
    ]
    input_paths = [
        MODEL,
        SECTIONS,
        PINNED_GRAM,
        SHORT_GRAM,
        SHORT_COORDS,
        PUBLISHED,
        LOCAL_MUTATION,
        PRIME_LOCAL,
        PLANNER,
        HISTORICAL,
        PHYSICAL_Q8,
    ]
    result = {
        "schema": "elkies-k3.r17-noncyclic-4a1-direct-fibration.v1",
        "status": "PASS_EXACT_NEW_NONCYCLIC_4A1_EQUATION_REVERSE_HOP_AND_MW13",
        "route": {
            "forward": "published R17 equation -> new 4A1/MW13 equation",
            "reverse": "target-free rootless marked-U selection -> published R17 equation",
            "old_fibre_degree_both_directions": 2,
            "reverse_target_frame_supplied": False,
            "reverse_selected_frame": "published R17 rootless frame",
        },
        "marked_u": {
            "basis_convention": "(F,F+O) followed by positive frame coordinates",
            "target_U_in_U_plus_short_R17": rows(new_u_short),
            "target_U_in_U_plus_pinned_R17": rows(new_u_pinned),
            "target_fibre_D": list(map(int, fibre)),
            "target_zero_O_prime": list(map(int, new_zero)),
            "cross_pairing_A": local["r17_example"]["relative_U"]["cross_pairing_A"],
            "primitive_split_transport_determinant": int(transport.det()),
            "nef": True,
            "physical_zero": True,
        },
        "prime_local_bridge": {
            "saturated_discriminant_group": "Z/4 + Z/8",
            "saturated_discriminant_invariants": [4, 8],
            "saturation_index": 1,
            "glue_order": 32,
            "maximal": True,
            "marked_graph_count": 32,
            "rootless_to_4A1_graph_multiplicity": 4,
        },
        "trace_section": {
            "short_R17_coordinates": list(map(int, trace_short)),
            "pinned_R17_coordinates": list(map(int, trace_pinned)),
            "published_R17_coordinates": list(map(int, trace_published)),
            "h_coefficients_low_to_high": polynomial_text(h),
            "Nx_coefficients_low_to_high": polynomial_text(Nx),
            "Ny_coefficients_low_to_high": polynomial_text(Ny),
            "M0_coefficients_low_to_high": polynomial_text(M0),
            "height": 12,
            "zero_intersection": 4,
        },
        "riemann_roch": {
            "method": "universal marked degree-two nonzero-trace chord compiler",
            "padding_k": -1,
            "trace_zero_intersection_c": 4,
            "degree_bounds": {"a": 7, "b": 1},
            "raw_coefficient_count": 10,
            "congruence_rank": int(rr_matrix.rank()),
            "kernel_dimension": int(rr_kernel.nrows()),
            "kernel_rows_a0_through_a7_b0_b1": [
                [rational_text(value) for value in row] for row in rr_kernel.rows()
            ],
            "pencil_coordinate": "u=L1/L0",
        },
        "genus_one_model": {
            "equation": "W^2=q(t,u)",
            "q_coefficients_in_t_low_to_high": polynomial_over_function_field_record(quartic),
            "radical_square_factor_coefficients_in_t_low_to_high": polynomial_over_function_field_record(square_factor),
            "radical_denominator_square_root_coefficients_in_t_low_to_high": polynomial_over_function_field_record(denominator_sqrt),
            "marked_zero": {
                "source_short_R17_coordinates": list(map(int, zero_short)),
                "source_pinned_R17_coordinates": list(map(int, zero_pinned)),
                "source_published_R17_coordinates": list(map(int, zero_published)),
                "old_base_map_u_of_t": rational_function_record(zero_base_map),
                "publication_base_map_s_of_t": rational_function_record(
                    Kt(evaluate_rational(s_of_u, zero_base_map))
                ),
                "t0": rational_function_record(t0),
                "W0": rational_function_record(v0),
                "maps_to_child_point_at_infinity": True,
            },
            "quartic_squarefree_over_Q_of_u": True,
            "radical_identity_verified": True,
        },
        "weierstrass_model": {
            "coordinate": "s",
            "equation": "Y^2=X^3+27*A0(s)*X+54*B0(s)",
            "A0_primitive_coefficients_low_to_high": polynomial_text(A_primitive),
            "B0_primitive_coefficients_low_to_high": polynomial_text(B_primitive),
            "A_coefficients_low_to_high": polynomial_text(A_publication),
            "B_coefficients_low_to_high": polynomial_text(B_publication),
            "degrees_A_B_affine_Delta": [8, 12, 22],
            "projective_discriminant_degree": 24,
            "discriminant_coefficients_low_to_high": polynomial_text(
                Delta_publication
            ),
            "finite_discriminant_factorization": publication_delta_factor_records,
            "fibre_configuration": "4 I2 + 16 I1",
            "root_system": "4A1",
            "rational_I2_locations": [
                "0",
                "1",
                rational_text(fourth_i2_s),
                "infinity",
            ],
            "infinity_orders_c4_c6_Delta": [0, 0, 2],
            "normalization": {
                "compiler_coordinate": "u",
                "compiler_gauge": rational_function_record(Ku(gauge)),
                "compiler_I2_locations_u_descending": [
                    rational_text(value) for value in i2_roots_u
                ],
                "s_of_u": rational_function_record(s_of_u),
                "u_of_s": rational_function_record(u_of_s),
                "formula": (
                    "s=(u-a)(b-c)/((u-c)(b-a)); "
                    "(a,b,c) map to (0,1,infinity)"
                ),
                "fourth_cross_ratio": rational_text(fourth_i2_s),
                "weighted_content_A": rational_text(A_content),
                "weighted_content_B": rational_text(B_content),
                "twist_parameter": rational_text(twist_parameter),
                "twist_square_class": "3",
                "coordinate_scale_square_root_of_twist_over_3": rational_text(
                    publication_scale
                ),
                "raw_to_publication": (
                    "u=u(s), X=den(s)^4*X_raw/q^2, "
                    "Y=den(s)^6*Y_raw/q^3"
                ),
            },
            "pointed_quartic_map_verified": True,
        },
        "frame_certificate": {
            "frame_gram": rows(frame),
            "determinant": int(frame.det()),
            "signed_roots": 8,
            "simple_root_coordinates": rows(root_coordinates),
            "simple_root_gram": rows(root_gram),
            "root_span_saturation_index": 1,
            "integral_isometry_to_local_mutation_frame": rows(local_isometry),
            "isometric_to_historical_current_4A1": False,
            "isometric_to_physical_q8_orbit376_4A1": False,
        },
        "mordell_weil": {
            "field": "QQ(s)",
            "arithmetic_rank": 13,
            "geometric_rank": 13,
            "torsion_order": 1,
            "status": "PASS_EXACT_SATURATED_ARITHMETIC_MW13_BASIS",
            "section_count": 13,
            "root_plus_section_coordinate_determinant": int(full_frame_generators.det()),
            "height_gram": [
                [rational_text(value) for value in row] for row in height_gram.rows()
            ],
            "height_gram_determinant": "237/4",
            "sections": section_records,
        },
        "reverse_birational_hop": {
            "target_selection": "root_rank=0, ADE=rootless; no target frame Gram supplied",
            "selected_marking": "published R17 (F_old,O_old)",
            "inverse_pointed_quartic_formula": {
                "z": "(4*W0^2*(x_g+c)-d^2)/(2*W0*y_g)",
                "W": "(x_g*z^2-d*z)/(2*W0)-W0",
            },
            "inverse_chord_formula": {
                "x_old": "(W*square_factor/denominator_sqrt+m^2-x(trace))/2",
                "y_old": "m*(x_old-x(trace))-y(trace)",
            },
            "generic_inverse_recovers_t_and_W": True,
            "generic_inverse_recovers_pencil_coordinate_u": True,
            "publication_base_change_u_of_s": rational_function_record(u_of_s),
            "generic_inverse_recovers_publication_coordinate_s": True,
            "literal_published_R17_equation_verified": True,
            "published_R17_equation": published["published_equation"]["form"],
        },
        "inputs": {relative(path): digest(path) for path in input_paths},
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "pari_version": ".".join(map(str, pari.version())),
            "required_features": [
                "Sage exact function fields",
                "Sage exact quadratic quotients",
                "PARI qfminim",
                "Sage integral quadratic-form isometry",
            ],
        },
        "reproducing_command": "sage -python elkies-k3/scripts/compile_r17_noncyclic_4a1_qq.sage --check",
        "proof_boundary": (
            "This exact replay compiles the marked degree-two pencil over QQ, points it at "
            "the certified physical zero, proves the fibre configuration 4I2+16I1, identifies "
            "the new non-historical 4A1 frame, and gives a saturated arithmetic MW13 basis. "
            "The generic inverse returns the literal published R17 equation through the "
            "target-free rootless marked-U control. It does not classify the J1 automorphism "
            "orbit or search for rank-jumping specializations."
        ),
    }

    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError(f"stored artifact differs from replay: {output}")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "R17NONCYCLIC4A1|h0=2|fibres=4I2+16I1|bridge=Z/4+Z/8|"
        "mw=13-saturated|reverse=target-free-R17|status={}|output={}".format(
            result["status"], relative(output)
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
