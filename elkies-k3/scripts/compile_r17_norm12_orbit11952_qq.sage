#!/usr/bin/env sage-python
"""Compile the orbit-11952 two-neighbor directly on the published R17 model.

For the norm-twelve section ``P_w`` the new fibre is

    D = O + P_w - F = (3,2,w).

Brandhorst--Elkies' degree-two-neighbor description realizes the two sections
of ``O_X(D)`` as

    a(t) + b(t) (y + y(P_w))/(x - x(P_w)),

with ``deg(a)<=7``, ``deg(b)<=1`` and one congruence modulo ``h(t)^2``.
This script solves that exact Riemann--Roch kernel, forms the resulting quartic
genus-one model, takes its classical Jacobian, and certifies the new elliptic
frame.  The section-recovery block also gives a saturated rank-17 basis on the
new equation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    Conic,
    EllipticCurve,
    PolynomialRing,
    QQ,
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
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
TARGET = ROOT / "artifacts/generated-results/elkies-2026-published-r17-target.json"
SPLITTING = ROOT / "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
CLASSIFICATION = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-isotropic-frame-classification-v1.json"
ALTERNATE = ROOT / "artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json"
BISECTIONS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
OUTPUT_103B2 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit103b2-direct-fibration-v1.json"
OUTPUT_08F72 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit08f72-direct-fibration-v1.json"
OUTPUT_08AB4 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit08ab4-direct-fibration-v1.json"
FRAME_103B2 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-isotropic-frame-v1.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def load_matrix(path: Path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def matrix_digest(value) -> str:
    payload = json.dumps(rows(value), separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def polynomial_text(poly) -> list[str]:
    if not poly:
        return ["0"]
    return [rational_text(poly[i]) for i in range(poly.degree() + 1)]


def rational_function_record(value):
    value = value.parent()(value)
    return {
        "numerator_coefficients_low_to_high": polynomial_text(value.numerator()),
        "denominator_coefficients_low_to_high": polynomial_text(value.denominator()),
    }


def reconstruct_basis(R, A, B, section_data):
    points = []
    for expected_index, record in enumerate(section_data["sections"]):
        assert int(record["basis_index"]) == expected_index
        x_coordinate = R([QQ(value) for value in record["x_coefficients_low_to_high"]])
        if expected_index == 0:
            y_coordinate = R([QQ(value) for value in record["y_coefficients_low_to_high"]])
        else:
            reference_x, reference_y = points[int(record["chord"]["reference_basis_index"])]
            slope = R([QQ(value) for value in record["chord"]["slope_coefficients_low_to_high"]])
            y_coordinate = reference_y + slope * (x_coordinate - reference_x)
        assert y_coordinate**2 == x_coordinate**3 + A * x_coordinate + B
        points.append((x_coordinate, y_coordinate))
    return points


def find_record(payload, label):
    for record in payload["construction"]["records"]:
        if record["label"] == label:
            return record
    raise KeyError(label)


def exact_square_root(poly):
    if not poly.is_square():
        raise ArithmeticError("expected an exact polynomial square")
    return poly.sqrt()


def evaluate_polynomial(poly, value):
    """Evaluate a univariate polynomial by Horner in a compatible field."""

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
        raise ArithmeticError("failed to invert a section's Mobius base map")
    return inverse


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-label",
        choices=(
            "norm12-orbit-11952",
            "norm12-orbit-103b2",
            "norm12-orbit-08f72",
            "norm12-orbit-08ab4",
        ),
        default="norm12-orbit-11952",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    default_outputs = {
        "norm12-orbit-11952": OUTPUT,
        "norm12-orbit-103b2": OUTPUT_103B2,
        "norm12-orbit-08f72": OUTPUT_08F72,
        "norm12-orbit-08ab4": OUTPUT_08AB4,
    }
    output = args.output or default_outputs[args.source_label]

    model = json.loads(MODEL.read_text())
    section_data = json.loads(SECTIONS.read_text())
    splitting = json.loads(SPLITTING.read_text())
    classification = json.loads(CLASSIFICATION.read_text())
    alternate_payload = json.loads(ALTERNATE.read_text())
    target_payload = json.loads(TARGET.read_text())
    source_record = find_record(splitting, args.source_label)

    Rt = PolynomialRing(QQ, "t")
    Kt = Rt.fraction_field()
    t = Rt.gen()
    Aold = Rt([QQ(value) for value in model["A_coefficients_low_to_high"]])
    Bold = Rt([QQ(value) for value in model["B_coefficients_low_to_high"]])
    Delta_old = -16 * (4 * Aold**3 + 27 * Bold**2)
    basis_xy = reconstruct_basis(Rt, Aold, Bold, section_data)
    Eold = EllipticCurve(Kt, [Aold, Bold])
    published_basis = [Eold(Kt(x), Kt(y)) for x, y in basis_xy]

    trace_data = source_record["trace_section"]
    h = Rt([QQ(value) for value in trace_data["h_coefficients_low_to_high"]])
    Nx = Rt([QQ(value) for value in trace_data["Nx_coefficients_low_to_high"]])
    Ny = Rt([QQ(value) for value in trace_data["Ny_coefficients_low_to_high"]])
    M0 = Rt([QQ(value) for value in trace_data["M0_coefficients_low_to_high"]])
    xP, yP = Kt(Nx / h**2), Kt(Ny / h**3)
    assert yP**2 == xP**3 + Aold * xP + Bold
    assert (M0 * Nx + Ny) % h**2 == 0

    # Proposition 2.17 of Brandhorst--Elkies, specialized to P.O=4 and
    # D=O+P-F.  Eight congruence coefficients constrain ten unknowns.
    columns = []
    for degree in range(8):
        columns.append((t**degree * Nx) % h**2)
    for degree in range(2):
        columns.append((-t**degree * Ny) % h**2)
    rr_matrix = matrix(QQ, 8, 10, lambda i, j: columns[j][i])
    rr_kernel = rr_matrix.right_kernel_matrix()
    if rr_matrix.rank() != 8 or rr_kernel.nrows() != 2:
        raise ArithmeticError("norm-twelve Riemann--Roch kernel is not two-dimensional")
    kernel_rows = [vector(QQ, row) for row in rr_kernel.rows()]
    ab = []
    for row in kernel_rows:
        a = Rt(list(row[:8]))
        b = Rt(list(row[8:]))
        assert (a * Nx - b * Ny) % h**2 == 0
        ab.append((a, b))
    (a0, b0), (a1, b1) = ab

    Ru = PolynomialRing(QQ, "u")
    Ku = Ru.fraction_field()
    u = Ru.gen()
    Stu = PolynomialRing(Ku, "t")
    Ftu = Stu.fraction_field()
    tt = Stu.gen()
    lift_t = lambda value: Stu([Ku(c) for c in Rt(value)])
    hh, NNx, NNy = map(lift_t, (h, Nx, Ny))
    AAold = lift_t(Aold)
    xxP, yyP = Ftu(NNx / hh**2), Ftu(NNy / hh**3)
    aa0, bb0, aa1, bb1 = map(lift_t, (a0, b0, a1, b1))
    numerator_m = aa1 - u * aa0
    denominator_m = u * bb0 - bb1
    slope_m = Ftu(numerator_m / (denominator_m * hh))
    radical = slope_m**4 - 6 * xxP * slope_m**2 - 8 * yyP * slope_m - 3 * xxP**2 - 4 * AAold
    radical_numerator = Stu(radical.numerator())
    radical_denominator = Stu(radical.denominator())
    square_factor = radical_numerator.gcd(radical_numerator.derivative()).monic()
    quartic, remainder = radical_numerator.quo_rem(square_factor**2)
    if remainder or quartic.degree() != 4 or quartic.gcd(quartic.derivative()).degree() != 0:
        raise ArithmeticError("residual chord did not produce a smooth quartic")
    denominator_sqrt = exact_square_root(radical_denominator)
    assert radical == Ftu(quartic * (square_factor / denominator_sqrt) ** 2)

    # The old zero is a shared zero.  At the pole of the chord slope it gives
    # the distinguished rational point on the quartic.
    if denominator_m.degree() != 1:
        raise ArithmeticError("the old zero did not give a linear quartic point")
    t0 = Ku(-denominator_m[0] / denominator_m[1])
    normalization = Ftu(denominator_sqrt / denominator_m**2)
    if Stu(normalization.numerator()).degree() or Stu(normalization.denominator()).degree():
        raise ArithmeticError("the chord radical normalization unexpectedly depends on t")
    normalization_u = Ku(Stu(normalization.numerator())[0] / Stu(normalization.denominator())[0])
    v0 = Ku(
        normalization_u * numerator_m(t0) ** 2
        / (hh(t0) ** 2 * square_factor(t0))
    )
    quartic_at_zero = Ku(quartic(t0))
    if v0**2 != quartic_at_zero:
        raise ArithmeticError("the shared old zero is not the pinned rational quartic point")
    assert v0**2 == quartic(t0)

    # Classical binary-quartic invariants, in the short normalization used by
    # the pointed quartic formulas below.
    e, d, c, b, a = [Ku(quartic[i]) for i in range(5)]
    invariant_I = 12 * a * e - 3 * b * d + c**2
    invariant_J = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
    Araw = Ku(-27 * invariant_I)
    Braw = Ku(-27 * invariant_J)
    factors_A = list(Ru(Araw.denominator()).factor())
    factors_B = list(Ru(Braw.denominator()).factor())
    if len(factors_A) != 1 or factors_A[0][1] != 8 or len(factors_B) != 1 or factors_B[0][1] != 12:
        raise ArithmeticError("unexpected Jacobian denominator profile")
    ell = factors_A[0][0].monic()
    if factors_B[0][0].monic() != ell:
        raise ArithmeticError("A and B do not have the same Jacobian pole")
    gauge = ell**2
    Achild = Ru(Araw * gauge**4)
    Bchild = Ru(Braw * gauge**6)
    Delta_child = Ru(-16 * (4 * Achild**3 + 27 * Bchild**2))
    if Achild.degree() != 8 or Bchild.degree() != 12 or Delta_child.degree() != 24:
        raise ArithmeticError("child equation lost the K3 degree profile")
    if Delta_child.gcd(Delta_child.derivative()).degree() or Achild.gcd(Delta_child).degree():
        raise ArithmeticError("child discriminant is not squarefree with nodal finite fibres")
    delta_factors = list(Delta_child.factor())
    if delta_factors != [(Delta_child / Delta_child.leading_coefficient(), 1)]:
        # Sage's factorization preserves the unit separately; check the actual
        # irreducible support rather than its normalization.
        if len(delta_factors) != 1 or delta_factors[0][1] != 1 or delta_factors[0][0].degree() != 24:
            raise ArithmeticError("child discriminant is not irreducible of degree 24")

    # Verify the pointed-quartic birational map and its compatibility with the
    # invariant model.  Put z=t-t0 and q(t0+z)=e+d*z+c*z^2+b*z^3+a*z^4.
    Sz = PolynomialRing(Ku, "z")
    z = Sz.gen()
    shifted = Sz(quartic(t0 + z))
    ee, dd, cc, bbb, aaa = [Ku(shifted[i]) for i in range(5)]
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
    Agen = 81 * (-c4g / 48)
    Bgen = 729 * (-c6g / 864)
    # The explicit X=9(x+b2/12), Y=27(y+(a1*x+a3)/2) normalization.
    assert Ku(Agen) == Araw and Ku(Bgen) == Braw

    pinned = load_matrix(PINNED)
    hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
    ns = block_diagonal_matrix(hyperbolic, -pinned)
    w = vector(ZZ, source_record["pinned_rank17_w"])
    fibre = vector(ZZ, [3, 2] + list(w))
    old_zero = vector(ZZ, [-1, 1] + [0] * 17)
    mate = fibre + old_zero
    complement = matrix(ZZ, [list(fibre * ns), list(mate * ns)]).right_kernel_matrix()
    transport = matrix(ZZ, [list(fibre), list(mate)] + rows(complement))
    frame = -(complement * ns * complement.transpose())
    alternate = matrix(ZZ, alternate_payload["rootless_frame"])
    is_alternate_target = args.source_label != "norm12-orbit-103b2"
    expected_frame = alternate if is_alternate_target else pinned
    rejected_frame = pinned if is_alternate_target else alternate
    if abs(transport.det()) != 1 or frame.det() != 948 or int(pari(frame).qfminim(2)[0]):
        raise ArithmeticError("new frame failed its primitive rootless certificate")
    isometry = pari(frame).qfisom(pari(expected_frame))
    if isometry == 0:
        raise ArithmeticError("new frame is not integrally isometric to the expected rootless frame")
    isometry = matrix(ZZ, isometry)
    if isometry * expected_frame * isometry.transpose() != frame:
        if isometry.transpose() * frame * isometry == expected_frame:
            isometry = matrix(ZZ, isometry.inverse().transpose())
        elif isometry.transpose() * expected_frame * isometry == frame:
            isometry = isometry.transpose()
        elif isometry * frame * isometry.transpose() == expected_frame:
            isometry = matrix(ZZ, isometry.inverse())
        else:
            raise ArithmeticError("PARI returned an unrecognized qfisom orientation")
    assert isometry * expected_frame * isometry.transpose() == frame
    if pari(frame).qfisom(pari(rejected_frame)) != 0:
        raise ArithmeticError("new frame unexpectedly matches the other rootless J2 class")

    classification_record = next(
        item for item in classification["classification"]["records"]
        if item["label"] == args.source_label
    )
    expected_frame_class = "alternate-Q80" if is_alternate_target else "published-R17"
    assert classification_record["frame_class"] == expected_frame_class
    assert classification_record["frame_gram_sha256"] == matrix_digest(frame)

    def point_on_child(t_section, x_section, y_section):
        """Transport one old-curve point on a D-fibre to the child Jacobian."""

        slope_section = evaluate_rational(slope_m, t_section)
        xP_section = evaluate_rational(xxP, t_section)
        yP_section = evaluate_rational(yyP, t_section)
        if y_section + yP_section != slope_section * (x_section - xP_section):
            raise ArithmeticError("source point does not lie on the compiled D-fibre chord")
        radical_root = 2 * x_section - (slope_section**2 - xP_section)
        W_section = radical_root * evaluate_polynomial(
            denominator_sqrt, t_section
        ) / evaluate_polynomial(square_factor, t_section)
        if W_section**2 != evaluate_polynomial(quartic, t_section):
            raise ArithmeticError("source point did not land on the quartic")
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
            raise ArithmeticError("transported point failed the child Weierstrass equation")
        return Ku(X_section), Ku(Y_section), Ku(W_section)

    def map_old_section(point):
        Xold, Yold = Kt(point[0]), Kt(point[1])
        L0 = a0 * (Xold * h**2 - Nx) + b0 * (Yold * h**3 + Ny)
        L1 = a1 * (Xold * h**2 - Nx) + b1 * (Yold * h**3 + Ny)
        old_base_map = Kt(L1 / L0)
        t_section = Ku(invert_mobius(old_base_map, u))
        if evaluate_rational(old_base_map, t_section) != u:
            raise ArithmeticError("old section base map did not invert to u")
        x_section = Ku(evaluate_rational(Xold, t_section))
        y_section = Ku(evaluate_rational(Yold, t_section))
        Xnew, Ynew, Wnew = point_on_child(t_section, x_section, y_section)
        return old_base_map, t_section, Xnew, Ynew, Wnew

    # Sixteen honest old sections together with one rational old bisection
    # form a unimodular basis in the new frame.  The bisection supplies the
    # missing index-two glue class that the degree-one old sections alone miss.
    selected_old_vectors_by_label = {
        "norm12-orbit-11952": [
        [0, -1, 0, 0, 0, 0, 0, -1, 1, 1, 0, -1, 1, 0, -1, 0, 1],
        [0, 0, 3, 1, 1, 2, -2, 0, 2, 0, 2, 0, 1, -2, 1, -2, -4],
        [-1, -1, 3, 0, -1, 2, -1, 1, 2, 1, 1, 0, 1, -1, 2, -2, -3],
        [0, -1, 1, 0, 0, 0, 0, 1, 2, 1, 0, -1, 1, -1, 1, -1, -1],
        [0, 0, 1, 0, 1, 1, -1, -1, 2, 1, 1, -1, 1, -1, 1, -2, -2],
        [0, 0, 2, 0, 0, 1, -2, 0, 2, 0, 1, 0, 1, -1, 1, -2, -2],
        [0, 0, 2, 1, 1, 2, -1, 0, 1, 0, 2, 0, 0, -1, 1, -2, -4],
        [0, 0, 2, 0, 1, 1, -1, -1, 1, 0, 1, 0, 1, -1, 1, -2, -2],
        [-1, -1, 2, 0, 0, 1, -1, 1, 2, 1, 1, 0, 1, -1, 1, -2, -2],
        [-1, 0, 3, 0, -1, 1, -2, 0, 2, 0, 1, 1, 1, -1, 2, -2, -2],
        [0, 1, 2, 1, 1, 1, -2, -2, 2, -1, 2, 0, 1, -1, 1, -2, -3],
        [0, 0, 2, 0, 1, 0, -1, 0, 2, 0, 1, 0, 1, -1, 1, -2, -2],
        [0, 0, 1, -1, 0, 1, -1, 0, 1, 1, 0, 0, 1, -1, 1, -1, -1],
        [0, 0, 0, 1, 0, 0, 0, -1, 2, 0, 1, 0, 1, -1, 0, -1, -2],
        [-1, 0, 2, -1, 1, 2, -1, 0, 1, 1, 1, 0, 1, -1, 1, -1, -2],
        [-1, 1, 1, -1, 1, 1, -1, -1, 1, 1, 1, -1, 1, -1, 1, 0, -1],
        ],
        "norm12-orbit-103b2": [
            [1, 0, -1, 1, 0, -1, 0, -1, 1, 0, 0, 0, 0, 0, -1, 0, 0],
            [0, 0, 1, -1, -1, 0, -1, 0, 0, 0, -1, 1, 1, 0, 0, 0, 1],
            [1, 0, 1, 1, 1, 2, -1, 0, 0, 0, 1, 0, 0, -1, 0, -1, -3],
            [-1, 0, 1, 0, -1, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 1, 0],
            [0, 0, 0, 1, 1, 1, 0, -2, 0, 0, 1, 0, 0, 0, -1, 0, -1],
            [0, 0, 1, 1, 0, 2, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, -2],
            [1, 0, 0, 2, 1, 0, 0, -1, 1, 0, 1, 0, 0, -1, -1, 0, -2],
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
            [0, -1, 1, 0, -1, 1, 0, 1, 0, 1, 0, 0, 1, -1, 0, 0, -1],
            [-1, -1, 1, -1, -1, 2, 0, 1, -1, 2, 0, 0, 0, 0, 0, 1, 0],
            [1, 0, 1, 1, 0, 1, -1, 0, 1, 0, 1, 1, 0, -1, 0, -1, -3],
            [0, -1, -1, 0, -1, 0, 1, 1, -1, 1, -1, 0, 0, 0, -1, 1, 1],
            [0, 0, 1, 1, 0, 1, -1, 0, 0, 0, 1, 0, -1, 0, 0, 0, -1],
            [1, -1, 0, 1, -1, -1, 0, 1, 1, 0, -1, 1, 1, -1, -1, 0, 0],
            [0, -1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0],
        ],
        # These sixteen norm-at-most-eight old sections are exactly the first
        # deterministic unimodular basis found in the degree-one shell for
        # orbit-08f72.  Together with orbit-04eb3 below their transported
        # child-frame rows have determinant one.
        "norm12-orbit-08f72": [
            [0, -1, 2, 0, -2, 1, 0, 3, 0, 0, 0, 1, -1, 0, 2, -1, -2],
            [0, -1, 0, -1, -1, 0, 1, 2, 0, 1, -1, 0, 0, 0, 1, 0, 0],
            [0, -1, 1, 0, -2, 1, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, -1],
            [0, 0, 2, 0, -1, 2, -1, 2, -1, -1, 0, 1, -1, 0, 2, -1, -2],
            [0, -1, 1, -1, -2, 0, 1, 3, 0, 0, -2, 1, 0, 0, 2, 0, 0],
            [0, 0, 2, -1, -2, 1, -1, 2, 0, -1, -1, 2, 0, 0, 2, -1, -1],
            [-1, -1, 2, -1, -2, 1, 0, 3, 0, 0, -1, 2, 0, 0, 2, -1, -1],
            [0, -1, 1, -1, -3, 0, 0, 3, 0, 0, -2, 2, 0, 0, 2, -1, 0],
            [0, -1, 0, -1, -2, -1, 1, 3, 0, 0, -2, 1, 0, 0, 1, 0, 1],
            [0, -1, 2, 0, -3, 1, 0, 3, 0, -1, -1, 2, 0, 0, 2, -1, -2],
            [0, -1, 0, 0, -1, 1, 1, 2, 0, 0, -1, 0, 0, 0, 1, 0, -1],
            [1, -1, 0, 0, -2, -1, 1, 3, 0, -1, -2, 1, 0, 0, 1, 0, 0],
            [1, 0, 0, 0, -1, 0, 0, 1, 0, -1, -1, 1, 0, 0, 1, -1, -1],
            [0, 0, 0, 0, -2, 0, 0, 1, 0, -1, -1, 1, 0, 0, 1, 0, 0],
            [-1, 0, 1, -1, -2, 0, 0, 2, 0, 0, -1, 1, 0, 0, 2, 0, 0],
            [0, -1, 1, 0, -2, 0, 0, 2, 0, -1, -1, 2, -1, 1, 1, -1, 0],
        ],
        # A deterministic degree-one-shell basis for orbit-08ab4.  These
        # sixteen old sections together with orbit-1ebca below have child
        # frame coordinate determinant one.
        "norm12-orbit-08ab4": [
            [0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, -1, 0, 0, 0, 1, 0],
            [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, -1],
            [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -1, 0, -1, 0, 1, 0],
            [0, 0, 0, 0, 1, 0, 1, 0, -1, 0, 0, -1, 0, 0, 0, 1, 0],
            [0, 0, 1, 0, 1, 0, 0, -1, 0, 0, 1, -1, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 1, 0, -1, -1, 0, 1, 0, 0, 0, 0, 0, -1],
            [0, -1, 0, 0, -1, 0, 1, 1, 0, 1, 0, -1, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, -1, 0, -1, 0, 0, -1],
            [0, 0, 1, 0, 1, 1, 0, 0, -1, 0, 1, 0, -1, 0, 0, 0, -1],
            [0, -1, 0, 0, 0, 0, 1, 0, -1, 1, 0, -1, 0, 0, -1, 1, 1],
            [0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, -1, -1, 0, 0, 1, -1],
            [1, -1, 0, 1, 0, -1, 1, 0, 1, 0, 0, 0, 0, 0, 0, -1, -1],
            [1, 0, -1, 0, 1, 0, 1, 0, 0, 1, 0, -1, 0, -1, 0, 0, -1],
            [1, 0, 1, 1, 0, 0, 0, 0, -1, -1, 1, 0, -1, 0, 0, 0, -1],
            [0, -1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, -1, 0, -1, -1],
            [0, -1, 1, 0, 1, 1, 0, 0, 0, 1, 1, -1, 0, 0, 0, -1, -1],
        ],
    }
    selected_old_vectors = selected_old_vectors_by_label[args.source_label]
    basis_change = matrix(ZZ, target_payload["pinned_identification"]["basis_change_matrix"])
    transport_inverse = transport.inverse()
    section_records = []
    new_mw_rows = []
    for index, entries in enumerate(selected_old_vectors):
        old_mw = vector(ZZ, entries)
        old_height = int(old_mw * pinned * old_mw)
        old_class = vector(ZZ, [(old_height - 2) // 2, 1] + list(old_mw))
        new_class = old_class * transport_inverse
        if new_class[1] != 1 or any(value not in ZZ for value in new_class):
            raise ArithmeticError("old section did not become an integral new section")
        new_mw = vector(ZZ, new_class[2:])
        published_mw = old_mw * basis_change.transpose()
        old_point = sum(
            (coefficient * point for coefficient, point in zip(published_mw, published_basis)),
            Eold(0),
        )
        old_base_map, t_section, Xnew, Ynew, Wnew = map_old_section(old_point)
        assert int(new_mw * frame * new_mw) in (4, 6, 8)
        new_mw_rows.append(new_mw)
        section_records.append(
            {
                "basis_index": index,
                "source": "published-R17 old section",
                "source_pinned_MW_coordinates": list(map(int, old_mw)),
                "source_published_MW_coordinates": list(map(int, published_mw)),
                "source_old_height": old_height,
                "new_frame_coordinates": list(map(int, new_mw)),
                "new_height": int(new_mw * frame * new_mw),
                "old_base_map_u_of_t": rational_function_record(old_base_map),
                "quartic_t_coordinate": rational_function_record(t_section),
                "X": rational_function_record(Xnew),
                "Y": rational_function_record(Ynew),
                "equation_verified": True,
            }
        )

    full_bisections = json.loads(BISECTIONS.read_text())
    glue_labels_by_source = {
        "norm12-orbit-11952": ["orbit-0adf9"],
        "norm12-orbit-103b2": ["orbit-1d5f2", "orbit-0abc2"],
        "norm12-orbit-08f72": ["orbit-04eb3"],
        "norm12-orbit-08ab4": ["orbit-1ebca"],
    }
    glue_labels = glue_labels_by_source[args.source_label]
    for glue_offset, glue_label in enumerate(glue_labels):
        glue_record = next(
            record for record in full_bisections["bisections"]
            if record["label"] == glue_label
        )
        glue_w = vector(ZZ, glue_record["pinned_rank17_w"])
        glue_class = vector(ZZ, [2, 2] + list(glue_w))
        glue_new_class = glue_class * transport_inverse
        if glue_new_class[1] != 1 or any(value not in ZZ for value in glue_new_class):
            raise ArithmeticError("glue bisection did not become an integral child section")
        glue_new_mw = vector(ZZ, glue_new_class[2:])

        glue_q = Rt([QQ(value) for value in glue_record["residual_chord"]["q_coefficients"]])
        Rconic = PolynomialRing(QQ, names=("T", "S", "Z"))
        Tconic, Sconic, Zconic = Rconic.gens()
        conic = Conic(
            Sconic**2
            - (glue_q[2] * Tconic**2 + glue_q[1] * Tconic * Zconic + glue_q[0] * Zconic**2)
        )
        has_point, conic_point = conic.has_rational_point(point=True)
        if not has_point:
            raise ArithmeticError("glue bisection cover has no rational point")
        parametrization, _ = conic.parametrization(point=conic_point)
        parameter_polys = parametrization.defining_polynomials()
        Rr = PolynomialRing(QQ, "r")
        Kr = Rr.fraction_field()
        r = Rr.gen()
        Tparam, Sparam, Zparam = [Rr(poly(r, 1)) for poly in parameter_polys]
        t_of_r = Kr(Tparam / Zparam)
        s_of_r = Kr(Sparam / Zparam)
        assert s_of_r**2 == evaluate_polynomial(glue_q, t_of_r)
        lifted = glue_record["lifted_section"]
        gx0 = Rt([QQ(value) for value in lifted["x0_coefficients"]])
        gx1 = Rt([QQ(value) for value in lifted["x1_coefficients"]])
        gy0 = Rt([QQ(value) for value in lifted["y0_coefficients"]])
        gy1 = Rt([QQ(value) for value in lifted["y1_coefficients"]])
        x_of_r = evaluate_polynomial(gx0, t_of_r) + evaluate_polynomial(gx1, t_of_r) * s_of_r
        y_of_r = evaluate_polynomial(gy0, t_of_r) + evaluate_polynomial(gy1, t_of_r) * s_of_r
        assert y_of_r**2 == evaluate_polynomial(Aold, t_of_r) * x_of_r + x_of_r**3 + evaluate_polynomial(Bold, t_of_r)
        L0r = evaluate_polynomial(a0, t_of_r) * (
            x_of_r * evaluate_polynomial(h, t_of_r) ** 2 - evaluate_polynomial(Nx, t_of_r)
        ) + evaluate_polynomial(b0, t_of_r) * (
            y_of_r * evaluate_polynomial(h, t_of_r) ** 3 + evaluate_polynomial(Ny, t_of_r)
        )
        L1r = evaluate_polynomial(a1, t_of_r) * (
            x_of_r * evaluate_polynomial(h, t_of_r) ** 2 - evaluate_polynomial(Nx, t_of_r)
        ) + evaluate_polynomial(b1, t_of_r) * (
            y_of_r * evaluate_polynomial(h, t_of_r) ** 3 + evaluate_polynomial(Ny, t_of_r)
        )
        u_of_r = Kr(L1r / L0r)
        r_of_u = Ku(invert_mobius(u_of_r, u))
        glue_t = Ku(evaluate_rational(t_of_r, r_of_u))
        glue_x = Ku(evaluate_rational(x_of_r, r_of_u))
        glue_y = Ku(evaluate_rational(y_of_r, r_of_u))
        glue_X, glue_Y, glue_W = point_on_child(glue_t, glue_x, glue_y)
        new_mw_rows.append(glue_new_mw)
        section_records.append(
            {
                "basis_index": len(selected_old_vectors) + glue_offset,
                "source": f"published-R17 rational bisection {glue_label}",
                "source_curve_class_in_U_plus_R17_minus": list(map(int, glue_class)),
                "source_cover_point_T_S_Z": [rational_text(value) for value in conic_point],
                "new_frame_coordinates": list(map(int, glue_new_mw)),
                "new_height": int(glue_new_mw * frame * glue_new_mw),
                "cover_base_map_u_of_r": rational_function_record(u_of_r),
                "quartic_t_coordinate": rational_function_record(glue_t),
                "X": rational_function_record(glue_X),
                "Y": rational_function_record(glue_Y),
                "equation_verified": True,
            }
        )

    section_coordinate_matrix = matrix(ZZ, new_mw_rows)
    if abs(section_coordinate_matrix.det()) != 1:
        raise ArithmeticError("the recovered child sections do not form a saturated basis")
    section_height_gram = section_coordinate_matrix * frame * section_coordinate_matrix.transpose()
    if section_height_gram.det() != 948 or int(pari(section_height_gram).qfminim(2)[0]):
        raise ArithmeticError("recovered section height Gram failed the rootless frame checks")

    frame_certificate = {
        "transport_rows_D_D_plus_O_complement": rows(transport),
        "transport_determinant": int(transport.det()),
        "frame_gram": rows(frame),
        "frame_gram_sha256": matrix_digest(frame),
        "determinant": int(frame.det()),
        "roots_of_norm_two": 0,
    }
    if is_alternate_target:
        frame_certificate.update({
            "integral_isometry_to_alternate_Q80": rows(isometry),
            "isometric_to_alternate_Q80": True,
            "isometric_to_published_R17": False,
        })
    else:
        frame_certificate.update({
            "integral_isometry_to_published_R17": rows(isometry),
            "isometric_to_published_R17": True,
            "isometric_to_alternate_Q80": False,
        })
    input_paths = [MODEL, SECTIONS, PINNED, TARGET, SPLITTING, CLASSIFICATION, ALTERNATE, BISECTIONS]
    if not is_alternate_target:
        input_paths.append(FRAME_103B2)
    result = {
        "schema": {
            "norm12-orbit-11952": "elkies-k3.r17-norm12-orbit11952-direct-fibration.v1",
            "norm12-orbit-103b2": "elkies-k3.r17-norm12-orbit103b2-direct-fibration.v1",
            "norm12-orbit-08f72": "elkies-k3.r17-norm12-orbit08f72-direct-fibration.v1",
            "norm12-orbit-08ab4": "elkies-k3.r17-norm12-orbit08ab4-direct-fibration.v1",
        }[args.source_label],
        "status": "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS",
        "divisor": {
            "label": args.source_label,
            "pinned_trace_vector_w": list(map(int, w)),
            "class_D_in_U_plus_R17_minus": list(map(int, fibre)),
            "identity": "D=O_old+P_w-F_old=(3,2,w)",
            "old_fibre_degree": 2,
            "old_zero_degree": 1,
            "shared_zero": True,
        },
        "riemann_roch": {
            "method": "Brandhorst-Elkies degree-two neighbor, Proposition 2.17",
            "section_form": "a(t)+b(t)*(y+y(P_w))/(x-x(P_w))",
            "degree_bounds": {"a": 7, "b": 1},
            "congruence": "a*Nx-b*Ny == 0 mod h^2",
            "constraint_matrix_rank": int(rr_matrix.rank()),
            "kernel_dimension": int(rr_kernel.nrows()),
            "kernel_rows_a0_through_a7_b0_b1": [
                [rational_text(value) for value in row] for row in rr_kernel.rows()
            ],
            "pencil_coordinate": "u=L1/L0",
        },
        "genus_one_model": {
            "equation": "W^2=q(t,u)",
            "q_coefficients_in_t_low_to_high": [
                rational_function_record(Ku(quartic[i])) for i in range(5)
            ],
            "distinguished_point_from_old_zero": {
                "t0": rational_function_record(t0),
                "W0": rational_function_record(v0),
                "identity_verified": True,
            },
            "chord_slope": "m=(a1-u*a0)/((u*b0-b1)*h)",
            "radical_square_factor_degree": int(square_factor.degree()),
            "radical_identity_verified": True,
            "quartic_squarefree_over_Q_of_u": True,
        },
        "weierstrass_model": {
            "coordinate": "u",
            "equation": "Y^2=X^3+A(u)*X+B(u)",
            "A_coefficients_low_to_high": polynomial_text(Achild),
            "B_coefficients_low_to_high": polynomial_text(Bchild),
            "gauge": rational_function_record(Ku(gauge)),
            "degrees_A_B_Delta": [int(Achild.degree()), int(Bchild.degree()), int(Delta_child.degree())],
            "discriminant_coefficients_low_to_high": polynomial_text(Delta_child),
            "discriminant_irreducible_over_Q": True,
            "fibre_configuration": "24 I1",
            "infinity_orders_c4_c6_Delta": [0, 0, 0],
            "pointed_quartic_map_verified": True,
        },
        "frame_certificate": frame_certificate,
        "sections": {
            "status": "PASS_EXACT_SATURATED_RANK17_BASIS",
            "rank": 17,
            "basis_source_profile": {
                "old_sections": len(selected_old_vectors),
                "old_rational_bisections": len(glue_labels),
            },
            "glue_source": glue_labels[0] if len(glue_labels) == 1 else glue_labels,
            "coordinate_matrix_in_compiled_frame": rows(section_coordinate_matrix),
            "coordinate_matrix_determinant": int(section_coordinate_matrix.det()),
            "height_gram": rows(section_height_gram),
            "height_gram_determinant": int(section_height_gram.det()),
            "roots_of_norm_two": 0,
            "records": section_records,
        },
        "inputs": {
            relative(path): digest(path)
            for path in input_paths
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "pari_version": ".".join(map(str, pari.version())),
            "required_features": ["Sage exact function fields", "PARI qfminim", "PARI qfisom"],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/compile_r17_norm12_orbit11952_qq.sage"
            + (
                ""
                if args.source_label == "norm12-orbit-11952"
                else f" --source-label {args.source_label}"
            )
        ),
        "proof_boundary": (
            (
                "This exact replay constructs H^0(X,O(D)), the quartic pencil, its pointed Jacobian, "
                "the squarefree degree-24 discriminant, and the primitive rootless alternate-Q80 frame. "
                "It also transports sixteen old sections and the orbit-0adf9 rational bisection to "
                "an explicit saturated rank-17 section basis on the new equation."
            )
            if args.source_label == "norm12-orbit-11952"
            else (
                "This exact replay constructs H^0(X,O(D)), the quartic pencil, its pointed Jacobian, "
                "the squarefree degree-24 discriminant, and the primitive rootless alternate-Q80 frame. "
                "It also transports sixteen old sections and the orbit-04eb3 rational bisection to "
                "an explicit saturated rank-17 section basis on the new equation."
            )
            if args.source_label == "norm12-orbit-08f72"
            else (
                "This exact replay constructs H^0(X,O(D)), the quartic pencil, its pointed Jacobian, "
                "the squarefree degree-24 discriminant, and the primitive rootless alternate-Q80 frame. "
                "It also transports sixteen old sections and the orbit-1ebca rational bisection to "
                "an explicit saturated rank-17 section basis on the new equation."
            )
            if args.source_label == "norm12-orbit-08ab4"
            else (
                "This exact replay constructs H^0(X,O(D)), the quartic pencil, its pointed Jacobian, "
                "the squarefree degree-24 discriminant, and the primitive rootless published-R17 frame. "
                "It transports fifteen old sections and both degree-one old rational bisections to "
                "an explicit saturated rank-17 section basis on the hidden equation."
            )
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored norm-twelve direct-fibration artifact differs from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "R17NORM12DIRECT|label={}|h0=2|quartic=4|A=8|B=12|Delta=24|fibres=24I1|"
        "frame={}|roots=0|sections=17-saturated|output={}".format(
            args.source_label, expected_frame_class, relative(output)
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
