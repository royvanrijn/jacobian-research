#!/usr/bin/env sage -python
"""Construct the exact degree-three S5 trace on the A11 equation.

Restrict the certified D12-to-A11 pencil to the exact shell curve S5.  Its
three conjugate A11 points are represented in the cubic algebra QQ(T)[V]/(m).
The one-dimensional kernel of their evaluations on

    L(4 O) = <1, x, x^2, y>

cuts out one residual point.  The negative residual point is the Abel--Jacobi
trace.  Thus the calculation uses only univariate rational-function arithmetic
and 3-by-4 linear algebra; it uses no section ansatz, Groebner basis, or Hensel
lift.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--marked-reduction",
    type=Path,
    action="append",
    default=None,
    help="marked modular zero artifact; repeat to require several good primes",
)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-a11-s5-trace-section-qq.json",
)
args = parser.parse_args()

A11 = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
PARENT = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
CANDIDATES = LOCAL / "q24-orbit42-exact-section-candidates-qq.json"
ZERO = LOCAL / "q24-orbit42-rational-zero-pole-sections-qq.json"
MODEL = LOCAL / "q24-orbit42-zero-pole-seeds-mod-43.json"
POINTED = LOCAL / "q24-a11-pointed-opposite-section-qq.json"
MARKED_PATHS = tuple(
    path.resolve() for path in (
        args.marked_reduction
        or [
            LOCAL / "q24-a11-pinned-zero-section-mod100003.json",
        ]
    )
)
INPUTS = (A11, PARENT, CANDIDATES, ZERO, MODEL, POINTED) + MARKED_PATHS
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

a11 = json.loads(A11.read_text())
parent = json.loads(PARENT.read_text())
candidates = json.loads(CANDIDATES.read_text())
zero = json.loads(ZERO.read_text())
model = json.loads(MODEL.read_text())["exact_model"]
pointed = json.loads(POINTED.read_text())
marked_rows = [json.loads(path.read_text()) for path in MARKED_PATHS]
assert a11["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert parent["status"] == "PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR"
assert candidates["status"] == "PASS_EXACT_Q42_ORBIT42_SECTION_CANDIDATES_QQ"
assert zero["status"] == "PASS_EXACT_Q42_RATIONAL_ZERO_POLE_SECTIONS_QQ"
assert pointed["status"] == "PASS_EXACT_A11_POINTED_OPPOSITE_SECTION_QQ"
assert all(
    row["status"] == "PASS_Q24_A11_PINNED_ZERO_SECTION_RECONSTRUCTION_MODP"
    for row in marked_rows
)

UQ = PolynomialRing(QQ, "u")
KU = UQ.fraction_field()
VQ = PolynomialRing(QQ, "V")
V = VQ.gen()
KV = VQ.fraction_field()
TQ = PolynomialRing(QQ, "T")
T = TQ.gen()
KT = TQ.fraction_field()


def evaluate_polynomial(poly, argument):
    return sum(poly[index] * argument**index for index in range(poly.degree() + 1))


def evaluate_u(value, argument):
    value = KU(value)
    return KV(evaluate_polynomial(VQ(value.numerator()), argument)) / KV(
        evaluate_polynomial(VQ(value.denominator()), argument)
    )


def evaluate_T(value, argument):
    """Compose rational functions with one final normalization.

    Homogeneous evaluation avoids the repeated polynomial gcds incurred by a
    term-by-term sum in the rational-function field.
    """

    value = KT(value)
    argument = KV(argument)
    N = VQ(argument.numerator())
    D = VQ(argument.denominator())

    def homogeneous(poly):
        poly = TQ(poly)
        degree = int(poly.degree())
        if degree < 0:
            return VQ.zero(), 0
        answer = VQ.zero()
        Npower = VQ.one()
        Dpowers = [VQ.one()]
        for unused in range(degree):
            Dpowers.append(Dpowers[-1] * D)
        for index, coefficient in enumerate(poly.list()):
            answer += QQ(coefficient) * Npower * Dpowers[degree - index]
            Npower *= N
        return answer, degree

    numerator, numerator_degree = homogeneous(value.numerator())
    denominator, denominator_degree = homogeneous(value.denominator())
    return KV(numerator * D**denominator_degree) / KV(
        denominator * D**numerator_degree
    )


def rational_square_root(value, ring):
    value = ring.fraction_field()(value)
    numerator = ring(value.numerator())
    denominator = ring(value.denominator())
    answer = value.parent().one()
    for polynomial, direction in ((numerator, 1), (denominator, -1)):
        if not polynomial:
            return value.parent().zero()
        leading = QQ(polynomial.leading_coefficient())
        if not leading.is_square():
            raise ArithmeticError("rational-function square has nonsquare leading coefficient")
        root = value.parent()(leading.sqrt())
        for factor, multiplicity in (polynomial / leading).factor():
            if int(multiplicity) % 2:
                raise ArithmeticError("restricted quartic value is not a square")
            root *= value.parent()(factor.monic()) ** (int(multiplicity) // 2)
        answer = answer * root if direction == 1 else answer / root
    if answer**2 != value:
        raise ArithmeticError("rational square-root reconstruction failed")
    return answer


def power_root(poly, exponent):
    poly = TQ(poly)
    leading = QQ(poly.leading_coefficient())
    if exponent % 2 == 0 and not leading.is_square():
        raise ArithmeticError("denominator leading coefficient is not a square")
    answer = TQ.one()
    for factor, multiplicity in (poly / leading).factor():
        if int(multiplicity) % exponent:
            raise ArithmeticError("section denominator is not a perfect power")
        answer *= factor.monic() ** (int(multiplicity) // exponent)
    return answer.monic()


def normalized_section(point, A, B):
    x, y = map(KT, point.xy())
    Zx = power_root(x.denominator(), 2)
    Zy = power_root(y.denominator(), 3)
    if Zx != Zy:
        raise ArithmeticError("x and y denominators give different section Z")
    Z = Zx
    X = TQ(x * Z**2)
    Y = TQ(y * Z**3)
    if Y**2 != X**3 + A * X * Z**4 + B * Z**6:
        raise ArithmeticError("normalized trace missed the exact A11 equation")
    return X, Y, Z


u_of_V = KV(str(a11["coordinate_change"]["u_of_V"]))
x_scale = KV(str(a11["coordinate_change"]["x_scale"]))
y_scale = KV(str(a11["coordinate_change"]["y_scale"]))
A12 = VQ([QQ(value) for value in parent["child"]["minimal_A_coefficients_low_to_high"]])
B12 = VQ([QQ(value) for value in parent["child"]["minimal_B_coefficients_low_to_high"]])

selected = candidates["candidates"][a11["selected_section"]["candidate_index"]]
X_parent_u = UQ([QQ(value) for value in selected["X_coefficients_low_to_high"]])
Y_parent_u = UQ([QQ(value) for value in selected["Y_coefficients_low_to_high"]])
Z_parent_u = UQ([QQ(value) for value in selected["Z_coefficients_low_to_high"]])
x_parent = x_scale * evaluate_u(KU(X_parent_u) / KU(Z_parent_u**2), u_of_V)
y_parent = y_scale * evaluate_u(KU(Y_parent_u) / KU(Z_parent_u**3), u_of_V)
if y_parent**2 != x_parent**3 + KV(A12) * x_parent + KV(B12):
    raise ArithmeticError("selected parent section missed D12")

denominator = VQ(x_parent.denominator()).monic()
Z_parent = VQ.one()
for factor, multiplicity in denominator.factor():
    if int(multiplicity) % 2:
        raise ArithmeticError("selected parent denominator is not a square")
    Z_parent *= factor.monic() ** (int(multiplicity) // 2)
X_parent = VQ(x_parent * Z_parent**2)
Y_parent = VQ(y_parent * Z_parent**3)
alpha_collision = QQ(model["I8star_root"])
collision_modulus = Z_parent**2
X_inverse = X_parent.inverse_mod(collision_modulus)
rr_pairs = []
for BB in (VQ.one(), V):
    AA = VQ((BB * Y_parent * X_inverse) % collision_modulus)
    AA -= AA(alpha_collision) / Z_parent(alpha_collision) ** 2 * Z_parent**2
    rr_pairs.append((KV(AA) / KV(Z_parent**2), KV(BB) / KV(Z_parent)))
(a0, b0), (a1, b1) = rr_pairs

s5 = zero["sections"][5]
x_s5_u = UQ([QQ(value) for value in s5["x_coefficients_low_to_high"]])
y_s5_u = UQ([QQ(value) for value in s5["y_coefficients_low_to_high"]])
x_s5 = x_scale * evaluate_u(KU(x_s5_u), u_of_V)
y_s5 = y_scale * evaluate_u(KU(y_s5_u), u_of_V)
if y_s5**2 != x_s5**3 + KV(A12) * x_s5 + KV(B12):
    raise ArithmeticError("exact S5 shell missed D12")
chord = (y_s5 + y_parent) / (x_s5 - x_parent)
Tmap = KV((a1 + b1 * chord) / (a0 + b0 * chord))
if max(Tmap.numerator().degree(), Tmap.denominator().degree()) != 3:
    raise ArithmeticError("S5 does not have new-base degree three")

quartic_coefficients = [KT(TQ(value)) for value in a11["quartic"]["coefficients_in_T_low_to_high"]]
quartic_on_s5 = sum(evaluate_T(value, Tmap) * V**index for index, value in enumerate(quartic_coefficients))
print("A11S5TRACE|stage=restricted_quartic|status=PASS", flush=True)
Wplus = rational_square_root(quartic_on_s5, VQ)
print("A11S5TRACE|stage=exact_ordinate|status=PASS", flush=True)


def exact_orientation_from_marked_reduction(marked):
    p = ZZ(marked["prime"])
    F = GF(p)
    R = PolynomialRing(F, "V")
    K = R.fraction_field()

    def reduce_rational(value):
        value = KV(value)
        numerator = R(
            [F(item.numerator()) / F(item.denominator()) for item in value.numerator().list()]
        )
        denominator = R(
            [F(item.numerator()) / F(item.denominator()) for item in value.denominator().list()]
        )
        return K(numerator) / K(denominator)

    def reduce_T_rational(value):
        value = KT(value)
        numerator = R(
            [F(item.numerator()) / F(item.denominator()) for item in value.numerator().list()]
        )
        denominator = R(
            [F(item.numerator()) / F(item.denominator()) for item in value.denominator().list()]
        )
        return K(numerator) / K(denominator)

    reduced_quartic = reduce_rational(quartic_on_s5)
    reduced_exact_root = reduce_rational(Wplus)
    modular_sampler_root = reduced_quartic.sqrt()
    if reduced_exact_root == modular_sampler_root:
        root_orientation = 1
    elif reduced_exact_root == -modular_sampler_root:
        root_orientation = -1
    else:
        raise ArithmeticError("exact S5 ordinate reduction disagrees with the modular square roots")
    modular_branch = int(marked["selected_orientation"]["S5_branch"])
    trace_signs = marked["selected_orientation"].get("trace_signs", {})
    equation_sign = int(trace_signs.get("S5", 1))
    reduced_anchor = reduce_T_rational(q_anchor)
    modular_anchor = reduce_T_rational(q_squared).sqrt()
    if reduced_anchor == modular_anchor:
        anchor_orientation = 1
    elif reduced_anchor == -modular_anchor:
        anchor_orientation = -1
    else:
        raise ArithmeticError("exact pointed anchor reduction disagrees with modular roots")
    return modular_branch * root_orientation, equation_sign, anchor_orientation

i8 = next(row for row in parent["child"]["finite_fibres"] if row["kodaira"] == "I8*")
i8_factor = VQ(str(i8["factor"]))
alpha = -QQ(i8_factor[0]) / QQ(i8_factor[1])
q_squared = sum(value * KT(alpha) ** index for index, value in enumerate(quartic_coefficients))
q_anchor = rational_square_root(q_squared, TQ)
e, d0, c0, b0q, a0q = quartic_coefficients
anchor_a = a0q
anchor_b = b0q + 4 * alpha * a0q
anchor_c = c0 + 3 * alpha * b0q + 6 * alpha**2 * a0q
anchor_d = d0 + 2 * alpha * c0 + 3 * alpha**2 * b0q + 4 * alpha**3 * a0q

marked_orientations = tuple(
    exact_orientation_from_marked_reduction(marked) for marked in marked_rows
)
selected_quartic_sign, selected_equation_sign, selected_anchor_orientation = (
    marked_orientations[0]
)
origin_translation_coefficient = 0 if selected_anchor_orientation == 1 else -3
print(
    f"A11S5TRACE|stage=orientation|quartic_sign={selected_quartic_sign}|"
    f"equation_sign={selected_equation_sign}|anchor_sign={selected_anchor_orientation}|"
    f"origin_translation={origin_translation_coefficient}Q|"
    f"raw_prime_orientations={marked_orientations}|"
    f"primes={tuple(int(row['prime']) for row in marked_rows)}|"
    "status=PASS",
    flush=True,
)

A11poly = TQ([QQ(value) for value in a11["child"]["minimal_A_coefficients_low_to_high"]])
B11poly = TQ([QQ(value) for value in a11["child"]["minimal_B_coefficients_low_to_high"]])
E11 = EllipticCurve(KT, [0, 0, 0, KT(A11poly), KT(B11poly)])


def exact_section_point(record):
    X = TQ([QQ(value) for value in record["X_coefficients_low_to_high"]])
    Y = TQ([QQ(value) for value in record["Y_coefficients_low_to_high"]])
    Z = TQ([QQ(value) for value in record["Z_coefficients_low_to_high"]])
    return E11(KT(X) / KT(Z**2), KT(Y) / KT(Z**3))


Qpoint = exact_section_point(pointed["section"])


SV = PolynomialRing(KT, "v")
v = SV.gen()


def lift_V(poly):
    return SV([KT(QQ(value)) for value in VQ(poly).list()])


modulus = lift_V(Tmap.numerator()) - T * lift_V(Tmap.denominator())
modulus = modulus.monic()
if modulus.degree() != 3 or modulus.gcd(modulus.derivative()).degree() != 0:
    raise ArithmeticError("S5 cubic algebra is not generically etale")


def algebra_inverse(value):
    """Invert in the cubic algebra by a fraction-free 3-by-3 adjugate."""

    value = SV(value) % modulus
    multiplication = matrix(
        KT,
        3,
        3,
        lambda row, column: SV((value * v**column) % modulus)[row],
    )
    common_denominator = TQ.one()
    for entry in multiplication.list():
        common_denominator = common_denominator.lcm(TQ(entry.denominator()))
    cleared = matrix(
        TQ,
        3,
        3,
        [TQ(entry * common_denominator) for entry in multiplication.list()],
    )
    determinant = TQ(cleared.det())
    if not determinant:
        raise ArithmeticError("element is not invertible in the cubic algebra")
    c00, c01, c02 = cleared[0]
    c10, c11, c12 = cleared[1]
    c20, c21, c22 = cleared[2]
    first_inverse_column = (
        c11 * c22 - c12 * c21,
        c12 * c20 - c10 * c22,
        c10 * c21 - c11 * c20,
    )
    print(
        "A11S5TRACE|stage=fraction_free_inverse|common_den_degree={}|"
        "det_degree={}|numerator_degrees={}|det_max_bits={}|numerator_max_bits={}|"
        "status=BEGIN_NORMALIZATION".format(
            common_denominator.degree(),
            determinant.degree(),
            tuple(item.degree() for item in first_inverse_column),
            max(abs(coefficient.numerator()).nbits() for coefficient in determinant.list()),
            max(
                abs(coefficient.numerator()).nbits()
                for item in first_inverse_column
                for coefficient in item.list()
            ),
        ),
        flush=True,
    )
    diagnostic_field = GF(100003)
    diagnostic_ring = PolynomialRing(diagnostic_field, "t")

    def diagnostic_reduction(poly):
        return diagnostic_ring(
            [
                diagnostic_field(coefficient.numerator())
                / diagnostic_field(coefficient.denominator())
                for coefficient in poly.list()
            ]
        )

    modular_determinant = diagnostic_reduction(determinant)
    modular_gcds = tuple(
        diagnostic_reduction(numerator).gcd(modular_determinant).monic()
        for numerator in first_inverse_column
    )
    modular_gcd_degrees = tuple(item.degree() for item in modular_gcds)
    modular_common_denominator = diagnostic_reduction(common_denominator).monic()
    print(
        f"A11S5TRACE|stage=fraction_free_inverse_modular_gcd|"
        f"degrees={modular_gcd_degrees}|"
        f"equals_common_denominator={tuple(item == modular_common_denominator for item in modular_gcds)}|"
        f"gcds={tuple(tuple(map(int, item.list())) for item in modular_gcds)}|"
        f"common_denominator={tuple(map(int, modular_common_denominator.list()))}|status=PASS",
        flush=True,
    )

    reduced_determinant, remainder = determinant.quo_rem(common_denominator)
    if remainder:
        raise ArithmeticError("cleared cubic determinant lost its common denominator factor")
    reduced_numerators = []
    for numerator in first_inverse_column:
        reduced_numerator, remainder = numerator.quo_rem(common_denominator)
        if remainder:
            raise ArithmeticError("cubic adjugate lost its common denominator factor")
        reduced_numerators.append(reduced_numerator)

    def normalized_coefficient(reduced_numerator):
        return KT(common_denominator * reduced_numerator) / KT(reduced_determinant)

    answer = SV(
        [normalized_coefficient(numerator) for numerator in reduced_numerators]
    )
    if (value * answer) % modulus != 1:
        raise ArithmeticError("fraction-free cubic inverse check failed")
    return answer


def residue(value):
    value = KV(value)
    # Reduce first: the stored ordinate may have a large ambient V-degree,
    # while its class in the cubic algebra has degree at most two.
    numerator = lift_V(value.numerator()) % modulus
    denominator = lift_V(value.denominator()) % modulus
    print(
        "A11S5TRACE|stage=ordinate_residue|num_degree={}|den_degree={}|"
        "num_coeff_degrees={}|den_coeff_degrees={}|status=BEGIN_INVERSE".format(
            numerator.degree(),
            denominator.degree(),
            tuple(
                (coefficient.numerator().degree(), coefficient.denominator().degree())
                for coefficient in numerator.list()
            ),
            tuple(
                (coefficient.numerator().degree(), coefficient.denominator().degree())
                for coefficient in denominator.list()
            ),
        ),
        flush=True,
    )
    return SV((numerator * algebra_inverse(denominator)) % modulus)


def quotient_inverse(value):
    value = SV(value) % modulus
    return algebra_inverse(value)


def pointed_residue(Wmap):
    """Evaluate the fixed pointed-quartic map directly in the cubic algebra."""

    w = residue(Wmap)
    u = SV(v - KT(alpha))
    q = KT(q_anchor)
    aa = KT(anchor_a)
    bb = KT(anchor_b)
    cc = KT(anchor_c)
    dd = KT(anchor_d)
    u2 = (u**2) % modulus
    u3 = (u2 * u) % modulus
    xg = SV(
        ((2 * q * (w + q) + dd * u) % modulus) * quotient_inverse(u2)
        % modulus
    )
    yg = SV(
        (
            (
                4 * q**2 * (w + q)
                + 2 * q * (dd * u + cc * u2)
                - dd**2 * u2 / (2 * q)
            )
            % modulus
        )
        * quotient_inverse(u3)
        % modulus
    )
    aa1 = dd / q
    aa2 = cc - dd**2 / (4 * q**2)
    aa3 = 2 * q * bb
    b2 = aa1**2 + 4 * aa2
    xres = SV((9 * (xg + b2 / 12)) % modulus)
    yres = SV((27 * (yg + (aa1 * xg + aa3) / 2)) % modulus)
    print("A11S5TRACE|stage=pointed_cubic_residue|status=PASS", flush=True)
    return xres, yres


def trace_via_L4O(xres, yres):
    print("A11S5TRACE|stage=L4O_evaluation|status=BEGIN", flush=True)
    columns = [SV.one(), xres, (xres**2) % modulus, yres]
    evaluation = matrix(
        KT, 3, 4, lambda row, column: columns[column][row]
    )
    kernel = evaluation.right_kernel().basis_matrix()
    print("A11S5TRACE|stage=L4O_kernel|status=PASS", flush=True)
    if kernel.nrows() != 1:
        raise ArithmeticError("L(4O) evaluation kernel is not one-dimensional")
    relation = kernel[0]

    RX = PolynomialRing(KT, "Xtrace")
    Xtrace = RX.gen()
    afun = sum(relation[index] * Xtrace**index for index in range(3))
    bfun = KT(relation[3])
    if not bfun:
        raise ArithmeticError("L(4O) relation has no y term")
    intersection = afun**2 - (Xtrace**3 + KT(A11poly) * Xtrace + KT(B11poly)) * bfun**2
    if intersection.degree() != 4:
        raise ArithmeticError("L(4O) residual intersection is not degree four")

    multiplication = matrix(
        KT,
        3,
        3,
        lambda row, column: SV((xres * v**column) % modulus)[row],
    )
    characteristic = RX(multiplication.charpoly())
    print("A11S5TRACE|stage=x_characteristic|status=PASS", flush=True)
    # All three conjugate x-values are intersections because the pointed map
    # is an exact quartic-to-Jacobian isomorphism and the L(4O) relation
    # vanishes in the cubic algebra.  Vieta therefore gives the fourth root
    # without a costly generic QQ(T)[X] quotient.
    total_intersection_root_sum = KT(-intersection[3] / intersection[4])
    conjugate_x_root_sum = KT(-characteristic[2] / characteristic[3])
    residual_x = KT(total_intersection_root_sum - conjugate_x_root_sum)
    residual_y = KT(-afun(residual_x) / bfun)
    residual = E11(residual_x, residual_y)
    trace = -residual
    return trace, {
        "evaluation_rank": int(evaluation.rank()),
        "kernel_dimension": int(kernel.nrows()),
        "x_characteristic_degree": int(characteristic.degree()),
        "intersection_degree": int(intersection.degree()),
        "residual_count_by_degree": 1,
        "exact_residual_x_by_vieta": True,
        "exact_residual_point_identity": True,
    }


xres, yres = pointed_residue(selected_quartic_sign * Wplus)
trace, certificate = trace_via_L4O(xres, yres)
print(
    f"A11S5TRACE|stage=L4O|quartic_sign={selected_quartic_sign}|status=PASS",
    flush=True,
)
signed_trace = (
    selected_equation_sign * trace
    + origin_translation_coefficient * Qpoint
)
section = normalized_section(signed_trace, A11poly, B11poly)


def reduction_matches(section, marked):
    p = ZZ(marked["prime"])
    F = GF(p)
    R = PolynomialRing(F, "T")
    K = R.fraction_field()
    A = R([F(value.numerator()) / F(value.denominator()) for value in A11poly.list()])
    B = R([F(value.numerator()) / F(value.denominator()) for value in B11poly.list()])
    E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])

    def reduce_exact(poly):
        return R([F(value.numerator()) / F(value.denominator()) for value in poly.list()])

    X, Y, Z = map(reduce_exact, section)
    point = E(K(X) / K(Z**2), K(Y) / K(Z**3))
    target = marked["selected_trace_sections"]["S5"]
    Xm = R(target["X_coefficients_low_to_high"])
    Ym = R(target["Y_coefficients_low_to_high"])
    Zm = R(target["Z_coefficients_low_to_high"])
    return point == E(K(Xm) / K(Zm**2), K(Ym) / K(Zm**3))


if not all(reduction_matches(section, marked) for marked in marked_rows):
    raise ArithmeticError("exact S5 trace missed at least one marked good-prime reduction")
quartic_sign = selected_quartic_sign
equation_sign = selected_equation_sign
X, Y, Z = section

payload = {
    "schema": "elkies-k3.h3-q24-a11-s5-trace-section-qq.v1",
    "status": "PASS_EXACT_MARKED_A11_S5_TRACE_SECTION_QQ",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "source": {
        "equation_shell_index": 5,
        "new_base_degree": 3,
        "cubic_modulus_degree": int(modulus.degree()),
        "selected_quartic_ordinate_sign": int(quartic_sign),
        "trace_equation_sign": int(equation_sign),
        "pointed_anchor_reduction_sign": int(selected_anchor_orientation),
        "origin_translation_coefficient_on_Q": int(origin_translation_coefficient),
    },
    "L4O_certificate": certificate,
    "section": {
        "X_coefficients_low_to_high": [str(value) for value in X.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y.list()],
        "Z_coefficients_low_to_high": [str(value) for value in Z.list()],
        "degrees_X_Y_Z": [int(X.degree()), int(Y.degree()), int(Z.degree())],
        "exact_weierstrass_identity": True,
    },
    "marked_reduction_primes": [int(row["prime"]) for row in marked_rows],
    "exact_reduction_matches_marked_trace_at_every_prime": True,
    "method": (
        "exact cubic algebra and the unique L(4O) relation; the cubic x-characteristic "
        "factor leaves one residual point whose negative is the trace"
    ),
    "large_Groebner_required": False,
    "proof_boundary": (
        "This exactly constructs and marks the A11 S5 trace. Combining it with the exact "
        "S7, S17, and Q sections to construct the pinned equation zero is the next step."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "A11S5TRACE|quartic_sign={}|equation_sign={}|degrees={}|primes={}|status={}".format(
        quartic_sign,
        equation_sign,
        tuple(payload["section"]["degrees_X_Y_Z"]),
        tuple(payload["marked_reduction_primes"]),
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
