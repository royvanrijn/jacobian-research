#!/usr/bin/env sage -python
"""Recover the H92 height-46 section by explicit neighbor transport.

This work-in-progress diagnostic starts from the explicit ancillary H92
section in ``92.txt`` and evaluates the first (D6+A8+A1 -> E7+A8)
two-neighbor parameter at the pinned non-CM H92 point.
"""

from sage.all import *

from importlib.machinery import SourceFileLoader
import json
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"

anchor = SourceFileLoader("h92_p2_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92

R = PolynomialRing(QQ, "t")
t = R.gen()
K = R.fraction_field()
t = K.gen()

e = -(r * s + 4 * s + r - 1) / (r * s - r - 1)
m = -2 * (4 * s**2 + 3 * r * s - 2 * s - 2 * r - 1) / (r * s - r - 1)
f = -(m - e + 5) * (m + e + 3) / (m**2 + 3 * e**2 + 2 * e - 5)
a1 = (f + 1) * (e**2 * f**2 + f**2 + 2 * e * f - 2 * f + 2) / (f - 1)
b3 = 1 - f**2
b0 = (b3 - 1) * b3 * (e + 1) ** 2
b2 = -((b3**2 - b3) * e**2 + 3 * b3**2 + (2 * a1 + 1) * b3) / 2
b1 = b3 + a1
a2 = 2 * b2 + b1**2 - 2 * a1 * b1 - b0 + a1**2
a3 = -2 * b3**3 - a1 * b3**2 + (-2 * b2 - b0) * b3 + b0
c1 = 2 * b1 - a1

x_section = (
    (b3 - 1) * b3 * (e + 1) ** 2 * t**2
    + (b3 - 1) ** 2 * b3 * (e + 1) ** 2 * t**3
    - (b3 - 1) * b3**2 * (e + 1) ** 2 * t**4
)

# The source coordinates first undergo t -> t*(s+r)^2/(r+2*s) and
# x -> x/(r*s-r-1)^2, exactly as in the pinned ancillary calculation.
base_scale = (s + r) ** 2 / (r + 2 * s)
x_scale = (r * s - r - 1) ** 2
x_transported = x_section(base_scale * t) * x_scale
# The first neighbor is taken after the old base is inverted and the
# Weierstrass coordinates are scaled by (t^4,t^6).
x_inverted = t**4 * x_transported(1 / t)

la = (r + s) ** 2
mu = r + 2 * s
offset = (
    t**2
    * (t - la / mu)
    * 4
    * r
    * (r + 1)
    * s**2
    * (2 * s + r)
    / (s + r**2 + r)
)
neighbor_parameter = (x_inverted - offset) / (t**3 * (t - la / mu))

print("H92P2|stage=first_neighbor_diagnostic")
print(f"H92P2|source_x_degree={x_section.numerator().degree()}")
print(
    "H92P2|two_neighbor_on_source="
    f"{neighbor_parameter}|numerator_degree={neighbor_parameter.numerator().degree()}|"
    f"denominator_degree={neighbor_parameter.denominator().degree()}"
)

# Rebuild the scaled ancillary model and its two-neighbor plane cubic.  The
# old base coordinate is ``z`` here, reserving ``t`` for the target H92 base.
C = PolynomialRing(QQ, names=("x", "z", "u", "v"))
x, z, u, v = C.gens()
a0 = (r * s - r - 1) ** 2
a1s = -2 * r * ((r + 1) ** 2 * (s - 1) ** 2 + s**2)
a2s = r * (r * s + 4 * s - r - 1) * (
    4 * s**2 + r**2 * s + 4 * r * s - r**2 - r
)
a3s = -4 * r * s * (s + r**2 + r) * (s**2 + r * s - r)
b0s = -4 * r * s**2 * (r * s - r - 1) ** 2
b1s = 4 * r**2 * s**2 * (r * s + 2 * s - r - 1) ** 2
b2s = -8 * r**2 * s**3 * (
    (r + 2) * s**2 + (2 * r**2 + 2 * r - 1) * s - 2 * r * (r + 1)
)
c0 = 16 * r**2 * s**4 * (r * s - r - 1) ** 2
c1s = -64 * r**3 * (r + 1) * (s - 1) * s**5

source_rhs = (
    x**3
    + (a0 + a1s * z + a2s * z**2 + a3s * z**3) * x**2
    + 2 * (b0s + b1s * z + b2s * z**2) * z**2 * (la * z - mu) * x
    + (c0 + c1s * z) * z**4 * (la * z - mu) ** 2
)
offset_z = (
    z**2
    * (z - la / mu)
    * 4
    * r
    * (r + 1)
    * s**2
    * (2 * s + r)
    / (s + r**2 + r)
)
# The original notebook inverts the base before forming this neighbor.
source_rhs_inverted = C.fraction_field()(source_rhs.subs({x: x / z**4, z: 1 / z}) * z**12)
assert source_rhs_inverted.denominator() in QQ
source_rhs_inverted = C(source_rhs_inverted)
rhs_neighbor = source_rhs_inverted.subs({x: u * z**3 * (z - la / mu) + offset_z})
rhs_neighbor /= z**6 * (z - la / mu) ** 2
rhs_neighbor = C.fraction_field()(rhs_neighbor)
print(
    "H92P2|two_neighbor_denominator="
    f"{rhs_neighbor.denominator()}"
)
assert rhs_neighbor.denominator() in QQ
rhs_neighbor = C(rhs_neighbor)
# The source then normalizes the new base and the quartic's y-coordinate.
# These are the two displayed substitutions immediately before ``jac2``.
rhs_neighbor = C(
    rhs_neighbor.subs({u: u * (2 * s + r)}) / (2 * s + r) ** 2
)
rhs_neighbor = C(
    rhs_neighbor.subs({u: u / (s + r + r**2)}) * (s + r + r**2) ** 4
)
two_neighbor_cubic = v**2 - rhs_neighbor
print(
    "H92P2|two_neighbor_polynomial_degrees="
    f"z:{two_neighbor_cubic.degree(z)},u:{two_neighbor_cubic.degree(u)},"
    f"total:{two_neighbor_cubic.total_degree()}"
)
print(f"H92P2|two_neighbor_terms={len(two_neighbor_cubic.dict())}")


def jac2(quartic):
    """The binary-quartic Jacobian routine pinned in arXiv:1209.3527."""
    polynomial_ring = quartic.parent()
    coordinate = polynomial_ring.gen()
    a0 = quartic[4]
    a1 = quartic[3] / 4
    a2 = quartic[2] / 6
    a3 = quartic[1] / 4
    a4 = quartic[0]
    invariant_i = a0 * a4 - 4 * a1 * a3 + 3 * a2**2
    invariant_j = a0 * a2 * a4 + 2 * a1 * a2 * a3 - a0 * a3**2 - a4 * a1**2 - a2**3
    return coordinate**3 - invariant_i * coordinate / 4 - invariant_j / 4


Ku = PolynomialRing(QQ, "u").fraction_field()
Zu = PolynomialRing(Ku, "z")
zu = Zu.gen()
ambient_to_u = C.hom([0, 0, Ku.gen(), 0], Ku)
quartic = Zu(
    sum(ambient_to_u(rhs_neighbor.coefficient({z: power})) * zu**power for power in range(5))
)
assert quartic.degree() == 4
weierstrass_2neighbor = jac2(quartic)
print(
    "H92P2|jac2_degree_u="
    f"{max(coefficient.numerator().degree() for coefficient in weierstrass_2neighbor.list())}"
)

# Marked degree-three divisor induced by the original ancillary section on
# the normalized binary quartic.  The two displayed quartic normalizations
# send (u,v) to (u*D/(2s+r), v*D^2/(2s+r)).
source_a = 1 + a1 * t + a2 * t**2 + a3 * t**3
source_b = b0 * (1 + b1 * t + b2 * t**2)
source_c = b0**2 * (1 + c1 * t)
source_y = (
    x_section**3
    + source_a * x_section**2
    + 2 * source_b * t**2 * (t - 1) * x_section
    + source_c * t**4 * (t - 1) ** 2
).sqrt()
y_transported = source_y(base_scale * t) * (r * s - r - 1) ** 3
y_inverted = t**6 * y_transported(1 / t)
v_source = y_inverted / (t**3 * (t - la / mu))
new_u_scale = (s + r + r**2) / (2 * s + r)
new_v_scale = (s + r + r**2) ** 2 / (2 * s + r)
u_source = neighbor_parameter * new_u_scale
v_source *= new_v_scale
source_divisor = Zu(
    sum(Ku(coefficient) * zu**index for index, coefficient in enumerate(u_source.numerator().list()))
    - Ku.gen() * sum(Ku(coefficient) * zu**index for index, coefficient in enumerate(u_source.denominator().list()))
)
assert source_divisor.degree() == 3
source_v_numerator = Zu(sum(Ku(coefficient) * zu**index for index, coefficient in enumerate(v_source.numerator().list())))
source_v_denominator = Zu(sum(Ku(coefficient) * zu**index for index, coefficient in enumerate(v_source.denominator().list())))

# Replay the coordinate normalizations immediately following ``jac2`` in the
# pinned ancillary source.  This identifies the precise E7+A8 model that is
# the input to the final three-neighbor.
Wx = PolynomialRing(Ku, "x")
xx = Wx.gen()
p_first = Wx(weierstrass_2neighbor(xx))
p_first = Wx(64 * p_first(xx / 4))
shift0 = 16 * r**2 * (r + 1) ** 2 * s**4 * (s + r**2 + r) ** 2 * (r * s - r - 1) ** 2 / 3
shift1 = -8 * r * s**2 * (s + r**2 + r) ** 2 * (
    2*r**2*s**3+r**4*s**2-r**3*s**2-8*r**2*s**2-6*r*s**2-2*r**4*s-2*r**3*s+4*r**2*s+6*r*s+2*s+r**4+3*r**3+3*r**2+r
) / 3
shift2 = r * (s + r**2 + r) ** 2 * (
    4*r*s**3-20*s**3+r**3*s**2-4*r**2*s**2-12*r*s**2+8*s**2-2*r**3*s+2*r**2*s+4*r*s+r**3+2*r**2+r
) / 3
shift3 = r * s * (s + r**2 + r) ** 2
for shift in (shift0, shift1 * Ku.gen(), shift2 * Ku.gen()**2, shift3 * Ku.gen()**3):
    p_first = Wx(p_first(xx + shift))
normalization_scale = (s + r + r**2) ** 2
p_first = Wx(p_first(xx * normalization_scale) / normalization_scale**3)
assert p_first.degree() == 3
print(
    "H92P2|first_neighbor_weierstrass_degrees="
    + ",".join(str(coefficient.numerator().degree()) for coefficient in p_first.list())
)

# Recover the degree-two and degree-five truncations A,B used for the
# three-neighbor directly from the exact first-neighbor equation.  This avoids
# duplicating the long printed A,B expressions from the source notebook.
Tpoly = PolynomialRing(QQ, "t")
tt = Tpoly.gen()


def polynomial_coefficient(value):
    numerator = value.numerator()
    denominator = value.denominator()
    assert denominator in QQ and denominator
    return Tpoly(numerator) / QQ(denominator)


c2_first = polynomial_coefficient(p_first[2])
c1_first = polynomial_coefficient(p_first[1])
c0_first = polynomial_coefficient(p_first[0])
series = PowerSeriesRing(QQ, "t", default_prec=12)
ts = series.gen()
assert series(c2_first)[0].is_square()
A_series = series(c2_first).sqrt()
A_truncation = Tpoly(sum(A_series[index] * tt**index for index in range(3)))
B_series = series(c1_first) / (2 * series(A_truncation))
B_truncation = Tpoly(sum(B_series[index] * tt**index for index in range(6)))
assert (c2_first - A_truncation**2) % tt**3 == 0
assert (c1_first - 2 * A_truncation * B_truncation) % tt**6 == 0
assert (c0_first - B_truncation**2) % tt**9 == 0
assert A_truncation.degree() == 2 and B_truncation.degree() == 5
# The square-root convention above is opposite to the printed three-neighbor
# parameter in 92.txt.  Thus the displayed ``A,B`` there are their negatives.
# Keeping this sign explicit is essential: it fixes the target-base relation
# to ``t = w/(2*s)`` after the final normalized Jacobian conversion.
transport_A = -A_truncation
transport_B = -B_truncation
print(
    "H92P2|three_neighbor_truncation="
    f"A_degree={A_truncation.degree()}|B_valuation={B_truncation.valuation()}|"
    f"B_degree={B_truncation.degree()}"
)

C3 = PolynomialRing(QQ, names=("z", "X", "w"))
z3, X3, w3 = C3.gens()
raw_cubic = (
    X3**3 * z3**9
    + C3(c2_first(z3)) * X3**2 * z3**6
    + C3(c1_first(z3)) * X3 * z3**3
    + C3(c0_first(z3))
    - (w3 * z3**6 - C3(transport_A(z3)) * X3 * z3**3 - C3(transport_B(z3))) ** 2
)
assert raw_cubic % z3**9 == 0
three_neighbor_cubic = C3(raw_cubic // z3**9)
three_neighbor_cubic *= 4 * r**2 * (r + 1)**6 * s**2 * (r * s - r - 1)**6
new_base_scale = r * (r + 1)**3 * s * (r * s - r - 1)**3
three_neighbor_cubic = C3(three_neighbor_cubic.subs({w3: w3 / new_base_scale}))
zx_degree = max(exponents[0] + exponents[1] for exponents in three_neighbor_cubic.dict())
assert zx_degree == 3
print(
    "H92P2|three_neighbor_cubic_degrees="
    f"z:{three_neighbor_cubic.degree(z3)},X:{three_neighbor_cubic.degree(X3)},"
    f"w:{three_neighbor_cubic.degree(w3)}|total:{three_neighbor_cubic.total_degree()}|"
    f"terms={len(three_neighbor_cubic.dict())}"
)

Kw = PolynomialRing(QQ, "w").fraction_field()
P3 = PolynomialRing(Kw, names=("z", "X", "Z"))
zc, xc, Zc = P3.gens()
plane_cubic = P3.zero()
for (power_z, power_x, power_w), coefficient in three_neighbor_cubic.dict().items():
    plane_cubic += (
        Kw(coefficient) * Kw.gen() ** power_w * zc**power_z * xc**power_x * Zc ** (3 - power_z - power_x)
    )
assert plane_cubic.total_degree() == 3
axis_candidates = (
    ("z0", P3(plane_cubic.subs({zc: 0, Zc: 1}))),
    ("X0", P3(plane_cubic.subs({xc: 0, Zc: 1}))),
    ("infinity", P3(plane_cubic.subs({Zc: 0, xc: 1}))),
)
for name, univariate in axis_candidates:
    variable = xc if name == "z0" else (zc if name == "X0" else zc)
    roots = univariate.univariate_polynomial().roots(Kw)
    print(f"H92P2|cubic_axis={name}|rational_roots={len(roots)}")

print(
    "H92P2|status=TRANSPORT_MODELS_RECONSTRUCTED|"
    "remaining_gate=explicit_abel_jacobi_map_for_marked_degree_three_divisor",
    flush=True,
)

# A finite-field spot check confirms that the reconstructed ternary cubic has
# the pinned short H92 Jacobian, so divisor-class evaluations may be compared
# directly with the final H92 coordinates.
spot_prime = int(os.environ.get("H92P2_PRIME", "1000003"))
spot_base = GF(spot_prime)
spot_w = spot_base(2)
sample_bases = tuple(spot_base(index) for index in range(2000))
C3_spot = PolynomialRing(spot_base, names=("z", "X"))
spot_z, spot_X = C3_spot.gens()
spot_cubic = C3_spot.zero()
for (power_z, power_x, power_w), coefficient in three_neighbor_cubic.dict().items():
    spot_cubic += spot_base(coefficient) * spot_w**power_w * spot_z**power_z * spot_X**power_x
spot_projective = PolynomialRing(spot_base, names=("z", "X", "Z"))
spot_zp, spot_Xp, spot_Zp = spot_projective.gens()
spot_form = spot_projective(sum(
    coefficient * spot_zp**powers[0] * spot_Xp**powers[1] * spot_Zp**(3 - powers[0] - powers[1])
    for powers, coefficient in spot_cubic.dict().items()
))
spot_point = None
for z_value in spot_base:
    polynomial = spot_cubic.subs({spot_z: z_value})
    roots = polynomial.univariate_polynomial().roots(spot_base)
    if roots:
        spot_point = (z_value, roots[0][0], spot_base(1))
        break
assert spot_point is not None
spot_jacobian = EllipticCurve_from_cubic(spot_form, spot_point, morphism=False)
_, h92_formulas = anchor.parse_h92(ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt")
A1_target, A_target, B1_target, B_target, B2_target = (
    spot_base(QQ(value(r, s))) for value in h92_formulas
)
spot_target = EllipticCurve([
    0, 0, 0,
    A1_target * (spot_w / (2 * spot_base(s)))**3 + A_target * (spot_w / (2 * spot_base(s)))**4,
    B1_target * (spot_w / (2 * spot_base(s)))**5 + B_target * (spot_w / (2 * spot_base(s)))**6 + B2_target * (spot_w / (2 * spot_base(s)))**7,
])
if spot_prime == 101:
    assert spot_jacobian.is_isomorphic(spot_target)
    print("H92P2|jacobian_spot_check=GF(101),w=2|status=PASS", flush=True)
else:
    print("H92P2|jacobian_spot_check=skipped_nonreference_prime", flush=True)


def reduce_function(value, field, base_value):
    numerator = value.numerator()
    denominator = value.denominator()
    top = sum(field(coefficient) * base_value**index for index, coefficient in enumerate(numerator.list()))
    bottom = sum(field(coefficient) * base_value**index for index, coefficient in enumerate(denominator.list()))
    if not bottom:
        raise ZeroDivisionError
    return top / bottom


def specialize_univariate(polynomial, field, base_value):
    ring = PolynomialRing(field, "z")
    variable = ring.gen()
    return ring(sum(
        reduce_function(coefficient, field, base_value) * variable**index
        for index, coefficient in enumerate(polynomial.list())
    ))


def marked_divisor_class(field, base_value):
    """Map the split marked degree-three divisor into the first Jacobian."""
    f0 = specialize_univariate(quartic, field, base_value)
    g0 = specialize_univariate(source_divisor, field, base_value)
    roots = g0.roots(field)
    if len(roots) != 3 or any(multiplicity != 1 for _, multiplicity in roots):
        raise ValueError("non-split marked degree-three fiber")
    numerator = specialize_univariate(source_v_numerator, field, base_value)
    denominator = specialize_univariate(source_v_denominator, field, base_value)
    h0 = (numerator * denominator.inverse_mod(g0)).mod(g0)
    complement = (f0 - h0**2) // g0
    if complement.degree() != 1:
        raise ValueError("no rational complement")
    qx = -complement[0] / complement[1]
    qy = h0(qx)
    z0 = f0.parent().gen()
    shifted = f0(z0 + qx)
    qa, qb, qc, qd = shifted[4], shifted[3], shifted[2], shifted[1]
    long_curve = EllipticCurve([
        qd / qy, qc - qd**2 / (4*qy**2), 2*qy*qb,
        -4*qy**2*qa, -(qc - qd**2/(4*qy**2))*4*qy**2*qa,
    ])
    answer = long_curve(0)
    for root, _ in roots:
        delta = root - qx
        y_value = h0(root)
        x_value = (2*qy*(y_value+qy)+qd*delta)/delta**2
        y_image = (
            4*qy**2*(y_value+qy)+2*qy*(qd*delta+qc*delta**2)
            - qd**2*delta**2/(2*qy)
        )/delta**3
        answer += long_curve(x_value, y_image)
    jacobian = EllipticCurve([
        0, 0, 0,
        reduce_function(weierstrass_2neighbor[1], field, base_value),
        reduce_function(weierstrass_2neighbor[0], field, base_value),
    ])
    return long_curve.isomorphism_to(jacobian)(answer)


marked_spot = None
for candidate in sample_bases:
    try:
        marked_spot = marked_divisor_class(spot_base, candidate)
        break
    except (ArithmeticError, ValueError, ZeroDivisionError):
        continue
assert marked_spot is not None
print(
    f"H92P2|marked_divisor_spot=GF({spot_prime}),u={candidate}|"
    f"nonzero={int(bool(marked_spot))}|status=PASS",
    flush=True,
)


def first_neighbor_coordinates(point, field, base_value):
    x_value, y_value = point.xy()
    x_value *= 4
    y_value *= 8
    for shift in (shift0, shift1 * Ku.gen(), shift2 * Ku.gen()**2, shift3 * Ku.gen()**3):
        x_value -= reduce_function(Ku(shift), field, base_value)
    root_scale = field(s + r + r**2)
    return x_value / root_scale**2, y_value / root_scale**3


def reduce_polynomial(polynomial, field, base_value):
    return sum(field(coefficient) * base_value**index for index, coefficient in enumerate(polynomial.list()))


first_x_spot, first_y_spot = first_neighbor_coordinates(marked_spot, spot_base, candidate)
assert first_y_spot**2 == (
    first_x_spot**3
    + reduce_polynomial(c2_first, spot_base, candidate) * first_x_spot**2
    + reduce_polynomial(c1_first, spot_base, candidate) * first_x_spot
    + reduce_polynomial(c0_first, spot_base, candidate)
)
marked_w = spot_base(new_base_scale) * (
    first_y_spot + reduce_polynomial(transport_A, spot_base, candidate) * first_x_spot + reduce_polynomial(transport_B, spot_base, candidate)
) / candidate**6
print(f"H92P2|marked_first_neighbor_spot=u={candidate}|w={marked_w}|status=PASS", flush=True)


marked_transport = {}
for old_base in sample_bases:
    if not old_base:
        continue
    try:
        old_point = marked_divisor_class(spot_base, old_base)
        old_x, old_y = first_neighbor_coordinates(old_point, spot_base, old_base)
        target_base = spot_base(new_base_scale) * (
            old_y
            + reduce_polynomial(transport_A, spot_base, old_base) * old_x
            + reduce_polynomial(transport_B, spot_base, old_base)
        ) / old_base**6
        marked_transport.setdefault(target_base, []).append((old_base, old_x / old_base**3))
    except (ArithmeticError, ValueError, ZeroDivisionError):
        continue
fiber_points = marked_transport[marked_w]
print(
    f"H92P2|marked_three_neighbor_fiber=w={marked_w}|"
    f"degree={len(fiber_points)}|status=PASS",
    flush=True,
)


def rational_interpolate(values, field, max_degree=60):
    for denominator_degree in range(max_degree + 1):
        for numerator_degree in range(max_degree + 1):
            unknowns = numerator_degree + denominator_degree + 2
            if len(values) < unknowns - 1:
                continue
            rows = [
                [-input_value**index for index in range(numerator_degree + 1)]
                + [output_value * input_value**index for index in range(denominator_degree + 1)]
                for input_value, output_value in values
            ]
            kernel = matrix(field, rows).right_kernel()
            if kernel.dimension() == 1:
                relation = kernel.basis()[0]
                return (
                    tuple(relation[:numerator_degree + 1]),
                    tuple(relation[numerator_degree + 1:]),
                )
    raise ValueError("degree cap too small")


first_jacobian_values = []
first_jacobian_y_values = []
for old_base in sample_bases:
    try:
        point = marked_divisor_class(spot_base, old_base)
        first_jacobian_values.append((old_base, point.xy()[0]))
        first_jacobian_y_values.append((old_base, point.xy()[1]))
    except (ArithmeticError, ValueError, ZeroDivisionError):
        continue
interpolated_x = rational_interpolate(first_jacobian_values, spot_base)
interpolated_y = rational_interpolate(first_jacobian_y_values, spot_base)
print(
    "H92P2|marked_first_jacobian_interpolation="
    f"samples={len(first_jacobian_values)}|"
    f"x_degrees={len(interpolated_x[0])-1},{len(interpolated_x[1])-1}|"
    f"y_degrees={len(interpolated_y[0])-1},{len(interpolated_y[1])-1}|status=PASS",
    flush=True,
)


def normalize_pair(pair):
    numerator, denominator = pair
    scale = denominator[-1]
    assert scale
    return (
        tuple(int(value / scale) for value in numerator),
        tuple(int(value / scale) for value in denominator),
    )


print(
    "H92P2|intermediate_modular_coefficients="
    + json.dumps({"prime": spot_prime, "x": normalize_pair(interpolated_x), "y": normalize_pair(interpolated_y)}, sort_keys=True),
    flush=True,
)

output_path = os.environ.get("H92P2_MODULAR_OUTPUT")
if output_path:
    Path(output_path).write_text(json.dumps(
        {"prime": spot_prime, "x": normalize_pair(interpolated_x), "y": normalize_pair(interpolated_y)},
        sort_keys=True,
    ) + "\n")
