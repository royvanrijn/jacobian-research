"""Check the mixed determinant and rank-two identities for the binary case."""

from sympy import Function, Matrix, Poly, cancel, expand, symbols


def jacobian(components, variables):
    return Matrix(components).jacobian(variables).det()


x, y, z, v = symbols("x y z v")
g, U, V = [Function(name)(x, y) for name in ("g", "U", "V")]
N = Matrix([x**2, x * y, y**2])
normal = Matrix([y**2, -2 * x * y, x**2])
tangent = Matrix([2 * x * U, y * U + x * V, 2 * y * V])
form = g * N
transverse = x * V - y * U
mixed = Matrix.hstack(form.diff(x), tangent.diff(y), tangent).det()
mixed += Matrix.hstack(tangent.diff(x), form.diff(y), tangent).det()
assert expand(mixed + 2 * transverse**2 * (2 * g + x * g.diff(x) + y * g.diff(y))) == 0
assert expand(normal.dot(tangent.diff(x)) + 2 * y * transverse) == 0
assert expand(normal.dot(tangent.diff(y)) - 2 * x * transverse) == 0

r, d, s = [Function(name)(x) for name in ("r", "d", "s")]
P, Q, S = [Function(name)(x, y) for name in ("P", "Q", "S")]
F = [x * z**2 + P * z + S, z**2 + Q, r * z + r.diff(x) * P + d]
det = Poly(jacobian(F, (x, y, z)), z)
r_prime = r.diff(x)
r_second = r.diff(x, 2)
d_prime = d.diff(x)
assert det.coeff_monomial(z**3) == 0
assert expand(
    det.coeff_monomial(z**2)
    - ((r - 2 * x * r_prime) * Q.diff(y)
       + 2 * (r_second * P + d_prime) * P.diff(y) + 2 * r_prime * S.diff(y))
) == 0

shift = x - r / (2 * r_prime)
S_from_z2 = shift * Q - r_second * P**2 / (2 * r_prime) - d_prime * P / r_prime + s
K = r_prime + r * r_second / r_prime
a = r.diff(x, 3) - 3 * r_second**2 / r_prime
primitive = (
    a * P**3 / 3
    + (d.diff(x, 2) - 3 * d_prime * r_second / r_prime) * P**2
    - (2 * d_prime**2 / r_prime + 2 * r_prime * s.diff(x)) * P
    - (K * P + r * d_prime / r_prime) * Q
)
assert cancel(
    det.coeff_monomial(z).subs(S, S_from_z2).doit() - primitive.diff(y)
) == 0

av, bv = [Function(name)(x) for name in ("av", "bv")]
quadratic_coordinates = [x * z**2 + v * z + bv * v**2, z**2 + av * v**2, r * z + r_prime * v + d]
constant_det = jacobian(quadratic_coordinates, (x, v, z)).subs(z, 0)
assert expand(
    Poly(constant_det, v).coeff_monomial(v**3)
    - (2 * r * (av * bv.diff(x) - bv * av.diff(x)) + r_prime * av.diff(x) - 2 * r_second * av)
) == 0

print("Verified the binary mixed determinants and rank-two coefficient identities.")

u = symbols("u")
ratio, leading, slope, offset, free = [
    Function(name)(x) for name in ("ratio", "leading", "slope", "offset", "free")
]
V = Function("V")(x, y)
R1 = slope * V + offset
R0 = slope.diff(x) * V**2 / (2 * ratio.diff(x)) + offset.diff(x) * V / ratio.diff(x) + free
F2 = ratio * u + V
F3 = leading * u**2 + R1 * u + R0
slice_det = Poly(jacobian([F2, F3], (x, y)), u)
assert expand(
    slice_det.coeff_monomial(u**2) - (ratio.diff(x) * slope - leading.diff(x)) * V.diff(y)
) == 0
assert cancel(slice_det.coeff_monomial(u)) == 0
assert cancel(
    slice_det.coeff_monomial(1)
    + (V**2 * (slope.diff(x) / (2 * ratio.diff(x))).diff(x)
       + V * (offset.diff(x) / ratio.diff(x)).diff(x) + free.diff(x)) * V.diff(y)
) == 0
print("Verified the final first-component fibre determinant.")
