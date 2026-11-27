"""Verify the determinant and fibre-coefficient identities in the square exclusion."""

from sympy import Function, Poly, cancel, expand, symbols


x, y, u = symbols("x y u")
L, M0, M1, J0, J1, J2, C4, H3, H20, H21, H10, H11, H12 = [
    Function(name)(x)
    for name in "L M0 M1 J0 J1 J2 C4 H3 H20 H21 H10 H11 H12".split()
]
G = L * y**2 + (M1 * u + M0) * y + J2 * u**2 + J1 * u + J0
H = C4 * y**4 + H3 * y**3 + (H20 + H21 * u) * y**2 + (H10 + H11 * u + H12 * u**2) * y
jacobian = Poly(G.diff(x) * H.diff(y) - G.diff(y) * H.diff(x), u, y)
assert expand(
    jacobian.coeff_monomial(u * y**4) - 4 * M1.diff(x) * C4 + M1 * C4.diff(x)
) == 0

p, b, f, f_prime, b_prime, k, m, d, q0, B2 = symbols(
    "p b f f_prime b_prime k m d q0 B2", nonzero=True
)
alpha = f * k / p
z = (u - d) / p - alpha * y
quadratic = Poly(b * z**2 + (q0 + k * y) * z + B2 * y**2, u, y)
cubic = Poly(z**3 + m * y**2 * z, u, y)
assert cancel(quadratic.coeff_monomial(y**2) - (b * alpha**2 - k * alpha + B2)) == 0
assert cancel(quadratic.coeff_monomial(u * y) - (k - 2 * b * alpha) / p) == 0
assert cancel(cubic.coeff_monomial(u * y**2) - (3 * alpha**2 + m) / p) == 0
assert cancel(cubic.coeff_monomial(u**3) - 1 / p**3) == 0

lam = b * f / p
theta = p * b_prime - 2 * b * f * b_prime
m_from_identity = k**2 * (1 - 2 * lam)**2 / p - 3 * f**2 * k**2 / p**2
B2_from_identity = -f_prime * k**2 / (2 * f * b_prime) - theta * m_from_identity / (3 * f * b_prime)
L_from_identity = b * alpha**2 - k * alpha + B2_from_identity
assert cancel(
    -3 * f * L_from_identity / k**2
    - ((1 - 2 * lam)**3 + 3 * lam**3 * p / b**2 + 3 * f_prime / (2 * b_prime))
) == 0

print("Verified the fibre determinant and final rational identity.")
