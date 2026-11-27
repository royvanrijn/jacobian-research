"""Replay coefficient identities for the square case on the planes b=s."""

from sympy import Function, Matrix, Poly, expand, symbols


def check_plane_coefficients():
    s, v = symbols("s v")
    gamma = symbols("gamma", nonzero=True)
    a, b0, linear_x, mixed, y_squared, d, e, j, k, m, n = symbols(
        "a b0 linear_x mixed y_squared d e j k m n"
    )
    c_x, c_y, c_0 = symbols("c_x c_y c_0")
    x = v
    y = (s - v**2 - a*v - b0) / gamma
    P2 = x*(a*x + gamma*y)
    Q2 = linear_x*x**2 + mixed*x*y + y_squared*y**2
    c = 2*x*y + c_x*x + c_y*y + c_0
    q = x**3*y + x*Q2 + y*P2 + d*y**2 + e*x*y + j*x**2 + k*x + m*y + n
    q -= 2*s*c/3

    assert expand(Poly(q, v).nth(5) - y_squared/gamma**2) == 0
    q = q.subs(y_squared, 0)
    assert expand(Poly(q, v).nth(4) + mixed/gamma - d/gamma**2) == 0
    q = q.subs(d, gamma*mixed)
    assert expand(Poly(q, v).nth(3).diff(s) - 1/(3*gamma)) == 0


def check_planar_determinant():
    v, z, b, q = symbols("v z b q")
    B, rho, C = [Function(name)(v) for name in ("B", "rho", "C")]
    plane = Matrix([b*z**2 + q*z + B, z**3 + rho*z + C])
    determinant = Poly(plane.jacobian([v, z]).det(), z)
    expected = (
        B.diff(v)*rho - q*C.diff(v),
        -2*b*C.diff(v) - q*rho.diff(v),
        3*B.diff(v) - 2*b*rho.diff(v),
    )
    assert all(expand(determinant.nth(i) - value) == 0 for i, value in enumerate(expected))


if __name__ == "__main__":
    check_plane_coefficients()
    check_planar_determinant()
    print("Verified the square-case coefficients on b=s and the planar determinant.")
