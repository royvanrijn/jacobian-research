"""Replay the Taylor and determinant identities in the square residual."""

from sympy import Function, Matrix, Poly, expand, symbols


def check_taylor_correction():
    u, w, a, b, epsilon = symbols("u w a b epsilon")
    veronese = Matrix([u**2, u*w, w**2])
    derivative = veronese.jacobian([u, w]) * Matrix([a, b])
    shift = {u: u + epsilon*a, w: w + epsilon*b}
    quadratic = Matrix([a**2, a*b, b**2])

    second_difference = veronese.subs(shift, simultaneous=True)
    second_difference -= veronese + epsilon*derivative + epsilon**2*quadratic
    assert all(expand(entry) == 0 for entry in second_difference)

    first_difference = derivative.subs(shift, simultaneous=True) - derivative
    assert all(expand(entry) == 0 for entry in first_difference - 2*epsilon*quadratic)


def check_determinants():
    s, v, z, delta = symbols("s v z delta")
    p, b = [Function(name)(s, v) for name in ("p", "b")]
    a, q, lower_b, c, rho, lower_c = [
        Function(name)(s, v) for name in ("A", "q", "B", "c", "rho", "C")
    ]
    field = Matrix([p*z + a, b*z**2 + q*z + lower_b, delta*z**3 + c*z**2 + rho*z + lower_c])
    determinant = Poly(field.jacobian([s, v, z]).det(), z)
    fifth = 3*delta*(p.diff(s)*b.diff(v) - p.diff(v)*b.diff(s))
    assert expand(determinant.nth(5) - fifth) == 0

    p, b = [Function(name)(s) for name in ("p", "b")]
    depressed = Matrix([p*z+a, b*z**2 + q*z + lower_b, delta*z**3 + rho*z + lower_c])
    determinant = Poly(depressed.jacobian([s, v, z]).det(), z)
    fourth = 3*delta*(p.diff(s)*q.diff(v) - b.diff(s)*a.diff(v))
    third = (
        3*delta*(p.diff(s)*lower_b.diff(v) + a.diff(s)*q.diff(v) - a.diff(v)*q.diff(s))
        + (p*b.diff(s) - 2*b*p.diff(s))*rho.diff(v)
    )
    assert expand(determinant.nth(4) - fourth) == 0
    assert expand(determinant.nth(3) - third) == 0


if __name__ == "__main__":
    check_taylor_correction()
    check_determinants()
    print("Verified the square residual Taylor and determinant identities.")
