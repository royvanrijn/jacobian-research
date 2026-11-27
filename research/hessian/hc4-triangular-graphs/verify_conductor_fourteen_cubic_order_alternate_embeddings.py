"""Replay the norm obstruction and Laurent chain rule for the cubic order."""

from sympy import Function, Poly, Rational, diff, expand, resultant, simplify, symbols


def check_norm_obstruction():
    z, U, Y, a, h, p, q = symbols("z U Y a h p q")
    cubic = z**3 + a*z + h - U
    candidate = q*z + p*(z**2 + 2*a/3)
    norm = Poly(resultant(cubic, Y-candidate, z), Y)
    assert norm.degree() == 3 and norm.LC() == 1 and norm.nth(2) == 0
    assert simplify(norm.nth(1) - (a*q**2 - 3*(U-h)*p*q - a**2*p**2/3)) == 0
    assert simplify(norm.nth(0) - (
        -(U-h)**2*p**3 + (U-h)*a*p**2*q - (U-h)*q**3
        - 2*a**3*p**3/27 - 2*a**2*p*q**2/3
    )) == 0

    p0, p1, q0, q1, q2 = symbols("p0 p1 q0 q1 q2")
    terms = {p: U**2 + p1*U + p0, q: q2*U**2 + q1*U + q0}
    A = Poly(expand(norm.nth(1).subs(terms, simultaneous=True)), U)
    B = Poly(expand(norm.nth(0).subs(terms, simultaneous=True)), U)
    assert A.nth(5) == -3*q2
    assert B.nth(7).subs(q2, 0) == 2*h - 3*p1
    top = {q2: 0, p1: 2*h/3}
    assert expand(A.nth(3).subs(top, simultaneous=True)
                  + 3*q0 - h*(-4*a**2+9*q1)/9) == 0
    top[q0] = h*(-4*a**2+9*q1)/27
    assert B.nth(1).subs(top, simultaneous=True).subs(h, 0) == 0

    assert expand(B.nth(5).subs(top, simultaneous=True)
                  - 2*h*(-4*a**3+9*a*q1+5*h**2+27*p0)/27) == 0
    top[p0] = (4*a**3-9*a*q1-5*h**2)/27
    assert expand(A.nth(1).subs(top, simultaneous=True) - 8*h**3*(q1-a**2)/27) == 0
    top[q1] = a**2
    assert expand(B.nth(3).subs(top, simultaneous=False) - 64*h**5/729) == 0


def check_chain_rule():
    z, y = symbols("z y")
    V = Function("V")(z, y)
    k = Function("k")(y)

    def jacobian(f):
        return diff(V, z)*diff(f, y) - diff(V, y)*diff(f, z)

    assert simplify(jacobian(V**Rational(3, 2))) == 0
    term = k*V**Rational(-1, 2)/2
    assert simplify(jacobian(term) - diff(k, y)*diff(V**Rational(1, 2), z)) == 0


if __name__ == "__main__":
    check_norm_obstruction()
    check_chain_rule()
    print("Verified the alternate cubic-order norm and Laurent identities.")
