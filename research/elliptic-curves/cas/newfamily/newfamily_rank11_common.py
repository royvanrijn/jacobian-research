from sage.all import *

def build_family(roots):
    """
    Six-root quartic construction.

    roots = six rational roots satisfying the quartic condition.

    Returns:
        T, x, P, S, Rquartic, I, J, A, B, Delta
    where the Jacobian is

        y^2 = x^3 + A(T) x + B(T)
            = x^3 - 27 I(T) x - 27 J(T).
    """

    roots = [QQ(r) for r in roots]

    RT = PolynomialRing(QQ, 'T')
    T = RT.gen()

    RX = PolynomialRing(RT, 'x')
    x = RX.gen()

    P = RX(1)

    for r in roots:
        P *= (x - T - r)
        P *= (x + T - r)

    assert P.degree() == 12
    assert P[12] == 1

    # ------------------------------------------------------------
    # Construct the monic degree-6 polynomial S whose square agrees
    # with P in degrees x^12,...,x^6.
    # ------------------------------------------------------------

    scoeff = {
        6: RT(1)
    }

    for degree in range(11, 5, -1):

        j = degree - 6

        known = RT(0)

        for i,si in scoeff.items():
            for k,sk in scoeff.items():
                if i+k == degree:
                    known += si*sk

        scoeff[j] = (P[degree] - known) / 2

    S = RX(0)

    for i,c in scoeff.items():
        S += c*x**i

    Rq = S**2 - P

    if Rq.degree() > 4:
        raise ArithmeticError(
            f"rootset does not give quartic: remainder degree {Rq.degree()}"
        )

    a = RT(Rq[4])
    b = RT(Rq[3])
    c = RT(Rq[2])
    d = RT(Rq[1])
    e = RT(Rq[0])

    # Binary quartic invariants.
    I = (
        12*a*e
        - 3*b*d
        + c**2
    )

    J = (
        72*a*c*e
        + 9*b*c*d
        - 27*a*d**2
        - 27*b**2*e
        - 2*c**3
    )

    A = -27*I
    B = -27*J

    Delta = -16 * (
        4*A**3
        + 27*B**2
    )

    return {
        "roots": tuple(roots),
        "RT": RT,
        "T": T,
        "RX": RX,
        "x": x,
        "P": P,
        "S": S,
        "Rq": Rq,
        "I": RT(I),
        "J": RT(J),
        "A": RT(A),
        "B": RT(B),
        "Delta": RT(Delta),
    }
