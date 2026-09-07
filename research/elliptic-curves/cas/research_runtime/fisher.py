"""Exact Fisher (2022), Theorem 3.1 / Remark 3.3 arithmetic.

Extracted from the repository's retained fixed-field replayer. No experiment,
field construction, BNF or point-search side effects occur on import.
https://antsmath.org/ANTSXV/papers/ANTS-XV_fisher.pdf
"""
from sage.all import PolynomialRing, QQ, ZZ, prod
R = PolynomialRing(QQ, "x")

def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def invariants(q):
    e, d, c, b, a = [q[i] for i in range(5)]
    return (12*a*e - 3*b*d + c*c,
            72*a*c*e - 27*a*d*d - 27*b*b*e + 9*b*c*d - 2*c**3)


def rational_unit(value, prime):
    value, prime = QQ(value), ZZ(prime)
    require(value != 0, "zero has no unit squareclass")
    valuation = value.valuation(prime)
    return valuation, value / prime**valuation


def residue(value, modulus):
    value = QQ(value)
    return int(value.numerator() * value.denominator().inverse_mod(modulus) % modulus)


def legendre_unit(value, prime):
    return 1 if pow(residue(value, prime), (int(prime)-1)//2, int(prime)) == 1 else -1


def local_square(value, prime):
    if value == 0:
        return True
    valuation, unit = rational_unit(value, prime)
    return valuation % 2 == 0 and (residue(unit, 8) == 1 if prime == 2
                                 else legendre_unit(unit, prime) == 1)


def hilbert_symbol(left, right, prime):
    """Elementary rational Hilbert formula; independent of PARI hilbert."""
    require(left != 0 and right != 0, "Hilbert arguments must be nonzero")
    if prime == "infinity":
        return -1 if left < 0 and right < 0 else 1
    va, ua = rational_unit(left, prime)
    vb, ub = rational_unit(right, prime)
    if prime == 2:
        a, b = residue(ua, 8), residue(ub, 8)
        exponent = (a-1)*(b-1)//4 + va*(b*b-1)//8 + vb*(a*a-1)//8
        return -1 if exponent % 2 else 1
    answer = -1 if (va*vb*((int(prime)-1)//2)) % 2 else 1
    if vb % 2:
        answer *= legendre_unit(ua, prime)
    if va % 2:
        answer *= legendre_unit(ub, prime)
    return answer


def primitive(q):
    denominator = q.denominator()
    integral = q * denominator
    content = ZZ(0)
    for coefficient in integral:
        content = content.gcd(ZZ(coefficient))
    require(content != 0, "zero pairing quadratic")
    answer = integral / content
    return -answer if answer.leading_coefficient() < 0 else answer


def fisher_gamma(quartics, square_root, I, J):
    """Check m^2=z1*z2*z3 and derive the primitive gamma of Theorem 3.1."""
    T = PolynomialRing(QQ, "phi")
    phi = T.gen()
    L = T.quotient(phi**3 - 3*I*phi + J, "phi_bar")
    phi = L.gen()
    S = PolynomialRing(L, "t")
    cubic = []
    for q in quartics:
        require(q.degree() == 4 and invariants(q) == (I, J), "quartic invariant mismatch")
        a, b, c = q[4], q[3], q[2]
        cubic.append((4*a*phi + 3*b*b - 8*a*c)/3)
    m = sum(QQ(c)*phi**i for i, c in enumerate(square_root))
    require(m != 0 and m*m == prod(cubic), "invalid cubic square-root identity")
    q = quartics[0]
    e, d, c, b, a = [q[i] for i in range(5)]
    hessian = R([3*d*d-8*c*e, 4*(c*d-6*b*e),
                 2*(2*c*c-24*a*e-3*b*d), 4*(b*c-6*a*d), 3*b*b-8*a*c])
    H = (4*phi*S(q.derivative(2)) + S(hessian.derivative(2)))/36
    H += QQ(2)/9*(I-phi*phi)
    require(H[2] == cubic[0], "H normalization mismatch")
    product = (m/cubic[0])*H
    return primitive(R([product[k].lift()[2] for k in range(3)]))

