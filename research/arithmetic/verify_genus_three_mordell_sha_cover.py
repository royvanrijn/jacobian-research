"""Exact characteristic-zero cover identities for HG8; SageMath 10.9.

Checks the family, every exceptional parameter via discriminants, map degrees
and signs. Selmer/Sha transfer and the use of MCT16 are written arguments.
"""
from sage.all import PolynomialRing, QQ

parameters = PolynomialRing(QQ, "a")
a = parameters.gen()
polynomials = PolynomialRing(parameters, "x")
x = polynomials.gen()
H = (x*x-1)*(1-a*(x*x-1)**3)
G = x*(x-1)*(1-a*(x-1)**3)
assert H.degree() == 8 and G.degree() == 5
assert H.discriminant() == 186624*a**5*(a+1)
assert G.discriminant() == -27*a*a*(a+1)**2
assert H(1) == H(-1) == 0
assert H.derivative()(1) == 2 and H.derivative()(-1) == -2

field = parameters.fraction_field()
functions = PolynomialRing(field, "x").fraction_field()
x = functions.gen()
H = functions(H)
U = 1/(x*x-1)
Y_coefficient = 1/(x*x-1)**2
assert H*Y_coefficient**2 == U**3-a
assert (x*x)*H == functions(G)(x*x)
assert max(U.numerator().degree(), U.denominator().degree()) == 2
assert max((x*x).numerator().degree(), (x*x).denominator().degree()) == 2
assert U(-x) == U and Y_coefficient(-x) == Y_coefficient
genus_two_ordinate_coefficient = x
assert genus_two_ordinate_coefficient(-x) == -genus_two_ordinate_coefficient
# The stated exception a=-1 is a real degeneration, rather than a discarded fixture.
singular = polynomials((polynomials.gen()**2-1)*(1+(polynomials.gen()**2-1)**3))
assert singular.gcd(singular.derivative()).degree() > 0
print("PASS: genus-three/genus-two discriminants and all excluded parameters;")
print("two degree-two quotient maps, elliptic equation Y^2=X^3-a, and signs.")
print("At a=2*R^2, R in Q*, both curves are smooth. No companion rank or Sha computed.")
