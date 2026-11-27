"""Exact symmetric-line, alternate-isogeny and paired-field checks for SL1–SL3."""
from sage.all import GF, QQ, ZZ, Matrix, PolynomialRing, pari

parameters = PolynomialRing(QQ, "r").fraction_field()
r = parameters.gen()
one_variable = PolynomialRing(parameters, "t")
t0 = one_variable.gen()
F0 = (1-4*r**3)*t0**6 + 12*r*(3*r**3-1)*t0**4 + 24*r*r*t0*t0 + 4
assert F0.discriminant() == 2**16*3**6*(4*r**3-1)*(8*r**3-1)**6
rational = one_variable.fraction_field()
t = rational(t0)
F = rational(F0)
quadratic_ring = PolynomialRing(rational, "z")
z = quadratic_ring.gen()
field = rational.extension(z*z-F, "y")
y = field.gen()
a, b = 1-4*r**3, 8*r**3-1
alpha1 = (t**3-6*r*r*t*t-2+y)/2
alpha2 = (t**3+6*r*r*t*t+2+y)/2
gamma = (t**3-6*r*r*t*t-2-y)*alpha2/2
eta = alpha1*alpha2
X1, Y1 = a*t*t-4*r, a*y
X2, Y2 = 4*(t**(-2)-r*r), 4*y/t**3
assert Y1*Y1 == X1**3+(6*r*r*X1+2*b)**2
assert Y2*Y2 == X2**3+(6*r*X2+4*b)**2
assert (Y1+6*r*r*X1+2*b)*gamma == (r*y+(6*r**3-1)*t*t+2*r)**3
assert (Y2+6*r*X2+4*b)*eta == (-2*(r*r*t*t-1))**3
Q, P, S = 4*r*r+2*r+1, 3*r*r+2*r+1, 6*r**3+6*r*r+4*r+1
Z = X2+4*Q
L = 2*(r+1)*Z-4*Q
assert Z == 4*(t**(-2)+P) and L == 8*(r+1)*t**(-2)+4*S
assert Y2*Y2 == Z**3-3*L*L
mirror = rational.extension(z*z+3*F, "rho")
rho = mirror.gen()
xi = -L/3+4*rho/(9*t**3)
conjugate = -L/3-4*rho/(9*t**3)
assert xi*conjugate == (Z/3)**3
assert xi+conjugate == -2*L/3

poly_ring = PolynomialRing(parameters, "U")
U = poly_ring.gen()
outer = U**3-Q*U+(2*r+1)*Q/3
assert outer.discriminant() == Q*Q*(2*r-1)**2
shanks_parameter = 6*(r+1)/(2*r-1)
scale = -(2*r-1)/3
shifted = U/scale + shanks_parameter/3
assert outer == scale**3*(shifted**3-shanks_parameter*shifted**2-(shanks_parameter+3)*shifted-1)
assert -36*(r*r+P)*6/6**3 == -Q
assert (144*(r+1)*r*r+72*S)/6**3 == (2*r+1)*Q/3

# An exact instance with a new cyclic outer field, of conductor 211.
ring = PolynomialRing(QQ, "T")
T = ring.gen()
p, h, r = ZZ(223), ZZ(18), ZZ(7)
assert p % 3 == 1 and h % 18 == 0
left, right = r*p-h, r*p+h
assert all(z.is_prime() for z in (left,p,right))
A, B, C = left**3,p**3,right**3
radicand = (A+C-B)**2-4*A*C
D = radicand.squarefree_part()
assert D == -99408409450866651 and D % 8 == 5 and D.valuation(3) == 1
traces = (A+B-C,B+C-A)
character_rows = []
for q, yq, root, exponents in [(37,15,26,(2,2)),(163,41,104,(1,0))]:
    finite = GF(q)
    yq, root = finite(yq),finite(root)
    assert yq*yq == finite(radicand)
    assert root**3 == 1 and root != 1
    values = [(finite(trace)+yq)/2 for trace in traces]
    assert tuple(value**((q-1)//3) for value in values) == tuple(root**e for e in exponents)
    character_rows.append(exponents)
assert Matrix(GF(3),character_rows).rank() == 2

Q, P, S = 4*r*r+2*r+1,3*r*r+2*r+1,6*r**3+6*r*r+4*r+1
outer_polynomial = T**3-Q*T+(2*r+1)*Q/3
assert outer_polynomial == T**3-211*T+1055
assert pari(outer_polynomial).nfdisc() == 211**2
first = pari(right**3*(T**3+1)-(left*T+p)**3)
second = pari(T**3-36*(h*h+P*p*p)*T+144*(r+1)*p*h*h+72*S*p**3)
assert first.nfdisc() == -D//3 and second.nfdisc() == 81*D
second_field = second.nfinit()
assert [(int(row[2]),int(row[3])) for row in second_field.idealprimedec(3)] == [(3,1)]
for q, expected in [(left,[(1,3)]), (p,[(1,1)]*3), (right,[(1,1)]*3)]:
    assert sorted((int(row[2]),int(row[3])) for row in second_field.idealprimedec(q)) == expected

print("Verified the general symmetric-line formulas and both cubic companions at r=7.")
