"""Exact HG15 coefficient, exceptional-case and rank-zero certificates.

Run with SageMath 10.9 from the repository root. This does not construct
a genus-three curve with three rational torsion directions.
"""
from sage.all import (QQ, PolynomialRing, EllipticCurve, matrix, vector)

A = PolynomialRing(QQ, names=("P", "t", "S", "d"))
P, t, S, d = A.gens()
K = A.fraction_field()
B = PolynomialRing(A, "X")
X = B.gen()
Q, R = (S+d)/2, (S-d)/2
p = ((X-t)**3+P*(X-t)-Q)*((X+t)**3+P*(X+t)-R)
g = X**2-t**2+2*P/3
r = g**3-p
expected = B([
    4*P*t**4-P*P*t*t/3+d*(t**3+P*t)+8*P**3/27-(S*S-d*d)/4,
    (3*t*t+P)*S,
    -4*P*t*t+3*d*t+P*P/3,
    S])
assert r == expected
assert p[4] == -3*t*t+2*P
disc = A(r.discriminant())
assert disc.degree(t) == 10
assert disc.coefficient({t: 10}) == 1024*P**4
print("Generic sextic/cubic identities and degree-ten discriminant passed.")

rx = A(r.derivative()(1))
rt = A(r.derivative(t)(1))
M = matrix(K, [[rx.derivative(z) for z in (S, d)],
               [rt.derivative(z) for z in (S, d)]])
H = K((P+3*(t-1)**2)*(P+3*(t+1)**2))
assert M.det() == H
sol = M.solve_right(vector(K, [-rx(S=0, d=0), -rt(S=0, d=0)]))
assert sol[0] == 2*P*(P+12*t*t)*(15*t*t-P-3)/(3*H)
assert sol[1] == (-2*t*P*(72*t**4+21*P*t*t-P*P+108*t*t-21*P-36)
                  /(3*H))
sub = A.hom([K(P), K(t), sol[0], sol[1]], K)
assert sub(rx) == sub(rt) == 0
assert sub(A(r(1))) == -P*(P+12*t*t)*(3*t*t-2*P-3)**3/(27*H)
assert g(1) == -(3*t*t-2*P-3)/3

rxp = rx.subs(P=-3*(t-1)**2)
rtp = rt.subs(P=-3*(t-1)**2)
rxm = rx.subs(P=-3*(t+1)**2)
rtm = rt.subs(P=-3*(t+1)**2)
assert rtp-rxp == -6*(t+1)*(t-1)**2*(3*t-1)**2
assert rtm+rxm == -6*(t-1)*(t+1)**2*(3*t+1)**2
assert rxp.subs(t=-1) == 6*(32-S-d)
assert rxm.subs(t=1) == 6*(32-S+d)
assert X**3-12*X-16 == (X-4)*(X+2)**2
assert g(1).subs({P: -QQ(4)/3, t: QQ(1)/3}) == 0
assert g(1).subs({P: -QQ(4)/3, t: -QQ(1)/3}) == 0
print("Nonsingular elimination and every singular-system case passed.")

# First variation of a cubic discriminant at a double and a simple root.
C = PolynomialRing(QQ, names=("a", "b", "h0", "h1", "h2", "h3", "e"))
aa, bb, hh0, hh1, hh2, hh3, ee = C.gens()
CX = PolynomialRing(C, "x")
xx = CX.gen()
perturbation = hh0+hh1*xx+hh2*xx*xx+hh3*xx**3
variation = C((aa*xx**2*(xx-bb)+ee*perturbation).discriminant())
assert variation.derivative(ee).subs({ee: 0}) == 4*aa**3*bb**3*hh0

# Newton's sum-of-cubes calculation in the collision pairs.
C = PolynomialRing(QQ, names=("xi", "tau"))
xi, tau = C.gens()
sum1, sum2 = 3*tau-xi, -3*tau-xi
prod1, prod2 = (sum1**2-2*xi*xi)/2, (sum2**2-2*xi*xi)/2
cubes = 2*xi**3+sum1**3-3*prod1*sum1+sum2**3-3*prod2*sum2
assert cubes/3 == xi*(9*tau*tau-xi*xi)
assert xi*xi-4*(-xi*xi/2) == 3*xi*xi
print("Discriminant variation and collision-to-Fermat identities passed.")

E = EllipticCurve(QQ, [0, 1]).minimal_model()
assert E.a_invariants() == (0, 0, 0, 0, 1)
assert E.rank(proof=True, use_database=False) == 0
assert E.rank_bound(algorithm="mwrank") == 0
assert E.torsion_subgroup().invariants() == (6,)
points = {tuple(z) for z in E.torsion_points()}
assert points == {(QQ(0), QQ(1), QQ(0)),
                  (QQ(-1), QQ(0), QQ(1)),
                  (QQ(0), QQ(-1), QQ(1)),
                  (QQ(0), QQ(1), QQ(1)),
                  (QQ(2), QQ(-3), QQ(1)),
                  (QQ(2), QQ(3), QQ(1))}
assert {(x, y) for x, y, z in points if z and x and abs(y) <= 1} == {(-1, 0)}
print("Fermat cubic: exact rank zero and all six rational points passed.")
print("HG15 internal exact certificate passed; no independent review asserted.")
