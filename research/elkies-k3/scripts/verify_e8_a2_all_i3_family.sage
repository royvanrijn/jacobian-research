from sage.all import PolynomialRing, QQ


names = (
    "t", "lam", "q0", "q1", "q2", "a0", "a1", "b0", "b1", "c"
)
R = PolynomialRing(QQ, names)
t, lam, q0, q1, q2, a0, a1, b0, b1, c = R.gens()

P = t * (t - 1) * (t - lam)
q = q0 + q1 * t + q2 * t**2
a = a0 + a1 * t
b = b0 + b1 * t

A = -3 * q**2 + P * a
B = 2 * q**3 - P * q * a + P**2 * b
relation = 12 * q * b - a**2 - c * P
residual = 9 * q**2 * c + 4 * a**3 - 54 * q * a * b + 27 * P * b**2
F = 4 * A**3 + 27 * B**2

# This identity is valid before imposing relation=0.  On that exact
# coefficient locus, the discriminant acquires the three required cubic
# factors P^3.
assert F - P**3 * residual == 9 * P**2 * q**2 * relation
assert A.degree(t) == 4
assert B.degree(t) == 7
assert residual.degree(t) == 5

# Exact CM boundary point: P=t^2(t-1), q=a=0, b=t.
cm = {
    lam: 0,
    q0: 0,
    q1: 0,
    q2: 0,
    a0: 0,
    a1: 0,
    b0: 0,
    b1: 1,
    c: 0,
}
assert relation.subs(cm) == 0
assert A.subs(cm) == 0
assert B.subs(cm) == t**5 * (t - 1) ** 2

# A concrete non-isotrivial open-stratum member.  This also checks that the
# coefficient locus is nonempty away from the CM boundary.
sample = {
    lam: QQ(-1) / 2,
    q0: 1,
    q1: 1,
    q2: 1,
    a0: 1,
    a1: 2,
    b0: QQ(1) / 12,
    b1: QQ(1) / 6,
    c: 2,
}
assert relation.subs(sample) == 0
Pt = PolynomialRing(QQ, "t")
ts = Pt.gen()
As = Pt(A.subs(sample))
Bs = Pt(B.subs(sample))
Ps = Pt(P.subs(sample))
Rs = Pt(residual.subs(sample))
assert Ps == ts * (ts - 1) * (ts + QQ(1) / 2)
assert Rs.degree() == 5 and Rs.is_squarefree()
assert Ps.gcd(Rs) == 1
j_sample = QQ(6912) * As**3 / (4 * As**3 + 27 * Bs**2)
assert j_sample.derivative() != 0

print("E8A2I3|A=-3*q^2+P*a")
print("E8A2I3|B=2*q^3-P*q*a+P^2*b")
print("E8A2I3|coefficient_relation=12*q*b-a^2-c*P=0")
print("E8A2I3|Delta=-16*P^3*R5|R5={}".format(residual))
print("E8A2I3|generic_fibers=II*+I3+I3+I3+5I1")
print("E8A2I3|moduli_parameters=5|weierstrass_scaling=1|dimension=4")
print("E8A2I3|sample=lambda=-1/2|residual_squarefree=1|j_nonconstant=1")
print("E8A2I3|cm_endpoint=lambda=0,q=0,a=0,b=t,c=0|fibers=II*+II*+IV")
print("E8A2I3|status=PASS")
