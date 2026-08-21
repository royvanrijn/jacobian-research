from sage.all import PolynomialRing, QQ


R = PolynomialRing(QQ, ("t", "r", "lam"))
t, r, lam = R.gens()

D = t * (t - 1)
A = -3 * r**2 * D**2
B = D**2 * ((t - lam) ** 3 - 2 * r**3 * D)
residual = (t - lam) ** 3 - 4 * r**3 * D
Delta = -16 * (4 * A**3 + 27 * B**2)
expected = -432 * D**4 * (t - lam) ** 3 * residual
assert Delta == expected

# K3 orders at infinity in the O(8), O(12), O(24) coefficient bundles.
assert A.degree(t) == 4
assert B.degree(t) == 7
assert Delta.degree(t) == 14
assert 8 - A.degree(t) == 4
assert 12 - B.degree(t) == 5
assert 24 - Delta.degree(t) == 10

# The residual cubic is squarefree on this explicit open set.  Its roots are
# automatically disjoint from 0, 1, and lam when r*lam*(lam-1) != 0.
residual_discriminant = residual.discriminant(t)
expected_residual_discriminant = -16 * r**6 * (
    -16 * r**6
    + 16 * r**3 * lam**3
    - 24 * r**3 * lam**2
    - 24 * r**3 * lam
    + 27 * lam**4
    + 16 * r**3
    - 54 * lam**3
    + 27 * lam**2
)
assert residual_discriminant == expected_residual_discriminant

# The j-invariant is nonconstant on the generic open set.
K = R.fraction_field()
j = K(6912 * A**3) / K(4 * A**3 + 27 * B**2)
assert j.derivative(t) != 0

# At (r,lam)=(0,0), recover Utsumi No.1 in the affine t-chart.
cm_substitution = {r: 0, lam: 0}
assert A.subs(cm_substitution) == 0
assert B.subs(cm_substitution) == t**5 * (t - 1) ** 2

open_factor = (
    r
    * lam
    * (lam - 1)
    * (
        -16 * r**6
        + 16 * r**3 * lam**3
        - 24 * r**3 * lam**2
        - 24 * r**3 * lam
        + 27 * lam**4
        + 16 * r**3
        - 54 * lam**3
        + 27 * lam**2
    )
)

print("E8A2MIX|A={}".format(A))
print("E8A2MIX|B={}".format(B))
print("E8A2MIX|Delta={}".format(Delta.factor()))
print("E8A2MIX|generic_fibers=II*+IV+IV+I3+3I1|j_nonconstant=1")
print("E8A2MIX|open_factor={}".format(open_factor))
print("E8A2MIX|cm_endpoint=r=0,lam=0|fibers=II*+II*+IV")
print("E8A2MIX|status=PASS")
