from sage.all import PolynomialRing, QQ, cyclotomic_polynomial


R = PolynomialRing(QQ, ("t", "lam", "mu"))
t, lam, mu = R.gens()

# The previously proposed all-IV lift of the abstract E8 + A2^3 root frame.
B = (t * (t - 1) * (t - lam)) ** 2 * (t - mu)
Delta = -432 * B**2
expected_delta = (
    -432
    * t**4
    * (t - 1) ** 4
    * (t - lam) ** 4
    * (t - mu) ** 2
)
assert Delta == expected_delta
assert B.degree(t) == 7

# In the K3 short-Weierstrass bundle, deg(B)=7 means ord_infinity(B)=5;
# A=0 has infinite order and ord_infinity(Delta)=10.  Thus the generic
# configuration is II* + 3 IV + II and j=0.
ord_B_infinity = 12 - B.degree(t)
ord_Delta_infinity = 24 - Delta.degree(t)
assert ord_B_infinity == 5
assert ord_Delta_infinity == 10

root_rank = 8 + 3 * 2
rho_target = 19
mw_rank_target = rho_target - 2 - root_rank
assert root_rank == 14
assert mw_rank_target == 3

# Over Qbar(t), (x,y) -> (zeta_3*x,y) acts on the MW group and satisfies
# Phi_3(rho)=rho^2+rho+1=0.  Hence MW tensor Q is a vector space over the
# quadratic field Q(zeta_3), so its Q-dimension is even.
phi3 = cyclotomic_polynomial(3)
assert phi3.degree() == 2
assert phi3.is_irreducible()
assert mw_rank_target % phi3.degree() != 0

print("E8A2PARITY|family_delta={}".format(Delta.factor()))
print("E8A2PARITY|j=0|cm_automorphism_order=3|mw_rank_parity=even")
print(
    "E8A2PARITY|target_root_rank={}|target_rho={}|target_mw_rank={}".format(
        root_rank, rho_target, mw_rank_target
    )
)
print("E8A2PARITY|status=OBSTRUCTED_ALL_IV_LIFT")
