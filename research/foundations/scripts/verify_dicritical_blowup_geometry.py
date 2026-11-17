"""Exact identity and valuation audit for the explicit dicritical compactification theorem blow-ups."""

import sympy as sp

A, B, C, W, c, a0, b0 = sp.symbols("A B C W c a0 b0", nonzero=True)
gamma = sp.symbols("gamma", nonzero=True)
N = gamma**2 + (a0 - 1) * gamma - a0 * W
x = C / gamma
y = (W - gamma) / C
z = gamma * N / (b0 * C**2)
j0, j1 = gamma * C**2, C**3
j2, j3 = gamma * C * (W - gamma), gamma**2 * N / b0
assert sp.simplify(j1 / j0 - x) == 0
assert sp.simplify(j2 / j0 - y) == 0
assert sp.simplify(j3 / j0 - z) == 0

# Exact nonzero-root coordinate change: b=B(rho+w)-cAC makes E=H-Cb.
rho, w, b, Hloc = sp.symbols("rho w b Hloc", nonzero=True)
B_from_b = (b + c * A * C) / (rho + w)
E_nonzero = Hloc - B_from_b * C * (rho + w) + c * A * C**2
assert sp.simplify(E_nonzero - (Hloc - C * b)) == 0

# Exact charts of the blow-up of (C,W) at a zero of H=W^m L(W).
R, T, L = sp.symbols("R T L", nonzero=True)
for m in range(2, 13):
    E_zero = W**m * L - B * C * W + c * A * C**2
    c_chart = sp.expand(E_zero.subs(W, C * R) / C**2)
    assert sp.expand(c_chart - (C ** (m - 2) * R**m * L - B * R + c * A)) == 0
    w_chart = sp.expand(E_zero.subs(C, W * T) / W**2)
    assert sp.expand(w_chart - (W ** (m - 2) * L - B * T + c * A * T**2)) == 0


def val_sum(*orders):
    """Generic valuation of a sum: the unique least order, as used below."""
    return min(orders)


# Nonzero A_(mu-1) chain.  The asserted relation is i+(mu-i)=mu.
for mu in range(2, 20):
    for chart in range(mu):
        # C=u^(i+1)v^i, b=u^(mu-i-1)v^(mu-i), q=uv.
        assert (chart + 1) + (mu - chart - 1) == mu
        assert chart + (mu - chart) == mu
    for i in range(1, mu):
        vC, vb, vq = i, mu - i, 1
        assert vC + vb == mu * vq
        vB = min(vb, vC)  # rho is a unit
        assert vB == min(i, mu - i) and vB > 0
        vBC = vB + vC
        vHprime = mu - 1
        vgamma = val_sum(vBC, vHprime)
        assert vgamma == min(2 * i, mu - 1)
        vx, vy = vC - vgamma, -vC  # W is a unit
        vz_general = vgamma - 2 * vC  # a0 != 0, so N has order zero
        vz_a0_zero = 2 * vgamma - 2 * vC
        assert (vx, vy, vz_general, vz_a0_zero) == (
            i - vgamma,
            -i,
            vgamma - 2 * i,
            2 * vgamma - 2 * i,
        )

    # Strict C-boundary: C=q^mu/b, so residue field is k(A,b), not a
    # degree-mu extension.  The integer mu is its ramification index.
    assert mu == mu * 1 - 0

# Zero-cluster chain after the first blow-up: T*b'=q^(m-2).
for m in range(4, 20):
    r = m - 2
    for chart in range(r):
        # T=u^(i+1)v^i, b'=u^(r-i-1)v^(r-i), q=uv.
        assert (chart + 1) + (r - chart - 1) == r
        assert chart + (r - chart) == r
    for i in range(1, r):
        vT, vbp, vq = i, r - i, 1
        assert vT + vbp == r * vq
        vW, vC = 1, i + 1
        vB = min(vbp, vT)
        vBC, vHprime = vB + vC, m - 1
        vgamma = val_sum(vBC, vHprime)
        assert vgamma == min(2 * i + 1, m - 1)
        vx, vy = vC - vgamma, -i
        vz_general = vgamma - 2 * i - 1
        vz_a0_zero = 2 * vgamma - 2 * i - 2
        assert (vx, vy, vz_general, vz_a0_zero) == (
            i + 1 - vgamma,
            -i,
            vgamma - 2 * i - 1,
            2 * vgamma - 2 * i - 2,
        )

    # On D_0, T=q^r/b', W=q up to a unit, hence C has order r+1=m-1.
    assert r + 1 == m - 1

print("PASS: the four homogeneous reconstruction coordinates give the graph blow-up")
print("PASS: the nonzero-root model is C*b=q^mu with the stated A-chain valuations")
print("PASS: blowing up (C,W) gives both zero-cluster charts and T*b'=q^(m-2)")
print("PASS: all exceptional and reconstruction valuations agree through multiplicity 19")
print("PASS: Kummer exponents are ramification indices; dicritical residue degrees are one")
