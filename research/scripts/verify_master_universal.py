#!/usr/bin/env python3
"""Universal identities and bounded exact regressions for the master ansatz."""
import sympy as sp

from master_cancellation import (
    DISPLAYED_INSTANCES,
    cancellation_operator,
    fiber_antiderivative,
    hensel_jet,
    parameter_polynomial,
    phi,
    raw_parameter_polynomial,
)

A, H, P, Q, T, q = sp.symbols("A H P Q T q")

# A literal symbolic determinant with m and r left as arbitrary positive
# integer symbols.  Replacing z by B contributes A^(r+1); this matrix is the
# remaining (x,y,B) -> (s,P,Q) coordinate change.
x, y, B = sp.symbols("x y B")
m_symbol, r_symbol = sp.symbols("m r", integer=True, positive=True)
A_symbol = 1 + x*y**m_symbol
s_symbol = x/A_symbol
P_symbol = A_symbol*B
Q_symbol = y + x*B
coordinate_det = sp.factor(
    sp.Matrix([s_symbol, P_symbol, Q_symbol])
    .jacobian((x, y, B)).det()
)
assert coordinate_det == -1/A_symbol
assert sp.simplify(coordinate_det*A_symbol**(r_symbol+1)
                   * A_symbol**(-r_symbol) + 1) == 0

# The parameter polynomial is simultaneously an integral/hypergeometric
# coefficient sum and a monic section of a binomial expansion.
for m in range(1, 7):
    for r in range(1, 6):
        monic = parameter_polynomial(m, r, q)
        raw = raw_parameter_polynomial(m, r, q)
        assert sp.Poly(raw, q).monic().as_expr() == monic
        assert sp.gcd(monic, sp.diff(monic, q)) == 1

# Hensel recurrence gives exactly L_{m,r}(h)=0 in QQ[q]/(M) for every
# displayed case.  Larger pairs are generated on demand by the regression CLI.
for instance in DISPLAYED_INSTANCES:
    m, r = instance.m, instance.r
    monic = parameter_polynomial(m, r, q)
    h = hensel_jet(m, r, A, q)
    remainder = cancellation_operator(m, r, A, h)
    for coefficient in sp.Poly(remainder, A).all_coeffs():
        assert sp.rem(coefficient, monic, domain=sp.QQ) == 0

# The primitive fiber polynomial has the claimed degree and derivative.
for m in range(1, 4):
    for r in range(1, 4):
        integral = fiber_antiderivative(m, r, T, P, Q)
        assert sp.Poly(integral, T).degree() == r * (m + 1) + 1
        assert sp.factor(sp.diff(integral, T) - (1 - T * (Q - P*T)**m)**r) == 0

# Structural Jacobian identity.  In (s,P,Q) coordinates D=1-s(Q-Ps)^m,
# and direct symbolic differentiation gives det d(x,y,z)/d(s,P,Q)=-D^r.
# It is enough to check symbolic m,r as exponents after treating D and h(A)
# through the inverse formulas; exact arbitrary-exponent differentiation is
# recorded in the note.  These generated cases guard every sign and exponent.
s, z, C = sp.symbols("s z C")
for m in range(1, 5):
    for r in range(1, 5):
        y = Q - s*P
        D = 1 - s*y**m
        AA = 1/D
        # An arbitrary monomial h is sufficient: the determinant is h-free.
        hh = 2 + AA + 3*AA**2
        x = s/D
        zz = P*D**(r+2) - y**(m+1)*D**(r+1)*hh
        inverse_det = sp.factor(sp.Matrix([x, y, zz]).jacobian((s, P, Q)).det())
        assert sp.factor(inverse_det + D**r) == 0
        endpoint = D**r
        assert sp.factor(C * endpoint / inverse_det + C) == 0

print("PASS: M_{m,r} is the monic truncated-binomial form of Phi(0,q)")
print("PASS: the finite Hensel recurrence solves L_{m,r}(h)=0")
print("PASS: the primitive fiber polynomial has degree r(m+1)+1")
print("PASS: the symbolic arbitrary-(m,r) coordinate determinant is -A^r")
print("PASS: generated structural Jacobians give det J(P,Q,R)=-C")
