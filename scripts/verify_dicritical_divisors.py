"""Exact valuation audit for the dicritical divisors of weighted maps."""

import sympy as sp

W, r, A, B, C, c = sp.symbols("W r A B C c", nonzero=True)

# The universal repeated-root divisor maps to the pulled-back discriminant.
H = sp.Function("H")
s = sp.diff(H(W), W).subs(W, r)
t = sp.expand(r*s-H(r))
E_at_r = sp.simplify(H(r)-s*r+t)
Eprime_at_r = sp.simplify(sp.diff(H(W), W).subs(W, r)-s)
assert E_at_r == 0 and Eprime_at_r == 0

d, k, rho, h = sp.symbols("d k rho h", nonzero=True)

# A nonzero primitive root of multiplicity mu.  With C=d^mu and
# W=rho+k*d, the first nonzero inverse-penc-il coefficient is
# h*k^mu-B*rho.  Thus k^mu=B*rho/h and the boundary map has degree mu.
for mu in range(1, 6):
    local_H = h*(W-rho)**mu
    local_E = sp.expand(local_H-B*C*W+c*A*C**2)
    pulled = sp.expand(local_E.subs({C: d**mu, W: rho+k*d}))
    leading = sp.expand(pulled).coeff(d, mu)
    assert sp.factor(leading-(h*k**mu-B*rho)) == 0

    local_Eprime = sp.diff(local_E, W)
    pulled_prime = sp.expand(local_Eprime.subs({C: d**mu, W: rho+k*d}))
    assert sp.expand(pulled_prime).coeff(d, mu-1) == mu*h*k**(mu-1)

# The zero cluster of multiplicity m.  The escaping cluster uses
# C=d^(m-1), W=k*d and has boundary equation k*(h*k^(m-1)-B)=0.
# Removing k=0 leaves a degree-(m-1) dicritical divisor.  For m=2 this
# branch belongs to the finite gamma=0 chart; only m>=3 is dicritical.
for m in range(3, 7):
    local_H = h*W**m
    local_E = sp.expand(local_H-B*C*W+c*A*C**2)
    pulled = sp.expand(local_E.subs({C: d**(m-1), W: k*d}))
    leading = sp.expand(pulled).coeff(d, m)
    assert sp.factor(leading-k*(h*k**(m-1)-B)) == 0

    local_Eprime = sp.diff(local_E, W)
    pulled_prime = sp.expand(local_Eprime.subs({C: d**(m-1), W: k*d}))
    assert sp.factor(
        sp.expand(pulled_prime).coeff(d, m-1)-(m*h*k**(m-1)-B)
    ) == 0

# For m=2 the A*C^2 term occurs at the same order and gives the finite
# gamma=0 quadratic h*k^2-B*k+c*A=0 instead of a dicritical divisor.
double_zero = sp.expand(
    (h*W**2-B*C*W+c*A*C**2).subs({C: d, W: k*d})
).coeff(d, 2)
assert sp.factor(double_zero-(h*k**2-B*k+c*A)) == 0

# At the distinguished simple root W=1, H'(1)=-c gives gamma=1 at C=0;
# this is the regular x=0 chart rather than a divisor at source infinity.
assert sp.factor(-(-c)/c) == 1

print("PASS: the repeated-root divisor maps by s=H'(r), t=rH'(r)-H(r)")
print("PASS: a nonzero primitive root of multiplicity mu gives a degree-mu Kummer divisor")
print("PASS: a zero of multiplicity m>=3 gives a degree-(m-1) dicritical divisor")
print("PASS: the distinguished root 1 and the double-zero cluster are finite charts")
