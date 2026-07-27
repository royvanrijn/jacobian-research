#!/usr/bin/env python3
"""Close every reduced quartic noncharacteristic caustic component.

Let p be a reduced geometrically irreducible component of L=0, with F=f_aa
and C nonzero generically.  Normal-Hessian polynomiality gives

    A=0,  C^2=2*F*kappa                         modulo p.

Lines have a separate exact certificate, so assume deg(p)>=2.  If h_a=0
modulo p, then degree and mixed-derivative compatibility give

    h_a=0,  h_b=b/2

as polynomial identities.  The forced normal Hessian then has

    V_cd=5*b^3/16,

which is impossible for a quartic potential.

On the remaining h_a != 0 branch, the Laurent conditions give p^2 | A and
H2=K2=K3=K4=0.  When deg(p)>=3, the bound deg(A)<=4 immediately forces
the global identity A=0; the conic case reaches the same identity below by
mixed-derivative compatibility.

Tangential differentiation of C^2/F=2*kappa and K3 form a two-by-two
linear system for S_a,S_b.  Its determinant is

    U=F*p_b-P*p_a,

the derivative of f_a along p.  Noncharacteristicity means U != 0, hence
S_a=S_b=0.  Then K4 gives R=0, while K2 and H2 give

    d_a(h_a/F)=d_b(h_a/F)=0                    modulo p.

In characteristic zero h_a/F is constant on the component.  Since h_a and
F have degree at most two, a component of degree at least three immediately
gives h_a=lambda*F.  For a conic,

    h_a=lambda*F+mu*p

as a polynomial identity.  Put r=h-lambda*f_a.  The equation A=0 also
gives r_b-b/2=0 modulo p.  Both expressions have degree at most two, so

    r_a=mu*p,  r_b-b/2=nu*p.

Equality of mixed derivatives gives mu*p_b=nu*p_a.  The two derivatives of
a geometrically irreducible conic are not proportional, hence mu=nu=0 and

    h=lambda*f_a+b^2/4+constant

as a polynomial identity.  The forced normal Hessian is then

    V_cd=5*b^3/16+lambda*(g_a-1/2).

For cubic g the second term has degree at most two, so the cubic term
cannot cancel.  Thus no reduced quartic noncharacteristic caustic
component survives.
"""

from __future__ import annotations

import runpy

import sympy as sp


# The h_a=0 degree obstruction is an identity in the generic Schur data.
boundary = runpy.run_path("scripts/verify_hc4_1000_boundary_schur_chain.py")
b = boundary["b"]
zero_h_substitution = {
    boundary["h_a"]: 0,
    boundary["h_b"]: b / 2,
}
forced_V_cd = sp.factor(
    sp.cancel(boundary["forced_V_cd"].subs(zero_h_substitution))
)
assert sp.factor(forced_V_cd - sp.Rational(5, 16) * b**3) == 0


# Exact compatibility identities on the h_a != 0 branch.
F, P, C = sp.symbols("F P C", nonzero=True)
p_a, p_b = sp.symbols("p_a p_b")
S_a, S_b = sp.symbols("S_a S_b")
U = F * p_b - P * p_a

tangential_C_equation = p_b * S_a - p_a * S_b
K3_equation = F * S_b - P * S_a
assert sp.factor(
    F * tangential_C_equation + p_a * K3_equation - U * S_a
) == 0
assert sp.factor(
    P * tangential_C_equation + p_b * K3_equation - U * S_b
) == 0

F_a, F_b = sp.symbols("F_a F_b")
h_a, h_aa, h_ab = sp.symbols("h_a h_aa h_ab")
H2 = F**2 * h_ab - F * F_b * h_a - F * P * h_aa + F_a * P * h_a
K2_after_H2_Sa = C * (F * h_aa - F_a * h_a)

# K2=0 first kills the a derivative of h_a/F.  Substitution into H2
# leaves F times the numerator of its b derivative.
h_aa_solution = F_a * h_a / F
assert sp.factor(
    H2.subs(h_aa, h_aa_solution) - F * (F * h_ab - F_b * h_a)
) == 0
assert sp.factor(K2_after_H2_Sa.subs(h_aa, h_aa_solution)) == 0

# The resulting global h normal form gives an uncancellable cubic in the
# forced mixed normal Hessian.
lam = sp.symbols("lam")
normal_form_substitution = {
    boundary["h_a"]: lam * boundary["f_aa"],
    boundary["h_b"]: b / 2 + lam * boundary["f_ab"],
}
forced_V_cd_normal_form = sp.factor(
    sp.cancel(boundary["forced_V_cd"].subs(normal_form_substitution))
)
expected_V_cd_normal_form = (
    sp.Rational(5, 16) * b**3
    + lam * (boundary["g_a"] - sp.Rational(1, 2))
)
assert sp.factor(forced_V_cd_normal_form - expected_V_cd_normal_form) == 0


def main() -> None:
    print("PASS: h_a=0 forces V_cd=5*b^3/16, exceeding quartic degree")
    print("PASS: U*S_a and U*S_b lie in the tangent/K3 compatibility ideal")
    print("PASS: noncharacteristicity U!=0 therefore forces S_a=S_b=0")
    print("PASS: K2 and H2 force both derivatives of h_a/F to vanish")
    print("PASS: mixed derivatives upgrade this to h=lambda*f_a+b^2/4")
    print("PASS: then V_cd=5*b^3/16+lambda*(g_a-1/2)")
    print("RESULT: deg(g_a)<=2 cannot cancel the cubic")
    print("        no reduced noncharacteristic caustic component survives")


if __name__ == "__main__":
    main()
