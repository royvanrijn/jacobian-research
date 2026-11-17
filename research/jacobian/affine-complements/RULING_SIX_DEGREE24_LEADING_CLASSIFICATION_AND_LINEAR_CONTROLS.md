# Degree-24 leading coefficients have four exact configurations

PF2RULED177. 17 November 2026. Actual complex polynomial plane Keller maps,
within the complete marked quadratic source model.

At leading degree D=24, every hypothetical survivor has
m=2, deg h=11, one actual A support of slope five and intermediate degree
b=5. After scaling its noncentral A root to one, its leading coefficients
have the form

    M=constant*u^2, h=constant*(u-1)^10*(u-tau).

The leading primitive equation holds precisely when tau=0 or tau is one
of the three distinct roots of

    F(tau)=40tau^3-105tau^2+84tau-14.                     (1)

This classifies the LEADING equation, including for a quadratic remainder.
The subsequent [complete degree-24 exclusion](RULING_SIX_DEGREE24_CRITICAL_QUARTIC_POLE_BUDGET_EXCLUSION.md)
closes this entire class, using the additional quartic pole budget. The linear controls below
show that several finite necessary tests remain weaker than the complete
companion resonance of
[PF176](RULING_SIX_COMPANION_RESONANCE_AND_LINEAR_REMAINDER_CLOSURE.md),
which excludes the entire linear branch.

## 1. The complete integer ledger at degree 24

Use PF165--166 and PF173 notation

    b3=M*h^2, m=deg M>=2, H=deg h, D=m+2H,
    b=[C(u,w):C(P,Q)], n_A>=1, S=sum of actual A slopes.

The full retained-source ledger gives

    ceil((D+3m)/6)<=b<=2S-5n_A<=H-5n_A.                  (2)

Suppose D=24. Then m is even, and (2) gives m+5n_A<=8.
Thus m=2,n_A=1,H=11. Now b>=5 and S<=floor(H/2)=5 force

    S=L=5, b=5.                                        (3)

The central M-root has multiplicity at least two, so a degree-two M
has only that root. Hence M is a constant times u^2. The one actual
A root c is off M, and ord_c h=2L=10. There is exactly one remaining
linear factor in h. It cannot coincide with c, because the A order
is exactly ten, not eleven. Scale u to make c=1. Multiplicative
constants in M and h do not affect the logarithmic leading ODE.
Consequently

    h=constant*(u-1)^10*(u-tau), tau!=1.                 (4)

The case tau=0 is allowed for the leading equation. It is excluded
in the LINEAR remainder branch by PF175, since that branch requires
h(0)!=0. It is not removed for a quadratic remainder here.

## 2. Exact leading primitive equation for all remaining positions

PF173 gives a polynomial rho satisfying

    rho_u-(G_u/(3G))*rho=1, G=M^2*h.                    (5)

If tau=0, then G is a nonzero constant times u^5*(u-1)^10.
The explicit degree-five primitive in PF173 proves (5).

Suppose tau!=0,1. Then, up to a constant,

    G=u^4*(u-1)^10*(u-tau), deg G=15.

The polynomial rho has degree five and a simple zero at each of
0,1,tau. Thus write

    rho=D0*s, D0=u*(u-1)*(u-tau),
    s=c2*u^2+c1*u+c0.

After dividing the cleared equation by D0, equation (5) becomes

    3D0*s_u+[-6u^2+(8tau-1)u-tau]*s=3.                 (6)

Its constant, linear and quadratic coefficients uniquely give

    c0=-3/tau,
    c1=3(8tau-1)/(2tau^2),
    c2=-3(40tau^2-25tau+4)/(10tau^3).                   (7)

The degree-four term cancels. The sole remaining coefficient, at u^3,
is

    -3F(tau)/(5tau^3).                                  (8)

Therefore (5) holds if and only if F(tau)=0. Conversely (7) gives
the polynomial rho explicitly at every such root, so no existence
obligation remains in this coefficient classification.

The discriminant of F is -1428840, so its three roots are distinct.
They are different from 0,1 and 3/2: the corresponding values of F
are -14,5 and 43/4. All divisions below are therefore by nonzero
constants. Together with tau=0, these are the four configurations
claimed in the title, after the specified source scaling.

## 3. Uniform controls passing the earlier linear tests

Fix ANY of the three roots tau of F. Put

    v=u-1, delta=1-tau,
    h=v^10*(u-tau), M=u^2,
    e=(u-tau)*v^4*(3-2tau-u), R2=e*w+f.                 (9)

These are coefficient controls, not actual Keller companions.
Their exact identities are

    e^2/h=v-3delta+4delta^3/v^2
      =partial_u[v^2/2-3delta*v-4delta^3/v].              (10)

Thus the third-character condition of PF175 holds; its first character
is zero in this linear branch. At the A root, ord_1 e=4<L=5, so
the pole-support condition holds. At the other h-root,

    ord_tau e=1=ord_tau h.

Consequently the unit-K h-root order test is met whenever that test
applies there. At the centre,

    h(0)=-tau!=0, e(0)=tau*(2tau-3)!=0.                 (11)

The stronger central tests also hold. The leading ODE holds because
F(tau)=0, by Section 2. Hence no combination of just these finite
tests excludes the three configurations.

The controls can be inserted into an actual regular source-image
polynomial, without asserting a mate. On

    w^2=u+u^2*z

choose K=v^5, B0=1/(4delta) and

    H3=(h*w+K)*z+B0*w, P=H3^2+e*w+f.                   (12)

At the A root the leading cubic is

    delta*Y^3+Y^2+Y/(4delta)
      =delta*Y*(Y+1/(2delta))^2,

so the required square quartic-root face is retained. The central
restriction of H3 is nonconstant affine; P has degree two on the
central source line and at u=1, degree four at u=tau, and degree
six on the remaining whole fixed-source-coordinate lines. Global
submersion is not claimed by these reached-line checks.

## 4. The fifth character is proper and can still be repaired finitely

For a linear remainder and

    H3=a3*w^3+a2*w^2+a1*w+a0,

PF175 gives

    m5=(5/16)*res_w[(e*w+f)^3/H3].

Expanding through the needed inverse-cubic coefficient gives

    m5/M=(5/16)*(e^2/h)*(3f-eK/h),                      (13)

since a3=h/M and a2=K/M. This is an unrestricted identity for a
linear remainder. Its rational differential must have zero residues.

For the controls (9),(12), eK/h=2delta/v-1. If f is constant,
the rational function in parentheses with its e^2/h factor has
residue 6delta^2 at v=0. Hence the fifth character fails for EVERY
constant f in this family, despite all tests in Section 3 passing.

The failure can be repaired by the polynomial choice

    f=-v/(2delta).

The function S=(e^2/h)*(3f-eK/h) is then

    S=-3v^2/(2delta)+11v/2-5delta
                         +4delta^3/v^2-8delta^4/v^3,

with rational primitive

    -v^3/(2delta)+11v^2/4-5delta*v
                         -4delta^3/v+4delta^4/v^2.        (14)

Thus even the first, third and fifth characters, leading ODE,
central conditions and named support orders can pass together.
This construction is uniform in all three roots of (1), without
numerical approximation or a mate-degree cutoff.

The COMPLETE companion resonance of PF176 still excludes these
linear data. They are not controls satisfying that resonance, the
full odd hierarchy, global submersion or the Keller map condition.
These quadratic-remainder configurations are excluded by the subsequent
complete degree-24 theorem; the coefficient classification here supplies
one of its retained inputs.

## 5. Proof and verification boundary

The checker exhausts the integer ledger (2) at D=24, verifies the
entire parameter equation (6)--(8), its discriminant and excluded
parameter values, all uniform rational identities and source-image
faces, and the fifth-character failure and repair. The full leading
classification follows from the factorization and exact coefficient
comparison, not from sampling parameter values.

Replay:

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_degree24_leading.py

No target membership, target integrality, actual mate or full ruling-six
closure is inferred. Independent whole-proof replay, formal verification
and external human review are not claimed.
