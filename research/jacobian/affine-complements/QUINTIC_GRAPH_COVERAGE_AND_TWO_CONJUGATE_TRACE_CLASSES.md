# Quintic graph coverage and two conjugate residual trace classes

PF2RULED155. 17 November 2026.

For the v2 Meng--Yang parent, every degree-at-most-five polynomial graph
with constant nonzero pulled-back Hessian determinant has ZERO quintic
plane trace. Its quartic/cubic trace belongs to one of two explicitly
specified conjugate classes over Q(sqrt(34)). Consequently there is no
such graph with rational coefficients, or over any characteristic-zero
field not containing sqrt(34).

This is coverage of the entire degree-five GRAPH ansatz, not an expansion
of the excluded sparse architecture. The two complex residual classes are
not excluded. Global coverage of hypothetical minimal plane counterexamples
by these graph coordinates remains open, as do HC4/JC2. The result is
field-valued: it does not set these coefficients to zero on a nonreduced
trace scheme or provide a new uniform ideal nilpotence exponent.

## 1. Statement and exact remaining forms

Use the v2 parent Psi=A^2+13A+2B from the retained
[graph note](../../foundations/HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md),
in source variables (x,y,p,q), with r=R(x,y,p,q). Let K have characteristic
zero, deg R<=5 and det Hess(Psi|r=R) a nonzero constant. Then

    T=R|x0=D_s y^3p+h y^3+i y^2p+e yp^2+f p^3
                         +rho_s y^2q+T_le2(y,p,q),                 (1)
    rho_s=-123/20+s*3sqrt(34)/40,
    D_s=-47/10+s*3sqrt(34)/5,   s in {+1,-1}.

All displayed lower parameters and all ten coefficients of T_le2 are
unrestricted by this statement; later necessary equations still apply.
Equation (1) is necessary, not sufficient or an existence assertion.
In particular [y^5]T=0 and T_4=D_s y^3p. If sqrt(34) is absent from K,
there is no such graph.

Degree five here means degree of the graph R. It is not the degree-five
potential convention in the external preprint audit. The conclusion is
about unmarked graphs; the marked collision route was closed separately
in all degrees.

## 2. Complete leading-trace reduction and plane equations

The retained field-valued HC4MYG5J theorem gives EVERY candidate trace as

    T=kappa y^5+A y^4+D y^3p+h y^3+i y^2p+e yp^2+f p^3
            +rho y^2q+T_le2.                                      (2)

This keeps all 18 surviving trace parameters. Write

    a=[pq]T_le2, b=[yq]T_le2, c=[q^2]T_le2,
    W=2a+3, Rho=20rho+123, V=16rho-2D+89,
    P=8rho^2+99rho+279, Q=160rho^2+1968rho+6021.

The forced first normal coefficient is S=(F_plane(T)-Cdet)/64.
Here Cdet is the desired constant; it is distinct from the parameter c.
The degree-six and degree-five plane faces give

    V^2+2kappa*Rho=0,                                             (3)
    32Va+40a*kappa-8Ve-12Vrho-39V+60kappa+Q=0,                     (4)
    c(4V+5kappa)=0,  fV=0,                                      (5)
    2Vi-(16rho+99)A-(16V+20kappa)b=0.                            (6)

These equations are derived from the exact plane formula in
[PF153](HC4_TRACE_QUOTIENT_AND_UNIFORM_CONE_NILPOTENCE.md), without deleting
any of the ten degree-at-most-two coefficients.

Suppose for contradiction kappa!=0. The top potential is
(x^3R_5)^2. Its zero highest Hessian face and the retained low-dimensional
Gordan--Noether argument in
[HC4MYG5N](../../foundations/HC4_MENG_YANG_QUINTIC_SCHUR_FRONTIER.md)
give a constant direction (0,v_y,v_p,v_q) annihilating R_5. Restriction
to kappa*y^5 makes v_y=0. Therefore the forced quartic S_4 has a nonzero
constant kernel direction in the (p,q)-plane. This conditional argument
is valid over characteristic-zero fields; the older rationality hypothesis
was needed only to exclude kappa=0, not here.

The coefficients [p^4]S_4=9f^2 and [p^3q]S_4=0. If v_p!=0,
the p^3 coefficient in its
directional derivative forces f=0. If v_p=0, the direction is partial_q.
Its q*y^2 coefficient forces c=0, and then its y^3 coefficient gives
P=WV. If V=0, (3) gives Rho=0, but P(-123/20)=-2727/100!=0.
Thus V!=0, and (5) gives f=0 in this chart too. In EVERY kappa!=0
candidate we may consequently use f=0.

For f=0, write the actual forced quartic as

    S_4=y^2(Z p^2+c K p q+256c^2 q^2)
                      +y^3(Bp+Gq)+H y^4,                         (7)
    Z=6Da+9D+64a^2-32ae-16a*rho+69a
                         +4e^2+2e*rho-(45/2)e-21rho-69/4,
    K=12D+256a-64e-32rho+138,
    G=(P-WV)/2+16Ac+256bc-32ci.

B,H need not be restricted below. Equations (7) follow from the complete
trace (2), including all its lower terms.

## 3. Repair-invariant second-jet witnesses on the COMPLETE trace

Let w7,w6,w5 denote respectively the coefficients y^7, y^5p, y^4q in
[x]D_(T+xS). By HC4MYGJ2, an arbitrary x^2U with deg U<=3 changes this
first transverse determinant by -192U. All three witnesses must vanish
for every degree-five graph. Their extraction uses the exact finite
Hessian-jet recipe of [PF154](FIRST_JET_ADMISSIBLE_NONTERMINATION_AND_SECOND_JET_QUOTIENT.md).

The retained derivation keeps all 18 parameters in (2), with f=0 justified
above, and forms S using its complete forced degree-at-most-four part.
On c=0 the three witnesses depend ONLY on

    (kappa,D,e,rho,a).                                            (8)

In particular A,i,b,h, every other lower coefficient and Cdet cancel
identically. This is the coverage step missing from the earlier sparse
trace exclusion HC4MYG5K. The exact full expressions, including the c!=0
terms, are retained in the [exact witness data](sources/second-jet-coverage/witnesses.json). The
checker rederives them from the actual parent Hessian jets using exact
polynomial multiplication. The cancellation is not inferred from samples.

## 4. All kernel directions on c=0 are excluded

First let the kernel be partial_q. Equation (7) gives P=WV.
We already have V!=0 and hence Rho!=0 from (3).
On W!=0, P!=0, substitute

    V=P/W, kappa=-P^2/(2Rho W^2),
    e=(32Va+40a*kappa-12Vrho-39V+60kappa+Q)/(8V).

Normalize the three witnesses as

    E7=w7 W^2 Rho/(24P),
    E6=-w6 W^2 Rho P/24, E5=-w5/48.                              (9)

They are polynomials in a,rho and generate the unit ideal. This is the
same three-equation obstruction as in the old slice, now verified for
the full trace by (8). The [retained exact certificate](sources/second-jet-coverage/q_kernel_unit_certificate.json) contains polynomial
multipliers u7,u6,u5 with

    u7 E7+u6 E6+u5 E5=1.                                         (10)

The checker verifies (10) by polynomial addition, rather than relying on
a stored assertion that a Groebner basis is one.

On the missing chart W=0, P=0, the unreduced witness w5 has remainder
540*(25rho+141) modulo P. The resultant of P and 25rho+141 is -15552,
so this chart is empty too. No denominator boundary is omitted.

Now let v_p!=0 and c=0. The p*y^2 coefficient of the kernel condition
on (7) gives Z=0. If V=0, then Rho=0, D=-47/10 and (4) gives
kappa=153/(100W), with W!=0. The actual witness is

    w7=66096/(125W)!=0.                                         (11)

If V!=0, then Rho!=0. Substitute

    kappa=-V^2/(2Rho),
    e=(V(32a-12rho-39)+20kappa*W+Q)/(8V).

Take polynomial numerators of w7/(24V), w6, w5 and Z, up to nonzero
rational normalizations retained in the input file.
Call them F1,F2,F3,F4. Their denominators are powers of the units V,Rho.
The [explicit polynomial certificate](sources/second-jet-coverage/V_membership_certificate.json) gives

    h1 F1+h2 F2+h3 F3+h4 F4=V.                                 (12)

Thus their vanishing contradicts V!=0. The certificate is an exact identity
in Q[a,V,rho], independent of parameter specialization; its 222KB of
rational multipliers are replayed directly. Its discovery used tracked
polynomial division, not a bounded point search. Equations (10)--(12)
exclude every constant-kernel direction on c=0.

## 5. The c!=0 chart is excluded over the algebraic closure

Here (5) gives kappa=-4V/5, so V!=0. Equation (3) then yields

    V=8Rho/5, kappa=-32Rho/25, Rho!=0,
    e=-(1120rho^2+15888rho+55503)/(64Rho),
    i=(16rho+99)A/(2V).                                         (13)

A constant kernel in (p,q) forces the determinant of the quadratic
coefficient matrix in (7) to vanish. Since c!=0, this is
1024Z-K^2=0. Under (13) its exact value is

    1024Z-K^2=-(24/25)(4000rho^2+49520rho+157969).                 (14)

The witness w7 simultaneously simplifies, for ALL other parameters, to

    w7=(192Rho/25)(160rho^2+1968rho+5841).                        (15)

The two quadratics in (14)--(15) have resultant 2544276275200.
An explicit rational Bezout identity is

    (91/1129384-5rho/1552903)(4000rho^2+49520rho+157969)
     +(125rho/1552903-24945/12423224)
                             (160rho^2+1968rho+5841)=1.           (16)

Thus the chart is empty over EVERY characteristic-zero extension field,
including complex rho. Its negative real discriminant alone is not used
as a complex exclusion. Together with Section 4 this proves kappa=0.

## 6. Deriving and retaining the two residual classes

Set kappa=0 in the necessary field-valued equations. Equation (3) gives
V=0; then (4) gives Q=0. The resultant

    resultant(Q,16rho+99)=-7776

and (6) force A=0. Solving Q=0 gives exactly the two rho values in (1),
and V=0 fixes D=(16rho+89)/2, giving D_s. This proves the coverage theorem
and the rational-coefficient degree-five exclusion.

The zero-kappa case need not have its top kernel in (p,q): its y-component
may be nonzero. Nothing above applies that dual-kernel restriction there.
Both conjugate residual classes are retained, including arbitrary f and
all degree-at-most-two coefficients. Their second-jet and later equations
still require analysis.

This result constrains FIELD-VALUED candidate graphs. PF153's uniform
cone-deviation J^4096 remains valid on actual trace quotients, but the
new field statement does not certify a power of kappa in the full
nonreduced ideal. Nor does it assert that PF154's A_2 ALONE is empty:
the full polynomial graph's necessary top-cone kernel was used.
No polynomial reduced body, increasing tail-power bound or termination
certificate on the remaining scheme is supplied. Preserve nilpotent
thickenings until an actual ideal-power certificate removes them.

The next concrete graph gate is the two zero-quintic-trace classes (1),
with their full normal-kernel constraints and actual multiplication.
Global [PF151 coverage](NORMAL_FORM_COVERAGE_AND_COMMON_LEADING_DIVISORS.md)
by graph/cyclic models and generic first-derivative integrality remain
open. Full finite normalization, codifferent and every ramified boundary
prime are retained, as are the 3D integrality controls and PF142's mixed
minimal retirement. No HC4/JC2 proof or actual counterexample is claimed.
The exact checker passes 61 controls. Independent whole-proof replay,
formal verification and external review remain unclaimed.
