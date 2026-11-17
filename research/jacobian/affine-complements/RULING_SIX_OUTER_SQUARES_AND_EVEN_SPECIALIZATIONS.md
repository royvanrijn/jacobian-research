# Odd companions force square outer faces and even finite specializations

PF2RULED168. 17 November 2026. Actual complex polynomial plane Keller maps.

At every off-M zero of the ruling-six leading coefficient, the OUTERMOST
Newton initial polynomial is a square. This uses the actual globally
reduced odd companion and the coupled leading coefficients, not only
the one-coordinate Newton polygon.

Consequently finite specialization degrees three and five are excluded.
Only degrees two and four remain. If all four escaping sheets have
the same outer slope, the quartic square-completion coefficient has
EXACTLY the same zero order as h. This includes every A pole divisor.

The [global divisibility theorem](RULING_SIX_GLOBAL_QUARTIC_DIVISIBILITY_AND_CUBIC_REMAINDER.md)
uses these conditions to reduce the regular square remainder to ruling three.
These are unrestricted necessary conditions on the remaining A class.
They do not prove that either lower-ruling square-completion piece is
a target function, or close JC2.

## 1. Actual polynomial class and every prime retained

Keep PF165's complete quadratic source-image class:

    w^2=N+M*z, M=u*A, N=u*Sigma^2,
    u=t^2, w=t*(Sigma+A*x), z=2Sigma*x+A*x^2,
    P=b3*z^3+w*(a2*z^2+a1*z+a0)+b2*z^2+b1*z+b0,
    b3=M*h^2, a2=2h*eta, eta in C[u].                    (1)

The coefficient eta is A2/2 in PF165. Set

    r4=b2-N*h^2-eta^2,
    T=h*w*z+eta*z,
    P=T^2+r4*z^2+w*(a1*z+a0)+b1*z+b0.                   (2)

All coefficient degrees are unrestricted. The actual fields are
F=C(P,Q), E=C(u,w), L=C(x,t), with full normalizations S_E and S
of B=C[P,Q]. They remain the ORIGINAL rings with all their primes.

Choose an actual globally reduced companion, which exists by the
retained elementary target reductions. PF165 proves

    s=deg_w Q*=3k, k>=3 ODD,
    p6=(h/M)^2, q_s=beta*(h/M)^k, beta in C*.              (3)

The resulting P,Q are still coordinates of the same actual target
ring B. The source volume is

    dP wedge dQ=du wedge dw/M.                            (4)

Let c be an arbitrary zero of h off M, v=u-c and j=ord_c h>0.
The polynomials P*,Q* in w have coefficients in C[[v]].
Write d=deg_w P*(c,w). PF52's whole-line criterion gives

    2<=d<=5.                                             (5)

On that whole line, w is affine-linear in an original source
coordinate. The same degree bound applies to every target coordinate
produced by a legitimate elementary target change.

## 2. The outer face is square

Give v weight -1 and w positive rational weight ell. For the
OUTERMOST escaping edge, define

    ell=max_(i<6,p_i!=0) (2j-ord_c p_i)/(6-i)>0,
    H=6ell-2j,
    in(P*)=v^-H I(v^ell*w), deg I=6.                     (6)

Because the specialization has degree d>=2 with unit coefficient,
H>=d*ell>=2ell. At least one lower term lies on this edge by
definition; I is not a monomial.

If K is the weighted degree of the actual companion, write

    in(Q*)=v^-K J(v^ell*w), deg J<=s.

Its top w term, from (3), has weight

    s*ell-kj=(k/2)*H,

so K>0 and K>=(k/2)H. The actual Jacobian has weight zero.
The possible leading Jacobian has weight H+K+1-ell>0, since
H>=2ell and k>=3. It must vanish. Direct differentiation gives

    H*I*J'-K*I'*J=0.                                    (7)

Comparing polynomial degrees in (7), or its logarithmic residues
at infinity, yields

    6K=H*deg J<=H*s=3kH.

Together with the opposite inequality above this forces

    K=(k/2)H, deg J=s.                                   (8)

At every root of I of multiplicity m, equation (7) gives
ord_root(J)=k*m/2. Since k is ODD and this is an integer,
every root multiplicity m is even. Over C,

    I=C*U(Y)^2, C in C*, deg U=3.                        (9)

This proof allows rational slopes and every coefficient cancellation.
It does not assume a simple root or bound a companion degree.

The same algebra applies after a p-independent Laurent/Puiseux
translation ONLY if the full degree-six coefficient is still on
the chosen face and the leading Jacobian weight remains positive.
No claim propagates (9) to a refined degree-two face automatically.

## 3. Escape mass rules out degrees five and three

The specialization P*(c,w)-p has exactly d finite generic roots.
This can also be seen by the all-order monic factor recursion
over C[p][[v]] retained in PF162. The remaining 6-d roots escape,
counting each normalized u-index and each residue conjugate.

For the outer initial I, the number of outer escaping roots,
with that multiplicity, is

    6-ord_(Y=0) I.

It is positive, since the edge contains another degree, and is
EVEN by (9). It cannot exceed 6-d.

For d=5 only one root escapes. This contradicts the positive
even outer mass immediately.

For d=3 there are three escaping roots. The outer face consumes
exactly two, leaving exactly ONE root at a smaller pole slope.
The separation by pole slope is independent of p, so that remaining
root has normalized u-index one and residue field C(p). It is
not a third conjugate of the outer double cluster.

PF162's degree-half specialization estimate makes EVERY coframe
end holomorphic when d>=3. The remaining single root therefore
has a finite actual Q residue.

Here is the exact graph obstruction. On that root, let lambda
be the first Laurent coefficient outside C after its finite
C-coefficient truncation. The PF159 leading-coefficient lemma
gives H0(lambda)=p for H0 in C[Y]. Since the whole Laurent
branch has residue C(p) and u-index one, lambda lies in C(p).
The degree of the polynomial map H0 is then one. Thus lambda
is affine in p.

The same Taylor-leading argument gives the finite companion
residue as a polynomial G(lambda), or zero. It is therefore
q0(p) in C[p]. The original t lifts over c are units, while
the source x is polar. The full-normalization closed-image
argument of PF165 gives the WHOLE graph

    Gamma={Q=q0(P)}

as an actual nonproper-value component. The retained published
Chau component Remark excludes it.

Therefore

    d in {2,4}                                           (10)

at EVERY off-M root of h. No root multiplicity, position or
number restriction has been imposed.

This proof retains the single leftover root; discarding it
after examining only the highest double cluster would conceal
the contradiction.

## 4. Exact r4 order when all four roots lie on one outer face

Suppose the outer face consumes four roots. Then d=2, and its
zero root has multiplicity two. Equation (9) has the form

    I=(aY^3+bY^2+c1Y)^2, a*c1!=0.                        (11)

Choose its square-root sign so that a is the leading coefficient
of h/M.

Since the degree-two coefficient of P* is a unit, H=2ell.
The degree-six coefficient gives 2j=4ell, so

    j=2ell.

The degree-five and degree-four coefficients in (11) are
2ab and b^2+2a*c1. From (1),(2) the actual coefficients are

    p5=2h*eta/M^2,
    p4=(r4+eta^2-2N*h^2)/M^2.                            (12)

It follows that eta/M has order at least ell, and its coefficient
at ell is b (zero is allowed if that exponent is absent).
The coefficient of p4 at order 2ell=j is b^2+2a*c1.
The contribution 2N*h^2/M^2 has order at least 2j>j.
Therefore

    ord_c r4=j=ord_c h,
    leading(r4/M^2)=2a*c1!=0.                            (13)

This remains true when b^2+2a*c1=0 and the displayed p4 face
coefficient itself vanishes. The difference in (12), not a
non-cancellation assumption on p4, proves (13).

Every quartic A pole face is a case of (11), with
U(Y)=aY(Y-c0)^2 and ell=L. Thus (13) holds at EVERY actual
pole divisor, with exact order 2L.

At any off-M h-root where ord r4<ord h, (13) consequently forces
the outer escaping mass to be exactly TWO, not four. Degree
three is excluded by Section3, so the specialization there has
degree two or four. This pinpoints the remaining two-cluster
analysis needed for global divisibility; it does not assume it.

## 5. Scope and controls

The graph step imports exactly the Chau component application
already retained and checked in PF161/PF165. Polynomial leading
coefficients and the complete specialization estimate are the
retained PF159/PF162 inputs. No new external classification is
imported.

The checker verifies (7), the exact coupled weight equality,
the multiplicity parity, the COMPLETE escape-mass possibilities
for degrees two through five, and the full coefficient identity
(12),(13), including a cancelled degree-four coefficient and
half-integral outer slopes. The checker passes 75 exact controls. The uniform proof is the
actual weighted Jacobian and the complete graph argument, not
bounded coefficient testing.

Replay with

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_outer_square.py

Only degrees three and five are closed here. The A class, both
regular lower-ruling pieces outside F, full derivative integrality
and normal-tail termination remain unresolved. Independent
whole-proof review, formal verification and external human
review remain unclaimed. HC4/JC2 is not proved.
