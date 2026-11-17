# Regular Bezout maps, polynomial slices and weighted Haar tails

15 November 2026. Internal written proofs with exact rational controls.
No novelty, independent specialist replay, formal verification or external
review is claimed. See the [checker](verify_exceptional_bezout_maps.py),
[frozen regular maps](exceptional_bezout_maps.json),
[receipt](exceptional_bezout_verification.json), and
[source audit](exceptional_bezout_source_audit.json).

The [quadratic model](QUADRATIC_TRACE_MODELS_AND_HAAR_RESIDUES.md) now has
an explicit regular map in the opposite direction. This yields a global
polynomial-slice description of the entire admissible algebra and an
exact weighted pullback of every ordinary SL_2 pure and marked identity.
It supplies a family with a fixed, certified infinite marked tail once
unweighted purity can be established. The canonical member fails purity
already in its first moment on both exceptional quotients. Other members
are not decided. Neither exceptional MZ question is closed.

## 1. A canonical regular completion by a binary Bezout identity

Use the forms f(z,w),g(z,w), degrees n,m and c=Jac(f,g)/(nm) from the
quadratic model. Thus (n,m,r)=(8,12,18) or(12,20,30), with r=n+m-2.
Put N=n+m-1=r+1. The forms f,g have no common projective zero.

The following map of binary-form spaces is an isomorphism:

    (A,B) -> f A+g B,
    Sym^(m-1)(C^2) direct_sum Sym^(n-1)(C^2)
                                  -> Sym^N(C^2).             (1)

It is injective by coprimality and the strict degree bounds, and the
dimensions agree. Apply its inverse, coefficientwise in x,y, to
(xw-yz)^N. There are unique polynomials A,B with

    f A+g B=(xw-yz)^N,
    deg_top A=deg_top B=N,
    deg_bottom A=m-1, deg_bottom B=n-1.                       (2)

Because f,g and the determinant are right Gamma invariant, uniqueness
proves that A,B are invariant as well. Define

    U=B, V=-A, F=f(z,w), G=g(z,w), C=c(z,w).                  (3)

On the source SL_2/Gamma, writing its ring as A_Gamma, we have

    U G-V F=1.                                               (4)

This gives a regular map Phi_0 to SL_2, and a regular map sigma_0 to
the quadratic ambient Z by also sending eta to C. Both are over the
common punctured base U_Gamma. No poles or unproved regularization
occur. The frozen artifact contains all32/36 terms of U,V for2O and
all56/118 terms for2I. The exact Sylvester systems have sizes20 and32;
this is one specified construction, not a degree or support sweep.

All regular maps to Z over this fixed base are precisely

    sigma_b^*X=U+F b, sigma_b^*Y=V+G b,
    sigma_b^*(F,G,eta)=(F,G,C), b in A_Gamma.                 (5)

Indeed the difference of two completions in (4) is (Fb,Gb), since
F,G generate the unit ideal by (4).

## 2. Exact Jacobians and the sharp minimum fibre degree

Let D=z partial_x+w partial_y on the source. Differentiating (2)
gives (DU)g=(DV)f in the polynomial ring before determinant1 is
imposed. Coprimality and bidegrees therefore give a binary form H
depending only on x,y, homogeneous of degree r, such that

    DU=F H, DV=G H.                                         (6)

It is Gamma invariant. In both models the degree-r invariant space
is one-dimensional, spanned by c. Thus H=kappa*c(x,y). The exact
construction gives

    kappa_O=3591/2, kappa_I=1117116/3125.                     (7)

For the general map (5), put

    H_b=H+D b, W_b=C H_b.                                   (8)

Then D(U+Fb)=F H_b and D(V+Gb)=G H_b. With the invariant source
volume Omega_Gamma and the standard target volume Omega_*,

    Phi_b^*Omega_* = nm W_b Omega_Gamma,
    sigma_b^*Omega_Z = nm H_b Omega_Gamma,                    (9)

where Omega_Z=p^*Omega_*/eta is the regular nowhere vanishing form
from the [branch theorem](BRANCH_CLEARING_AND_HAAR_AUTOMORPHISM_RIGIDITY.md).
For example on F!=0, differentiate the target slice (U+Fb)/F and
use dF wedge dG=nm C dz wedge dw. This proves (9) on a dense chart;
all forms in (9) are regular, so the identities hold globally.

The canonical U,V have D-depth exactly N: by (6),

    D^N U=F*kappa*r!*C, D^N V=G*kappa*r!*C.                  (10)

In fact **every regular map (5) has maximum coordinate depth at
least N**, so the canonical map attains a sharp global minimum.
If b has depth>N, the added terms cannot cancel their top depth
against U,V. If b has depth<=N and reduced both coordinates to
depth<N, (10) would force

    D^N b=-kappa*r!*C.

But the [all-order higher-image theorem](QUADRATIC_TRACE_MODELS_AND_HAAR_RESIDUES.md)
gives D^N(A_Gamma) intersect ker D equal to the bottom invariant
forms of degree>=N. The nonzero right side has degree r=N-1,
a contradiction. Thus the lower bound19/31 is independent of the
degree of the freely chosen b, not evidence from a bounded search.

## 3. The full admissible algebra has a global polynomial slice

Retain R,delta,eta and the exact image S of A_Gamma in R. Add one
variable t and the locally nilpotent derivation

    Delta=delta+eta partial_t on R[t].                       (11)

For a in A_Gamma, define its invariant extension

    i(a)=sum_(j>=0) (-t)^j/(j!*eta^j) *delta^j(Psi^*a).      (12)

This is a finite polynomial in R[t], by exact admissibility. It is
killed by Delta, and evaluation t=0 recovers Psi^*a. Conversely,
the coefficient recurrence for any Delta-invariant polynomial forces
exactly (12). Hence

    i(A_Gamma)=ker Delta.                                   (13)

There is now an explicit **global slice**

    s=i(U)*Y-i(V)*X, Delta s=i(U)G-i(V)F=1.                  (14)

The elementary finite Taylor expansion for a locally nilpotent
derivation with a slice proves

    R[t]=(ker Delta)[s] isomorphic to A_Gamma[s].             (15)

In particular the original quotient and the quadratic ambient variety
have isomorphic cylinders. This alone supplies no product-functional
transfer. The canonical period of Omega_Z has T(1)=0, whereas the
normalized original Haar has h(1)=1; L=T(eta*-) is a different weighted
period. No compatibility of auxiliary-variable functionals is asserted.

The retraction rho=exp(-s Delta) has

    rho X=i(U), rho Y=i(V), rho t=i(nu),
    nu=V*a_f-U*a_g,
    a_f=(x f_z+y f_w)/n, a_g=(x g_z+y g_w)/m.                 (16)

To check the last formula, use

    i(a_f)=eta X-Ft, i(a_g)=eta Y-Gt,
    i(nu)=t-eta s.

Since rho is onto ker Delta, (16) proves the explicit six-generator
statement

    A_Gamma=C[F,G,C,U,V,nu].                                 (17)

These generators also have a useful differential relation. Directly
from (4),(6) and a_f G-a_g F=C,

    D nu=C H-1, D(nu+C b)=W_b-1.                            (17a)

Thus the Jacobian weight differs from1 by an explicit unipotent
derivative. It cannot be dropped inside a nonconstant phase moment:
integration by parts differentiates the phase as well.

It also gives a complete-intersection presentation in six variables:

    A_Gamma isomorphic to
    C[X,Y,F,G,eta,t]/(XG-YF-1, eta^2-K, s).                  (18)

Here the polynomial s is specified by the finite formula (12),(14),
not expanded in the retained arithmetic receipt. The first two
equations define R[t], a domain; s is a polynomial coordinate in
(15). They are a regular sequence, proving (18). The universal slice
and presentation proof is written; the small checker does not certify
every expansion implicit in (12). No cancellation conjecture is
assumed in this argument.

This provides a global finite construction for all exact admissibility
cancellations, rather than replacing (1) of the quadratic model by
an unconstrained ambient expression. It does not classify which of
these expressions have vanishing pure moments.

## 4. Exact weighted transfer for every target pullback marker

Let T(a_0+eta a_1)=h_*(a_1) and L(a_0+eta a_1)=h_*(a_0).
For every b in A_Gamma and every a in R,

    h_Gamma(H_b sigma_b^*a)=T(a),
    h_Gamma(W_b sigma_b^*a)=L(a).                            (19)

Proof. Both sigma_b and Psi are over the identity of U_Gamma;
the additive torsor projections are homotopy equivalences. Thus
Psi sigma_b induces the identity on H^3. Pulling back
Psi^*Omega_Gamma=eta Omega_Z/(nm) through sigma_b and using (9)
gives

    h_Gamma(W_b)=1.                                         (20)

The normalized period T has T(eta)=1. Pull its defining top-period
classes back through sigma_b and apply (9),(20), obtaining the
first identity in (19). Apply it to eta*a for the second. This
argument retains the branch locus; it is not an etale trace argument.

In particular, for arbitrary target regular polynomials P,Q and every
j>=0,

    h_Gamma(W_b (Q o Phi_b)(P o Phi_b)^j)=h_*(Q P^j).        (21)

All source functions in this identity are regular and invariant.
It is valid for all b and all orders, not just the canonical map.
The markers in (21) are Q o Phi_b; those in (19) are sigma_b^*a.
These identities alone do not cover every original marker of A_Gamma.
The [Haar weight-invariance theorem](HAAR_WEIGHT_INVARIANCE_AND_ORIGINAL_MARKERS.md)
removes a nonzero weight from qualitative safety only AFTER all original
weighted marked tails have been controlled. It does not remove W_b
from individual identities or transfer all-positive purity.
Also H_b is never identically zero: otherwise (20) would fail.
Every Phi_b is dominant, but its Jacobian always vanishes on C=0.
No choice of b can make W_b the constant1, or turn these fixed-base
maps into etale maps. Dropping W_b in (21) has no justification.

For the canonical map, (7),(20) give the exact normalization

    h_Gamma(c(x,y)c(z,w))=2/3591 or3125/1117116.              (22)

These norms and a nonconstant second-moment transfer are independently
computed algebraically by the checker within this same internal run.
They are not an independent specialist replay.

## 5. A family with the entire marked tail already supplied

Use the pinned Long Hopf pair on the target SL_2:

    u=-2XY, v=2FG, tau=XG+YF, uv+tau^2=1,
    P_H=(1+u)(v-(2+u)tau^2), Q_H=u.

Its retained all-order input is

    h_*(P_H^j)=0,
    h_*(Q_H P_H^j)=c_j=4^j(j!)^2/(2j+1)!>0, j>=1.

Define regular source functions

    P_b=P_H o Phi_b, M_b=W_b (Q_H o Phi_b).                  (23)

Then for **every** b, the exact identities (21) supply

    h_Gamma(W_b P_b^j)=0,
    h_Gamma(M_b P_b^j)=c_j>0 for every j>=1.                 (24)

In particular **every P_b is outside the original nullcone**, since
the fixed multiplier M_b has a nonzero infinite tail. This family
already avoids the universally cleared unstable phases. The remaining
purity condition is substantive, not a missing choice of multiplier.

Therefore any b for which the **unweighted** pure moments
h_Gamma(P_b^j) all vanish would immediately give an actual exceptional
MZ failure, with its fixed regular marker already specified in (23).
This is a sufficient construction family, not an equivalence between
existence of such a b and failure in the entire coordinate algebra.

The canonical b=0 does not solve that remaining condition. Its torus
homogeneity makes the weight-zero part of P_0 equal to1-3tau_0^2.
The exact calculations give

    h_2O(tau_0^2)=40224721/118915225,
    h_2O(P_0)=-1758938/118915225 !=0,

    h_2I(tau_0^2)=694662221/2181764186,
    h_2I(P_0)=97777523/2181764186 !=0.                        (25)

In both cases h_Gamma(W_0 tau_0^2)=1/3 exactly, as required by
(21); the ordinary first moment in (25) exhibits why that weight
cannot silently be omitted. No search over b, support or degree was
run. The all-degree new target is to force the unweighted pure
identity in (23), or prove an obstruction for this entire family,
while retaining the broader pure-null classification problem in S.
