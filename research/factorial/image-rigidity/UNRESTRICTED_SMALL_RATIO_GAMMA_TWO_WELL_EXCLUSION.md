# The small-ratio sector without a lower bound on a/d

Separate scratch proof variant for whole adversarial review. The already
reviewed restricted proof SMALL_SUBCRITICAL_RATIO_GAMMA_TWO_WELL_EXCLUSION.md
remains unchanged. No canonical source or status is changed here.

For integers a>=1,d>=2 write D=a+d and

    A_m(x)=sum_(j=0)^m binom(m,j)*(am+1)_(dj)*x^(m-j),
    c=a/(d-1)^2,       J=d^2/D.

**Theorem UnrestrictedSmallRatioUniformEventualExclusion.** There exist
unevaluated epsilon_0>0 and N_0 such that for every integer a>=1,d>=2
with c<=epsilon_0 and every N>=N_0, the ACTUAL pair A_N,A_(N-1) is
coprime over Q. The same cutoff works for every nonzero coefficient pair and
every characteristic-zero field.

There is NO restriction d<=3a. Unlike the restricted proof, this variant
uses the retained OLD sparse theorem as an indispensable logical premise.
It does not assert an unconditional fixed-annulus bound for every individual
root at every order. The root capture used here is conditional on the order
constraint forced by a hypothetical consecutive common root.

## 1. The old theorem and the order constraint at a bad pair

The actual retained statement, lower_exponent_initial_and_boundary.tex,
thm:sparse, is:

    a>=1,d>=2,N>=2 and d^2>=4(a+d)log(8N)
    imply that the first N factorial moment polynomials are individually
    squarefree and pairwise coprime.                              (U1)

Therefore a hypothetical bad terminal N necessarily has

    J<4log(8N).                                                  (U2)

We will exclude every sequence of bad pairs with J->infinity,N->infinity.
There is no requirement that a tend to infinity. We use (U2) to control
the ratio of the two Gaussian widths, not to replace actual critical
actions by model values.

Put R=D/a. Always R>1. If R>=2 then d>=D/2, whence J>=D/4 and R<=D<=4J,
using a>=1. Consequently, whenever J>=1,

    R<=4J,      log R=O(log J).                                 (U3)

At m=N or N-1 in a sequence satisfying (U2), this implies

    log R=o(m),     sqrt(R*m)*exp(-delta*m)=O(1/m)                (U4)

for EVERY fixed delta>0. The implicit eventual bound in (U4) is uniform:
J<=4log(8N), N-1<=m<=N, and the exponential in m dominates the displayed
fixed polynomial/logarithmic factors.

## 2. Exact two-Gamma normalization and the identric seam

For ell>=1 define on S>0

    f_ell(S)=(S/ell)^ell*exp(ell-S).

For the actual integer exponents these are entire functions of S. Set

    W=exp(D log D-a log a-d),
    R_0=D^D exp(-D)=W*a^a exp(-a),
    S_0=W^(1/d).

Exactly, without any asymptotic replacement,

    S^a(Wq+S^d)exp(-S)/R_0=q f_a(S)+f_D(S),
    I_m(Wq)=R_0^m J_m(q),
    J_m(q)=integral_0^infinity(q f_a+f_D)^m dS,
    I_m(w)=(am)!*A_m(m^d w)/m^(Dm+1).                            (U5)

The integer-order integrand is entire and all displayed scalar factors
are nonzero. Our finite contour deformations use only integer orders.
No complex-order identity across a phase zero is needed.

The seam is the identric mean:

    log S_0=(1/d)integral_a^D log u du,
    sqrt(aD)<=S_0<=(a+D)/2,   S_0>=D/e,
    f_D/f_a=(S/S_0)^d.                                         (U6)

The first two inequalities follow from the chord and midpoint bounds for
the concave function log. For the third, the average of the increasing
function log over [a,D] is at least its average over [0,D], which equals
log D-1. In particular the left and right halves have their appropriate
single-well dominance: f_D<=f_a on [0,S_0] and f_a<=f_D on [S_0,infinity).

Both seam heights have a uniform exponential loss, even when R is unbounded:

    f_D(S_0)<=exp(-J/8),    f_a(S_0)<=exp(-J/288).                (U7)

Indeed D-S_0>=d/2, and the Gamma rate-function integral gives
-log f_D(S_0)>= (D-S_0)^2/(2D)>=J/8. If R<=4, then
S_0-a>=d/3 and -log f_a(S_0)>=J/18. If R>=4, use e<3,
S_0>=D/e>D/3 and a<=D/4 to get S_0-a>D/12. Thus
-log f_a(S_0)>D/288>=J/288. No delicate best constant is required.

The center-to-center cross losses are

    K_L=log(W/a^d)=D log(D/a)-d,
    K_R=log(D^d/W)=d-a log(D/a),
    K_L>=J/2, K_R>=J/2,
    K_L/d=(R/(R-1))log R-1>=log R-1.                            (U8)

The first two lower bounds follow by integrating (D-u)/u and (u-a)/u
over [a,D], and replacing the denominator by D.

## 3. Complete real Gamma tails, including a bounded left shape

The elementary inequalities

    log(1+x)-x<=-x^2/2                 (-1<x<=0),
    log(1+x)-x<=-x^2/[2(1+x)]         (x>=0)

give

    integral_0^infinity f_ell(S)dS<=6sqrt(ell),                  (U9)

for every ell>=1. For example on 0<S<=2ell use v=(S-ell)/sqrt(ell)
and bound by exp(-v^2/4); for S>=2ell the logarithmic derivative is at
most -1/2, giving a remaining tail at most2.

Fix one v_0 in (0,1/4). Monotonicity on either side of ell and the same
inequalities give a uniform strict loss, for every ell>=1:

    sup_(|S-ell|>=v_0 sqrt(ell)) f_ell(S)
       <=k_0:=exp[-v_0^2/(2(1+v_0))]<1.                       (U10)

This statement does not require ell>=R^2 for a large local radius.
The endpoints ell+-v_0 sqrt(ell) lie on the positive real axis.

A useful stronger LEFT-half L1 bound avoids even the width ratio.
On [0,S_0], the logarithmic derivative of f_D is at least
D/S_0-1>0. Hence

    integral_0^S_0 f_D(S)dS
       <=S_0/(D-S_0)*f_D(S_0)
       <=(2D/d)*exp(-J/8)
       <=4sqrt(a)*exp(-J/8)       (J>=1),                      (U11)

since D/(d sqrt(a))=1/sqrt(a)+sqrt(a)/d<=1+1/sqrt(J)<=2.
This estimate controls a long interval even when a is bounded or D/a
is unbounded; the interval length is not silently treated as uniform.

## 4. Uniform analytic Morse families without assuming a->infinity

For b in [0,1] define the removable analytic family

    H_b(v)=[Log(1+bv)-bv]/b^2,
    H_0(v)=-v^2/2.                                             (U12)

It is jointly holomorphic for |v|<1/2 and b in a small complex
neighborhood of [0,1]. It satisfies H_b'(0)=0,H_b''(0)=-1 with
all required higher derivatives uniformly bounded. The natural Gamma
logarithm at ell is H_(1/sqrt(ell)), because S=ell+sqrt(ell)v.
Thus bounded integer a is part of a compact nondegenerate family,
not a spurious Gaussian-limit assumption.

On every fixed local disc |v|<=v_1<1/2, the cross ratios tend uniformly
to zero as J->infinity:

    |f_D/f_a|(a+sqrt(a)v)
       <=exp[-K_L+2v_1*d/sqrt(a)] ->0,
    |f_a/f_D|(D+sqrt(D)v)
       <=exp[-K_R+2v_1*sqrt(J)] ->0.                           (U13)

For the second inequality use (U8). For the first, if R stays bounded,
d/sqrt(a)=sqrt(JR) and (U8) dominates its O(sqrt(J)) term. If R tends
to infinity, (U8) gives K_L/d>=log R-1, whereas
2v_1/sqrt(a)<=2v_1; this again tends to infinity after subtraction.
Every sequence with J->infinity has a subsequence in one of these cases,
which proves the asserted uniform convergence. Derivative convergence
follows from Cauchy estimates on a slightly larger fixed disc.

For q in a fixed annulus away from0, the exact local phases are

    log q+H_(1/sqrt(a))(v)+Log(1+q^(-1)f_D/f_a),
    H_(1/sqrt(D))(v)+Log(1+q f_a/f_D).                           (U14)

The uniform analytic implicit-function and Morse constructions give
critical points S_L=a+sqrt(a)v_L, S_R=D+sqrt(D)v_R, with v_i->0,
and actual normalized critical values lambda_i such that

    lambda_L/q->1, lambda_R->1,
    Re(S_R-S_L)/d->1, |S_R-S_L|<=K d.                           (U15)

The last statement uses sqrt(a)/d<=1/sqrt(J) and
sqrt(D)/d=1/sqrt(J). The normalized Morse Jacobians c_L,c_R tend to1,
because the critical curvatures tend uniformly to -1. Their q-derivatives
are uniformly bounded on a slightly smaller fixed annulus. No global
complex log q is needed: choose finitely many local q-charts when required.

The same local assertions hold for the right well when |q|<=1/2,
including q=0, and for the left well in the reciprocal variable z=1/q
when |z|<=1/2. The corresponding perturbing cross ratio is multiplied
by a bounded parameter, with no division by that parameter.

The natural local intervals lie in the appropriate seam halves: from
S_0>=sqrt(aD) one obtains

    (S_0-a)/sqrt(a)>=sqrt(D)-sqrt(a)
       =d/(sqrt(D)+sqrt(a))>=sqrt(J)/2,
    (D-S_0)/sqrt(D)>=sqrt(J)/2.                             (U16)

Thus the fixed v_0 intervals are disjoint and strictly inside their
respective halves once J is large.

## 5. Conditional capture of ALL common-pair roots into a fixed annulus

For |q|<=1/2 the right local Morse packet has actual lambda_R->1.
Outside its fixed local interval, split the real ray at S_0.
On its left half the modulus of q f_a+f_D is at most1/2+exp(-J/8).
On its right half but outside the right interval it is at most
k_0+(1/2)exp(-J/288). Thus the whole outside contour is uniformly
below the selected critical modulus by a fixed margin. Its full L1 norm
is at most9sqrt(D), by (U9). The short local Morse connectors retain
the same strict loss by the compact local family (U12)--(U14).
The standard symmetric local Gaussian calculation consequently gives a
nonzero leading term and a uniform relative O(1/m) error as J,m->infinity.
So J_m has no roots with |q|<=1/2.

For |q|>=2 divide by q and put z=1/q. The factor is f_a+z f_D.
Its left local Morse packet has critical value tending to1. Outside
that interval, the left half has modulus at most
k_0+(1/2)exp(-J/8), and the right half at most1/2+exp(-J/288).
Again there is a uniform strict height gap. The full L1 norm is now

    <=6sqrt(a)+3sqrt(D).

After division by the leading sqrt(a/m) Gaussian, its nonlocal relative
bound is C sqrt(m)*(1+sqrt(R))*exp(-delta m). The local Gaussian error
is uniformly O(1/m) by the compact Morse family. For m=N,N-1 of a
HYPOTHETICAL BAD PAIR, (U2)--(U4) absorb the width factor and make the
whole nonlocal error O(1/m). The actual critical powers are retained.

It follows that both actual common-pair arguments satisfy

    1/2<|q_N|,|q_n|<2,    q_m=x/(m^d W), n=N-1,              (U17)

along any bad sequence with J,N->infinity. This is conditional capture
at the two orders constrained by (U2). We do NOT claim (U17) for every
individual root at arbitrary small m and unbounded D/a.

## 6. Two complete exact partial periods in the annulus

Split the actual positive ray at the support-fixed S_0:

    P_L=integral_0^S_0(q f_a+f_D)^m dS,
    P_R=integral_S_0^infinity(q f_a+f_D)^m dS,
    J_m=P_L+P_R.                                               (U18)

On the annulus 1/3<=|q|<=3, outside the fixed LEFT interval the modulus
is at most |q|k_0+exp(-J/8), compared with the actual left peak |q|(1+o(1)).
Outside the RIGHT interval it is at most k_0+3exp(-J/288), compared
with the right peak1+o(1). Both are fixed strict relative gaps.
The complete L1 bounds are, by (U9)--(U11),

    integral_0^S_0(|q| f_a+f_D)dS<=22sqrt(a),
    integral_S_0^infinity(|q| f_a+f_D)dS<=24sqrt(D).             (U19)

Finite local Morse deformations and the entire integer integrand now give
TWO ACTUAL complete Gaussian packets, uniformly as J,m->infinity:

    P_L=sqrt(a)*lambda_L^m*sqrt(2pi/m)*c_L*(1+R_L),
    P_R=sqrt(D)*lambda_R^m*sqrt(2pi/m)*c_R*(1+R_R),
    |R_i|<=K/m.                                                (U20)

Here the packet-relative estimates themselves need no restriction on R:
the left cross L1 was bounded in its OWN width in (U11), and the right
width is at least the left one. In a symmetric Morse interval the odd
first Jacobian term integrates to zero; the quadratic Taylor remainder
is O(1/m) relative to the Gaussian. The short local connectors and every
whole remaining real piece have the just-proved fixed exponential loss.
The Gaussian Jacobians and their q-derivatives are uniformly controlled
by (U12)--(U15). No derivative of an unspecified remainder will be used.

Set locally

    D_*=Log(lambda_R/lambda_L),
    b=Log(sqrt(R)*c_R/c_L).

Then

    D_*=-Log q+o(1),
    |b|<=K+(1/2)log R,       |b_q|<=K.                          (U21)

The factor sqrt(R) is support-fixed, so it does not enlarge b_q. At any
actual zero in the inner annulus, (U18)--(U20) give

    Re[mD_*+b]=O(1/m).                                         (U22)

Only endpoint values of these uniform relative errors will be used.

## 7. Endpoint capture precedes interpolation; actual actions stay exact

By (U17), both own-order common-root arguments are in the fixed annulus.
Equations (U21)--(U22) and log R=o(N) give

    |q_N|->1, |q_n|->1,
    Lambda=(N/n)^d=|q_n|/|q_N|->1,       d/N->0.                (U23)

The use of (U2) here is essential for this proof variant. After a subsequence
both endpoints and their joining radial path lie in one fixed q-chart.
For m=N-theta, q_theta=(N/m)^d q_N, define

    E(theta)=mD_*(q_theta)+b(q_theta).

Stationarity of Phi=aLog S+Log(w+S^d)-S gives exactly

    w Phi_w=(D-S)/d,
    d q(D_*)_q=-(S_R-S_L)=-DeltaS,
    E'(theta)=-D_*-DeltaS+(d q_theta/m)b_q.                     (U24)

All leading actions are ACTUAL. Noninteger m here interpolates only
these critical actions and leading amplitudes; no period identity is
asserted at a noninteger order.

The upper endpoint zero gives

    |Re D_*(q_N)|<=K(1+log R)/N.

Along the path the second identity in (U24) adds at most Kd/N. Since
log R=log(1+d/a)<=d/a<=d and a>=1, (U15),(U21),(U24) imply

    Re E'(theta)<=-d/2                                      (U25)

eventually, uniformly throughout the adjacent path. But both integer
endpoint values have real part O(1/N) by (U22), a contradiction.
This excludes every bad sequence with J,N->infinity.

At no point were lambda_L^m or lambda_R^m replaced by their model powers.
It is allowed that m exp(-K_L) or m exp(-K_R) tends to infinity.

## 8. Uniform small-sector quantifiers and fields

If c=a/(d-1)^2 tends to zero, a>=1 implies d>=1+1/sqrt(c), and

    J=d^2/(a+d)>=1/(c+1/d)>=1/(c+sqrt(c))->infinity.             (U26)

If no epsilon_0,N_0 existed, choose counterexamples with c<=1/j and
N>=j. They satisfy J,N->infinity, and the old sparse premise forces
(U2), contradicting (U25). This proves the stated uniform eventual
small-ratio theorem WITHOUT a lower bound on a/d or a finite lower-one
certificate.

The rational pair's absence of a common complex zero is equivalent to
coprimality over Q. A rational Bezout identity persists over every
characteristic-zero field. Exactly

    Lambda((X^a(A+B X^d))^m)=B^m(am)! A_m(A/B)

when B!=0; a zero coefficient leaves a monomial with nonzero factorial
moments. The cutoff is therefore independent of the coefficients and
the characteristic-zero field. No claim about a fixed nonclosed field
is substituted for complex root absence.

## 9. Scope and dependence compared with the restricted proof

The old sparse theorem (U1) is a genuine indispensable premise of THIS
argument's conditional root capture and width control. The elementary
Gamma L1 bounds, compact family (U12), seam split and stationary envelope
are all proved here. No 392 relative-gap theorem, 407 theorem, 411 theorem,
sharp 404 lower annulus, root-density limit, or old a=0/a=1 finite-prefix
computation is a premise. This is a classical Laplace method applied to
an exact decomposition, with no historical-priority claim.

The reviewed restricted proof remains a separate, stronger-in-dependence
route when d<=3a: it requires neither (U1) nor an order/width constraint.
This extension removes that parameter restriction by importing (U1),
while explicitly refusing an unrestricted individual-root annulus claim.
Both statements are eventual and non-effective. This note supplies no
all-order Strong Factorial result or global finite bad-triple assertion.
