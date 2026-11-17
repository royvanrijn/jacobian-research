# Weighted torus covers and all zero-weight triples

16 November 2026. Internal written proofs with retained named inputs.
Independent replay, formal verification and external review are not claimed.
Both COMPLETE exceptional MZ cases remain **OPEN**; review is deferred.

The symmetric restriction is removed: every original phase with at most
three torus weights, one of them zero, satisfies MZ on BOTH exceptional
quotients, in every degree and over every commutative C-algebra. For
opposite charged weights the exact eventual-purity criterion is nilpotence
of the neutral part and the primitive balanced product. All polynomial
recompositions are classified over C. The proof uses a weighted torus
curve whose entire unweighted moment functional is preserved.

This is a milestone in the FULL exceptional programme. Three-weight
supports WITHOUT zero, larger supports and the full corrected-component
implication remain unresolved. Section7 states the coverage/preservation
lemma still missing; it is not replaced by a claim of full closure.

See the [checker](verify_haar_weighted_conic.py),
[receipt](haar_weighted_conic_verification.json),
[source audit](haar_weighted_conic_source_audit.json) and
[ownership amendment](exceptional_mz_ownership_20261116.json).

## 1. The weighted curve realizes EVERY original polynomial moment

Write P=p+U+V with weights 0,a,-b, a,b>0. Put

    g=gcd(a,b), r=a/g, s=b/g, M=r+s, gcd(r,s)=1,
    q=U^s V^r.

Weight projection gives, for every m,

    h(P^m)=sum_(M*l<=m) m!/[(m-M*l)!(s*l)!(r*l)!]
                                      *h(p^(m-M*l)q^l).        (1)

On the original exceptional algebra the joint neutral functional on p,q
is exactly the normalized compact sphere period in C[Q2]^O or C[Q2]^I.
For q!=0 use the normalization of the weighted torus threefold

    W: A^s B^r=q(x), S=p(x)+A+B,
    Theta=Omega wedge dA/(r A)
          =Omega wedge dS/(r A-s B).                   (2)

The second equality follows by differentiating the relation along the
fibre; all base differentials wedge to zero with Omega.

Here is an actual compact-cycle construction, including the root issue.
Normalize the finite cover t^r=q of the affine sphere, keeping ALL
components with their total generic degree r. Perturb its compact sphere
cycle M0 within its homology class to be transverse to the smooth reduced
q=0 divisor and avoid that divisor's finitely many singular/intersection
points. The full inverse image M' is a compact oriented branched cycle
of degree r. It lies in the smooth locus of the normalized cover: near
a smooth divisor with multiplicity k, normalization splits the cyclic
cover into smooth ramified charts. A resolution is an isomorphism there.

Over M' times |z|=1 use A=z^r, B=t z^-s. This maps into the normalized
W; locally a finite free action on the nonzero z coordinate describes
the map, so its image is in the smooth locus along the chosen cycle.
The actual pullback of Theta is Omega wedge dz/z. Divide its pushed
cycle by r and normalize the circle by 1/(2pi i). Its mass is1 and

    CT_z(p+z^r+t z^-s)^m
      =sum_l m!/[(m-M*l)!(s*l)!(r*l)!] p^(m-M*l)t^(r*l).

Since t^r=q, (1) is recovered on ALL orders. Thus for every polynomial R,

    h(R(P))=integral_C R(S) Theta/(2pi i).              (3)

C is a compact rational cycle; multiplying it by r gives an integral
cycle of nonzero mass, sufficient for the projective-period theorem.
Perturbing M0 preserves every joint period since all polynomial multiples
of the holomorphic top form are closed. No weighted-Haar, marker,
subsequence, scalar norm or primitive-purity reflection is used.

## 2. Finite coefficients: the original form on a normalized monomial model

Compactify the base sphere and resolve the rational map [1:p^M:q] to P2.
Resolve all boundary divisors, divisors of the relevant functions and
the volume divisor to SNC. Over points where p,q are finite, Omega is
holomorphic by the retained neutral-height volume theorem (q is
nonconstant in the application). Resolve q=0 on this smooth base. Locally
q is a unit times x^i y^j. Absorb its unit into A by an analytic root;
regularity of an algebraic meromorphic form can be checked analytically.

The normalization of A^s B^r=x^i y^j is toric. On its dense torus,
Theta is, up to a nonzero constant and a regular coefficient,

    x y * dlog x wedge dlog y wedge dlog A.            (4)

On a smooth toric refinement its order at a primitive ray v is
v(x)+v(y)-1. The cone satisfies nonnegative valuations of A,B,x,y and
s v(A)+r v(B)=i v(x)+j v(y). A nonzero ray cannot have v(x)=v(y)=0,
since then both other valuations vanish. Thus v(x)+v(y)>=1 and (4)
has no pole. Regular coefficients of Omega only add zeros.

For precision, a smooth refinement can be constructed as follows.
Triangulate the cone. In any nonsmooth simplicial cone choose a nonzero
lattice point of its half-open fundamental parallelepiped and subdivide
at its primitive ray. Its coefficients are all <1; every replacement
cone has strictly smaller lattice multiplicity. If the point is on a
face, the same coefficient calculation applies in all affected cones.
The number of cones of maximal multiplicity decreases, then the maximal
multiplicity decreases; iteration terminates. Refinement with unchanged
support is proper by the toric valuative criterion. This supplies the
local resolution used in (4). The basic cone, character-divisor and
smooth-chart construction is the classical one in
[Cox's minicourse, LectureI Sections4--5](https://dacox.people.amherst.edu/lectures/toric.pdf).
The ray-volume calculation and application here are written arguments;
no unread classification of singularities is imported.

## 3. A stronger infinity estimate from the ACTUAL charged degrees

Let D be the sphere degree of q. Its leading zero multiplicity m at
each point of the infinity conic satisfies m<=D/3 for O and <=D/6 for I.
Let delta_U,delta_V be the minimum matrix degrees of U,V. The associated
graded ring of C[SL2] is C[x,y,z,w]/(xw-yz), a domain, so minimum
degrees ADD under multiplication. A neutral polynomial of minimum
matrix degree 2D has sphere degree D: torus projection gives equal row
degrees, and every row-balanced monomial factors into quadratic row
products; the converse is the quadratic sphere presentation. Hence

    D=(s*delta_U+r*delta_V)/2,
    D/M>=4 on 2O, D/M>=6 on 2I.                      (5)

The last inequalities use the retained minimum nonconstant invariant
matrix degrees8/12. They hold without a homogeneity assumption on U,V.

For any boundary valuation E on the projective sphere put u=ord_E(B0)>0
for its boundary equation, h=ord_E(H_q), and let L be its log discrepancy.
The original sphere-volume order is k=L-1-2u. If q has pole beta=Du-h>0,
then

    k+1+beta/M=L-2u+(Du-h)/M>0.                       (6)

Indeed if L>=2u this is immediate. Otherwise 1<=L/u<=2. At a leading
root, the elementary threshold bound and m<=D/3 give h<=DL/3. Writing
c=D/M>=4 and t=L/u, the right side is at least

    u[t-2+c(1-t/3)] >= u*c/3>0.

The expression decreases on [1,2], so its minimum is at2. Outside the
leading roots h=0 and the inequality is stronger. Thus (6) holds at
EVERY boundary divisorial valuation, not just on one selected resolution.
Every finite-q valuation still has k>=0 by the retained finite-volume
theorem.

For p nonconstant retain the earlier all-valuation half-pole estimate

    p has pole beta_p>0 => k+ceil(beta_p/2)>=0,         (7)

as proved in the quadratic-suspension theorem. Estimates (6) and (7)
serve different branches of the higher-degree cover; replacing (6) by
an arbitrary fractional-height assertion would be unjustified.

## 4. All infinity branches of the higher-degree cover

The equation A^s(S-p-A)^r=q is monic up to sign in A of degree M.
Its normalization over the compactified base times A1_S is therefore
finite and projective over A1_S. We show its ORIGINAL Theta extends
on a smooth projective resolution at every finite S boundary point.
The resolved map [1:p^M:q] gives the following charts. Further base
resolutions below preserve (6),(7), which hold at ALL valuations.

If q dominates and p^M/q vanishes at the point, adjoin h0=q^(1/M).
The functions p/h0 and S/h0 are regular on the normalized base cover;
at infinity they tend0. Set alpha=A/h0, beta=B/h0. Then
alpha^s beta^r=1 and alpha+beta=(S-p)/h0. At their limiting roots
alpha+beta=0, r alpha-s beta=M alpha!=0. The finite root cover is
etale there and Theta is a unit times (Omega/h0) wedge dS.
This form is holomorphic on a smooth resolution of the root cover:
over a base divisor with ramification e its exact order is

    e(k+1+beta_q/M)-1>=0,                              (8)

by (6) and integrality. Every divisor on the finite base cover restricts
to a base divisorial valuation after a birational base resolution.
Thus this check includes all its exceptional valuations. Finite further
etale root choices cannot introduce a pole. If p^M/q is a unit instead,
the p-dominant chart below applies. Points where neither function has
a pole were treated in Section2.

In the p-dominant infinity chart, f=1/p is regular and vanishes, and
lambda=q/p^M is regular. Set alpha=A/p, beta=B/p and T=Sf-1. Then

    alpha+beta=T, alpha^s beta^r=lambda,
    Theta=Omega f wedge dS/(r alpha-s beta).           (9)

Whenever the denominator is a unit, (7) makes Omega f holomorphic;
its pullback to the normalized finite cover remains holomorphic. This
also covers roots alpha=0 or beta=0: their other coordinate is a unit
since T=-1 at f=0, and the denominator in (9) is nonzero.

There is only ONE possible nonzero critical value. Put

    c=s^s*r^r/M^M, alpha_c=sT/M, beta_c=rT/M.

The derivative of alpha^s(T-alpha)^r is
alpha^(s-1)(T-alpha)^(r-1)(sT-M alpha). At alpha_c its second derivative
is -M alpha_c^(s-1)beta_c^(r-1)!=0. Therefore, near this critical root,
polynomial division and taking a local unit square root give a coordinate
v with

    v^2=lambda-c T^M,
    r alpha-s beta=unit*v.                            (10)

This is an actual local double-branch calculation, not semisimplicity
of the whole degree-M cover. Set eta=lambda-c(-1)^M. Then

    v^2=eta+f L(S,f),
    partial_S L(S,0)=-c M(-1)^(M-1)!=0.               (11)

Principalize (f,eta) on the BASE. If eta=f h, the cover is
v^2=f(h+L); its second branch divisor is transverse to all components
of f=0 by (11). Square removal gives volume exponents
k+ceil(ord f/2)>=0, exactly (7).

If f=eta h, the cover is v^2=eta(1+h L). At a boundary point where h
vanishes the second factor is a unit. Where h is a unit, its S derivative
is a unit and it is transverse to the boundary. Resolve eta,f and the
volume divisor to SNC; this is legitimate after the chart principalization.
The volume exponents after square removal are

    k+ord f-floor(ord eta/2)
        >=k+ceil(ord f/2)>=0,                         (12)

because ord eta<=ord f componentwise in this chart. The retained
explicit SNC-double-cover blowups preserve the residue form, so (12)
controls the corners as well as generic points. Local unit factors in
(10) are regular and nonzero and do not affect pole orders.

These charts cover every finite-S infinity point. Local holomorphic
extension transfers to a global smooth projective resolution by common
resolutions and the divisorial criterion for regular top forms. Choose
it to retain the compact cycle in its smooth open locus. The retained
projective ORIGINAL-volume period theorem and (3) now exclude EVERY
nonzero polynomial H(P) whenever U,V are both nonzero.

## 5. Complete field law and every coefficient ring

For both exceptional algebras, U,V nonzero imply q=U^s V^r is
nonconstant: a constant nonzero product would make nonzero-weight
functions units. The preceding proof therefore applies. If UV=0, one
charged function is zero and weight projection reduces all H(P) moments
to those of H(p). Complete neutral closure gives the exact field law

    H(P) eventually pure iff
       H=0, or [UV=0, p=c constant, H(c)=0].           (13)

All surviving phases are one-sided or zero and safe against EVERY
original marker. In particular P itself is eventually pure iff p=0
and UV=0. The symmetric law is retained as the case r=s=1.

Over ANY commutative C-algebra R, the exact criterion is

    P eventually pure iff p AND q=U^s V^r are nilpotent.        (14)

Necessity is the same finite Qbar coefficient-field embedding at every
prime, followed by finite nilpotent-coefficient lifting, as in the
symmetric theorem. No embedding of the whole R or reducedness is
assumed. Conversely, if p^nu=0 and q^kappa=0, (1) has no surviving
term once m>nu-1+M(kappa-1).

For a fixed original marker Q of maximum absolute weight B_Q, separate
the neutral count c<=nu-1. The remaining charged part has the exact
two-weight bound already proved: if its total count is at least
M*kappa+ceil(B_Q/min(a,b)), weight zero forces a factor q^kappa.
Thus a sufficient all-original-marker cutoff is

    m>=nu-1+M*kappa+ceil(B_Q/min(a,b))
       =>h_R(QP^m)=0.                                (15)

The earlier symmetric cutoff is sharper in its special case and is not
replaced. Nilpotence characterizes EVENTUAL purity, not the first moment.
No uniform numerical nilpotence exponent is computed. The field all-H
classification is not copied verbatim to arbitrary rings.

If the two charged weights have the SAME sign, their sums cannot
contribute weight zero. Eventual purity is then exactly nilpotence of
p by complete neutral closure and residue lifting. If p^nu=0, every
marker vanishes when m>nu-1+B_Q/min(|a|,|b|). Together with (14),
this proves MZ for EVERY support contained in {0,a,b}, in all degrees
and over every R. Zero or missing charged parts are included. Haar
transport gives the same result for every automorphic image and every
maximal left torus; a nonlinear automorphism need not retain the same
matrix-degree bound.

## 6. Consequence for the full residual problem

A full exceptional failure must, after EVERY algebraic automorphism
and for EVERY maximal left torus, either have at least four weights or
have THREE weights with NO zero weight. The latter case is still open;
this is not an unconditional four-weight lower bound. Larger supports,
general mixed Bezout parameters and distinct Hopf nonlinear composites
are not closed by this weighted-curve argument.

The positive and negative full targets remain unchanged: prove every
actual exceptional eventually-pure phase has ZERO corrected local
component/torus-valuation vector, or produce an actual invariant pure
phase with a fixed original bad marker. The Gaussian thread's general
formal-symbol and rational-slope results are retained inputs; this
exceptional thread does not duplicate its constructive Bernstein work.

## 7. The exact full-class coverage and preservation lemma still missing

The new proof covers a zero-weight term plus TWO charged terms. It
depends on their weighted fibre map having one nonzero SIMPLE critical
value. A general actual phase has a finite Laurent expansion on the
original torus torsor, with coefficients that are sections of associated
line bundles. A rational trivialization can introduce poles that the
original compact Haar cycle cannot avoid. Several charged terms can
also have multiple or colliding critical branches, whose local volume
orders are not controlled by the calculation in (10)--(12).

A sufficient FULL exclusion/safety bridge would prove the following
dichotomy for EVERY actual exceptional phase f: either certify its
ALL-original-marker eventual safety, or construct a compact nonzero-mass
period on a smooth projective map with ORIGINAL holomorphic volume
whose integral of R(S) equals h(R(f)) for EVERY polynomial R. The second
branch excludes eventual purity; the first must retain the actual markers.
Neither coverage of this dichotomy nor its moment-preservation step is
proved for the residual multi-term class. Requiring the proper-model
branch for ALL phases is false: the known pure safe one-weight phases
already prohibit it.

A more flexible positive bridge could instead preserve the FULL actual
unit extension/marked-period module before quotienting by polynomials
and force its corrected local orders zero under unweighted purity.
It must retain literal1 AND actual global relations; matching only the
scalar unit Cauchy function is insufficient. On a pure phase that scalar
function can already be rational even when original marked tails fail.
No such full-module bridge is supplied here. These are concrete missing
coverage/preservation obligations, not equivalent criteria claimed as
a solution or an external-review gate. Both FULL MZ cases remain OPEN.
