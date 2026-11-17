# Every Wright Jacobian pair must mix grading weights

14 November 2026. PF2GRADE1 is an all-height obstruction on the uniquely
marked ring of PF2WRIGHT1. No nonconstant homogeneous element of
A=C[Y,XY,X^2Y,X^3Y+X] has a Jacobian mate in A. Consequently no polynomial
target presentation of a candidate Keller pair can have a homogeneous
coordinate, even after adding a constant. The unrestricted mixed-weight
problem, finiteness and actual ruled compatibility remain open.

The core proof is elementary and uses no degree cutoff. Only the final
corollary on a grading-invariant target algebra imports plane C*-action
linearization, with its exact primary source pinned below.

## The full weight spaces

Write Z=X^2Y, T=XY, D=X^3Y+X. The grading field

Delta=-X*d/dX+2Y*d/dY

has weights2,1,0,-1 on Y,T,Z,D. Its source-volume divergence is one.
The relations T^2=YZ, TD=Z(Z+1), YD=T(Z+1) show, for n>=1,m>=0,

\[
 \begin{aligned}
 A_0&=\mathbb C[Z],\\
 A_{-n}&=D^n\mathbb C[Z],\\
 A_{2m}&=Y^m\mathbb C[Z],\\
 A_{2m+1}&=T Y^m\mathbb C[Z].
 \end{aligned}\tag{1}
\]

These are equalities of full weight spaces, not a bounded monomial list.
For instance a monomial Y^b*T^a*D^d of weight -n has
 d=a+2b+n, hence equals D^n*(TD)^a*(YD^2)^b times a polynomial in Z.
For positive weight, reduce T^2=YZ and cancel D against T or Y to
obtain the displayed positive generator. The source expressions in
C[X,Y] show the remaining coefficient space is exactly C[Z].

## No homogeneous constant-Jacobian pair

If nonzero homogeneous P,Q have weights r,s, their Jacobian has weight
r+s-1. A nonzero constant therefore requires r+s=1. After swapping,
put r=-n, s=n+1, n>=0. All possibilities in(1) have the form

\[
 P=X^n f(Z),\qquad Q=X^eY^l g(Z),
 \qquad 2l-e=n+1,
\]

where e is0 or1, l>=1, f=(Z+1)^n p(Z), and p,g are nonzero polynomials.
For n=0 this means P=p(Z), Q=XY*g(Z).
Direct differentiation gives the universal identity

\[
 \boxed{J(P,Q)=Z^{l-1}
  [nlfg+nZfg'+(n+1)Zf'g].}\tag{2}
\]

If n>=1, put d=deg(f)>=n and c=deg(g)>=0. The bracket has degree d+c
and leading coefficient equal to the nonzero product of the two leading
coefficients times

nl+nc+(n+1)d>0.

Thus the Jacobian is a nonconstant polynomial in Z. In characteristic
zero the positive integer factor cannot cancel, over arbitrary complex
coefficients. If n=0, (2) is Z*p'(Z)*g(Z), either zero or divisible by Z.
Neither case permits a nonzero constant. This exhausts all integer weights.

## One homogeneous coordinate already suffices

Suppose P is homogeneous and Q=sum Q_s is its finite weight decomposition.
The equation J(P,Q)=c!=0 in weight zero forces

J(P,Q_(1-r))=c.

Every Q_s belongs to A, so this contradicts the homogeneous-pair result.
The same holds for P+c0, since a target translation removes the constant.
It holds in every polynomial target presentation: target automorphisms
preserve the ring and multiply the source Jacobian by a nonzero constant.
It is not a prohibition on arbitrary homogeneous elements of the target
algebra that are not polynomial coordinates.

## Both extreme pairs commute and the weight range straddles one

Center P,Q at the source origin and denote their minimum and maximum
nonzero grading weights by r_-,r_+ and s_-,s_+. A hypothetical Keller pair
must satisfy the strict inequalities

\[
 \boxed{r_-+s_-<1<r_++s_+.}\tag{3}
\]

If the first sum were greater than one, there would be no weight-zero
Jacobian term. If it equalled one, the only weight-zero term would be the
Jacobian of the two minimum components, again an impossible homogeneous
mate. The maximum case is the same argument. Moreover the unique lowest
and highest weight Jacobian terms must vanish, so

J(P_(r_-),Q_(s_-))=J(P_(r_+),Q_(s_+))=0.

These extreme homogeneous pairs are algebraically dependent. They are
necessary constraints on mixed-weight candidates, not an exclusion of all
such candidates.

## No polynomial Keller target algebra is invariant under this grading

Let B=C[P,Q] be a polynomial subalgebra of A with J(P,Q) a nonzero
constant. Suppose Delta(B) is contained in B. Polynomial interpolation
in the finitely many eigenvalues of each element isolates every weight
component inside B. Thus B is graded, and the induced action on
Spec(B)=A2 is an algebraic C*-action.

Import Gutwirth's plane linearization theorem, in the independently proved
form of Flenner--Zaidenberg,
[On the uniqueness of C*-actions on affine surfaces, arXiv:math/0406239v2,
Corollary4.6](https://arxiv.org/html/math/0406239v2).
It supplies polynomial plane coordinates which are homogeneous for the
action. Their pullbacks are homogeneous elements of A and still have
nonzero constant source Jacobian, a contradiction. Therefore

\[
 \boxed{\Delta(\mathbb C[P,Q])\not\subseteq\mathbb C[P,Q].}\tag{4}
\]

Only this corollary uses the external linearization theorem. The
[source manifest](sources/grading-linearization-manifest.json) pins that
primary result; its full proof is not independently replayed here.
Finiteness and degree four are not needed for(1)-(4).

## Relation to the B-pole filtration

Use f0,f1,f2,f3 from PF2WRIGHT1, each of B-pole at most one.
The same grading field acts by the exact identities

\[
 \begin{aligned}
 \Delta(f0)&=2f0+2f2^2-4f1f3,\\
 \Delta(f1)&=f1-2f3+4f2f3,\\
 \Delta(f2)&=-2f3^2,\\
 \Delta(f3)&=-f3.
 \end{aligned}\tag{5}
\]

Section generation and the Leibniz rule give Delta(V_m) subset V_(m+1)
for every m. This permits a finite calculation in each pole layer,
but does not imply the target algebra is grading-stable. In fact(4)
proves stability is impossible for a Keller pair.

## Verification boundary and next gate

The [producer](compute_wright_grading.py),
[exact certificates](wright_grading_certificates.json) and
[checker](verify_wright_grading.py) verify(2),(5), bounded weight-module
controls, degree arithmetic and the pinned external source. Weight-space
coverage, the all-height degree argument and the stability corollary are
written/imported proofs, without formal verification, independent whole-
proof implementation or external review.

The separate [trace audit](TRACE_GRADING_BOUNDARY.md) explains why field
trace has constant divergence but does not automatically preserve the
branch curve. Consequently no branch Euler theorem is inferred merely
from this grading. The remaining compatible quartic map must be genuinely
mixed in grading weights and have non-invariant target algebra; the
first open pole pairs remain(6,8),(7,8).

The subsequent [all-height pole-three theorem](POLE_THREE_ALL_HEIGHTS.md)
excludes every(3,n) presentation, independently of the2022 degree bound.
The current first open reduced pairs are (6,8),(7,8).

The adopted frontier(6,8),(7,8) incorporates the subsequent
[PF2RULED12 full(5,8) exclusion](POLE_FIVE_EIGHT_EXCLUSION.md); it is
not an additional exclusion proved by this earlier note.
