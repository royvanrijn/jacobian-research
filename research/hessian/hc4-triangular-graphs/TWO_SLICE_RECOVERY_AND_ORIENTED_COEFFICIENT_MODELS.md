# Two-slice global recovery and oriented coefficient models

15 November 2026. Claims **HC4MYJV0/1/2**. The full quotient/torsor
certificate now recovers global source coordinates. A coordinate inverse
does not have to be supplied separately once the exact two-slice identities
below are proved. Mere quotient type, annihilator support or local periods
still do not supply those identities.

More generally, the reduced plinth and missing-line geometry recover a
global polynomial plane bundle over the actual invariant base. The remaining
obstruction is an oriented generating-pair problem for the coefficient ideal.
This gives an all-degree reduction, not a constructed nonstandard action.

## 1. Statements

Throughout R=C[X1,X2,X3,X4], D is a nonzero locally nilpotent derivation,
and its actual kernel is A=B[h], B=C[b1,b2], with b1,b2,h algebraically
independent. Equality is of embedded algebras, not abstract isomorphism.

**JV0, general recovery.** Suppose I=(F,G) in B is radical with N>=2
points, pl(D)=IA, V(IR) has codimension at least two in Spec R, and the
quotient (b,h):Spec R -> Spec A reaches exactly one point (b_i,c_i) on
each line {b_i} times A1_h. Then B is a coordinate subalgebra:

    R=B[U,V], h=d0(b)+U f(b)+V g(b), (f,g)=I.             (1)

Both the coordinate system and its polynomial inverse exist globally;
D is a B-linear translation in these coordinates. No constant-Jacobian
condition on F,G is needed for this recovery theorem. The plinth,
codimension and missing-image assumptions are substantive global hypotheses.

Choose p,q in R with Dp=G,Dq=-F, and set a=Fp+Gq in A. The class
of a modulo IA is independent of those choices. Its components are

    a(b_i,h)=lambda_i(h-c_i), lambda_i in C*.            (2)

With arbitrary polynomial interpolants lambda,d of lambda_i,-lambda_i*c_i,
the entire source has the global presentation

    R=B[h,P,Q]/(FP+GQ-lambda(b)h-d(b)).                  (3)

The row (F,G,-lambda) is unimodular. A free basis of its kernel gives
the coordinates in (1). If mu in C* is the determinant of a compatible
three-column completion, their coefficient Jacobians satisfy

    J_b(f,g)(b_i)=-mu J_b(F,G)(b_i)/lambda_i.            (4)

Thus the point residue slopes lambda_i/J_b(F,G)(b_i), not just the
abstract plinth ideal, retain the unsquared coefficient orientation.

**JV1, elementary exact certificate.** Suppose instead that F,G in B
are coprime and there are global p,q with Dp=G,Dq=-F for which

    Fp+Gq = ell h+d(b) modulo IA, ell in C*.             (5)

Then there are invariant corrections to p,q producing P,Q with
FP+GQ=ell h+d, and R=B[P,Q]. This conclusion needs neither JV0's
image assumptions nor the projective-module theorem used there.

If J_b(F,G)=j in C*, (5) therefore produces an actual unit-Hessian
source member, with determinant ell^-4 j^2. Conversely, for h with
N>=2 finite Morse critical points, every unit-Hessian source member
supplies these data through Tang and JT2. Hence this kernel/two-slice/
affine-socle certificate is equivalent to source-orbit unit Hessian.
It replaces a separately supplied full coordinate inverse by global
identities that prove and construct that inverse.

For the JU setting with Keller generators, equality of the full socle
class and its common affine residues is sufficient: it is exactly (5),
by the injective Cech map A/IA -> A_FG/(A_F+A_G). Matching only the
annihilator or only an abstract quotient is not this certificate.

**JV2, oriented pairs and universal models.** Let f=(f1,f2) generate
a radical ideal I of N>=2 points in B. Two ordered generating pairs
f and k of I lie in the same GL2(B) orbit if and only if

    J_b(k)(b_i)/J_b(f)(b_i)

is one common nonzero constant. This is an unsquared condition. For
h=d0+U f1+V f2 in actual full coordinates, its whole B-preserving
polynomial source orbit contains a unit-Hessian member exactly when
the GL2(B) orbit of f contains a globally Keller pair. Even arbitrary
nonlinear fibre automorphisms add no further coefficient pairs to this
test. Thus nonconstant J_b(f)(b_i) excludes the whole relative orbit.

For an explicit all-N uniform-orientation model put

    C_N(z)=product_{i=1}^N (z-i),
    a_N(z) C_N'(z)=1 modulo C_N(z), deg a_N<N,
    T_N=a_N(z)x P+C_N(z)Q.                              (6)

For each N>=2 the following are equivalent: T_N has a unit-Hessian
member in its B=C[x,z]-preserving source orbit; it has one via a
B-affine source change; and some plane Keller map has a fibre of
exactly N points. No N is enumerated or asserted to occur. Existence
for at least one N is equivalent to failure of JC2, hence of HC4.
These auxiliary models are not asserted to be realizable by the original
affine graph families, or to be equivalent to K2,h1 or S_rho.

## 2. Retained and imported inputs

Retain [JT](PRIMITIVE_QUOTIENTS_PLINTH_AND_GLOBAL_ACTION_CERTIFICATES.md),
[JU](QUOTIENT_HOLES_TORSOR_RESIDUES_AND_LOCAL_PERIOD_BARRIER.md),
[JS](GLOBAL_ADDITIVE_ACTIONS_AND_THE_K2_TRANSLATION_BARRIER.md),
[HR](RELATIVE_CRITICAL_GRAPH_CLASSIFICATION.md),
[HG](GLOBAL_VALUE_ORBITS.md),
[JK](SUSPENDED_MULTIPLIER_ORBITS_AND_BASE_PRESERVING_CLOSURE.md), and
[Tang](../tang-hc4-audit/STRUCTURAL_REVIEW.md). The elementary slice
theorem and JU Cech injection are reused without replaying their checkers.

JV0 additionally uses these established algebraic results, with their
hypotheses checked below:

- Quillen--Suslin: finite projective modules over C[b1,b2] are free.
  [Dicks, Theorem 6 and proof](https://arxiv.org/html/2112.02272v1)
  reproduces the theorem and identifies Quillen 1976 and Suslin 1976.
- A hypersurface in a polynomial ring is Cohen--Macaulay, by the
  regular-sequence quotient result in [Stacks 00N7](https://stacks.math.columbia.edu/tag/00N7).
- The normality criterion [Stacks 031S](https://stacks.math.columbia.edu/tag/031S)
  applies to a Cohen--Macaulay domain regular in codimension one.
- A normal Noetherian domain is the intersection of its height-one
  localizations, [Stacks 031T](https://stacks.math.columbia.edu/tag/031T).
- For the final existential bridge, the quasi-finite normalization
  factorization in [Stacks 02LR](https://stacks.math.columbia.edu/tag/02LR)
  makes a separated birational quasi-finite map to a normal target
  an open immersion.

These are explicit proved imports, not new claims or an unproved
cancellation assumption. JV1 and the oriented-pair argument do not use
Quillen--Suslin. No action classification is imported.

## 3. From the two slices to the entire source

Since F,G belong to the plinth, choose p,q as above. Then a=Fp+Gq
is invariant. On the two opens the slice theorem gives exactly

    R_F=A_F[q], R_G=A_G[p],
    t_G=p/G, t_F=-q/F, t_G-t_F=a/(FG).                  (7)

First prove JV1. By (5), write a-ell h-d=F alpha+G beta with
alpha,beta in A, and set P=p-alpha,Q=q-beta. These are global
polynomials and have the same D derivatives. Now FP+GQ=ell h+d.
The four elements b1,b2,P,Q are algebraically independent: over Frac B,
their invariant H=FP+GQ is transcendental and a slice is independent
of H. Put R0=B[P,Q]. Equation (7) gives

    R_F=(R0)_F, R_G=(R0)_G.

All rings lie in the same fraction field. Because F,G are coprime
in the polynomial ring R0, its two localizations intersect in R0.
Thus R0 subset R subset (R0)_F intersect (R0)_G=R0. This proves
global generation, not just agreement of generic fibres.

The inverse is constructive. Express each native X_i in the two slice
charts and substitute h=(FP+GQ-d)/ell. This yields A_i/F^m and B_i/G^n
with numerators polynomial in b,P,Q. Equality in the common fraction
field gives G^n A_i=F^m B_i. Coprimality forces exact divisibility of
A_i by F^m; its quotient is the required polynomial inverse coordinate.
The two inverse identities follow from generation and algebraic
independence. The mixed Hessian block of h in these actual coordinates
is ell^-1 D_b(F,G), and its vertical block is zero. Its determinant
is ell^-4 j^2. This proves JV1, including the stated converse via JT2.

## 4. General plinth-hole recovery

For JV0 use codimension at least two in the polynomial source to obtain
R=R_F intersect R_G. A function regular on the two opens has no possible
prime denominator: no source prime divisor contains both F and G.

Let xi=[a/(FG)] in H^1(Spec A minus V(IA),O). Then

    pl(D)=ann_A(xi).                                    (8)

One direction follows by subtracting f t_F,f t_G from a global local
slice with derivative f. Conversely, if f xi=0, invariant corrections
glue f t_F,f t_G to a function in R_F intersect R_G=R with derivative
f. This proves (8). The JU injection, applied over B[h], gives
ann_A(xi)=IA:a. Since pl(D)=IA and A/IA=product_i C[h], every
a_i(h)=a(b_i,h) is a nonzero polynomial.

Form the integral hypersurface

    S=B[h,P,Q]/(FP+GQ-a(b,h)).                          (9)

Its equation is primitive and linear in P,Q over the UFD B[h], so
irreducible. Its localizations S_F,S_G map isomorphically onto R_F,R_G,
identifying P,Q with p,q. Hence S embeds in R with the same fraction
field. The locus V(F,G) in S is a finite union of planes indexed by
the roots of the nonzero a_i; its dimension is at most two.

Outside that locus S is smooth. At b_i and a root h=r of a_i, the
P,Q derivatives vanish, the h derivative is -a_i'(r), and the two
base derivatives are an invertible linear system in P,Q. Invertibility
holds because F,G generate the reduced point ideal locally. A singular
point can therefore occur only at a repeated root, with a unique P,Q.
The singular locus is finite. S is a four-dimensional hypersurface
domain, hence Cohen--Macaulay and regular in codimension one. The
imported normality criterion proves S normal. Its removed locus has
codimension at least two, so the intersection theorem gives

    S=S_F intersect S_G, hence R=S.                     (10)

The native polynomial source is smooth. Thus none of the indicated
singular points occurs: every a_i is squarefree. Formula (10) also
identifies the reached points on each quotient line with the roots
of a_i. The one-point image hypothesis now forces (2).

Choose lambda,d by the Chinese remainder theorem. Subtract invariant
corrections from P,Q as in Section 3 to replace a by lambda h+d.
Since lambda is nonzero at every point of V(I), the ideal
(F,G,lambda) is B. Let r=(F,G,-lambda) and choose a column e with
r e=1. Its kernel K is a finite projective B-module of rank two:
B^3=B e direct-sum K. Quillen--Suslin gives a free basis v1,v2 of K.
The matrix M=(e,v1,v2) is invertible over B, with determinant mu in C*.
Every solution of r(P,Q,h)^T=d is uniquely

    (P,Q,h)^T=d e+U v1+V v2.                           (11)

The inverse matrix gives polynomial inverse formulas for U,V. Together
with (10) this proves a full source coordinate system and inverse.

D acts on (P,Q,h) by (G,-F,0), so it is a B-linear translation in
U,V. Its coefficient gcd is one: otherwise that factor would divide
both F and G. Its kernel is the polynomial algebra generated by its
linear invariant. Since it is B[h], one-variable generator comparison
forces h to be affine-linear in U,V with constant nonzero scale.
The ideal of its two coefficients is I: compare the image ideals of
D in the two global presentations and contract to B. This proves (1).

For (4), take f=(v1)_3,g=(v2)_3 in (11). At each b_i these vanish.
Differentiating r v_j=0 at b_i gives
lambda_i d(v_j)_3=(v_j)_1 dF+(v_j)_2 dG. The third component of
v1 cross v2=mu r^T equals -mu lambda_i. Taking the two determinants
gives (4). This also checks its unsquared sign and common factor.

## 5. Oriented generating pairs and every relative source degree

Write k=M f for some polynomial two-by-two matrix, since both pairs
generate I. At each reduced zero,

    J_b(k)=det(M) J_b(f).

Necessity of a common nonzero ratio under GL2(B) follows. For sufficiency
suppose this ratio is u in C*. Then u-det M belongs to I=(k1,k2).
Choose r,s with u-det M=s k1-r k2. The matrix

    M'=M+(r,s)^T(-f2,f1)

still maps f to k, and its determinant is det M+s k1-r k2=u.
It is therefore invertible over B. This proves the full pair criterion
in every degree; no stable-range assertion is being assumed.

Now take h=d0+U f1+V f2 in actual coordinates. Its relative critical
scheme is N disjoint affine planes over the base points, with zero
vertical Hessian along them. Its N full critical points are Morse,
one above each zero, because J_b(f)(b_i) is nonzero. Under any
B-preserving source change the same relative scheme and zero vertical
block persist, by relative Jacobian transformation and congruence.

If the changed function g has unit Hessian, its gradient is not an
automorphism since N>=2. Apply the retained Tang plane W and the HR
three-rank calculation. Rank two projection of W to the base would
make the relative scheme one A2 graph, contradicting its N components.
Rank one, with zero vertical Hessian along the relative scheme, forces
the coefficient C_U to vanish on C=0 and then identically, making
the Tang pair triangular and invertible. This is again impossible.
Rank zero remains: W is the fibre plane and g is affine-linear in
the two fibre variables, with an everywhere Keller coefficient pair.

Any polynomial fibre automorphism producing an affine-linear g has
coefficient pair given by its linear part at U=V=0. That matrix is
in GL2(B), since its fibre Jacobian is a nonzero constant. A base
automorphism composes the pairs by a polynomial base automorphism
with constant Jacobian and does not enlarge the existence criterion.
Conversely a GL2(B) coefficient change is realized by its transpose
on the fibre variables and has a full polynomial inverse. Its full
Hessian is a unit precisely when its new coefficient Jacobian is.
This proves the whole-relative-orbit criterion in JV2.

In particular J_b(f)(b_i) must be the same unsquared constant. If the
native full critical Hessians agree, their squares already agree;
mixed signs still fail this stronger criterion. This explains the
old h0,h1 base obstruction: f=(2x,z^2-1) has Jacobians +4 and -4.
It says nothing about source changes that replace the actual base.

## 6. Explicit uniform-orientation all-N models

The inverse a_N in (6) exists since C_N is squarefree. The pair
f_N=(a_N(z)x,C_N(z)) generates (x,C_N), and its Jacobian is
a_N C_N', equal to one at every zero. T_N has exactly N Morse
critical points, all of value zero and full Hessian determinant one.
Its native full Hessian determinant is (a_N C_N')^2, generally not
constant. No global Hessian condition is inferred from its point values.

Every ordered finite set of N distinct plane points can be sent to
(0,1),...,(0,N) by an actual polynomial automorphism. Here is a direct
construction. Choose a determinant-one linear frame with separating
second coordinates v_i; shear the first coordinate to make it zero
at the points. Then three interpolating shears successively send
(0,v_i) to (i,v_i), to (i,i), and to (0,i). Each shear has a polynomial
inverse and determinant one. This is a finite interpolation argument
for arbitrary N, not a point-set classification assumption.

If a Keller pair has exactly those N zeros after this source change,
it generates (x,C_N) by reducedness. Its point Jacobian is one common
nonzero constant. The oriented-pair lemma puts it in the GL2(B) orbit
of f_N. The associated transpose fibre map therefore produces a
unit-Hessian member of T_N's relative source orbit. The converse follows
from the relative criterion: its Keller coefficient pair has precisely
N zeros. The construction even gives the B-affine equivalence in JV2.

A nonautomorphic plane Keller map is generically finite of degree
at least two. Indeed degree one gives a birational etale morphism into
normal A2, hence an open immersion; the absence of nonconstant units
in its polynomial source rules out a missing divisor, and normal
codimension-two extension rules out a remaining finite complement.
Equivalently one may use the classical birational Keller theorem.
Over a nonempty finite-etale target open its fibre therefore has some
N>=2 points. Conversely a fibre with N>=2 forbids injectivity. This
proves the existential JC2 bridge; individual fibre-cardinality questions
remain open. The standard birational/quasi-finite open-immersion fact
used in this last bridge is recorded separately in the source context.

No h_N model is asserted to lie in the original affine graph frontier.
For that frontier the four old tests,K2 and S_rho remain the actual
models. The result explains why uniform oriented residues lead toward
the plane Keller issue, while the full source frame is recoverable from
the specified global action data. It does not classify those actions.

## 7. Assurance and next boundary

The checker controls the exact correction,localization,unimodular-row,
orientation and model identities. The all-degree conclusions use the
written proofs and the explicit algebraic imports, not bounded searches.
No necessary nonstandard action for K2 or S_rho has been constructed.
No coefficient ideal with N>=2 is proved to admit Keller generators.
No endpoint or unrestricted HC4/JC2 problem is solved. Independent
replay,formal verification and external review remain absent.
