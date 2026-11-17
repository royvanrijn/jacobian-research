# Spectral denominators for rank one projective Weyl carriers

15 November 2026. **QPD1.** The scalar spectral-denominator theorem extends
to endomorphism algebras of rank-one finite projective Weyl modules.
Every rational operator with finite joint adjoint orbit has a denominator
polynomial in the commuting images. Thus rational individual mates exist
exactly when the global spectral translation classes are torsion.

Rational commutants are automatically fractions of actual commutants.
Together with scalar idealizer units, this extends QTB1's birational,
Galois-factor and radical tests to nonfree projective carriers. A scalar
free frame and finite generation of the centralizer are unnecessary.
This does not turn spectral torsion into vanishing, remove every
rational denominator from a mate, or construct a finite Long carrier.

The scalar theorem is established machinery:
[Bavula, Lemma3.7](https://arxiv.org/html/math/0401263).
The argument below transports its denominator-ideal mechanism through
a rank-one Morita equivalence. The unit fact is an immediate corollary
of the retained QHI1 symbol theorem, not a new independent unit theorem.

## The actual endomorphism algebra

Let L have characteristic zero, A=A_n(L), D=Frac A, and let P be a
nonzero finite projective left A-module of generic rank one. Put
\[
                         E=\operatorname{End}_A(P)^{op}.  \tag{1}
\]
A fractional-ideal realization identifies E with right multipliers
of a nonzero left ideal I in D. Since A is simple, a nonzero finite
projective module is a generator: its trace ideal is nonzero and
therefore A. Hence P gives a Morita equivalence between left E-modules
and left A-modules.

Choose n commuting algebraically independent F_i in E and set
\[
 P_F=L[F_1,\ldots,F_n],\quad S=P_F\setminus\{0\},\quad
                         B=\operatorname{Cent}_E(F_1,\ldots,F_n). \tag{2}
\]
The notation P_F is the position-image polynomial algebra, distinct
from the projective module P. No finite-generation assumption on B
is imposed. We assume B commutative when invoking QTB1's field
conclusions; QSC1 supplies this for Long deformations in the stated
faithful formal frame.

## Scalar units and scalar algebraic constants

QHI1 proves, in a normalized fractional realization, gr E is a
subalgebra of L[x,p] and E in Bernstein degree at most zero is L.
Degree is additive on nonzero division-ring products. If q is a
unit of E, both q and its inverse have nonnegative degree and their
degrees sum to zero. Thus
\[
                                  E^\times=L^\times.     \tag{3}
\]
This requires only finite generation, torsion-freeness and rank one
of the underlying A-module; projectivity is used later for Morita
denominator clearing. Conjugating a fractional realization preserves
degree and the conclusion.

There are also no nontrivial algebraic L-constants anywhere in D.
If 0!=q in D is algebraic over L, a minimal polynomial and additivity
of degree exclude both positive and negative degree. In degree zero,
its rational leading symbol is algebraic over L in L(x,p), so it is
some c in L. Then q-c has negative degree unless zero, and is still
algebraic over L, a contradiction. Therefore
\[
                     q\in D,\ q\text{ algebraic over }L
                                           \ \Longrightarrow\ q\in L. \tag{4}
\]
This proof uses only the Ore symbol calculus already established in
[QHI1](HAMILTONIAN_IDEALIZER_CRITERION.md).

## The module dimension bound needed for clearing

For a finite A-module M, the Joseph--Rentschler endomorphism bound and
the Weyl second filter inequality give
\[
 \operatorname{Cdim}\operatorname{End}_A(M)
           \le\operatorname{Kdim}_A M
           \le\operatorname{GKdim}_A M-n.                 \tag{5}
\]
Here Cdim is the maximum number of commuting algebraically independent
elements; the second inequality uses filter dimension one for A.
The external inputs are [Bavula, Theorems1.2,1.4 and3.3](https://arxiv.org/html/math/0401263).
They concern the Weyl algebra A and its finite modules, not an assumed
finite-type associated graded algebra for E.

If a finite A-module M has generic rank zero, every one of its finitely
many generators is killed by a nonzero element of A. It is a quotient
of finitely many A/I_j with I_j nonzero. Their associated graded
quotients have growth dimension at most 2n-1. Hence
\[
 M\ne0,\quad\operatorname{rank}_A M=0
                  \ \Longrightarrow\
                  \operatorname{Cdim}\operatorname{End}_A(M)\le n-1.
                                                               \tag{6}
\]
Also Kdim_A P=n: P contains a copy of A and is a finite projective
summand of a finite free module. Equation (5)'s first inequality
therefore gives Cdim E<=n. In particular if B is commutative,
H=Frac B is algebraic over Q=Frac P_F: a transcendental b in B
would provide n+1 independent commuting elements in E.

## Clear a finite adjoint orbit with one spectral polynomial

Let q in D have finite-dimensional joint adjoint orbit
\[
              V=\operatorname{span}_L
                  \{\operatorname{ad}(F)^\alpha q:\alpha\ge0\}. \tag{7}
\]
There is a nonzero polynomial s in P_F with s*q in E.

To prove this, consider the left ideal
\[
                       J=\{a\in E:aV\subset E\}.          \tag{8}
\]
It is nonzero. One may use the classical quotient division ring of
the rank-one Morita order E, which is D. Alternatively, realize P as
a left ideal I in A. Then I is contained in its right idealizer E.
For each finite v in a basis of V, the left ideal
{a in A:a*v in I} is nonzero by Ore common denominators; intersect
these finitely many left ideals with I. Any nonzero element of the
intersection lies in J. Nonzero left ideals of the Ore domain A
have nonzero finite intersections.

Because V is adjoint-stable,
\[
                   (aF_i)v=a[F_i,v]+(av)F_i\in E
                    \quad(a\in J,\ v\in V),               \tag{9}
\]
so J*F_i is contained in J. Thus P_F lies in the right normalizer
N_E(J). If J=E, take s=1. Otherwise Morita equivalence identifies
End_E(E/J) with End_A(P/PJ). The module P/PJ has generic rank zero:
for 0!=a in J, right multiplication by a is injective on rank-one P
and P*a is contained in PJ. It is finite as an A-module.

If P_F intersect J were zero, the polynomial algebra P_F would
inject into N_E(J)/J=End_E(E/J)^op. That contradicts (6). Therefore
\[
                 0\ne s(F)\in J\cap P_F,\qquad s(F)q\in E. \tag{10}
\]
No degree bound for s or algorithmic stopping bound is asserted.
We clear only finite joint adjoint orbits; no Ore condition for S
inside the whole algebra E is assumed.

## Rational centralizers are spectrally saturated

A rational commutant has one-dimensional adjoint orbit. Its cleared
numerator in (10) still commutes with every F_i. Since S is central
in B, we obtain
\[
                   \operatorname{Cent}_{D}(F)=S^{-1}B.    \tag{11}
\]
If B is commutative, the left side is a field, so this is Frac B.
Consequently a rational commuting operator in D automatically has
the membership required by QTB1: one does not need to exhibit a
ratio of actual commutants separately. This conclusion requires
the complete n-tuple in the rank-n projective Weyl order; it is
not a statement about arbitrary partial commuting sets or algebras.

Changing a fractional realization by P->P*g conjugates E, F and q
together. Equations (3), (10) and (11) are invariant under that
actual rational frame change. No scalar free frame is introduced.

## Spectral torsion is exactly rational individual admission

Fix a nonzero scalar eta in L. For each i let alpha_i be the global
Koszul translation class represented by eta*e_i in the complex with
d0(T)=([F_j,T])_j. The polynomial algebra P_F acts by left
multiplication, commuting with all the adjoint differentials.
Its annihilator is
\[
 I_i=\{s\in L[t]:\ [F_j,U]=\eta\delta_{ij}s(F)
                                      \text{ for some }U\in E\}. \tag{12}
\]
Independence of the F_i ensures nonzero s gives nonzero s(F).

Then
\[
              I_i\ne0
       \quad\Longleftrightarrow\quad
       \exists T\in D:\ [F_j,T]=\eta\delta_{ij}\quad
                                      \text{for every }j. \tag{13}
\]
The forward implication uses T=s(F)^-1*U. For the converse, the
joint orbit of such T is contained in span{T,1}. Equation (10)
provides s(F)T in E, which witnesses s*alpha_i=0.
If several rational individual mates exist, the product of their
spectral denominators clears them all individually.

This is a fixed linear equation theorem. Torsion does not mean
alpha_i=0 and need not supply commuting individual mates. Full
rational flatness still has its curvature obstruction; finite full
admission needs the actual zero classes and actual curvature
primitive of QSC1.

For the rank-one coefficient-curve algebra D_C=K(C)[a;eta*delta],
choose a smooth affine open X on which the coefficients of S and
the nonzero vector field delta are regular and delta is nowhere
zero. Then S is in D(X), Frac D(X)=Frac D_C, and the published
scalar Lemma3.7 clears any rational mate by a polynomial in S.
For every positive-order S in D_C, its universal spectral class
therefore satisfies
\[
              [\eta]\in D_C/[S,D_C]\text{ is torsion over }K[u]
                       \quad\Longleftrightarrow\quad
                       [S,T]=\eta\text{ has a rational mate}. \tag{14}
\]
The vector field trivializes the tangent line bundle on that open;
all choices involve only finitely many removed zeros and poles.
The same denominator-ideal proof works over characteristic-zero
constant fields. This supplies the converse left unused in QST1.
Scalar spectral localization still differs from arbitrary localization
at moving coefficients of the universal (v,u)-module; no such
module-localization identification is asserted.

## QTB1 applies before a projective free frame

Assume B is commutative and actual individual mates T_i lie in E.
Their derivations preserve B, commute there, and extend the standard
coordinate partials to H. Algebraicity of H/Q follows from (5).
If c in H is horizontal, differentiate its monic minimal polynomial
over Q. Minimality makes every coefficient horizontal, hence in L.
Thus c is algebraic over L, and (4) gives c in L.

All QTB1 hypotheses are now automatic:
\[
                 H/Q\text{ algebraic},\quad
                 H^{\nabla}=L,\quad P_F\cap B^\times=L^\times.
                                                               \tag{15}
\]
This uses neither a scalar free frame nor finite generation of B.
Accordingly a nontrivial finite Galois factor of the full rational
centralizer (11) forbids simultaneous actual individual mates.
Birationality forces B=P_F and gives finite full completion in E.
For a regular Long corner, QET1 then supplies integral descent and
the prescribed classical frame; QMC1 supplies the separate matrix
consequence once its full finite-left graph conditions hold.

The radical version also extends. Under a finite constant-field
extension, P remains finite projective and rank one: extend an ideal
realization inside the scalar Weyl algebra. Finite projective Hom
commutes with base change, so E extends to the endomorphism algebra
of that extended module. The units, dimension, clearing and constant
arguments remain valid. Centralizer commutativity survives as the
kernel of the finite list of adjoint maps. Thus adjoining roots of
unity and descending as in QTB1 proves
\[
 q\in\operatorname{Cent}_D(F),\quad q^k\in Q,\quad k\ge2
                                    \quad\Longrightarrow\quad q\in Q
                                                               \tag{16}
\]
whenever actual individual mates exist. A new rational radical is
therefore an obstruction even when its representative is not in E.

For generic corners this is an endomorphism-algebra statement, not
a scalar A_n realization. A partial module lacking rank one,
projectivity or a justified division realization does not acquire
these conclusions. General matrix-valued positions on a higher-rank
left module remain outside the rank-one theorem.

## A nonfree control with torsion but no actual mate

Use QHI1's normalized charge-one module
M=A*u+A*(v+u^-1), [v,u]=1. Its idealizer contains u^2,u^3,H0=u*v
and L0=v^2-2*u^-2. Its graded algebra has no degree-one symbols,
and degree at most zero is scalar.

Take F=u^2. The rational commutant theta=u and rational mate
T=-(2u)^-1*v satisfy
\[
       F\theta=u^3\in E,\qquad FT=-H0/2\in E,\qquad
                              [F,-H0/2]=F.               \tag{17}
\]
Neither theta nor T is in E. The rational centralizer is L(u),
and the actual centralizer is the nonnormal cusp ring
\[
                              B=L[u^2,u^3].              \tag{18}
\]
Indeed an element of E commuting with u^2 is Laurent in u.
Subtract its polynomial terms of degrees at least two and its
constant, all in E. A remaining linear term is forbidden by
gr E; a remaining negative part is forbidden by nonnegative degree.
This proves (18) in every degree. The fraction field L(u) is
quadratic Galois over L(F), so the field obstruction excludes every
actual mate. QHI1 also gives that exclusion independently.

If xi is the spectral variable, the exact annihilator is
I_F=(xi): equation (17) kills xi*[1], whereas [1] is nonzero.
Polynomial division then excludes an annihilator with nonzero
constant term. Thus this genuine nonfree projective control
distinguishes torsion from zero, and rational from actual admission.

The radical has an entirely actual cleared certificate:
(u^3)^2=F*(u^2)^2. This is the same quadratic extension as theta^2=F,
with theta outside E but automatically in Frac B. No free frame
or algebraic root extraction is assumed.

## Rank and completion boundaries

The rank-one unit conclusion does not hold for arbitrary matrix
algebras: I+x*e12 in Mat_2(A) is a nonconstant unit with inverse
I-x*e12. A formal unit in a completed algebra likewise need not
be an actual unit of E. The coefficient-curve ring D(X) may also
have nonconstant units; (14) uses denominator clearing, not the
rank-one Weyl unit theorem.

Spectral clearing has no uniform denominator-degree bound here.
A finite prefix of annihilator tests proves no nontorsion theorem.
The result does not force a rational mate when the spectral class
is nontorsion, make a PDO mate rational, or identify general moving
coefficient poles with spectral localization. It closes the rational
mate/annihilator equivalence for the specified complete systems.

## Evidence and next construction gate

The standalone checker verifies retained QHI1/QPF1/QTB1 inputs,
charge-one multiplier memberships and exact denominator identities,
the cusp semigroup controls, rational versus actual symbols,
fractional conjugation transport, scalar spectral examples and
the higher-rank matrix-unit boundary. The universal Morita,
dimension, denominator clearing, constants, centralizer equality,
torsion equivalence and projective radical descent use written
proofs with the stated external dimension theorems. No independent
theorem replay, formal verification, external review or novelty
priority is asserted.

For a new rank-one Long carrier, inspect the full rational
centralizer in its Weyl division realization. Rational Galois-factor
or radical certificates now automatically descend to its actual
centralizer field, without finite generation or a scalar free frame.
For positive construction, rational individual mates certify spectral
torsion; actual vanishing and finite flat completion are still
required. No new finite Long pair, carrier, full lift or
all-representative nonlifting is supplied.
