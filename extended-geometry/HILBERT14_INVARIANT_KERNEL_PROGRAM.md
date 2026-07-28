# Hilbert--14 invariant-kernel program

## Outcome of the first experiment

The normalized quadratic--cubic factorization slice gives an exact positive
control, not a counterexample.  If

\[
 R=k[X_{2,3}]
\]

and \(D_{10},D_7\) are the commuting primitive locally nilpotent
derivations constructed in
[the factorization-slice audit](QUADRATIC_CUBIC_FACTORIZATION_SLICE.md#exact-kernels-and-the-hilbert--14-test),
then

\[
 \ker D_{10}\cap\ker D_7=k[K,H,V]\simeq k^{[3]},
\]

\[
 \ker D_7=k[K,H,V,s]\simeq k^{[4]},\qquad
 \ker D_{22}=k[K,H,V,W]\simeq k^{[4]}.
\]

The remaining displayed action \(D_{10}\) has a global slice, so its kernel
is finitely generated too.  Thus saturating the most natural source/factor
gauge shears on this fivefold terminates after the two regularized
denominator classes

\[
 G={J-4U\over a},\qquad
 H={K^2+8J-16U\over a^2}.
\]

This is useful negative information: the first place to seek a
non-finitely-generated Keller-adjacent invariant ring is not the primitive
\((2,3)\) action itself.

The next \((2,4)\) experiment now terminates as well.  For its natural
commuting primitive pair,

\[
 \ker D_{10}\cap\ker D_0
 =k[a,p,U,N,M,C,Q,S].
\]

The generic quotient \(M^2+4UN^2=256a^4\) initially collapses two boundary
coordinates, but the three regular fractions \(C,Q,S\) restore them.  The
last relation is linear in \(S\) modulo \(a\), which prevents an infinite
Hensel ladder.  See the
[quadratic--quartic audit](QUADRATIC_QUARTIC_HILBERT14_SLICE.md).
The noncommuting third Euclidean action closes as well:

\[
 \ker D_{10}\cap\ker D_0\cap\ker D_2
 =k[a,U,N,M,C],
\]

a threefold whose \(a=0\) quotient is the cusp
\(M^2+64U^3=0\).

The multiboundary comparison is now exact too.  On
\[
 k[s^2,s^3,t^2,t^3,X,Y,U,V]
\]
the two commuting cusp LNDs have a non-finitely-generated common kernel:
the classes \(s^2t^2(X+sY)^m(U+tV)^n\) escape every finite bidegree
rectangle modulo \((s^4,t^4)\).  The
[multiboundary control](MULTIBOUNDARY_HILBERT14_CONTROL.md) also proves that
the apparent leading pair in every tangent-normalized factorization slice is
disjoint, since \(1=ad+bp\in(a,p)\).

## 1. The saturation ledger

Let \(R\) be a finitely generated domain, let \(\pi\in R\) define an
irreducible boundary divisor, and suppose a localized action has an
explicit invariant algebra \(B_\pi\subset R_\pi\).  The global invariant
candidate is the intersection

\[
 \mathcal S_\pi(B,R)=R\cap B_\pi.
\]

Record its pole filtration

\[
 C_n=\{b\in B:\ b/\pi^n\in R\},\qquad
 \mathcal S_\pi(B,R)=\sum_{n\geq0} C_n\pi^{-n}.       \tag{1.1}
\]

This gives a reproducible fork.

1. **Terminating certificate.**  Exhibit finitely many regularized
   fractions, prove that they generate after localization, and prove that
   reduction modulo \(\pi\) is injective.  A minimal-pole argument then
   eliminates every further denominator.
2. **Non-terminating certificate.**  Produce \(F_n\in C_n\pi^{-n}\) of
   unbounded auxiliary degree and a quotient in which products of positive
   filtration vanish.  Then a finite set of generators has bounded degree
   but some \(F_n\) escapes.

The \((2,3)\) computation realizes the first branch.  Equations (59)--(65)
of the slice audit give both the localized generators and the boundary
injectivity certificate.

## 2. Exact non-terminating control

Maubach's example supplies the comparison that a search script must be able
to distinguish.  Put

\[
 A=k[T^2,T^3],\qquad
 S=A[X,Y,Z]\big/
 \left(Z^2-T^4(T^2X+T^3Y)^2-1\right)
\]

and

\[
 D=T^3{\partial\over\partial X}
   -T^2{\partial\over\partial Y},\qquad P=X+TY.
                                                               \tag{2.1}
\]

The derivation is locally nilpotent and preserves the relation.  In the
normalization obtained by adjoining \(T\),

\[
 \ker D=k[T,Z,P].
\]

Its intersection with \(S\) is not finitely generated.  The explicit ladder

\[
 F_n=T^2P^n\in S\cap k[T,Z,P]\qquad(n\geq0)           \tag{2.2}
\]

is regular because every coefficient \(T^{2+j}\) in its binomial expansion
belongs to \(k[T^2,T^3]\).

Here is the uniform obstruction.  Every invariant of positive \(P\)-degree
lies in the conductor ideal \((T^2,T^3)S\).  Modulo \(T^4\), products of two
such invariants vanish.  If finitely many generators had maximum
\((X,Y)\)-degree \(d\), then all surviving linear terms would have degree at
most \(d\), whereas \(F_{d+1}\) has the surviving term
\(T^2X^{d+1}\).  This is impossible.

The local checker
[`verify_hilbert14_saturation_ladders.py`](../scripts/verify_hilbert14_saturation_ladders.py)
replays (2.1), the cusp-semigroup membership of (2.2), and the
modulo-\(T^4\) escape through a configurable finite degree.  The arbitrary
\(d\) argument above is the proof; the bounded replay is only a regression
control.  The construction and proof are from Maubach,
[Infinitely generated Derksen and ML invariant](https://www.math.ru.nl/~maubach/Papers/InfGenDerInv_v12.pdf),
Proposition 2.3.

## 3. Three Keller-derived sources of actions

### 3.1 Gauge shears

Factor operations and Tschirnhausen changes give honest unipotent actions
on coefficient rings.  They are the safest starting point because local
nilpotence is triangular before imposing the equations.  The current
inventory contains:

- the primitive \(D_{10},D_7,D_{22}\) actions on \(X_{2,3}\);
- the vertical factor shear with polynomial slice in the
  [normalized linear--quadratic model](../verified/NORMALIZED_FACTORIZATION_MODEL.md);
- upper and lower binary-variable shears in
  [cubic gauge straightening](../cancellation/CUBIC_GAUGE_STRAIGHTENING.md).

For every new factorization slice, the first experiment should form all
commuting primitive gauge derivations, compute the generic joint quotient,
and run the pole filtration (1.1) at every omitted leading-coefficient
divisor.

### 3.2 Hamiltonian lifts

For a plane Keller pair \((P,Q)\) with \(\{P,Q\}=1\), the Hamiltonian field

\[
 \delta_P=\{P,-\}
\]

satisfies \(\delta_P(P)=0\) and \(\delta_P(Q)=1\), up to the bracket sign
convention.  It is a derivation and has a slice on \(k[P,Q]\), but it must
**not** be called locally nilpotent on \(k[x,y]\) without proof.  In fact,
local nilpotence there would force a polynomial \(\mathbb G_a\)-action with
a global slice and would essentially solve the coordinate problem for that
pair.  This makes Hamiltonian local nilpotence a sharp screen rather than a
free source of examples.

The computational task is instead to extend \(\partial/\partial Q\) from
\(k[P,Q]\) to the inverse-equation field, compute its pole orders along the
normal boundary, and test whether multiplication by a kernel element makes
the lift both regular and locally nilpotent.  A persistent pole is useful
obstruction data even when no action results.

### 3.3 Translation symmetries of inverse equations

For a normalized inverse equation

\[
 g(T;P,Q)=0,
\]

target translation in \(Q\) lifts formally by

\[
 \partial_Q(T)=-{\partial_Qg\over\partial_Tg}.        \tag{3.1}
\]

Equation (3.1) exposes the discriminant denominator explicitly.  Normalize
the inverse-equation algebra, evaluate its valuation at every missing
boundary prime, and saturate the numerator derivation by the smallest
invariant boundary factors.  Only after checking preservation of the
integral algebra and local nilpotence does this become a
\(\mathbb G_a\)-action.  The coefficient ideals \(C_n\) in (1.1) are the
right data to retain even when it does not.

## 4. Dimension and boundary screens

Two screens prevent low-yield searches.

- A single LND kernel on a normal affine threefold is finitely generated by
  the standard transcendence-degree-two finiteness theorem.  A
  non-finite kernel experiment therefore needs a sufficiently large
  ambient ring, a nonnormal base such as the cusp control, or an
  intersection/Derksen construction.
- A global slice gives \(R=(\ker D)[s]\), hence finite generation of \(R\)
  immediately gives finite generation of \(\ker D\).  Slice actions remain
  valuable for intersections, but not as single-kernel counterexamples.

For boundary compactifications there is already a genuine non-finite
algebra in the surrounding apparatus: when there are at least eight
labelled finite branch values, the universal receiver is
\(\overline M_{0,b+2}\) with \(b+2\geq10\), whose Cox ring is not finitely
generated.  The precise scope and the reason this does not automatically
descend to the Keller graph are recorded in
[the wonderful-pullback audit](BRANCH_GRAPH_WONDERFUL_PULLBACK.md#hilbert--14-boundary-consequence).

## 5. Next Keller-attached experiment

The multiboundary algebraic control is complete, but the obvious
factorization leading divisors cannot meet.  The next construction must
therefore leave the tangent-normalized family.  The smallest concrete
options are:

1. add a second reconstruction variable to the direct two-boundary
   suspension;
2. retain the third divisor forced by its one-variable Jacobian and run a
   three-index conductor ledger; or
3. lift gauge/translation derivations to one of the completed
   two-boundary charts already produced by the plane log-boundary compiler.

For each option, reduce the putative invariant algebra modulo the squares of
the conductor ideals.  A finite rectangle bound proves termination; a
bidegree escape as in the two-cusp control proves non-finite generation.
A bounded Gröbner search without one of those certificates remains an
experiment, not a Hilbert--14 result.
