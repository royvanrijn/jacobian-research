# Canonical boundary-incidence invariants

This note extracts the C04--C24 boundary comparison from its coordinates and
turns the entire ramification-labelled incidence geometry of the target
boundary into a stable polynomial left--right invariant.  No boundary
component has to be individually distinguishable.  When two components are
uniquely marked by the finite normalization, their intersection is a
canonical subobject of this stronger invariant.

Throughout, `k` is a characteristic-zero field and

\[
 F:X=\mathbb A_k^d\longrightarrow Y=\mathbb A_k^d
\]

is a dominant Keller map.  Put `K=k(Y)` and `L=k(X)`, and let

\[
 \pi_F:\overline X_F=\operatorname{Norm}_Y(L)\longrightarrow Y
                                                               \tag{1}
\]

be the finite normalization.  The factorization of `F` through (1) identifies
`X` with an open subscheme of `bar(X)_F`; write
`partial_F=bar(X)_F minus X`.  These facts and the invariance of the resulting
boundary valuations are proved in
[the resolvent--ramification note](RESOLVENT_RAMIFICATION_SIGNATURE.md).

## 1. The labelled boundary-incidence object

Let `D_F` be the set of prime divisors `Z` of `Y` dominated by a prime divisor
`E` of `partial_F`.  Attach to `Z` the arithmetic boundary profile

\[
 \Lambda_F(Z)=
 \left(\{(e(E/Z),f(E/Z)):E\subset\partial_F,\ E\mapsto Z\},
       \ell_Z(F)\right),                                      \tag{2}
\]

where the braces denote a multiset and

\[
 \ell_Z(F)=\sum_{E\mapsto Z}e(E/Z)f(E/Z).                     \tag{3}
\]

After base change to an algebraic closure, each arithmetic prime and profile
is replaced by its geometric constituents.  Intrinsic marking rules may use
either the arithmetic profile or this geometric refinement, provided the
selected arithmetic target divisor is unique.

One may refine (2) by the boundary valuation functionals on `k[X]`.  Nothing
below depends on using only the coarse profile displayed in (2).  The set
`D_F` is finite: its members are images of the finitely many divisorial
irreducible components of the closed boundary in the noetherian finite
normalization.

For every nonempty subset `S subset D_F`, define

\[
 I_S(F)=\left(\bigcap_{Z\in S}Z\right)_{\mathrm{red}}.          \tag{4}
\]

If `S subset T`, there is a canonical closed immersion
`I_T(F) -> I_S(F)`.  The **labelled boundary-incidence object**
`mathfrak I(F)` consists of:

- the finite set `D_F`, with each vertex `Z` labelled by `Lambda_F(Z)` and
  any chosen intrinsic valuation refinements;
- every reduced multiple intersection `I_S(F)`; and
- all the closed immersions induced by inclusions among the indexing sets.

This is an unlabelled finite diagram except for the intrinsic finite-cover
labels: components with identical profiles may be permuted.  Thus the whole
object is defined even when no individual divisor can be singled out.

### Marked subobjects

An **intrinsic marking rule** is an isomorphism-invariant predicate on this
finite-cover data which selects exactly one member of `D_F`.  Thus the rule is
not allowed to mention a coordinate equation, polynomial degree, or a chosen
resolvent.  Examples are “the unique divisor with one boundary prime of
ramification index two” and “the unique divisor all of whose boundary primes
are unramified and whose positive sheet loss is prescribed by the rest of the
signature.”

Fix two intrinsic marking rules `p` and `q` which select distinct divisors

\[
 Z_p(F),Z_q(F)\in D_F.                                        \tag{5}
\]

Define the **reduced marked boundary intersection**

\[
 I_{p,q}(F)=\bigl(Z_p(F)\cap Z_q(F)\bigr)_{\mathrm{red}}.       \tag{6}
\]

The data in (6) are target data.  The boundary primes upstairs are what make
the two target components in (5) canonical.

## 2. Invariance theorem

### Theorem 2.1 (stable boundary-incidence theorem)

Let `F,G:A^d_k -> A^d_k` be dominant Keller maps.  Then:

1. if `G=beta F alpha` for polynomial automorphisms `alpha,beta`, the target
   automorphism `beta` induces an isomorphism of labelled incidence diagrams
   \[
      \mathfrak I(F)\xrightarrow{\sim}\mathfrak I(G);          \tag{7}
   \]
2. for every `a>=0`,
   \[
      \mathfrak I(F\times\operatorname{id}_{\mathbb A^a})
      =\mathfrak I(F)\times\mathbb A^a,                        \tag{8}
   \]
   where the product on the right is taken at every scheme in the diagram;
3. consequently the affine-cylinder class of `mathfrak I(F)` is a polynomial
   stable left--right invariant of `F`.

**Proof.**  A source automorphism identifies the extensions `K subset L`,
their finite normalizations, the affine-source opens, and all boundary primes
over each target divisor.  It therefore identifies `D_F`, its labels, and
every intersection and incidence arrow.

A target automorphism transports `K subset L` and, by functoriality of
integral closure, transports the finite normalization and all data in
(2)--(3).  It bijects `D_F` with `D_G`, preserving all labels.  For every
nonempty `S subset D_F`, scheme-theoretic intersection commutes with `beta`,
as does passage to the reduction, and hence

\[
 \beta(I_S(F))=I_{\beta(S)}(G).                               \tag{9}
\]

These isomorphisms commute with all the closed immersions in the diagram,
which proves (7).

After adjoining identity coordinates the field extension is the purely
transcendental base change

\[
 L(t_1,\ldots,t_a)/K(t_1,\ldots,t_a).
\]

Its finite normalization, affine-source open, boundary primes, target
divisors, ramification indices, residue degrees, and sheet losses are the
products of the old objects with `A^a`.  Flat base change by the geometrically
reduced scheme `A^a` commutes with the scheme-theoretic intersections in (4)
and with reduction.  Every object and arrow is therefore the corresponding
old object or arrow times `A^a`.  This proves (8), and the two assertions give
the stable conclusion.  QED

### Corollary 2.2 (canonical marked intersections)

Suppose intrinsic rules `p_1,...,p_s` uniquely select distinct members of
`D_F`.  Then

\[
 \left(Z_{p_1}(F)\cap\cdots\cap Z_{p_s}(F)\right)_{\mathrm{red}}
                                                                    \tag{10}
\]

has a canonical stable isomorphism class.  In particular, (6) is a polynomial
stable left--right invariant whenever `p` and `q` exist.

**Proof.**  The marking rules canonically select one entry of the labelled
diagram and are respected by every diagram isomorphism in Theorem 2.1.
Equation (8) gives its behavior under stabilization.  QED

### Corollary 2.3 (usable stable obstructions)

Any invariant of reduced affine schemes which is unchanged by product with
affine space gives a stable left--right obstruction when applied to (10).  In
particular this holds for:

- connectedness, equivalently the existence of a nontrivial idempotent in
  the coordinate ring;
- the number of connected components after geometric base change;
- irreducibility and the number of irreducible components; and
- the unit group modulo constants.

**Proof.**  For a reduced affine `k`-algebra `A`, the idempotents and units of
`A[t_1,...,t_a]` are those of `A`; its minimal primes are the extensions of
the minimal primes of `A`.  The same statements hold after base change to an
algebraic closure.  Apply Corollary 2.2.  QED

The qualifier “intrinsic” is essential.  Intersecting two coordinate
hyperplanes seen in a particular inverse resolvent does not define an
invariant until their generic points have been identified with uniquely
marked target divisors in the canonical finite normalization.

## 3. C04 versus C24

Let `m>1`.  For `C24_(m,1)` and a generic weighted C04 seed of inverse degree
`m+2`, the intrinsic signatures uniquely mark:

- the target divisor with one boundary prime of `(e,f)=(2,1)` and sheet loss
  two; and
- the other target divisor, whose boundary primes are unramified and whose
  total sheet loss is `m-1`.

The uniform calculations in
[the boundary-obstruction note](BOUNDARY_INTERSECTION_OBSTRUCTION.md) give

\[
 I_{p,q}(C24_{m,1})\simeq\mathbb A^1\sqcup\mathbb G_m,
 \qquad
 I_{p,q}(C04_{m+2})\simeq\mathbb A^1.                         \tag{11}
\]

The first coordinate ring has a nontrivial idempotent and the second does
not.  Corollary 2.3 therefore yields the following coordinate-free form of
the former comparison theorem.

### Corollary 3.1

For every `m>1`, no `C24_(m,1)` map, including an arbitrary allowed tail, is
polynomially left--right equivalent, even after adjoining identity
coordinates, to a generic weighted C04 map of inverse degree `m+2`.

This is an application of Theorem 2.1, not part of its hypotheses.  The full
labelled incidence object is available for every nonproper
characteristic-zero Keller map; unique marking is needed only to name a
particular intersection without retaining the rest of the diagram.
