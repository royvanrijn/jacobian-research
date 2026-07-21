# Canonical boundary-incidence invariants

This note extracts the C04--C24 boundary comparison from its coordinates and
turns the incidence geometry of the target boundary into a stable polynomial
left--right invariant.  It also records the exact hypotheses: quasi-finiteness
is essential for the canonical open immersion, while reduction of the target
intersections is optional rather than necessary.

The target diagrams defined here are the target quotients of the reduced and
scheme-theoretic levels in
[the boundary-invariant ladder](BOUNDARY_INVARIANT_LADDER.md), whose formal
level adds completed neighborhoods, finite-stratum conductors, differents,
and boundary valuation filtrations.

## 1. Correct general setup

Let `k` be a field.  Let `X=Spec(A)` and `Y=Spec(B)` be normal integral affine
`k`-varieties, and let

\[
 F:X\longrightarrow Y
\]

be dominant and quasi-finite.  Put `K=k(Y)` and `L=k(X)`.  Then `L/K` is
finite.  Since `Y` is excellent, its normalization in `L` is finite:

\[
 \pi_F:\overline X_F=\operatorname{Norm}_Y(L)\longrightarrow Y. \tag{1}
\]

These hypotheses will be called the **admissible boundary setup**.  They hold
for every polynomial Keller map between affine spaces: a nonzero constant
Jacobian makes the morphism étale, hence quasi-finite, and affine spaces are
normal and excellent.  The invariant theorem itself does not use
characteristic zero once the admissible setup is assumed.

Logically, the incidence theorem only needs the normalization (1) to be
finite and the induced factor `j_F` to be an open immersion, with a normal
target in codimension one so that the displayed DVR labels are defined.  The
normal affine quasi-finite finite-type package above is a transparent,
checkable sufficient hypothesis set.  Excellence may be replaced by the
explicit finiteness of normalization; it is not otherwise used.

### Proposition 1.1 (canonical affine-source open)

In the admissible boundary setup, `F` factors uniquely through an open
immersion

\[
 j_F:X\hookrightarrow\overline X_F.                            \tag{2}
\]

**Proof.**  Let `bar(B)` be the integral closure of `B` in `L`.  If
`u in bar(B)`, then `u` is integral over `B` and hence over `A`, because
`B subset A`.  Normality of `A` gives `u in A`.  Thus
`B subset bar(B) subset A`, which defines the factorization (2).

The map `j_F` is of finite type and birational.  It is quasi-finite because
each of its fibers is contained in a fiber of the quasi-finite map `F`.
Zariski's Main Theorem now says that a quasi-finite separated birational map
to the normal scheme `bar(X)_F` is an open immersion.  QED

Write `partial_F` for the reduced closed complement of (2).  This pair is
canonical.  No projective compactification, blow-up sequence, or primitive
resolvent is part of its definition.

### Proposition 1.2 (purity and noncontraction of the boundary)

In an admissible boundary setup:

1. every irreducible component of `partial_F` has codimension one in
   `bar(X)_F`; and
2. every such component maps onto a prime divisor of `Y`.

In particular, no boundary divisor of the canonical finite normalization can
have image of codimension greater than one.

**Proof.**  We first use a general affine-open purity argument.  Let
`V=Spec(C)` be a normal noetherian integral affine scheme and let `U subset V`
be a dense affine open.  Suppose that an irreducible component of `V-U`, with
generic prime `p`, has height at least two.  After replacing `V` by a
principal neighborhood `D(g)` of `p`, chosen to avoid the other finitely many
irreducible components of `V-U`, the complement of `U intersect D(g)` in
`D(g)` has codimension at least two.  The intersection `U intersect D(g)` is
affine because it is a principal open of the affine scheme `U`.  Normality
and the codimension-two extension theorem give

\[
 \Gamma(U\cap D(g),\mathcal O)=\Gamma(D(g),\mathcal O)=C_g.
\]

The open immersion between these two affine schemes is therefore an
isomorphism, contradicting the presence of `p`.  Thus `V-U` is pure of
codimension one.  Apply this with `V=bar(X)_F` and `U=X`; the former is affine
because it is finite over affine `Y`, and both schemes are normal.

For the image assertion, write `Y=Spec(B)`, `bar(X)_F=Spec(bar(B))`, and let
`q subset bar(B)` be the height-one prime of a boundary component.  Put
`p=q intersect B`.  The prime `p` is nonzero: after localizing at
`B-{0}`, the finite domain `bar(B) tensor_B k(Y)` is the field `k(X)`, whose
only prime is zero.  Because `B` is normal, going down holds for the integral
domain extension `B subset bar(B)`.  A chain of primes below `p` therefore
lifts to a chain below `q`, so `ht(p)<=ht(q)=1`.  Hence `ht(p)=1`.  Since the
finite map is closed, the image is exactly `V(p)`, a prime divisor of `Y`.
QED

### Divisorial boundary images and labels

Let `D_F` be the finite set of prime divisors `Z` of `Y` which are images of
prime divisors `E` contained in `partial_F`.  Finiteness follows because the
divisorial `E` are codimension-one irreducible components of a closed subset
of the noetherian scheme `bar(X)_F`.  Proposition 1.2 shows both that these
are all irreducible components of the boundary and that every one has a
prime-divisor image.

At the generic points of `E` and `Z`, normality gives an extension of DVRs.
Define

\[
 \Lambda_F(Z)=
 \left(\{(e(E/Z),f(E/Z)):E\subset\partial_F,\ E\mapsto Z\},
       \ell_Z^\partial(F)\right),                              \tag{3}
\]

where the braces are a multiset,

\[
 f(E/Z)=[\kappa(E):\kappa(Z)],\qquad
 \ell_Z^\partial(F)=\sum_{E\subset\partial_F,\ E\mapsto Z}
 e(E/Z)f(E/Z).                                                \tag{4}
\]

The superscript in (4) matters: this is only the boundary contribution to the
degree over `Z`.  Primes over `Z` whose generic points lie in the affine open
`X` are not included.  Even over an algebraically closed field, the residue
degree `f(E/Z)` can exceed one; it is the generic degree of `E -> Z`, not
merely the size of a constant-field orbit.

The construction is arithmetic over the fixed field `k`.  If geometric data
are desired, repeat it after base change and normalization.  One must not
replace residue degrees by orbit sizes without separately proving that the
geometric maps have degree one.

Optional refinements, such as boundary valuation functionals on `k[X]`, may
be added only when they are specified functorially under source and target
isomorphisms.  Coordinate pole vectors and raw resolvent factors are not
eligible labels by themselves.

## 2. Scheme-theoretic and reduced incidence diagrams

Each `Z in D_F` is given its canonical reduced prime-divisor structure.  For
every nonempty subset `S subset D_F`, define the scheme-theoretic intersection

\[
 J_S(F)=\bigcap_{Z\in S}^{\mathrm{sch}} Z\subset Y             \tag{5}
\]

by the sum of the prime ideals of the members of `S`, and define

\[
 I_S(F)=J_S(F)_{\mathrm{red}}.                                \tag{6}
\]

For `S subset T`, there are canonical closed immersions
`J_T(F) -> J_S(F)` and `I_T(F) -> I_S(F)`.  Let `mathfrak J(F)` be the finite
diagram consisting of `D_F`, the labels (3), all `J_S(F)`, and all incidence
arrows.  Let `mathfrak I(F)` be its termwise reduction.

Thus `mathfrak J(F)` retains tangency and nilpotent intersection structure;
`mathfrak I(F)` retains only reduced incidence geometry.  Components with
identical intrinsic labels may be permuted, so neither diagram requires a
component to be individually distinguishable.

### Marked subobjects

An **intrinsic marking rule** is a predicate on the entire functorial finite-
cover data which selects exactly one member of `D_F` and is preserved by
isomorphisms of admissible boundary setups.  It may not mention a coordinate
equation, ordinary polynomial degree, or a chosen resolvent.

If rules `p_1,...,p_s` select distinct divisors, write

\[
 J_{p_1,...,p_s}(F)=
 Z_{p_1}(F)\cap^{\mathrm{sch}}\cdots\cap^{\mathrm{sch}}Z_{p_s}(F),
 \qquad
 I_{p_1,...,p_s}(F)=J_{p_1,...,p_s}(F)_{\mathrm{red}}.         \tag{7}
\]

## 3. Corrected invariance theorem

### Theorem 3.1 (stable boundary-incidence theorem)

Let `F:X -> Y` and `G:X' -> Y'` be admissible boundary setups.

1. If `G=beta F alpha` for isomorphisms `alpha:X' -> X` and
   `beta:Y -> Y'`, then `beta` induces isomorphisms of labelled diagrams
   \[
      \mathfrak J(F)\simeq\mathfrak J(G),\qquad
      \mathfrak I(F)\simeq\mathfrak I(G).                     \tag{8}
   \]
2. For every `a>=0`, stabilization gives termwise equalities
   \[
      \mathfrak J(F\times\operatorname{id}_{\mathbb A^a})
      =\mathfrak J(F)\times\mathbb A^a,
      \qquad
      \mathfrak I(F\times\operatorname{id}_{\mathbb A^a})
      =\mathfrak I(F)\times\mathbb A^a.                       \tag{9}
   \]
3. Consequently the affine-cylinder classes of both diagrams are stable
   left--right invariants.  This includes stable equivalences in which
   different numbers of identity coordinates are added to maps of different
   original dimensions, provided the stabilized source and target dimensions
   agree.

**Proof.**  A source isomorphism identifies `L/K`, its integral closure, and
the distinguished open copy of the source.  A target isomorphism transports
the base field inclusion and its integral closure.  Functoriality and
uniqueness of integral closure therefore identify the finite normal models,
their boundary primes, the DVR extensions, and the labels (3)--(4).  The
target isomorphism carries prime ideals to prime ideals, commutes with their
sums, and commutes with radicals.  Hence it transports both the full
scheme-theoretic intersections (5) and their reductions (6), compatibly with
all incidence arrows.  This proves (8).

For stabilization, the standard polynomial-extension theorem for integral
closure gives

\[
 \operatorname{Norm}_{Y\times\mathbb A^a}
       (L(t_1,\ldots,t_a))
 =\overline X_F\times\mathbb A^a.                             \tag{10}
\]

The distinguished open is `X times A^a`, so every divisorial boundary image
is exactly `E times A^a` over `Z times A^a`; no new divisorial component is
created from a higher-codimension part of the boundary.  The corresponding
DVR is the Gauss extension of the old valuation.  Its ramification index is
unchanged, and

\[
 [\kappa(E)(t_1,\ldots,t_a):\kappa(Z)(t_1,\ldots,t_a)]
 =[\kappa(E):\kappa(Z)],                                     \tag{11}
\]

so the residue degree and boundary contribution are unchanged.

If the prime ideal of `Z` is `p_Z`, then the stabilized ideal is
`p_Z k[Y,t_1,...,t_a]`.  Sums of these extended ideals are the extensions of
their sums, proving the first equality in (9).  Moreover, for every ideal
`I subset k[Y]`,

\[
 \sqrt{I k[Y,t_1,\ldots,t_a]}
 =\sqrt I\,k[Y,t_1,\ldots,t_a],                               \tag{12}
\]

because a polynomial ring over a reduced ring is reduced.  This proves the
reduced equality and hence (9).  Combining (8)--(9) proves the stable claim.
QED

### Corollary 3.2 (canonical marked intersections)

Every marked scheme-theoretic intersection `J_(p_1,...,p_s)(F)` and every
marked reduced intersection `I_(p_1,...,p_s)(F)` has a canonical stable
isomorphism class.  Any invariant unchanged by affine cylinders gives a
stable left--right obstruction.  For reduced affine intersections this
includes connectedness, geometric connected-component count,
irreducibility, irreducible-component count, and the unit group modulo
constants.  For the full scheme-theoretic intersections it also includes
reducedness and, in the noetherian case, the nilpotency index of the
nilradical.

**Proof.**  The marking rules select entries functorially from the diagrams in
Theorem 3.1.  For a reduced ring `C`, the idempotents and units of
`C[t_1,...,t_a]` are those of `C`, and its minimal primes are the extensions
of the minimal primes of `C`.  After geometric base change, apply the same
argument to the geometric reduction.  For any ring `R`,
`Nil(R[t_1,...,t_a])=Nil(R)[t_1,...,t_a]`; if the nilradical is nilpotent,
its least nilpotency exponent is therefore unchanged.  QED

## 4. Adversarial hypothesis audit

The following weakenings or shortcuts are invalid.

1. **Generically finite is not enough.**  The map
   `(x,y) -> (x,xy)` on `A^2` is generically birational but has the line
   `x=0` as a fiber over the origin.  Its target normalization is `A^2`, and
   the map to that normalization is not an open immersion.  Quasi-finiteness,
   automatic for Keller maps, is the missing hypothesis.
2. **Normality of the source is used.**  It is what places the integral
   closure of `B` inside `A` and hence produces the factorization (2).
3. **Finiteness of normalization is used.**  Excellence (or at least the
   Nagata property needed here) prevents the normalization from being
   nonfinite.  Finite-type schemes over a field satisfy this automatically.
4. **Finiteness prevents divisorial contraction.**  Proposition 1.2 is not a
   generic property of compactifications: for example, the exceptional
   divisor of a blow-up may map to a codimension-two center.  It holds here
   because the canonical normalization map is finite.  Moreover, affineness
   and normality make the reduced boundary pure of codimension one.  If one
   drops those hypotheses, a theory retaining all boundary components and
   their images is required.
5. **A displayed resolvent is not an exhaustion proof.**  To mark a target
   divisor intrinsically, every divisorial boundary prime in the canonical
   normalization must be accounted for and the proposed profile must be
   unique among all of them.
6. **Raw discriminants have the wrong scheme structure in general.**  They
   may contain powers from ramification, leading-coefficient degeneration, or
   a nonprimitive generator.  Vertices of `mathfrak J(F)` are reduced prime
   divisors.  Their scheme-theoretic intersection uses the sum of those prime
   ideals, not an unsaturated resultant.
7. **Reduction loses information but does not create invariance.**  The full
   intersections `J_S(F)` are already canonical and stable.  Passing to
   `I_S(F)` discards tangency multiplicities and embedded nilpotents; it is
   justified when only support, connectedness, or idempotents have actually
   been computed.  For example, the reduced divisors `y=0` and `y-x^2=0` in
   `A^2` meet scheme-theoretically in `k[x]/(x^2)`, while their reduced
   intersection is one point.  Stabilization preserves the former as
   `k[x,t_1,...,t_a]/(x^2)`; it does not erase the nilpotent.
8. **Geometric residue degree is not automatically one.**  Algebraic closure
   of the constant field does not force a generically finite map
   `E -> Z` to be birational.

## 5. C04 versus C24

Let `m>1`.  The C04--C24 application needs two family-specific inputs beyond
Theorem 3.1:

1. the normalization calculations exhaust all divisorial boundary images;
2. among that exhaustive list, the profiles uniquely mark the divisor with
   one boundary prime of `(e,f)=(2,1)` and boundary contribution two, and the
   divisor whose boundary primes are unramified with total contribution
   `m-1`.

The normalization and two-chart calculations in the resolvent--ramification
note supply the branches, while the
[boundary-exhaustion certificate](BOUNDARY_EXHAUSTION_CERTIFICATE.md),
Theorem 5.1, uses support containment and exact `sum ef` degree accounting to
prove that these are all components of the canonical normalization boundary
and all of their target images.  The uniform calculation in
[the boundary-obstruction note](BOUNDARY_INTERSECTION_OBSTRUCTION.md) gives
the **reduced** marked intersections

\[
 I_{p,q}(C24_{m,1})\simeq\mathbb A^1\sqcup\mathbb G_m,
 \qquad
 I_{p,q}(C04_{m+2})\simeq\mathbb A^1.                         \tag{13}
\]

The first coordinate ring has a nontrivial idempotent and the second does
not.  Corollary 3.2 therefore proves that no `C24_(m,1)` map, including an
arbitrary allowed tail, is polynomially left--right equivalent, even after
stabilization, to a generic weighted C04 map of inverse degree `m+2`.

The strengthened calculation in the boundary-obstruction note also
identifies the full scheme-theoretic intersections.  C24 has a nilradical of
index `m(m+1)`, while the weighted intersection is reduced.  The older
idempotent obstruction uses only (13), but the nilpotent-sensitive obstruction
is independently stable by Corollary 3.2.

The target diagrams `mathfrak J(F)` and `mathfrak I(F)` still forget how the
boundary divisors and their intersections sit in the finite normalization.
Their canonical bipartite enhancement, including upstairs intersections,
finite residue covers, scheme images, and the relative different, is defined
in
[UPSTAIRS_DOWNSTAIRS_BOUNDARY_INCIDENCE.md](UPSTAIRS_DOWNSTAIRS_BOUNDARY_INCIDENCE.md).
