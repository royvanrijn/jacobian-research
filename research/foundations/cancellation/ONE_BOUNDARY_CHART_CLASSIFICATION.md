# One-boundary chart extraction and the missing monotonicity theorem

This note attacks the upstream chart problem for controlled-boundary Keller
suspensions.  It does not revisit the cancellation reconstruction after a
reciprocal link has been supplied.  Its conclusions are:

1. there is a precise, intrinsic height-one criterion which produces the
   common principal open and hence the sign `+1` or `-1` of a one-boundary
   chart;
2. in the negative case that criterion already produces the reciprocal
   equation with exponent exactly one;
3. factoriality plus an exact relative Picard-rank count eliminates hidden
   exceptional divisors on a graph model;
4. effectivity of a polynomial target chart eliminates a hidden target
   ledger once reciprocal cancellation has made the source ledger zero;
5. divisor minimality, smooth boundary, factorial graph completion, and a
   bounded relative Picard rank do **not** by themselves force the proposed
   weighted/reciprocal dichotomy;
6. the extra missing property is boundary monotonicity in the positive
   branch and residue-lattice saturation in the tangential directions.

The countermodels below are elementary.  They show that the phrase
"divisor-minimal" must not be used as a substitute for these properties.
They also isolate a corrected conditional chart theorem which is strong
enough to feed the marked plane-core theorem and the already-complete
reciprocal reconstruction theorem.

Work over a characteristic-zero field `k`.  All rings are integral,
finitely generated `k`-algebras.

## 1. Divisorial pole support

Let `R` and `S` be normal domains with a fixed identification

\[
 \operatorname{Frac}R=K=\operatorname{Frac}S.          \tag{1.1}
\]

For a height-one prime `P` of `R`, write `v_P` for its normalized
valuation.  Define the forward divisorial pole support by

\[
 \operatorname{Pole}_R(S)
 =\{P\in(\operatorname{Spec}R)^{(1)}:
       v_P(f)<0\text{ for some }f\in S\}.              \tag{1.2}
\]

It is enough to test a finite set of algebra generators of `S`.  Define
`Pole_S(R)` symmetrically.  These sets depend only on the two affine models
inside `K`, not on a chosen rational formula or graph resolution.

If `a in R` and `d in S` are prime elements, call `(R,a;S,d)` a
**height-one saturated one-boundary link** when

\[
 S\subset R[a^{-1}],\qquad R\subset S[d^{-1}],         \tag{1.3}
\]

and both boundary equations become units in the opposite localization:

\[
 d\in R[a^{-1}]^*,\qquad a\in S[d^{-1}]^*.             \tag{1.4}
\]

Equivalently, all poles of the two affine coordinate rings and all zeros
and poles of the opposite boundary equation are supported on the declared
height-one boundary.  Conditions (1.3)--(1.4) are the useful precise
replacement for "one-boundary" at the affine height-one level.

### Proposition 1.1 -- intrinsic principal-open extraction

If `(R,a;S,d)` is a height-one saturated one-boundary link, then the field
identification induces an isomorphism

\[
 \boxed{S[d^{-1}]\simeq R[a^{-1}].}                    \tag{1.5}
\]

If, in addition, `R` and `S` are factorial and
`R^*=S^*=k^*`, then

\[
 \boxed{d=c a^\epsilon,\qquad
        c\in k^*,\quad\epsilon\in\{1,-1\}.}           \tag{1.6}
\]

Thus a negative link is intrinsically reciprocal, with exponent exactly
one.  No choice of compactification or local parameter is involved.

#### Proof

The first inclusion in (1.3), together with the first assertion of (1.4),
extends to a homomorphism

\[
 S[d^{-1}]\longrightarrow R[a^{-1}].                  \tag{1.7}
\]

The other two assertions give the reverse homomorphism.  Both are induced
by the identity of `K`, so they are inverse.

In a factorial domain with only scalar units,

\[
 R[a^{-1}]^*/k^*=\mathbb Z[a],qquad
 S[d^{-1}]^*/k^*=\mathbb Z[d].                         \tag{1.8}
\]

The isomorphism (1.5) induces an automorphism of these rank-one lattices.
Every automorphism of `Z` is multiplication by `+1` or `-1`, which proves
(1.6).  QED

Normality makes the criterion directly valuation-theoretic.  Indeed, the
Krull intersection formula gives

\[
 R[a^{-1}]
 =\bigcap_{P\in(\operatorname{Spec}R)^{(1)},\ P\ne(a)}R_P
 \subset K.                                            \tag{1.9}
\]

Consequently `S\subset R[a^{-1}]` is equivalent to saying that the forward
pole support is contained in `V(a)`.  What must additionally be checked is
that `d` has no zero away from `V(a)`.  Pole support alone does not imply
that the rational map lands in the opposite principal open.

### Corollary 1.2 -- the reciprocal height-one link

Suppose the source and core affine charts of a suspension are factorial
affine spaces with scalar units, and their maximal height-one common chart
is saturated in the sense of (1.3)--(1.4).  If its boundary orientation is
negative, then

\[
 X\setminus V(A)\simeq Z\setminus V(D),qquad
 D=cA^{-1}.                                            \tag{1.10}
\]

After rescaling `D`, this is exactly the reciprocal height-one link required
at the start of the cancellation bridge.  In particular, the reciprocal
exponent is not a numerical assumption and is not obtained from the
Jacobian ledger: it is forced by primitivity of the rank-one unit lattice.

This proves one substantial part of the supply problem.  The remaining
reciprocal data are tangential: the valuations of `P` and `Q-sP`, and the
noncontraction of the boundary under the invariant-plane map.

## 2. What divisor minimality should mean

There are three different kinds of minimality in the current suspension
problem.  They should not be conflated.

1. **Affine height-one minimality:** (1.3)--(1.4) has one prime boundary on
   each affine model.  This produces the common principal open and its sign.
2. **Graph-divisor minimality:** on a chosen proper graph model, every prime
   exceptional divisor is on the recorded list.
3. **Ledger minimality:** the relative-Jacobian identity has no coefficient
   at an unrecorded divisorial valuation, including valuations whose centers
   have codimension at least two on both original affine charts.

The first property does not imply the second or third.  A Rees valuation
over a codimension-two center is invisible in the complement of a prime
affine divisor.  Conversely, a component at infinity on a projective graph
can be irrelevant to the affine principal-open isomorphism but still carry
a nonzero Jacobian coefficient.

The following definition is therefore the one used in this note.

### Definition 2.1 -- divisor-minimal suspension square

Fix projective normal completions of all four affine charts, normal graph
models `Gamma_X` and `Gamma_Y` for `alpha` and `beta`, and a normal model
dominating the pullback of `Gamma_Y` to the source function field.  A
one-boundary suspension square is **divisor-minimal on these models** when:

1. every affine link is height-one saturated;
2. every prime exceptional divisor of each graph projection is recorded;
3. the logarithmic chain-rule identity is checked at every prime divisor of
   the relevant models, after pulling target valuations to the source, and
   its support outside the recorded boundary is empty;
4. deleting any recorded boundary prime destroys one of properties 1--3.

This is model-dependent unless one quantifies over all higher models.  The
model-independent version is the same statement for the relative Jacobian
`b`-divisors.  For the present chart problem, a fixed factorial graph model
with an exact relative Picard-rank audit is sufficient.

## 3. A relative Picard-rank audit for hidden valuations

The suggested factorial-graph hypothesis has a precise use, but it is not
to straighten the chart.

### Proposition 3.1 -- exceptional primes consume relative Picard rank

Let

\[
 f:\Gamma\longrightarrow X                             \tag{3.1}
\]

be a projective birational morphism of normal `Q`-factorial varieties.  The
classes of the prime `f`-exceptional divisors form a basis of
`N^1(Gamma/X)`.  In particular,

\[
 \boxed{\#\{\text{prime }f\text{-exceptional divisors}\}
       =\rho(\Gamma/X).}                               \tag{3.2}
\]

#### Proof

Let `M` be a `Q`-Cartier divisor on `Gamma`.  Since `X` is `Q`-factorial,
`f_*M` is `Q`-Cartier.  Hence

\[
 M-f^*f_*M                                              \tag{3.3}
\]

is an exceptional `Q`-divisor, proving that exceptional prime classes span
the relative Neron--Severi space.

If an exceptional `Q`-divisor `E` is numerically trivial over `X`, then both
`E` and `-E` are `f`-nef.  The negativity lemma applied in the two
directions gives `E<=0` and `-E<=0`; hence `E=0`.  The exceptional prime
classes are therefore independent.  QED

### Corollary 3.2 -- exact graph audit

If a normalized graph completion is factorial, and a projection has
relative Picard rank `rho`, then a displayed list of `rho` distinct
exceptional prime divisors is complete.  Any additional exceptional
valuation would contradict (3.2).

This is exactly where factoriality and bounded relative Picard rank help.
A mere upper bound is useful only after the displayed divisors already
consume it.  Smoothness of the boundary controls local parameters and
residue maps, but does not count exceptional primes.

Corollary 3.2 does not audit nonexceptional components of the projective
boundary.  Those must be checked by the principal divisors of the boundary
equations and by the ledger on `Gamma`.

## 4. Eliminating a hidden target ledger

Consider the logarithmic chain rule for a Keller suspension:

\[
 R_\alpha+r\alpha^*E=F^*R_\beta.                       \tag{4.1}
\]

The following elementary effectivity observation is useful in the
reciprocal branch.

### Proposition 4.1 -- zero reciprocal source ledger forces zero polynomial target ledger

Let all affine charts be affine spaces.  Assume:

1. `F` is dominant and has constant nonzero Jacobian;
2. `beta` is a polynomial morphism, so `R_beta` is effective;
3. on every divisorial valuation the reciprocal source terms cancel:
   \[
    R_\alpha+r\alpha^*E=0.                             \tag{4.2}
   \]

Then

\[
 \boxed{R_\beta=0.}                                   \tag{4.3}
\]

#### Proof

Equation (4.1) gives `F^*R_beta=0`.  If a prime divisor `C` occurred in the
effective divisor `R_beta`, choose its irreducible polynomial equation `c`.
Dominance makes `c\circ F` a nonconstant polynomial.  It therefore has a
nonzero effective principal divisor in the source affine space, contradicting
`F^*R_beta=0`.  QED

Thus a polynomial target chart cannot hide a positive ledger after exact
reciprocal cancellation.  Rational target charts are different: their
relative Jacobian divisors can have negative coefficients, and (4.1) alone
does not prevent cancellation among unrecorded target valuations.

If `beta` is also birational, (4.3) makes it a birational etale polynomial
map of affine space and hence an automorphism.  Indeed it is quasi-finite;
Zariski's Main Theorem makes it an open immersion, a missing divisor would
create a nonconstant unit on its affine-space source, and a codimension-two
complement is removed by normal Hartogs extension.  Thus such a target
chart can be absorbed by polynomial left equivalence.

In the positive weighted branch the target ledger is not zero and must not
be discarded.  The weighted identity

\[
 3V(x)+3V(\gamma)=3F^*V(C)                             \tag{4.4}
\]

is recorded data, not a hidden defect.

## 5. Three countermodels to an overstrong chart theorem

The following examples explain why the results above do not yet prove the
global dichotomy.

### 5.1 A same-orientation one-boundary suspension

Let

\[
 R=k[a,y,z],\qquad S=k[w,q,p],                          \tag{5.1}
\]

and, for `r>=1`, define the rational chart

\[
 \alpha(a,y,z)=(w,q,p)=(a^{-r}y,a,z).                  \tag{5.2}
\]

Its inverse is the polynomial map

\[
 a=q,\qquad y=q^rw,\qquad z=p.                         \tag{5.3}
\]

The principal opens `a!=0` and `q!=0` are isomorphic and the orientation is
positive: `q=a`.  The chart Jacobian is

\[
 \det D\alpha=-a^{-r}.                                \tag{5.4}
\]

Now take the coordinate-preserving plane core, stabilized by `p`,

\[
 \Phi(w,q,p)=(q,q^rw,p).                               \tag{5.5}
\]

Its Jacobian is `-q^r`, with sole reduced critical support `q=0`.  The
composition is

\[
 \Phi\circ\alpha(a,y,z)=(a,y,z).                       \tag{5.6}
\]

Hence it is a polynomial Keller map, indeed the identity.  The suspension
has one smooth boundary, one valuation, no target ledger, and exact
Jacobian cancellation, but its chart is neither a positive polynomial chart
nor a negative reciprocal link.  For `r>1` its `A^1` core is also one of the
higher-power cores not covered by the weighted simple-ramification theorem.

This example is divisor-minimal in the naïve sense of "one pole divisor."
It is excluded by requiring positive-orientation charts to be polynomial,
or by a nontriviality condition which forbids a vertical chart from being a
one-sided inverse factor of the core.  Neither exclusion follows from the
word minimal.

### 5.2 Positive orientation does not imply polynomial extension

Let

\[
 R=k[a,y,z],\qquad S=k[d,u,v]                           \tag{5.7}
\]

inside their common fraction field via

\[
 d=a,\qquad u=y+a^{-1},\qquad v=z.                     \tag{5.8}
\]

Then

\[
 S[d^{-1}]=R[a^{-1}],\qquad d=a,                       \tag{5.9}
\]

so the link has positive orientation and smooth factorial boundary on both
sides.  Nevertheless `u` has a pole at `a=0`, and neither the sign nor the
unit lattice makes the map polynomial.

This also survives a very small graph audit.  First compactify only the
`u`-line to `P^1`.  The rational function is represented on the affine
source by

\[
 [ay+1:a],                                             \tag{5.10}
\]

whose two entries have no common zero there.  Thus no affine graph divisor
detects the failure of polynomiality.  If both the `y`- and `u`-lines are
compactified to `P^1`, the homogeneous formula is

\[
 [Y_0:Y_1]\longmapsto[aY_0+Y_1:aY_1].                \tag{5.10a}
\]

It has the single boundary base point `a=Y_1=0`.  Blowing up that point
resolves the map and its inverse.  The resulting graph is smooth and
factorial, its boundary is simple normal crossing, and each graph projection
has relative Picard rank one.  Consequently even factorial graph
completion, smooth boundary, and the smallest nonzero relative Picard rank
do not force polynomiality in the chosen affine target chart.

### 5.3 The normal direction does not determine tangential valuations

Fix `r>=1` and `1<=ell<=r+1`.  Put

\[
 d=a^{-1},\qquad
 u=a^\ell y,\qquad
 v=a^{r+2-\ell}z.                                     \tag{5.11}
\]

Then

\[
 k[d,u,v,d^{-1}]=k[a,y,z,a^{-1}],                     \tag{5.12}
\]

the orientation is negative, and

\[
 \det{\partial(d,u,v)\over\partial(a,y,z)}=-a^r       \tag{5.13}
\]

for every `ell`.  However

\[
 (v_a(u),v_a(v))=(\ell,r+2-\ell)                      \tag{5.14}
\]

varies with `ell`.

Thus the reciprocal normal valuation `v_a(d)=-1` and the complete Jacobian
coefficient do not determine the tangential valuation vector.  In the
cancellation bridge, the conditions

\[
 v_A(Q-sP)=0,\qquad v_A(P)=1                           \tag{5.15}
\]

are genuinely residue-theoretic.  They cannot be recovered from the
rank-one unit lattice or the determinant `b`-divisor alone.

## 6. Boundary monotonicity

The obstruction in Sections 5.1--5.2 is exactly a negative valuation of an
affine target coordinate at a positive-oriented source boundary.

### Definition 6.1

A height-one saturated link `(R,a;S,d)` is **forward boundary-monotone** if

\[
 v_{(a)}(f)\ge0\qquad\text{for every }f\in S.           \tag{6.1}
\]

It is enough to check polynomial generators of `S`.  Define backward
boundary monotonicity symmetrically.

### Proposition 6.2 -- monotonicity is the polynomial-extension criterion

Assume `R` is normal and `S\subset R[a^{-1}]`.  If `a` is the only possible
forward pole divisor, then the rational chart
`Spec R \dashrightarrow Spec S` extends
to a morphism exactly when it is forward boundary-monotone.

#### Proof

If the chart extends, every element of `S` belongs to `R` and has
nonnegative valuation at `(a)`.  Conversely, elements of `S` already have
nonnegative valuation at every other height-one prime by the pole-support
hypothesis.  Condition (6.1) therefore puts them in every height-one DVR of
`R`.  Normality and the Krull intersection formula give `S\subset R`.  QED

This criterion is intrinsic, finite on a chosen coordinate set, and stable
under polynomial changes of source and target coordinates.  It is the exact
extra assertion needed to turn the positive sign into a polynomial chart.

### Theorem 6.3 -- signed one-boundary chart theorem under monotonicity

Let `R` and `S` be factorial polynomial rings with scalar units, and suppose
they form a height-one saturated one-boundary link.  Assume that every
positive-oriented link in the suspension square is forward
boundary-monotone.  Then every elementary link is of exactly one of the
following two types:

1. **positive polynomial type:** `d=ca`, and the rational chart extends to
   a polynomial morphism across `a=0`;
2. **negative reciprocal type:** `d=ca^(-1)`, and the two affine spaces are
   reciprocal completions of their common principal open.

If a factorial graph completion is also given and the displayed exceptional
primes consume the relevant relative Picard ranks, then there are no hidden
graph-divisorial valuations.  If the target chart in a negative branch is
polynomial and the reciprocal source ledger is zero, then it has no hidden
target ledger.

#### Proof

Proposition 1.1 gives the sign and the two equations for `d`.  In the
positive case Proposition 6.2 gives polynomial extension.  Corollary 3.2
and Proposition 4.1 give the last two assertions.  QED

This is a chart theorem, not yet the full suspension classification.  It
proves that under one additional natural hypothesis---boundary
monotonicity---the only normal-direction mechanisms are precisely the
polynomial and reciprocal mechanisms occurring in the repository.

## 7. Intrinsic primitive residue vectors

The normal primitive vector is now settled:

\[
 v_A(D)=+1\quad\text{or}\quad -1.                      \tag{7.1}
\]

Tangential vectors live on the normalization of the critical boundary and
require a separate saturation statement.

Let `E` be a smooth rational affine curve with at most two punctures and let
`bar E` be its smooth completion.

### Proposition 7.1 -- residue-lattice saturation

1. If `E\simeq\mathbb A^1` and a regular function `w|_E` induces a finite
   map of degree one to `A^1`, then `w|_E` is a polynomial coordinate.
2. If `E\simeq\mathbb G_m` and a boundary function
   `Y\in\mathcal O(E)^*` spans a nonzero saturated subgroup of
   `\mathcal O(E)^*/k^*`, then, after scaling and
   possibly interchanging the two punctures,
   \[
    (v_0(Y),v_\infty(Y))=(1,-1).                       \tag{7.2}
   \]
3. If another unit `s|_E` is nonconstant, then for a unique `m>=1`, after
   the same choice of orientation and a scalar rescaling,
   \[
    s|_E=Y^{-m}\quad\text{or}\quad Y^m.                \tag{7.3}
   \]

#### Proof

The first assertion is a finite birational map between normal affine
curves, hence an isomorphism.  For the second, the unit lattice modulo
scalars is `Z`.  Its only nonzero saturated subgroup is the whole lattice,
so `Y` is a generator and has the primitive vector (7.2), up to sign.  The
last assertion follows because every unit on `G_m` is a scalar monomial in
`Y`.  QED

This replaces a numerical valuation-vector assumption by an intrinsic
saturation condition.  It still does not prove that the primitive unit has
the required ambient lift

\[
 Y=Q-sP.                                               \tag{7.4}
\]

That lift is an embedding statement, not a consequence of divisor theory.
Once (7.4) is supplied and `s|_E=Y^(-m)`, the marked plane-core theorem gives

\[
 D=1-s(Q-sP)^m                                        \tag{7.5}
\]

and the reciprocal reconstruction theorem applies.  Similarly, the
degree-one statement in the `A^1` case supplies the coordinate marking used
by the weighted plane-core theorem.

In geometric degree three this last `A^1` marking is automatic for every
coordinate-preserving core: the Abhyankar--Moh degree reduction in
[`CUBIC_MARKING_EXTRACTION.md`](CUBIC_MARKING_EXTRACTION.md) straightens the
degree-two boundary projection by triangular shears preserving the core
coordinate.  The `G_m` residue and lifting problem remains; the same note
gives an infinite cubic toric defect atlas showing that puncture count alone
cannot settle it.

The exact remaining tangential problem is therefore:

> Prove that a nontrivial divisor-minimal Keller suspension makes the image
> of the distinguished ambient boundary functions saturated in the residue
> coordinate ring, and in the two-place case that its primitive generator
> has the affine-linear lift (7.4).

This is strictly narrower than classifying arbitrary birational charts of
`A^3`.

There is a second, transverse saturation condition on the source boundary.
Let `V` be the valuation ring of `v_A`, with maximal ideal `m_A`.

### Proposition 7.2 -- intrinsic transverse markings

Suppose `Y,P\in K^*` satisfy:

1. `Y\in V^*` and its residue is nonconstant on `V(A)`;
2. `P\in\mathfrak m_A` and its class in
   `\mathfrak m_A/\mathfrak m_A^2` is nonzero.

Then

\[
 \boxed{v_A(Y)=0,\qquad v_A(P)=1.}                    \tag{7.6}
\]

#### Proof

The first assertion says exactly that `Y` is a valuation-ring unit.  The
second says that `P` generates the conormal line of the smooth prime
boundary at its generic point, rather than lying in its square.  These are
respectively the valuation statements in (7.6).  QED

Proposition 7.2 is deliberately intrinsic: it asks that the displayed
functions give a nonconstant residue coordinate and a primitive conormal
coordinate.  It does not choose a uniformizer or prescribe a numerical
valuation vector in advance.  Example (5.11) shows that this transverse
saturation is additional information, not a consequence of the normal
reciprocal equation.

## 8. Conditional suspension dichotomy

Combining the preceding chart theorem with the existing downstream results
gives the following honest conditional statement.

### Theorem 8.1 -- one-boundary dichotomy with saturated boundary coordinates

Consider a coordinate-preserving, divisor-minimal one-boundary Keller
suspension whose normalized critical curve is smooth rational with at most
two punctures, and whose one-place core has reduced ramification.  Assume:

1. its affine links are height-one saturated;
2. positive-oriented links are boundary-monotone;
3. a factorial `Q`-factorial graph completion has exactly the displayed
   relative Picard ranks;
4. the critical-boundary coordinate map is residue-saturated in the sense
   of Proposition 7.1;
5. in the two-place case the primitive residue unit has the lift
   `Y=Q-sP`, the unit `s|_E` has negative exponent relative to this oriented
   generator, the pair `(Y,P)` is transversely saturated in the sense of
   Proposition 7.2, the invariant-plane boundary is not contracted, and
   the target chart is polynomial and birational.

Then, after polynomial left--right equivalence, its height-one chart and
plane core are of one of the following two types:

\[
\begin{array}{c|c|c}
 E&\text{chart}&\text{core}\\ \hline
 \mathbb A^1&\text{positive polynomial}&
 (w,q)\mapsto(q,wq-H(w))\\
 \mathbb G_m&\text{negative reciprocal}&
 (s,Q)\mapsto\left(Q,c\int_0^s
          \{1-t(Q-Pt)^m\}^r\,dt\right).
\end{array}                                             \tag{8.1}
\]

In the second row the suspension is the reciprocal cancellation
construction by the completed reciprocal-branch theorem.  There is no
unrecorded exceptional valuation or effective target ledger.  In the first
row (8.1) proves the polynomial weighted **core and chart mechanism**; the
identification of an arbitrary positive polynomial vertical chart with the
specific weighted formulas remains a separate positive-branch
straightening problem.

#### Proof

Theorem 6.3 supplies the signed polynomial/reciprocal chart dichotomy and
Corollary 3.2 removes hidden graph divisors.  Proposition 7.1 supplies the
primitive one- or two-place residue vector, while Proposition 7.2 supplies
the reciprocal height-one tangential markings.  The marked plane-core
theorem then gives the two displayed normal forms.  In the reciprocal row,
Proposition 4.1 removes an effective target ledger and the existing
reciprocal reconstruction theorem produces the cancellation slice and
finite jet.  No downstream cancellation reconstruction is repeated here.
QED

## 9. Exact use of divisor minimality

The proof audit is now short.

Divisor minimality is used to:

1. reduce the affine pole and zero support to the two declared prime
   boundaries, enabling Proposition 1.1;
2. make the boundary unit lattice rank one, so its automorphism is `+1` or
   `-1` and not a larger integral matrix;
3. compare the displayed exceptional primes with the full relative Picard
   rank in Corollary 3.2;
4. ensure that (4.1) has been tested on every divisorial valuation rather
   than only on the original affine charts.

It is **not** used, and is not sufficient, to:

1. make a positive-oriented rational chart polynomial;
2. exclude the same-orientation inverse suspension (5.2)--(5.6);
3. determine tangential valuations such as `v_A(P)`;
4. produce the primitive lift `Q-sP` on a `G_m` critical curve;
5. prove boundary noncontraction;
6. straighten an arbitrary positive polynomial chart to the explicit
   weighted vertical formulas.

Items 1--2 require boundary monotonicity or a stronger nontriviality
condition.  Items 3--5 are residue and embedding statements.  Item 6 is the
remaining weighted-chart theorem.

## 10. Next attack

The next useful proof should not classify all birational self-maps of
`A^3`.  It should prove one of the following consequences of the Keller
square.

1. **Positive monotonicity:** a positive-oriented elementary link with a
   generically nonbirational plane core has no negative boundary valuation
   on its affine target coordinate ring.
2. **Residue saturation:** the boundary restriction of the distinguished
   core coordinates contains a primitive generator of the `G_m` unit
   lattice.
3. **Linear lifting:** that generator can be chosen as `Q-sP` after a
   polynomial shear preserving the family parameter `P`.

The same-orientation example (5.2)--(5.6) says that any proof of item 1 must
use nontriviality of the core or of the finite map, not merely the ledger.
The family (5.11) says that any proof of items 2--3 must use residue
geometry, not only the normal valuation and its discrepancy.

These three statements are now the precise upstream chart frontier.
