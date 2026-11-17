# Mixed-sign gradings: literature reconciliation

## Status and source version

This is a literature reconciliation, not a new theorem or computation.  It
does not change any entry of
[`MATH_STATUS.json`](../MATH_STATUS.json).  Repository claims below are only
the claims already recorded there, cited by status identifier when useful.

The external source is T. Shaska,
[*Graded Keller maps and the Jacobian Conjecture*](https://arxiv.org/abs/2607.20210v2),
arXiv:2607.20210.  Version 1 was submitted on 22 July 2026; the current
version 2 was submitted on 25 July 2026.  This version distinction matters:
v2 adds the general \((1,-r,-s)\) quotient exponent, the bounded parameter
schemes, higher-dimensional extensions, and the explicit classification
boundary.  References below use the numbering in v2.  The
[22 July version](https://arxiv.org/abs/2607.20210v1) is retained only to
identify the paper named in the original comparison request.

The short conclusion is:

> Shaska's classification problem concerns Keller maps **carrying a chosen
> algebraic \(\mathbb G_m\)-equivariance**, organized by the signature of
> its weights; some rows are complete and other rows are reductions or open.
> The repository's stable-left--right results ask a different, stronger
> equivalence question on narrower explicit families: which normalization
> and boundary data can be reconstructed from the map after arbitrary
> polynomial source and target changes and after adjoining identity
> variables?

Neither framework subsumes the other.

## 1. Dictionary

Let

\[
 G:\mathbb A^n_{\mathbf w}\longrightarrow\mathbb A^n_{\mathbf q}
\]

be equivariant for diagonal source and target actions.  Shaska calls this a
*graded map*.  The repository usually calls it *torus-equivariant* or
*weight-homogeneous* and sometimes calls \(\mathbf q\) the vector of target
degrees.  For a Keller map, Shaska's Corollary 2.3 proves that
\(\mathbf q\) is a permutation of \(\mathbf w\).

The terminology translates as follows.

| Shaska | Repository usage |
|---|---|
| elliptic | every weight is nonzero and all have one sign |
| parabolic | at least one weight is zero |
| hyperbolic | positive and negative weights both occur; usually “mixed-sign” |
| quotient | the categorical quotient `Spec` of the weight-zero invariant ring |
| graded equivalence | source and target automorphisms preserving the selected grading |
| stable polynomial left--right equivalence | arbitrary polynomial left--right equivalence after adjoining identity coordinates |

The last two rows are not synonyms.  Shaska's set
\(\mathcal K(n,\mathbf w)\) fixes a weight multiset and quotients by graded
source and target automorphisms.  A stable polynomial left--right
equivalence need not preserve any displayed grading or quotient
coordinates.

There are also three different uses of toric or stacky language which must
not be conflated.

1. Shaska's acting \(\mathbb G_m\) is a symmetry group of the affine map.
2. A repository phrase such as “the intrinsic two-torus” can mean that a
   normalized boundary stratum is the variety \(\mathbb G_m^2\).  It does
   not assert a two-torus of source--target symmetries.
3. Shaska's “stacky line” is the stabilizer-\(\mu _2\) stratum of the
   chosen target action.  The repository's marked admissible-cover stack
   and its quotient-stack stabilizers arise from finite normalization or
   moduli data.  They are related in viewpoint, but are not the same stack.

## 2. The foundational grading and quotient

For the foundational map

\[
\begin{aligned}
a&=(1+xy)^3z+y^2(1+xy)(4+3xy),\\
b&=y+3x(1+xy)^2z+3xy^2(4+3xy),\\
c&=2x-3x^2y-x^3z,
\end{aligned}
\]

Shaska uses

\[
\mathbf w=(1,-1,-2),\qquad \mathbf q=(-2,-1,1).
\]

Thus

\[
\sigma_t(x,y,z)=(tx,t^{-1}y,t^{-2}z),\qquad
\tau_t(a,b,c)=(t^{-2}a,t^{-1}b,tc).
\]

The source weights are forced up to common scaling by the monomials of
\(c\).  The target weights are their reversal, as the nonzero
anti-diagonal linear part requires.  This is exactly the grading used in
the repository's
[plane-quotient audit](../plane-jc/FOUNDATIONAL_QUOTIENT_AUDIT.md) and
[torus filtered module](TORUS_FILTERED_LR_MODULE.md).

Shaska's quotient coordinates are

\[
u=xy,\qquad v=x^2z,\qquad P=bc,\qquad Q=ac^2.
\]

Writing

\[
L=2-3u-v,
\]

the third component is \(c=xL\).  Theorem 7.1 gives

\[
\operatorname{Jac}_{u,v}(P,Q)=2L^2
\]

for the foundational determinant \(-2\).  The repository's
plane-quotient audit writes the target generators in the opposite order,

\[
(P_{\rm repo},Q_{\rm repo})=(ac^2,bc),
\]

and consequently obtains

\[
\operatorname{Jac}_{u,v}(P_{\rm repo},Q_{\rm repo})=-2L^2.
\]

These are the same identity; the sign is only the transposition of the
target quotient coordinates.  Both calculations also identify
\(V(L)\) as a line on the source quotient contracted to the quotient
origin.  The descended plane map is therefore not Keller, so it does not
conflict with Shaska's dimension-two theorem.

The repository's
[weighted invariant-coordinate reduction](WEIGHTED_INVARIANT_JACOBIAN_REDUCTION.md)
proves, for weights \((1,-1,-k)\),

\[
\operatorname{Jac}_{u,v}(bc,ac^k)
   =-\Lambda^k\det JG.                                \tag{2.1}
\]

Shaska v2, Proposition 8.2 and Remark 8.4, contains the wider formula for
\((1,-r,-s)\),

\[
\det JG
 =\pm\frac{\operatorname{Jac}_{u,v}
 (B\Lambda^r,A\Lambda^s)}
 {\Lambda^{r+s-1}},                                   \tag{2.2}
\]

and the analogous higher-dimensional exponent
\(\sum r_i-1\).  Therefore (2.1) should now be described as the
\(r=1,s=k\) specialization, together with the repository's explicit
family-level factorization for every admissible weighted seed.  It is not
a generalization of Shaska v2.

## 3. What Shaska proves

The following is the external theorem boundary relevant here.

1. **Elliptic case, every dimension.**  Theorem 3.1 proves over
   \(\mathbb C\) that an elliptically graded Keller map is a graded
   polynomial automorphism.  Properness of the weighted action is the
   decisive input.
2. **Every signature in dimension two.**  Theorem 3.4 treats arbitrary
   algebraic \(\mathbb G_m\)-actions on source and target.  After
   polynomial linearization, the elliptic, parabolic, and hyperbolic cases
   become triangular, affine, and linear, respectively.
3. **The foundational quotient.**  Sections 4--7 recover the cubic inverse
   equation, degree three, the contracted line, the discriminant and
   escaping sheets, the rational thin image, the stabilizer-\(\mu _2\)
   line, and the order-two quotient-Jacobian identity.
4. **One-positive-weight hyperbolic maps.**  Sections 8--9 write maps of
   weights \((1,-r_2,\ldots,-r_n)\) as invariant-polynomial data and reduce
   their Keller condition to one multilinear equation.  With degree bounds,
   counterexamples are rational points of an explicit multiprojective
   scheme.
5. **Classification boundary.**  Section 10 completely removes the
   elliptic case and every signature in dimensions at most two.  It reduces
   a parabolic case with elliptic nonzero part to the ordinary Jacobian
   conjecture in the number of zero-weight variables.  Hyperbolic actions
   with one positive weight equal to one are reduced to the parameter
   schemes.  One positive weight at least two, and more than one positive
   weight, remain open in general.

Thus “classification by signature” is an organizing classification with
complete and reduced rows, not a closed normal-form theorem for every
hyperbolic action.

## 4. Repository results that genuinely strengthen the comparison

The word *strengthen* here always means “a stronger equivalence or
intrinsicness conclusion under narrower family hypotheses,” not a broader
classification by weights.

### 4.1 Stable normalization functoriality (`S1`)

[Stable normalization functoriality](../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md)
starts with a dominant quasi-finite map \(F:U\to Y\), not with a chosen
action.  It reconstructs

\[
\mathcal B(F)=
(\overline X_F\to Y,\ U\hookrightarrow\overline X_F,\ \partial_F)
\]

from the map itself.  It proves transport under arbitrary polynomial
left--right isomorphism and polynomial stabilization of:

- boundary primes and their target images;
- ramification and residue degrees;
- scheme intersections and nilpotency indices;
- relative-differential Fitting ideals;
- node pairings and conductors; and
- affine-versus-boundary status of divisorial branches.

This is stronger than equivariance for purposes of stable classification:
it survives changes which do not preserve the displayed torus action.
It is also conditional in a precise sense: a named root, quotient
coordinate, or compactification point is not invariant until a
family-specific theorem characterizes it uniquely from \(\mathcal B(F)\).

### 4.2 Weighted intrinsic boundary (`WB1`, `D1`, `F2`)

On the boundary-clean weighted-seed locus, with zero an exact double root
and all nonzero roots simple, the
[decorated-normalization theorem](DECORATED_NORMALIZATION_INVARIANT.md)
reconstructs from the map:

- exactly two ordered target boundary images
  \((Z_\Delta,Z_0)\);
- the ramified discriminant normalization after deleting \(Z_0\);
- the full Fitting divisor \(\operatorname{div}(H'')\), including
  multiplicities;
- the second-boundary center; and
- on the stated clean open, the distinguished affine root sheet.

This strengthens the coarse quotient picture because the quotient origin
collapses information.  For example, the target plane \(c=0\) and several
orbit types can have the same coarse quotient image, whereas the
Zariski--Main package retains which normalized branches lie in the affine
source and which are boundary branches.

The hypotheses are narrower than Shaska's: these statements concern the
repository's weighted seed construction on certified clean loci, not every
hyperbolically graded map.

The resulting marked Hessian divisor has stable-moduli dimension \(N-3\)
(`D1`), and adding the affine root sheet generically recovers the normalized
seed exactly (`F2`).  This is a classification by **stable polynomial
left--right class** inside the weighted family, not by chosen weight
signature.  Shaska's parameter schemes and graded-automorphism quotient do
not make this stable identification.

### 4.3 Stable separation from other construction skeletons (`RQG2`)

The
[quadratic-versus-weighted stable-separation theorem](../verified/QUADRATIC_WEIGHTED_STABLE_SEPARATION.md)
proves that no boundary-clean weighted map of degree \(N\geq4\) is stably
polynomially left--right equivalent to an admissible quadratic-gauge map.
After deleting the intrinsically ordered second boundary, their normalized
ramified strata are

\[
\mathbb A^1\times\mathbb G_m
\quad\hbox{and}\quad
\mathbb G_m^2,
\]

with unit ranks one and two.  These unit ranks survive affine
stabilization.

This has no direct analogue in Shaska's classification: it separates two
left--right construction classes even when a visible quotient or inverse
polynomial presentation makes them look similar.  It does not classify all
mixed-sign graded maps.

### 4.4 Orbit-wide exclusion of grading (`LTS2`)

For the repository's geometric-degree-four quadratic-gauge example, the
[intrinsic algebraic-torus theorem](../cancellation/NO_ALGEBRAIC_TORUS_EQUIVARIANCE.md)
computes a decorated ramified-normalization automorphism group \(\mu _5\)
and zero infinitesimal automorphism algebra.  It then proves:

\[
\text{no polynomial left--right representative of that map is
equivariant for a nontrivial algebraic torus.}
\]

This is strictly stronger than checking that one displayed formula has no
linear grading: it excludes nonlinear polynomial conjugates as well.  It is
not an extension of Shaska's positive-weight or plane theorem, because it
concerns a different three-dimensional, geometric-degree-four,
noninvertible map.

There is an essential stable caveat.  After adjoining an identity variable,
the new variable carries a tautological torus action.  The theorem therefore
does **not** claim literal torus-freeness after stabilization.  Its stable
conclusion is only that every connected torus acting on the pulled-back
decoration is vertical over the intrinsic \(\mathbb G_m^2\) stratum; no
splitting of that vertical action is proved.

### 4.5 Bounded three-variable equivariant classification

The
[degree-four equivariant classification](../verified/DEGREE_FOUR_EQUIVARIANT_CLASSIFICATION.md)
proves over characteristic-zero fields that every coordinate-degree-at-most-
four Keller map on \(\mathbb A^3\) equivariant for a nontrivial **linear**
\(\mathbb G_m\)-action is a tame automorphism.  It exhausts the bounded
weight-support coefficient schemes, including all mixed-sign and zero-weight
patterns.

This is a concrete bounded strengthening within dimension three.  Shaska
v2 constructs the general parameter schemes and proves several emptiness
rows, but explicitly leaves many bounded hyperbolic strata open.  The
repository theorem is limited to degree at most four in coordinates after
linearization; a nonlinear linearizing conjugacy can raise degree.

## 5. The stacky overlap, stated narrowly

Shaska's stabilizer-\(\mu _2\) line is

\[
\{b=c=0,\ a\ne0\}.
\]

For a rational target \((\alpha,0,0)\), the two points on the contracted
locus satisfy

\[
x^2=-\frac1{4\alpha};
\]

they are rational exactly when \(-\alpha\) is a square.  This is a
stabilizer-stratum refinement invisible on the coarse quotient.

The repository's intrinsic boundary theory makes a parallel methodological
move: it retains finite normalized stratum maps, inertia, Fitting divisors,
conductors, and selected affine sheets rather than only their coarse target
images.  The
[intrinsic-selector attack](INTRINSIC_SELECTOR_ATTACK.md) shows why this
extra structure matters: a coarse \(\mu _3\) symmetry of a symmetric
quintic incidence cover is not a stable polynomial left--right
self-equivalence, because it exchanges the unique affine sheet with
pole-two boundary sheets.

This is strong conceptual overlap, but not an identification of the two
stacky strata.  In particular:

- Shaska's \(\mu _2\) is isotropy of the chosen affine torus action;
- the selector attack's \(\mu _3\) is a finite symmetry of a normalized
  incidence cover; and
- the quadratic-gauge \(\mu _5\) is a decorated-stratum or moduli
  stabilizer.

Their orders and their functorial origins are different.

## 6. Reconciled claims

The repository should use the following formulations.

1. The foundational map is Shaska-hyperbolic, with the unique primitive
   source weights \((1,-1,-2)\) up to sign and target permutation.
2. The repository quotient formula agrees exactly with Shaska's after
   swapping the two target invariant generators.
3. The exponent-\(k\) repository formula for \((1,-1,-k)\) is a
   specialization of Shaska v2's exponent \(r+s-1\), not a broader external
   result.
4. Shaska's elliptic theorem and dimension-two theorem are broader than any
   corresponding repository calculation and should be cited as the
   literature results.
5. `S1`, `WB1`, `D1`, and `F2` strengthen the **invariance and stable
   classification** of the repository's clean weighted families; they do
   not strengthen Shaska's global signature classification.
6. `RQG2` separates weighted and quadratic-gauge stable classes by
   intrinsic boundary unit rank, a question outside a fixed chosen grading.
7. `LTS2` upgrades “this formula is not visibly graded” to “no polynomial
   left--right representative is torus-equivariant” for the repository's
   quartic example, with the stated stabilization caveat.
8. The phrase “stacky stratum” must always specify whether it means action
   isotropy, an incidence/admissible-cover stack, or a moduli stabilizer.

## 7. Continuations and open angles

### 7.1 A simplification that closes immediately: graded formal triviality

The repository's formal source-triviality theorem has an equivariant
strengthening which is now recorded in
[the canonical note](FORMAL_ORBIT_TRIVIALITY.md#equivariant-corollary).
If \(F\) and an Artin deformation \(\mathcal F\) are equivariant for the
same source and target actions, the unique polynomial source automorphism
\(\alpha\) with

\[
\mathcal F=F\circ\alpha
\]

commutes with the source action.  The proof conjugates \(\alpha\) by the
universal group action and invokes uniqueness over the nilpotent extension
of the group coordinate ring.

For Shaska's fixed-weight parameter schemes this has a sharp consequence:

> the Artin deformation functor at every graded Keller point, modulo the
> full graded polynomial source-automorphism functor, is a point.

This does **not** prove that
\(\mathcal Y^\circ_{1,2,(4,3,1)}\) is one
\(\mathcal G\)-orbit.  It shows why tangent ranks and nonreduced local
directions cannot answer that question.  Any obstruction must be one of:

- failure of the formal trivializer to algebraize over a reduced base;
- growth beyond a prescribed degree or support filtration;
- passage between different global orbit components; or
- a global normalization-boundary invariant.

This replaces an unrestricted infinitesimal search by a filtered
algebraization problem.

### 7.2 The foundational parameter scheme versus the repository slice

For Shaska's
\((r,s,\mathbf d)=(1,2,(4,3,1))\), the coefficient spaces are

\[
\dim V_A=13,\qquad \dim V_B=9,\qquad\dim V_\Lambda=3,
\]

with

\[
V_A=(u^2,v)_{\le4},\qquad
V_B=(u,v)_{\le3},\qquad
V_\Lambda=k[u,v]_{\le1}.
\]

The repository's
[foundational weighted coefficient scheme](FOUNDATIONAL_WEIGHTED_COEFFICIENT_SCHEME.md)
is a closed linear slice of this space.  Its \(A\)-support omits

\[
v^2,\ uv^2,\ v^3,\ u^2v^2,\ uv^3,\ v^4,
\]

and its \(B\)-support omits

\[
v^2,\ uv^2,\ v^3.
\]

It also fixes the anti-diagonal linear part and the two effective diagonal
gauges.  On that slice the Keller scheme is the single nonreduced point

\[
k[\varepsilon]/(\varepsilon^2).
\]

The equivariant formal-triviality corollary explains the nilpotent tangent:
it is a graded source-orbit direction in the full polynomial automorphism
functor, even though it is not removed by the small diagonal gauge.

The exact next calculation is therefore not another tangent matrix.  It is
the saturation of the full \(13+9+3\) coefficient scheme by the graph of
the liftable graded automorphism action, followed by a test for more than one
reduced orbit.  A practical intermediate target is to cover the
\(([B],[\Lambda])\)-base of Shaska's Proposition 8.12 by rank-minor charts
and compare every surviving chart with the foundational orbit.  The main
obstruction is degree growth: the unique formal graded trivializer need not
belong to the finite-dimensional subgroup preserving the bounds
\((4,3,1)\).

### 7.3 The smallest unresolved signature has a sharper first target

Shaska asks whether
\(\mathcal K(3,(1,-1,-1))\) is empty and identifies
\((2,1,d_\Lambda)\) and its transpose as the first open invariant-degree
rows.  For these weights,

\[
G=(x^{-1}A(u,v),x^{-1}B(u,v),x\Lambda(u,v)),
\]

and an invariant monomial of degree \(d\) lifts to source degree
\(2d-1\) in \(A,B\) and \(2d+1\) in \(x\Lambda\).

The repository's
[degree-four equivariant classification](../verified/DEGREE_FOUR_EQUIVARIANT_CLASSIFICATION.md)
therefore already removes the first subrow

\[
(2,1,1)\quad\text{and}\quad(1,2,1):
\]

their coordinate degree is at most three, and every such Keller map is
tame.  The smallest genuinely new finite problem is

\[
\boxed{(1,-1,-1),\qquad
\mathbf d=(2,1,2)\ \text{or}\ (1,2,2),}
\]

whose coordinate degree is five.  This is a better next elimination target
than an unbounded search over the signature.  Dependency-graph acyclicity
should remove triangular supports before forming coefficient ideals.

### 7.4 Recover the grading from the intrinsic boundary

Existence of a nontrivial algebraic torus action is preserved by ordinary
polynomial left--right conjugacy, but the displayed weight vector is not
automatically an intrinsic label.  The foundational grading is unique in
its displayed coordinates; what is not yet proved is an orbit-wide
statement such as

\[
\operatorname{Aut}^{\sharp}(F)^\circ\simeq\mathbb G_m
\]

with this \(\mathbb G_m\) equal to the foundational grading up to conjugacy.

This is the natural bridge between Shaska's Remark 8.15 and the repository's
boundary method.  The first deliverable is to compute the action of
\(\operatorname{Aut}^{\sharp}(F)\) on the full intrinsic normalization
package—ordered target boundary images, affine/boundary sheets, Fitting
divisor, and conductor—and prove that its connected kernel is exactly the
known grading torus.  The component group would then control both uniqueness
of the grading and the field-of-moduli obstruction.

The quartic theorem `LTS2` is the opposite endpoint: there the connected
torus is killed entirely.  The foundational cubic should be the positive
calibration where exactly one torus survives.

### 7.5 Relate the quotient stack to the Zariski--Main boundary

The coarse quotient forgets precisely the data needed by the intrinsic
boundary theory.  A direct comparison should use the stack morphism

\[
[\mathbb A^3_{\rm source}/\mathbb G_m]
\longrightarrow
[\mathbb A^3_{\rm target}/\mathbb G_m]
\]

together with the finite normalization of the original affine map.  For the
foundational example, the first exact goals are:

1. identify the inertia-\(\mu _2\) line and its inverse image in the
   normalized cover;
2. show how the coarse quotient origin separates into the intrinsic
   \(Z_\Delta\) and \(Z_0\) incidence data before coarse contraction; and
3. compare Shaska's Kummer square class with the residue field and inertia
   labels of the corresponding normalized boundary branches.

The obstruction is functorial: Zariski--Main normalization is formed from
the function-field extension and affine open immersion, while a quotient
stack is formed from the chosen action.  No general commutation theorem
between these constructions is currently supplied.

### 7.6 Make the quotient exponent coordinate-free

Shaska's exponent

\[
m-1=\sum_{i\ge2}r_i-1
\]

is obtained from explicit quotient coordinates.  It should instead be the
coefficient of the contracted divisor in a relative canonical or
determinant-of-cotangent identity on the quotient stack.  A coordinate-free
statement would:

- explain the exponent without choosing invariant generators;
- show exactly which part is preserved when the grading is intrinsically
  recovered;
- extend naturally to singular quotients using reflexive canonical modules;
  and
- replace the scalar factor \(\Lambda\) by a Fitting or determinant ideal
  when the weight-one module has rank greater than one.

The first test is to rederive
\(\operatorname{Jac}(P,Q)=\kappa\Lambda^{r+s-1}\)
as an equality of sections of canonical modules on the quotient, then take
the valuation along \(V(\Lambda)\).  The main obstruction is that outside
the cyclic weight-one case the contracted locus need not be Cartier.

### 7.7 The two omitted hyperbolic shapes

The coordinate-free formulation points to concrete models for Shaska's two
open shapes.

1. For weights \((2,-1,-1)\), the invariant ring is the quadric cone
   \[
   k[xy^2,xyz,xz^2]\simeq k[p,q,r]/(q^2-pr).
   \]
   The graded pieces are reflexive modules, not free modules over a
   polynomial ring.  The immediate task is to express the Keller condition
   as a section of the cone's reflexive canonical module and compute its
   divisor class.  The repository's Cox and boundary-lattice methods are
   relevant here.  The likely obstruction to a scalar \(\Lambda\) formula
   is the nontrivial class group at the vertex.
2. For the first genuinely multi-positive case, represented after sign
   normalization by \((1,1,-1,-1)\), the weight-one piece has rank two over
   the invariants.  A scalar contracted factor should be replaced by the
   Fitting ideal of the map from the positive-weight module to its target
   counterpart.  The quotient is determinantal, and the natural boundary
   is multicomponent.  The repository's multiboundary ledgers are a closer
   model than the one-divisor formula.

These are structural programs, not bounded searches.  They should follow,
not precede, the canonical-module reformulation of the exponent.

### 7.8 The direct intrinsic-classification gap

The closest existing bridge is the
[controlled-divisor equivariant program](../cancellation/CONTROLLED_DIVISOR_EQUIVARIANT_CLASSIFICATION.md).
Its missing implication is

\[
\text{intrinsic one-boundary marked cover}
\Longrightarrow
\text{coordinate-preserving affine-linear quotient pencil}.
\]

Closing this would turn the repository's boundary data into an actual
normal-form theorem for a Shaska-hyperbolic stratum.  The sharp obstructions
are already isolated there: higher incidence coordinate, unsaturated
two-place marking, nontrivial target ledger, or an additional reconstruction
divisor.  This is the highest-value conceptual continuation after the
foundational single-orbit calculation.

### 7.9 Priority order

The recommended order is:

1. compute the full foundational graded orbit quotient
   \(\mathcal Y^\circ_{1,2,(4,3,1)}/\mathcal G\), using formal graded
   triviality to avoid uninformative tangent work;
2. run the finite coordinate-degree-five
   \((1,-1,-1)\), \((2,1,2)\) coefficient classification;
3. compute \(\operatorname{Aut}^{\sharp}(F)\) on the intrinsic foundational
   boundary package;
4. build the quotient-stack/Zariski--Main comparison for the foundational
   \(\mu _2\) line;
5. formulate the quotient exponent through canonical modules; and only then
6. attack the singular-quotient and multi-positive hyperbolic shapes.

The one-boundary incidence-extraction theorem can proceed in parallel with
items 1--4, because it uses the same foundational boundary package but not
the singular quotient generalization.

## References

- T. Shaska, *Graded Keller maps and the Jacobian Conjecture*,
  [arXiv:2607.20210v2](https://arxiv.org/abs/2607.20210v2), 25 July 2026.
- The repository theorem-status authority is
  [`MATH_STATUS.json`](../MATH_STATUS.json); the generated table is
  [`STATUS.md`](../STATUS.md).
