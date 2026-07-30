# Programme 1: logarithmic boundary types

> **Status.** This is a definition and a first classification census.  It
> introduces no new realization theorem.  Every row labelled `realized`,
> `obstructed`, or `open` below points to the canonical note containing that
> result.  In particular, a balanced determinant ledger or a package passing
> the stage-one compiler is not by itself a polynomial Keller realization.

The purpose of the programme is to classify the boundary mechanism rather
than a preferred coordinate formula.  The weighted and reciprocal families
therefore appear as the first two cells of one category.  Their formulas are
representatives of those cells, not the definitions of the cells.

Work over a characteristic-zero field \(k\).  Geometric statements about
punctures are made after base change to \(\bar k\).

## 1. The category

### 1.1 Objects

A **logarithmic boundary package**
\(\mathfrak B\) consists of the following data.

1. **Normalized critical curves.** A finite selected set \(I\), a reduced
   critical curve \(C_i\) for each \(i\in I\), and its finite normalization
   \[
   \nu_i:\widetilde C_i^\circ\longrightarrow C_i.
   \]
   Write \(\overline C_i\) for the smooth projective completion of
   \(\widetilde C_i^\circ\).

2. **Puncture divisors.** A reduced divisor
   \[
   P_i=\overline C_i\setminus\widetilde C_i^\circ.
   \]
   The package retains the divisor, not only its cardinality.  It also
   retains the puncture-divisor homomorphism
   \[
   \operatorname{div}_{P_i}:
   \mathcal O(\widetilde C_i^\circ)^*/k^*
   \longrightarrow
   \operatorname{Div}_{P_i}^0(\overline C_i).
   \tag{1.1}
   \]

3. **Source reconstruction valuations.** A finite marked set
   \(\mathcal R\subset k(X)^*\) of reconstruction functions, a finite set
   \(\mathcal V\) of divisorial valuations on a common normalized graph, and
   the matrix
   \[
   M_{\mathcal R}=(v_E(r))_{E\in\mathcal V,\ r\in\mathcal R}.
   \tag{1.2}
   \]
   When an affine chart is asserted, retain also its actual value semigroup,
   its pole bounds, and a polar-completeness certificate.  Replacing the
   actual semigroup by its saturation is not allowed.

4. **Ramification and different.** For every selected upstairs prime
   \(E\) over a target prime \(Z\), retain
   \[
   E\longrightarrow Z,\qquad
   e(E/Z),\qquad k(Z)\subset k(E),\qquad f(E/Z),
   \tag{1.3}
   \]
   the different ideal, and the completed height-one map.  In characteristic
   zero the generic DVR extension is tame and the different exponent is
   \(e-1\), but the residue map and the completed equation remain part of the
   object.

5. **Conductor.** Retain the full conductor square of each finite curve
   normalization:
   \[
   \begin{array}{ccc}
   \overline{\mathcal C}_i&
      \longrightarrow&\widetilde C_i^\circ\\
   \downarrow&&\downarrow\nu_i\\
   \mathcal C_i&\longrightarrow&C_i ,
   \end{array}
   \qquad
   \mathfrak c_i=
   \operatorname{Ann}_{\mathcal O_{C_i}}
   \bigl(\nu_{i*}\mathcal O_{\widetilde C_i^\circ}/
               \mathcal O_{C_i}\bigr).
   \tag{1.4}
   \]
   Thus a node remembers its two normalization branches and a cusp remembers
   its thick normalization point.  The conductor is not compressed to its
   length or \(\delta\)-invariant.

6. **Target boundary images.** Retain the target boundary primes and strata,
   every map \(E\to Z\), every map \(C_i\to T_i\) to a selected target
   stratum, and all equal-image and intersection incidences.  This separates
   two source primes with the same target image from two primes mapping to
   different target components.

7. **Affine-sheet markings.** Every divisorial sheet has the color
   \[
   \epsilon(E)=
   \begin{cases}
   0,&E\text{ meets the distinguished affine source},\\
   1,&E\text{ is a missing reconstruction boundary}.
   \end{cases}
   \tag{1.5}
   \]
   The package also retains any intrinsically selected unramified affine
   root sheet and its specialization to conductor strata.  This mark is
   essential: conductor, Fitting, and node-pairing data alone do not remove
   rerooting ambiguity.

8. **Determinant-ledger relation.** On one common normalized graph retain
   the rational square
   \[
   \begin{array}{ccc}
   X&\stackrel F\dashrightarrow&Y\\
   q_X\downarrow&&\downarrow q_Y\\
   Z&\stackrel\Phi\dashrightarrow&T
   \end{array}
   \]
   and, for every \(E\in\mathcal V\), the signed coefficient
   \[
   \ell_E=
   v_E(R_{q_X})+
   v_E(q_X^*R_\Phi)-
   v_E(R_F)-
   v_E(F^*R_{q_Y}).
   \tag{1.6}
   \]
   The support list must be complete: an unrecorded Rees valuation or a
   redundant spectator row is a defect even if the displayed coefficients
   sum to zero.

9. **Incidence strata.** For reducible selected loci retain the reduced,
   scheme-theoretic, and completed intersections of selected components.
   This final block records how the preceding eight blocks are coupled; it
   is not an additional numerical invariant.

The ninth block is forced by the first eight bullets in a single-component
normal-crossing example, but not for multiple selected critical components.
It is therefore written explicitly.

### 1.2 Compatibility axioms

The blocks are required to satisfy the following consistency conditions.

1. Restricting the source valuation matrix to a selected normalized curve
   gives its puncture-divisor rows.  A function declared regular on the
   reconstruction open has no unrecorded negative valuation.
2. The unit lattice on \(C_i\) is the conductor-descending sublattice of the
   puncture-unit lattice on \(\widetilde C_i^\circ\).  Thus conductor
   characters may change the index or rank, but they cannot be added after
   the puncture calculation as unrelated labels.
3. Every completed height-one equation gives the recorded ramification
   index and different, and its residue map is the restriction of the
   recorded target-boundary map.
4. Component intersections, equal-target incidences, and affine-sheet
   specializations commute with the normalization and conductor squares.
5. The valuation list is complete for every divisor occurring in the four
   terms of (1.6).  In the realized Keller subgroupoid, \(R_F=0\), and every
   ramified sheet is boundary-colored because the affine source is étale.

These axioms make the package a single diagram rather than a tuple of
unrelated numerical invariants.

### 1.3 Morphisms

Let \(\mathsf{LBT}_k\) be the groupoid whose objects are logarithmic boundary
packages and whose arrows are isomorphisms of the entire marked diagram.
An arrow consists of:

- a bijection of selected critical components, source primes, target primes,
  and recorded strata;
- isomorphisms of the normalized curves and target strata carrying
  \(P_i\), the conductor squares, and all completed local maps to their
  counterparts;
- an isomorphism of the common function-field diagrams carrying
  \(\mathcal R\), \(\mathcal V\), their valuation matrix, and actual value
  semigroups to the corresponding data;
- preservation of \(e\), residue extensions, differents, target images,
  affine colors, and affine-sheet marks; and
- equality of the transported ledger rows \(\ell_E\).

All arrows are invertible.  This is deliberate: the present programme
classifies types up to marked isomorphism.  Contractions and forgetful maps
belong to separate functors, not to the equivalence relation.

There are useful full subgroupoids

\[
\mathsf{LBT}^{\mathrm{bal}}_k
=\{\mathfrak B:\ell_E=0\text{ for every }E\},
\qquad
\mathsf{LBT}^{\mathrm{Kell}}_k
\subset\mathsf{LBT}^{\mathrm{bal}}_k,
\tag{1.7}
\]

where the second contains packages actually realized by a polynomial Keller
map and its canonical finite normalization.  Membership in
\(\mathsf{LBT}^{\mathrm{bal}}_k\) is a divisor calculation.  Membership in
\(\mathsf{LBT}^{\mathrm{Kell}}_k\) additionally requires polynomial
reconstruction, an affine-space source, and the correct affine-sheet open.

The automatic tame identity

\[
(e-1)+1-e=0
\tag{1.8}
\]

for a finite reduced log pair is not the ledger condition (1.6).  Equation
(1.8) explains codimension-one log crepancy; (1.6) is the source/core/target
Jacobian cancellation which must survive polynomial algebraization.

### 1.4 Cell coordinates

For a connected selected normalization define its coarse discrete signature

\[
\sigma(\mathfrak B)=
\left(
g,\ |P|,\ \operatorname{rank}\mathcal O(\widetilde C^\circ)^*/k^*,\
\delta,\ \operatorname{rank}M_{\mathcal R},\
|\!I\!|,\ \{(e,f,\operatorname{diff},\epsilon)\},\
\ell
\right).
\tag{1.9}
\]

A **classification cell** is a full subgroupoid cut out by a fixed discrete
signature together with fixed primitive valuation, conductor, target-image,
and affine-sheet incidence types.  The term “cell” does not assert that its
coarse moduli space is an affine cell.  Continuous seed or curve moduli may
remain inside it.

For a smooth rational normalization,

\[
\operatorname{rank}\mathcal O(\widetilde C^\circ)^*/k^*
=|P|-1.
\tag{1.10}
\]

Outside the smooth rational case, genus and conductor data are independent
and must not be inferred from the unit rank.

## 2. First classification cells

### Cell W: weighted, one puncture

This is the first cell.

\[
\overline C=\mathbb P^1_W,\qquad
P=\{\infty\},\qquad
\widetilde C^\circ=\mathbb A^1_W,\qquad
\mathcal O(\widetilde C^\circ)^*/k^*=0.
\tag{2.1}
\]

The primitive affine-sheet mark is the coordinate \(W\), not a boundary
unit.  In a coordinate-preserving plane core it gives

\[
(W,q)\longmapsto(q,Wq-H(W)),\qquad H'=h,
\tag{2.2}
\]

once the reduced critical section \(q=h(W)\), the preserved coordinate, and
the primitive conormal are marked.  This normal form is the categorical
definition of the weighted core cell; admissible seed polynomials \(H\) give
its polynomial Keller realizations.

For the weighted suspension the generic discriminant prime has

\[
(e,f,\operatorname{diff})=(2,1,1).
\tag{2.3}
\]

The distributed ledger has the characteristic selected row

\[
v_\gamma(J_\alpha)+v_\gamma(\gamma)
=2+1=3
=v_\gamma(J_\beta\circ F).
\tag{2.4}
\]

The target discriminant image can have cusps and nodes.  Those are conductor
strata inside Cell W, not new punctures of its normalization.  The exact
conductor is retained by (1.4).

**Status:** realized, in all admissible weighted degrees.  The marked
one-place core normal form and the suspension polynomiality theorem are in
[Log geometry of controlled-boundary suspensions](LOG_GEOMETRY_OF_SUSPENSIONS.md)
and the [weighted-seed theorem](../verified/WEIGHTED_SEED_THEOREM.md).

### Cell R: reciprocal, two punctures

This is the second cell.

\[
\overline C=\mathbb P^1_Y,\qquad
P=\{0,\infty\},\qquad
\widetilde C^\circ=\mathbb G_m,\qquad
\mathcal O(\widetilde C^\circ)^*/k^*\simeq\mathbb Z.
\tag{2.5}
\]

The primitive unit mark has divisor

\[
\operatorname{div}_P(Y)=(1,-1),
\tag{2.6}
\]

and the cancellation subcell additionally retains

\[
\operatorname{div}_P(s)=(-m,m),\qquad
D=1-s(q-ps)^m.
\tag{2.7}
\]

For controlled exponent \(r\), the finite discriminant ramification data is

\[
(e,f,\operatorname{diff})=(r+1,1,r).
\tag{2.8}
\]

At the reciprocal source boundary \(A=D^{-1}\), the reconstruction valuation
row is

\[
v_A(s,Y,P,B,D)=(-1,0,1,0,-1),
\tag{2.9}
\]

and the ledger is the zero--pole cancellation

\[
v_A(J_\alpha)+r\,v_A(D)=-r+r=0.
\tag{2.10}
\]

The completed polynomiality equations force the reciprocal integral core and
the global slice; they are part of realization, not consequences of the
puncture count.

**Status:** the primitive cancellation subcell is realized.  The exact
marked-link classification is in
[Log geometry of controlled-boundary suspensions](LOG_GEOMETRY_OF_SUSPENSIONS.md)
and the executable recognition procedure is in the
[reciprocal-link classifier](RECIPROCAL_LINK_CLASSIFIER.md).

There is at least one further realized two-puncture cell.  The
root-engineered quadratic-gauge family has the same coarse curve and
ramification signature, while its extracted chart mark is \(S^2\), not
\(S\).  It therefore lies outside Cell R even though its first seven
minimal-boundary predicates pass.  This is the first proof that the category
must retain reconstruction markings and cannot classify by
\((g,|P|,e,f)\) alone; see the
[operational invariant pipeline](MINIMAL_BOUNDARY_INVARIANT_PIPELINE.md).

### Three-puncture rational cells

Here

\[
\overline C=\mathbb P^1,\qquad
P=\{0,1,\infty\},\qquad
\Lambda_P\simeq\mathbb Z^2
\tag{2.11}
\]

with saturated basis

\[
\operatorname{div}(Y)=(1,0,-1),\qquad
\operatorname{div}(Y-1)=(0,1,-1).
\tag{2.12}
\]

Two cells are currently separated.

1. **Direct two-centre reciprocal cell.** Its controlled character is
   \[
   (a,b,-a-b),\qquad a,b\ge1,
   \tag{2.13}
   \]
   and its reciprocal determinant ledger is balanced for every \(a,b,r\).
   The first polynomiality moment nevertheless fails for every such triple.
   This entire cell is obstructed within the stated one-reconstruction,
   identity-target chart.

2. **Double-incidence cell.** It has two independently selected critical
   primes \(D_0,D_1\), multiplicities \(a,b\), and intersection curve
   \(\mathbb P^1\setminus\{0,1,\infty\}\).  Its polynomial core satisfies
   \[
   \operatorname{Jac}\Phi=-D_0^aD_1^b.
   \tag{2.14}
   \]
   Reciprocal cancellation on the boundary complement is exact, but that
   complement has two nonconstant units and is not affine space.  A
   polynomial affine-space completion remains open after the proved
   low-degree and rank-drop screens.

The obstruction and the surviving cell are respectively in the
[puncture-rank frontier](PUNCTURE_RANK_FRONTIER.md) and the
[three-puncture double-incidence core](../extended-geometry/THREE_PUNCTURE_DOUBLE_INCIDENCE_CORE.md).
There is no current theorem saying that these are all three-puncture cells.

## 3. Conductor refinements: nodal and cuspidal images

A singular critical image is classified by its normalization and conductor
square, not by the words “node” or “cusp” alone.

### Rational one-component normalization

The basic conductor models in \(B=k[t]\) are

\[
\begin{array}{c|c|c|c}
\text{image}&A\subset B&\mathfrak cB&
\text{normalization fibre}\\ \hline
\text{node}&k+t(t-1)k[t]&t(t-1)B&\{0,1\}\\
\text{cusp}&k[t^2,t^3]&t^2B&\operatorname{Spec}k[t]/(t^2).
\end{array}
\tag{3.1}
\]

Both have \(\delta=1\), but their finite conductor maps are different.

- The cuspidal conductor occurs in a realized Cell-W object: starting from
  \(k[u^2,u^3]\subset k[u]\) reconstructs the foundational cubic weighted
  map and its distributed ledger.  It is a conductor stratum of the weighted
  cell, not a new construction mechanism.
- For both the node and cusp, the separated chart obtained by inverting only
  the conductor and adjoining affine variables is obstructed: the required
  reconstruction pole contradicts polynomiality, and the localization has
  nonconstant units.
- The symmetric three-boundary Cox fill makes the reconstruction coordinates
  polynomial but leaves one dualizing-form pole.  Its normalized nodal fill
  is smooth but not affine space; its cuspidal fill is normal and singular.

These statements are proved in the
[conductor-first cusp realization](../extended-geometry/CONDUCTOR_FIRST_FOUNDATIONAL_CUSP_KELLER.md),
[one-chart obstruction](../extended-geometry/CONDUCTOR_FIRST_ONE_CHART_OBSTRUCTION.md),
and [three-boundary Cox-fill obstruction](../extended-geometry/CONDUCTOR_THREE_BOUNDARY_COX_FILL_OBSTRUCTION.md).
An ambient-coupled nodal Keller realization is not currently known.

### Nodal \(\mathbb G_m\) normalization

If the normalization is \(\mathbb G_m\) and \(1\) is identified with \(-1\),
then a unit \(ct^m\) descends precisely when \(m\) is even.  The conductor
therefore changes the unit lattice from \(\mathbb Z\) to the index-two
sublattice \(2\mathbb Z\) without changing its rank.  This nodal rational
package passes the stage-one conductor, unit, and adjunction gates when its
involution preserves the unordered conductor pair.  Passing is only a
feasibility result; affine-space realization is open.  See the
[boundary-package compiler](../extended-geometry/BOUNDARY_PACKAGE_COMPILER.md).

## 4. Genus-one normalization

The lowest smooth positive-genus cell presently retained is

\[
g(\overline C)=1,\qquad |P|=1,\qquad
\mathcal O(\widetilde C^\circ)^*/k^*=0,\qquad
\mathfrak c=\mathcal O_C.
\tag{4.1}
\]

It is not Cell W: both have scalar units and one puncture, but their
log-canonical degrees are respectively

\[
\deg(K_{\overline C}+P)=1\quad\text{and}\quad -1.
\tag{4.2}
\]

The first exact benchmark is a degree-three cover of \(\mathbb P^1\) with
six simple branch values.  Its paired transpositions generate \(S_3\), and
Riemann--Hurwitz gives

\[
2g-2=3(-2)+6=0.
\tag{4.3}
\]

A plane-cubic adjunction row and the one-puncture unit calculation also pass.
This defines a nonempty stage-one genus-one cell of
\(\mathsf{LBT}^{\mathrm{bal}}_k\) at the abstract-package level.  No root
equation, reconstruction open isomorphic to affine space, or polynomial
Keller realization is known.  Its status is therefore **open**, not
realized; see the
[boundary-package compiler](../extended-geometry/BOUNDARY_PACKAGE_COMPILER.md).

## 5. Multiple selected critical components

For \(|I|>1\), componentwise signatures do not determine the type.  The
category retains:

\[
\left\{
\text{component normalizations},\
\text{intersection strata},\
\text{target-image partition},\
\text{valuation-lattice span},\
\text{sheet colors}
\right\}.
\tag{5.1}
\]

The first low-complexity cells are:

| cell | selected components | coupling | current status |
|---|---:|---|---|
| double incidence | \(2\) | two primitive puncture characters meet on a smooth three-puncture curve | balanced polynomial core; affine Keller completion open |
| tetrahedral degree four | \(2\) index-two primes | same branch divisor, three-puncture selected curve | abstract stage one passes; exact symbolic replay exposes an affine-colored ramified prime |
| Fano degree seven | \(2\) index-two primes plus three unramified primes | one involution branch in the degree-seven action | abstract stage one passes; realization open |

The last two rows are feasibility fixtures, not constructions.  They are
included because they show why a selected set \(I\), target-image partition,
and affine colors must be categorical data.  Collapsing them to a single
“ramification index two” label loses the obstruction.

## 6. Census and frontier

The first census is:

| order | cell | \((g,|P|)\) | conductor | ledger | realization status |
|---:|---|---:|---|---|---|
| 1 | W: weighted | \((0,1)\) | variable node/cusp strata allowed | distributed positive \(2+1=3\) | realized |
| 2 | R: reciprocal cancellation | \((0,2)\) | zero on the smooth core | reciprocal \(-r+r=0\) | realized |
| 3 | QG: quadratic gauge | \((0,2)\) | zero on the smooth core | reciprocal, quadratic mark | realized; not Cell R |
| 4 | 3R: direct two-centre reciprocal | \((0,3)\) | zero | balanced | polynomiality obstructed in the complete stated chart |
| 5 | 3DI: double incidence | \((0,3)\) | zero on the selected intersection | two-factor balanced core | affine-space completion open |
| 6 | N/C: conductor refinements | normalization genus \(0\) | nonzero finite conductor map | separated ledger obstructed; weighted cusp realized | mixed by subcell |
| 7 | E1: elliptic, one puncture | \((1,1)\) | zero in the smooth benchmark | stage-one balanced | realization open |
| 8 | MC: multiple selected components | componentwise | component and intersection conductors | vector ledger | mixed; no completeness theorem |

This table is a classification of the low-complexity packages presently
supported by exact repository results.  It is not an exhaustiveness theorem
for all objects of \(\mathsf{LBT}_k\).

The smallest classification questions are now precise.

1. Does every realized smooth rational one-puncture package with primitive
   affine mark lie in Cell W, after marked polynomial left--right
   equivalence?
2. Which primitive or higher-incidence marks split the smooth rational
   two-puncture locus beyond Cells R and QG?
3. Can any three-puncture or multiple-component balanced core acquire an
   affine-space reconstruction open without collapsing its valuation lattice?
4. Can a genus-one stage-one package satisfy the actual value-semigroup,
   affine-sheet, and polynomial reconstruction gates?
5. Which conductor squares admit a distributed ledger genuinely different
   from a conductor stratum of Cell W?

Questions 1 and 2 refine the existing minimal-boundary gateway.  Questions
3--5 are the first new classification directions exposed by this category.

## 7. Exact regressions

No new computation is claimed here.  The cells are replayed by the existing
targeted commands:

```bash
.venv/bin/python scripts/verify_minimal_boundary_pipeline.py
.venv/bin/python scripts/verify_puncture_rank_frontier.py
.venv/bin/python scripts/verify_three_puncture_double_incidence.py
.venv/bin/python scripts/verify_conductor_first_one_chart_obstruction.py
.venv/bin/python scripts/verify_conductor_first_foundational_cusp_keller.py
.venv/bin/python scripts/verify_conductor_three_boundary_cox_fill.py
python3 scripts/verify_boundary_package_compiler.py
```

The category is assembled from the intrinsic finite-normalization layers in
[the decorated-normalization invariant](../extended-geometry/DECORATED_NORMALIZATION_INVARIANT.md),
the signed relative-canonical ledger in
[Log geometry of controlled-boundary suspensions](LOG_GEOMETRY_OF_SUSPENSIONS.md),
and the exact stage-one schema in the
[boundary-package compiler](../extended-geometry/BOUNDARY_PACKAGE_COMPILER.md).
