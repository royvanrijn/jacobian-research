# Tate--Jacobian tests from the explicit Keller maps

This note separates three notions which should not be conflated:

1. a two-sided inverse in a global Tate algebra;
2. an analytic inverse section on one target residue tube; and
3. a chosen algebraic inverse sheet on a maximal rigid domain.

For the foundational Keller map, the first fails at every odd prime, the
second has an exact `3/1/0` classification at every prime at least five, and
the centered inverse at the origin has exact isotropic radius one but does
not belong to the closed-unit Tate algebra.  The third question reduces
effectively to factorization of the inverse polynomial over an affinoid
algebra, but a complete branch-by-branch Berkovich-domain classification is
still open.

## 1. Relation with the 2025 Tate--Jacobian theorem

Hamada, Kato, and Komiya define the Tate--Jacobian conjecture for an
`I`-adically topologized ring and prove:

- if the topology is Hausdorff, `TJC(R,I,n)` is equivalent to
  `JC(R/I,n)`;
- under the stated characteristic-zero residue assumptions, the
  all-dimensional Tate and classical conjectures are equivalent; and
- `JC(C)` is equivalent to the assertion that, for all but finitely many
  primes `p`, every Keller map over `C_p` has a two-sided inverse in
  `C_p<X_1,...,X_n>`.

The last assertion concerns one global restricted power series inverse on
the closed unit polydisc.  It is much stronger than the existence of local
inverse branches.  See
[Hamada--Kato--Komiya, *A Tate algebra version of the Jacobian
Conjecture*](https://doi.org/10.1016/j.jpaa.2025.108129), especially
Theorem 0.4.

The [foundational map](../verified/FOUNDATIONAL_GEOMETRY.md) is defined over
`Z`, fixes the origin, and has determinant `-2`.  Its three integral source
points

\[
 (0,0,-1/4),\qquad (1,-3/2,13/2),\qquad(-1,3/2,13/2)
\]

have the common integral target `(-1/4,0,0)` at every odd prime.  If its
formal inverse at the origin were a two-sided element of
`C_p<A,B,C>^3`, evaluation of `G(F(X))=X` at these three points would give
three different values for `G(-1/4,0,0)`.  Hence:

> **Global Tate obstruction.** For every odd prime `p`, the formal inverse
> at the origin of the foundational map is not a two-sided inverse in
> `C_p<A,B,C>^3`.

Thus one fixed integral map witnesses failure of the global Tate conclusion
at every odd prime.  This is a direct test of the global formulation, not
merely a failure at a bad reduction prime.

## 2. Residue tubes and a general local lemma

Let `K=C_p`, with valuation ring `O`, maximal ideal `m`, and residue field
`k=overline(F_p)`.  For `bar y in k^n`, write

\[
 ]\bar y[=\operatorname{red}^{-1}(\bar y)\subset O^n.
\]

This is an open unit polydisc.  It is not itself the affinoid closed unit
polydisc: a function on the whole tube means a compatible analytic function
on every centered closed polydisc of radius `rho<1`.  A `Z_p` residue ball
`y_0+pZ_p^n`, by contrast, becomes a closed radius-`p^{-1}` polydisc after a
choice of lift and is described by an ordinary rescaled Tate algebra.

We use the following standard multivariate Hensel consequence.

> **Residue-tube lemma.** Let `F in O[X]^n` and suppose
> `det JF in O^*`.  For every
> `bar x in k^n` with `bar F(bar x)=bar y`, restriction gives an analytic
> isomorphism
> \[
>       F:]\bar x[\;\xrightarrow{\sim}\;]\bar y[.
> \]
> In particular, every geometric special-fiber preimage produces one
> integral analytic inverse branch on the complete target tube.

Indeed, choose lifts, apply Hensel to obtain an exact lifted point, and
translate it to the origin.  The linear term is in `GL_n(O)`, while every
higher term is strictly contracting on `m^n`.  Successive approximation is
unique and analytic in the target variables.  The same argument over each
closed radius `rho<1` gives the compatible Tate representatives.

This lemma classifies integral tube branches by the geometric special fiber.
It does not say that a branch descends to `Q_p`: descent requires the
corresponding special-fiber point to be `F_p`-rational.

## 3. Exact tube classification for the foundational cubic

For a target `(a,b,c)`, put

\[
 P_{a,b,c}(T)=cT^3-2T^2+bT-2a
\]

and

\[
 Q=27a^2c^2-18abc+16a+b^3c-b^2.
\]

The exact identities from
[Image and nonproperness](../verified/IMAGE_AND_NONPROPERNESS.md) are

\[
 \operatorname{Disc}_T(P)=-4Q
\]

and

\[
 \Gamma=V(3bc-4,\;12a-b^2).
\]

For a finite simple root `t`, the inverse branch is

\[
\begin{aligned}
 x&=\frac2{P'(t)},\\
 y&=t-\frac{P'(t)}2,\\
 z&=\frac54P'(t)^2-\frac32tP'(t)-\frac c8P'(t)^3.
\end{aligned}
\tag{3.1}
\]

The projective root at `T=infinity` is regular on the chart `s=1/T`:

\[
 c=2s-bs^2+2as^3,\qquad
 D=1-bs+3as^2,
\tag{3.2}
\]

\[
 x=\frac{s}{D},\qquad y=b-3as,
\tag{3.3}
\]

and `z` is a polynomial in `(a,b,s)` with value `a-4b^2` at `s=0`.

Let `p>=5`.  Reducing these formulas over `k=overline(F_p)` gives the exact
geometric fiber table

| reduced target | affine special-fiber points | integral analytic branches on the tube |
|---|---:|---:|
| `bar Q != 0` | 3 | 3 |
| `bar Q = 0`, `bar y notin bar Gamma` | 1 | 1 |
| `bar y in bar Gamma` | 0 | 0 |

The proof is elementary.  Off `Q=0`, the projective cubic has three simple
roots, with (3.2)--(3.3) supplying the root at infinity when `c=0`.  On
`Q=0` away from `Gamma`, it has one simple root and one repeated root; only
the simple root reconstructs.  On `Gamma`,

\[
 P_{a,b,c}(T)
 =c\left(T-\frac{2}{3c}\right)^3,
\quad
 (a,b)=\left(\frac4{27c^2},\frac4{3c}\right),
\]

so no root reconstructs.  Each residue point of `bar Gamma` lifts to this
characteristic-zero omitted curve by lifting `c in k^*`; consequently its
whole tube cannot have an inverse section.  The other rows follow from the
residue-tube lemma.

This answers the “thin set” question in the geometric `C_p` sense:

> **All but a certified thin set.** The globally noninvertible foundational
> map has an integral analytic inverse branch on every geometric integral
> residue tube except the codimension-two family `bar Gamma`; it has three
> branches off the boundary discriminant and one on its complement
> `V(bar Q)-bar Gamma`.

There is an important arithmetic qualification.  A branch over a
`Z_p`-residue ball exists over `Q_p` only when the corresponding root is
`F_p`-rational.  The generic inverse polynomial has `S_3` monodromy, so a
positive proportion of rational residue targets have no rational root.
Thus “all but a thin set” is false over `Q_p` without extending scalars to
`C_p`, even though it is true for geometric residue tubes.

## 4. Exact centered polydisc radii

For a chosen inverse germ `g` at a target `y_0`, define its isotropic section
radius to be the supremum of `r` such that `g` converges and
\(F\circ g=\mathrm{id}\) on every centered closed polydisc of radius
`rho<r`.  This definition records a possible nonattained boundary radius.

### 4.1 The formal inverse at the origin

At `(a,b,c)=(0,0,0)`, the source point is the projective-root branch `s=0`.
For

\[
 |a|,|b|,|c|<1,
\]

equation (3.2), written as

\[
 s=\frac{c+bs^2-2as^3}{2},
\]

is strictly contracting on `m`; moreover `D=1-bs+3as^2` is a unit.
Equations (3.2)--(3.3) therefore give the inverse germ on the entire open
unit polydisc.

On the other hand,

\[
 y_\Gamma=(4/27,4/3,1)
\]

lies on the omitted curve `Gamma` and lies in the centered closed unit
polydisc for every `p>=5`.  No inverse section can be evaluated there.
Therefore:

\[
 \boxed{R_{\mathrm{origin}}=1,\quad\text{not attained}.}
\]

Equivalently, the inverse is analytic on every smaller closed polydisc but
is not a member of the closed-unit Tate algebra.  Notice that the origin
itself lies on `Q=0`.  The raw discriminant therefore has distance zero
from the center, yet the surviving affine sheet has radius one.

### 4.2 The three branches at the rational collision

Put `y_* = (-1/4,0,0)`.  Since `Q(y_*)=-4`, every `p>=5` has three inverse
branches throughout

\[
 \max(|a+1/4|,|b|,|c|)<1.
\]

The same omitted point `y_Gamma` is at sup distance exactly one from `y_*`.
It prevents every branch from extending to the complete closed unit
polydisc.  Hence all three branches have the same exact centered radius:

\[
 \boxed{R_1=R_2=R_3=1,\quad\text{not attained}.}
\]

For the two finite-root branches, the one-parameter slice `b=c=0` makes the
boundary explicit:

\[
 t=\mathord\pm\frac12\sqrt{1-4(a+1/4)}.
\]

The collision of these two roots occurs at distance one.  The projective
branch is controlled by (3.2); its obstruction on the full polydisc is the
same omitted boundary.

### 4.3 What the boundary discriminant actually controls

The correct rule is sheet-sensitive.

- `V(Q)` is the locus where the full three-sheet affine cover ceases to be
  finite.  Its distance controls the largest centered polydisc on which all
  three sheets can split simultaneously.
- At a smooth point of `V(Q)`, two sheets escape but the simple-root sheet
  remains analytic.  That sheet can cross the discriminant.
- `Gamma` is the triple-root/omitted locus.  It obstructs every affine
  section and bounds the origin branch even though that branch starts on
  `V(Q)`.

Thus “radius equals distance to the discriminant” is correct only after one
specifies which discriminant or boundary stratum belongs to the chosen
sheet.

## 5. Effective criterion for the other explicit families

Let `D` be a closed target polydisc with affinoid algebra `A_D`.  On the
finite reconstruction chart, a branch of an explicit marked-root Keller map
over `D` is equivalent to a linear factor of its inverse polynomial over
`A_D`, together with invertibility of the reconstruction derivative.

For the weighted family this is

\[
 E_{A,B,C}(W)=H(W)-BCW+cAC^2,
\qquad E'(W)\in A_D^*,
\]

on `C!=0`.  For the quadratic-gauge family it is

\[
 E_{\Pi,B,C}(S)=0,\qquad
 \partial_SE_{\Pi,B,C}(S)\in A_D^*,
\]

on `Pi!=0`.  The special `C=0` and `Pi=0` charts must be retained rather
than discarded by denominator clearing.  Their exact algebraic fiber
tables are already proved in
[the weighted-seed theorem](../verified/WEIGHTED_SEED_THEOREM.md) and
[the quadratic-gauge nonproperness paper](../papers/quadratic-gauge-nonproperness/README.md).

At primes of good reduction, the residue-tube lemma transfers those tables
directly:

- each affine geometric special-fiber point yields one integral analytic
  branch on the target tube;
- repeated inverse roots with vanishing reconstruction derivative are lost
  boundary sheets; and
- a target tube containing a lifted omitted point admits no global section.

This is a finite algorithm: factor the reduced inverse polynomial, retain
the simple reconstructing roots in every affine/projective chart, and then
use Hensel lifting.  For a larger affinoid polydisc, replace residue
factorization by factorization in `A_D` and test the derivative in
`A_D^*`.

## 6. Sharp local formulations already decided

The foundational map gives the following exact outcomes.

1. “A local inverse at the origin of a Keller map lies in the closed-unit
   Tate algebra” is false at every odd prime.
2. “A Tate inverse on one residue tube, or on all geometric residue tubes
   outside codimension two, forces global invertibility” is false.
3. “A chosen inverse branch must stop at the boundary discriminant” is
   false: the surviving branch crosses `V(Q)-Gamma`.
4. “Every geometric target residue tube has an integral inverse branch” is
   false, sharply on `bar Gamma`.
5. “Every `Q_p` residue class outside an algebraic thin set has a
   `Q_p`-valued branch” is false because Frobenius fixed-point obstructions
   occur on a positive-density set.

## 7. Remaining research problems

The following are not proved here.

- Determine the maximal connected Berkovich or rigid domain of each
  individual cubic branch, not merely its maximal centered isotropic
  polydisc.
- Compute branch radii at every target stratum and at the bad primes `2`
  and `3`.
- For every weighted and quadratic-gauge seed, identify the exact
  no-section locus rather than only the nonproperness hypersurface and the
  available special-slice fiber tables.
- Express the maximal radius directly from the cluster picture/Newton
  polygon of the inverse polynomial, including leading-coefficient drops.
- Quantify separately the geometric `C_p` tube cover and the arithmetic
  `Q_p`-rational tube cover using the certified `S_N` Frobenius statistics.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_padic_inverse_branches.py
```

The checker verifies the determinant, both inverse charts, the discriminant
identity, the triple-root omitted curve, the integral collision, the radius
boundary witness, and finite-field regressions for the rational-root/fiber
formula at `p=5,7,11`.  The uniform residue-tube statements use the
multivariate Hensel lemma and the exact algebraic fiber classification above;
the finite-field samples are regressions, not their proof.
