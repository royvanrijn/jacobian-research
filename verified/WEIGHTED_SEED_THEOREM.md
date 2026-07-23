# Universal weighted-seed inverse and monodromy theorem

This note isolates the part of the weighted construction that is independent
of the special quartic boundary geometry. It supplies one theorem schema for
inverse geometry, monodromy, and finite-field statistics; exact affine images
and nonproperness strata require additional boundary hypotheses.

## Universal pencil theorem

Let `k` be a characteristic-zero field, let `H in k[W]` have degree `n>=2`,
and put

\[
E(W)=H(W)-sW+t.
\]

Then:

1. `E` is irreducible over `k(s,t)`.
2. Its repeated-root discriminant is an irreducible curve normalized by
   \[
   r\longmapsto (s,t)=\bigl(H'(r),rH'(r)-H(r)\bigr).
   \]
3. The normalization is generically one-to-one.
4. The geometric and arithmetic monodromy groups of the generic degree-`n`
   cover are both `S_n`.

### Connectedness

The polynomial `E` has degree one in `t`.  In a factorization in
`k[s,W,t]`, one factor has `t`-degree zero.  Writing the other as `B_0+B_1t`,
comparison of the coefficient of `t` gives `AB_1=1` in `k[s,W]`; hence both
`A` and `B_1` are units and the factorization is trivial.  Since the leading
coefficient of `E` as a polynomial in `W` is the nonzero constant leading
coefficient of `H`, `E` is primitive over `k[s,t]`.  Gauss's lemma therefore
gives irreducibility in `k(s,t)[W]`.  After base change to an algebraic
closure of `k`, the same argument applies verbatim.  The geometric cover is
connected, so its monodromy group is transitive on the `n` roots.

### Discriminant normalization

A repeated root `r` satisfies

\[
0=E'(r)=H'(r)-s,
\qquad
0=E(r)=H(r)-rH'(r)+t,
\]

which gives the displayed parameterization. Its coordinate degrees are
`n-1` and `n`. If its generic degree were `e`, then `e` would divide both pole
orders at infinity. Since `gcd(n-1,n)=1`, one has `e=1`. Thus the
parameterization is birational and the discriminant curve is irreducible.

For generic `r`, `H''(r)!=0`, so `r` is exactly a double root.  Birationality
of the normalization also says that the generic discriminant polynomial has
no second repeated root.  A transverse parameter has local form

\[
\epsilon+c(W-r)^2+\text{higher terms},\qquad c\ne0,
\]

and the corresponding quadratic tame local extension exchanges the two
nearby roots.  Generic discriminant inertia is therefore a transposition.

### Full symmetric monodromy

The proof and its precise universal quantifier are isolated in the
[standalone universal-monodromy note](UNIVERSAL_SYMMETRIC_MONODROMY.md).
For every fixed `H`, all but finitely many linear tilts `H(W)-s_0W` are Morse:
the normalization above is birational, so a generic vertical line has `n-1`
simple, pairwise distinct branch values.  The classical Morse-polynomial
theorem then gives `S_n` monodromy on that vertical line, forcing the ambient
two-parameter geometric group to be `S_n`; the arithmetic group contains it
and is therefore also `S_n`.

The standalone note cites Serre's standard branch-cycle theorem, treats
monomials, Chebyshev polynomials, decomposable polynomials, affine symmetries,
and repeated critical values explicitly, and records the former
inertia-and-purity proof as an alternative.

## Consequence for admissible weighted seeds

For the weighted construction, take `p=H'` and impose

\[
p(0)=0,
\quad p(1)=-c,
\quad H(0)=H(1)=0,
\quad c\ne0,
\quad p'(1)/c\ne-2.
\]

These conditions make the affine weighted lift polynomial.  The
[tangent-map core theorem](TANGENT_MAP_CORE.md) applies directly to the plane map

\[
 \Phi_H(W,\gamma)=
 \bigl(H'(W)+c\gamma,\,W(H'(W)+c\gamma)-H(W)\bigr),
 \qquad \det D\Phi_H=-c^2\gamma.
\]

and proves at once its inverse pencil, critical divisor, discriminant
normalization, Hessian Fitting divisor, and reconstruction pole.  It also
identifies the weighted threefold as the suspension

\[
 (BC,cAC^2,C)=(\Phi_H(W,\gamma),C).
\]

Its determinant conclusion is therefore `det DG_H=b_0c`; no separate
three-variable Jacobian calculation is needed here.

Its inverse pencil is

\[
E_{A,B,C}(W)=H(W)-BCW+cAC^2.
\]

On `C!=0`, differentiation gives

\[
\gamma=-E'(W)/c,
\qquad
x=-cC/E'(W).
\]

Thus simple roots reconstruct finite source points and repeated roots are
exactly the reconstruction poles. The universal theorem gives `S_n`
monodromy for every such degree-`n` inverse pencil, including canonical and
deformed seeds.

This statement is deliberately on `C!=0`.  Globally, normalize the finite
marked-root incidence and retain the open on which every reconstruction
coordinate is regular.  This includes the finite `x=0` and `gamma=0` charts
without misclassifying extra primitive-root branches at `C=0`.  The expanded
two-chart derivation is retained in the
[core-support archive](../archive/core-support/WEIGHTED_MARKED_ROOT_MODEL.md).

## What is universal and what remains conditional

The following data are now universal:

- connected generic inverse cover;
- irreducible repeated-root discriminant and its normalization;
- generic transposition ramification;
- geometric and arithmetic monodromy `S_n`;
- no nonidentity deck transformation of the generic inverse cover for
  `n>=3`;
- identification of repeated roots with reconstruction poles on `C!=0`.

Indeed the sheet stabilizer in the natural `S_n` action is `S_(n-1)`.  For
`n>=3` it is self-normalizing, so the deck group
`N_(S_n)(S_(n-1))/S_(n-1)` is trivial.  This removes generic target-fixed
sheet symmetries, but it does not by itself extend the affine root sheet
through collision strata.  That separate valuative problem is complete as
`H3`: the coarse affine mark extends uniquely over every DVR limit.

The following still depend on the seed and must not be inferred from
monodromy alone:

- the exact affine image and omitted locus;
- the full `C=0` fiber decomposition;
- self-intersections and higher singularities of the discriminant curve;
- whether extra zeros of `H` introduce additional boundary branches;
- radical decompositions of the nonproperness hypersurface.

Their generic behavior is proved by the projective-dual theorem in
`GENERIC_DISCRIMINANT_CURVE.md`: a uniform contact-incidence argument leaves
only ordinary cusps and bitangent nodes in every inverse degree, while special
seeds retain the seed-dependent higher collisions listed here. Exact
certificates through inverse degree ten provide computational regressions.

For the canonical family `H_d(W)=W^d(1-W)`, the only primitive zeros are the
distinguished points `0` and `1`. Its global geometry is now proved in
[the canonical-family image theorem](../extended-geometry/CANONICAL_FAMILY_IMAGE.md).
Deformed seeds should be treated afterward because their extra primitive
zeros create additional `C=0` branches.  The earlier one-extra-zero derivation
is retained in
[archive/geometry-support](../archive/geometry-support/DEFORMED_SEED_BOUNDARY.md).

## Finite-field interpretation

At primes of good reduction, factorization and rational-root statistics of
`E` are governed by Frobenius conjugacy classes in `S_n`. In particular, the
number of simple rational roots is the number of fixed points of the natural
permutation. This explains the histograms in the bounded seed scan.

The `S_n` theorem predicts the large-field distribution through function-field
Chebotarev; it does not by itself prove an exact formula for every finite
field. Bad characteristics, boundary fibers, and exact error terms remain a
separate arithmetic step.

The good-reduction asymptotic, its `S_n` fixed-point probabilities, and the
transfer from the pencil to the three-dimensional target are proved in
[the finite-field theorem](../extended-geometry/FINITE_FIELD_CHEBOTAREV.md).

### Clean-room reproduction

A second derivation uses dependency-free sparse-polynomial arithmetic and an
alternative vertical-line branch-cycle argument.  It rechecks polynomiality,
the constant Jacobian, marked-root reconstruction, normalization, and `S_n`
monodromy without importing the project weighted-model implementation.  Run

```bash
python3 scripts/audit_weighted_independent.py
```

The detailed former audit narrative is retained in
[archive/core-support](../archive/core-support/README.md).  The public command
sequence is in [REPRODUCE.md](../REPRODUCE.md).

Christopher D. Long's later GMC and Mathieu-conjecture papers are external
consequences motivated by JC(3), not reviews of this weighted-seed theorem and
not derivations from its inverse pencil.  Motivated by Long's Lagrange--Good
search architecture, a subsequent repository theorem now does turn every
nonconstant normalized inverse seed into a uniform four-real-Gaussian witness
family; this is a new internal bridge, not a reinterpretation of Long's direct
formulas.  See the [weighted Gaussian bridge](../extended-geometry/WEIGHTED_GAUSSIAN_BRIDGE.md)
and the exact provenance split in the
[external-consequences note](../extended-geometry/EXTERNAL_CONSEQUENCES_AND_PROVENANCE.md).
