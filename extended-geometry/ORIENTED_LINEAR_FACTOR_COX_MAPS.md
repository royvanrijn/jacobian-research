# All-arity oriented Cox maps and their dicritical divisors

The oriented-discriminant target absorbs the complete Jacobian unit for every
number of ordered linear factors.  A normalization tree determines which
collision divisors remain affine and which become dicritical.  The result is
an explicit family of polynomial constant-Jacobian morphisms on smooth
Cox--Kummer charts with arbitrarily many distinct dicritical boundary
valuations.

Work over a characteristic-zero field `k`.

## 1. Normalization tree without removing the other collisions

Fix `s>=3` and a bipartite spanning tree `G` chosen as in the
[all-arity ledger theorem](COX_LEDGER_LINEAR_FACTORS.md), with

\[
 d_s=|\det M_G|=
 \begin{cases}
 1,&s\text{ odd},\\
 2,&s\text{ even}.
 \end{cases}                                         \tag{1}
\]

For

\[
 L_i=u_iX+v_iY,\qquad
 P=\prod_iL_i=\sum_{j=0}^sp_jX^{s-j}Y^j,
\qquad r_{ij}=u_iv_j-v_iu_j,
\]

choose the irreducible hyperplane

\[
 m=p_0+p_s.
\]

Define

\[
 \overline Y_{s,G}=
 \{r_e=1\ (e\in E(G)),\quad m=1\}
 \subset\mathbb A^{2s}.                              \tag{2}
\]

Unlike the reciprocal Cox suspension, (2) does not invert the non-tree
resultants.  Their collision divisors are part of the affine source.

The selected boundary matrix consists of the tree-edge columns and `E`; it
has determinant `+-d_s`.  Thus (2) is the normalized Cox slice when `s` is
odd and its canonical `mu_2` Kummer slice when `s` is even.  It is a smooth
affine `s`-fold.

## 2. Oriented target and constant Jacobian

Let

\[
 \widetilde T_s=
 \{m(P)=1,\quad D^2=\operatorname{Disc}(P)\}.         \tag{3}
\]

Its smooth locus contains the generic simple-double-root divisor.  All
statements below may be restricted to any affine principal open of that
smooth locus meeting this generic divisor.

The binary-form discriminant identity is

\[
 \operatorname{Disc}(P)
 =\left(\prod_{i<j}r_{ij}\right)^2.                  \tag{4}
\]

On (2), the tree resultants equal one.  Hence the primitive orientation

\[
 D_G=\prod_{f\notin E(G)}r_f                         \tag{5}
\]

satisfies `D_G^2=Disc(P)`.  Define

\[
 \boxed{
 \Psi_{s,G}:\overline Y_{s,G}\longrightarrow
 \widetilde T_s,\qquad
 (L_1,\ldots,L_s)\longmapsto(P,D=D_G).
 }                                                    \tag{6}
\]

This morphism is polynomial across every non-tree collision divisor.

### Theorem 2.1

On the smooth target locus, (6) is etale and has constant residue Jacobian

\[
 \boxed{
 J_{\Psi_{s,G}}=\pm\frac{d_s}{2}.
 }                                                    \tag{7}
\]

Its generic degree is

\[
 \boxed{
 \deg\Psi_{s,G}=\frac{d_s\,s!}{2}.
 }                                                    \tag{8}
\]

### Proof

The universal ambient determinant gives, after taking the residue along
the tree normalization equations and `m=1`,

\[
 \mu_{s,G}^*(dp_0\wedge\cdots\wedge dp_{s-1})
 =
 \pm d_sD_G\,\omega_{\overline Y}.                   \tag{9}
\]

On `D!=0`, the hypersurface residue form on (3) is

\[
 \Omega_{\widetilde T}
 =
 \frac{dp_0\wedge\cdots\wedge dp_{s-1}}{2D}.          \tag{10}
\]

Pulling back and using `D=D_G` proves (7).  Both sides are regular top forms
on the smooth loci, so the equality extends across the generic divisor
`D=0`.

The unoriented multiplication cover has degree `d_s s!`.  Adjoining `D`
chooses one of the two Vandermonde orientations, giving (8).  QED

For odd `s`, scaling one target form by `2` normalizes the determinant to
one.  For even `s`, `d_s=2` and the displayed determinant is already `+-1`.

## 3. Collision allocation

At a generic point of `D=0`, the target polynomial has one double root and
`s-2` distinct simple roots.  A branch of the ordered cover records the pair
of factor positions occupied by the colliding roots.

There are two cases.

### Non-tree pair

If the colliding positions form a nonedge `f` of `G`, no normalization
equation degenerates.  The branch extends inside the affine source along

\[
 R_f=(r_f=0)\subset\overline Y_{s,G}.                 \tag{11}
\]

There are

\[
 \binom{s}{2}-(s-1)
 =\frac{(s-1)(s-2)}2
\]

such finite collision divisors.

### Tree pair

If the colliding positions form a tree edge `e`, the raw semi-invariant
which is required to equal one tends to zero.  Solving the normalization
weights forces at least one factor scale to have a pole.  The branch
therefore leaves every affine compact set and defines a boundary valuation

\[
 \mathcal D_e\longrightarrow(D=0).                  \tag{12}
\]

The oriented discriminant has order one along this valuation, so the map is
generically unramified there.

More explicitly, if `b_e` is the standard basis vector recording the
vanishing selected resultant, the factor-scaling valuation vector `v_e`
satisfies

\[
 M_G^{\mathsf T}v_e=-b_e.
\]

Since `M_G` is invertible over `Q`, distinct tree edges give distinct
columns of `-(M_G^{\mathsf T})^{-1}` and hence distinct boundary valuations.
Consequently:

### Theorem 3.1

The normalized inverse graph of (6) has at least

\[
 \boxed{s-1}
\]

distinct dicritical boundary primes over the generic oriented discriminant,
one for every edge of `G`.  For the unimodular odd-`s` chart these are
exactly the `s-1` boundary primes over that generic divisor.  The remaining
`(s-1)(s-2)/2` collision types extend as finite affine divisors.

The branch count is consistent with (8).  For every pair of positions, the
other `s-2` roots have `(s-2)!` orderings, with the Kummer index `d_s`
retained.  Thus

\[
 \binom{s}{2}d_s(s-2)!
 =\frac{d_s s!}{2}.                                  \tag{13}
\]

## 4. Three factors and affine descent

For `s=3`, a tree has two edges and one nonedge.  The nonedge collision stays
inside the affine source, while the two tree-edge branches are the two
dicritical divisors of the
[oriented cubic chart](ORIENTED_CUBIC_COX_CHART.md).

For the tangent hyperplane `p_1=1`, quotienting by the involution which
orders the two nonedge factors gives the
[affine descent](ORIENTED_CUBIC_AFFINE_DESCENT.md).  The quotient is the
foundational cubic Keller map on `A^3`, but it exchanges and merges the two
tree-edge dicritical divisors.

Thus affine descent exists in the first arity precisely at the cost of the
independent dicritical marking.

## 5. Finite-field limitation

The oriented target removes the weighted `C=0` excess in every arity.
However, the generic cover still comes from ordered factorizations.
Squarefree split targets have many rational lifts, while irreducible targets
with compatible discriminant orientation can have none.  Hence these maps
are not permutation reductions on fields containing the corresponding
split targets.

The [toric Cox family](TORIC_COX_PERMUTATION_MAPS.md) proves that permutation
reductions are compatible with Cox determinant ledgers in principle.  What
is missing here is exceptional, nonsymmetric monodromy.

## 6. Status of the Cox-ledger conjecture

For ordered linear factors, the Cox-ledger construction is now complete:

1. the full boundary rank and Smith obstruction are known;
2. minimal normalization trees exist in every arity;
3. the universal ambient determinant is exact;
4. separated and compressed reciprocal suspensions are explicit;
5. the oriented target removes all reciprocal coordinates;
6. the oriented maps have constant Jacobian; and
7. arbitrarily many distinct dicritical boundary primes are explicit.

Two stronger affine-space goals remain open:

1. retain at least two of these dicritical primes as distinct invariants
   after passing to an affine-space source and target; and
2. combine that geometry with exceptional monodromy giving permutation
   reductions.

For three factors, the first goal is now achieved on the source side but
not simultaneously on the target side.  The
[affine-source triple-root construction](AFFINE_SOURCE_TRIPLE_ROOT_COX_MAP.md)
has an `A^3` source completion and four distinct dicritical target
components, with constant residue Jacobian on smooth target charts.  Its
target is a five-branch hypersurface rather than affine three-space, and
its exact finite-field profile rules out permutation behavior.

The later
[Danielewski multi-dicritical family](DANIELEWSKI_MULTI_DICRITICAL_FAMILY.md)
achieves the arithmetic objective after changing the boundary model: it
has an `A^3` source, a smooth target, arbitrarily many geometric dicritical
divisors, and bijective reductions whenever the nonzero boundary roots
have no rational residue-field points.  For `P=x(x^2+1)`, this occurs at
every prime `p=3 mod 4`.  The remaining defect is exact: the target has
class `L^3+(deg(P)-1)L^2`, not `L^3`.

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/verify_oriented_linear_factor_cox_maps.py
```

The checker verifies the normalization indices, discriminant and residue
Jacobians in the first arities, the all-arity degree/component count, and
the exact division of collision pairs into finite and dicritical types.
