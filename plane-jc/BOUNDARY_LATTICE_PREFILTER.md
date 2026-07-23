# Boundary-lattice prefilter for plane Newton chains

> **Status and scope.**  This is a necessary-condition filter for a proposed
> translation of Newton-chain data into a complete smooth boundary.  It is not
> a new exclusion of an admissible chain.  In particular, the published
> `(8,28) -> (11/4,7)` chain has not yet been translated into the complete
> proximity data required by the filter.

The mechanical half of that translation is now implemented by the
[log-boundary compiler](LOG_BOUNDARY_COMPILER.md).  It accepts certified
positive local branch scales, constructs their regular toroidal blowups, and
feeds the resulting complete matrix to this checker.  It deliberately does
not infer those scales from corners.

Let `X` be a smooth projective rational surface and let

\[
 U=X\setminus(D_1\cup\cdots\cup D_r).
\]

Use a basis of the free lattice `Pic(X)` and put the class `[D_i]` in the
`i`-th column of an integral matrix

\[
 B:\mathbb Z^r\longrightarrow\operatorname{Pic}(X).
\]

Because `O(X)^*=k^*`, divisor localization is

\[
 0\longrightarrow \mathcal O(U)^*/k^*
 \longrightarrow\mathbb Z^r\xrightarrow{B}\operatorname{Pic}(X)
 \longrightarrow\operatorname{Pic}(U)\longrightarrow0.       \tag{1}
\]

Consequently Smith normal form simultaneously reports:

- the free rank of `ker(B)`, hence the rank of boundary-supported units;
- the free rank of `coker(B)`;
- the nontrivial Smith factors, hence the torsion in `Pic(U)`.

For `U=A^2`, both end invariants in (1) vanish, so

\[
 \boxed{r=\rho(X),\qquad \operatorname{SNF}(B)=I_r.}          \tag{2}
\]

The cokernel must be zero, not merely torsion-free.  For the Laurent chart
`U=G_m x A^1`, the correct condition is instead

\[
 \boxed{\operatorname{coker}(B)=0,\qquad\ker(B)=\mathbb Z.}   \tag{3}
\]

Thus a Laurent-stage matrix has one more boundary column than its rank.  It
must not be tested against (2).

## What the input must contain

An admissible-chain row or list of corners is not yet a boundary matrix.  A
candidate certificate must supply:

1. the affine chart being compactified and its expected unit rank;
2. the initial completion;
3. every irreducible boundary prime, not only the dicritical primes;
4. every blowup center and its one or two boundary proximity parents;
5. every final strict-transform class in a fixed Picard basis;
6. boundary valuations of the relevant coordinates and map components;
7. target-boundary pullback multiplicities, recorded separately from the
   primitive source-divisor classes.

Item 7 prevents a common matrix error: a ramification or pullback
multiplicity is not a reason to replace the primitive column `[D_i]` by
`m_i[D_i]`.

The GGHV definition gives corners, directions, homogeneous edge polynomials,
and bracket identities.  It motivates divisorial valuations but does not by
itself specify items 2--5.  Triangular Laurent transformations similarly
change local parameters; they are not, without further work, a global
proximity certificate.

## Why boundary blowups preserve the test

Start with `P^2`, whose affine-plane boundary has matrix `(1)`.  Suppose a
boundary point is blown up and let `epsilon_i` record which old boundary
components contain the center.  In the total-transform basis obtained by
adjoining the new exceptional class, the new boundary matrix is

\[
 B'=
 \begin{pmatrix}
 B&0\\
 -\epsilon&1
 \end{pmatrix}.                                      \tag{4}
\]

Hence `det(B')=det(B)`.  Any complete boundary constructed honestly by
successive boundary blowups of the line at infinity is automatically
unimodular.  The prefilter is therefore primarily an audit of completeness,
primitivity, and realizability of the chain-to-boundary transcription.  It is
not expected to see higher coefficient moduli once that transcription is
correct.

If all boundary classes are present, the full intersection matrix provides a
secondary equivalent determinant check.  In the basis `H,E_1,...,E_s`, let

\[
 J=\operatorname{diag}(1,-1,\ldots,-1),\qquad Q=B^tJB.
\]

For an affine-plane completion, `det(Q)=` $\pm1$.  A selected dicritical
submatrix does not have this property in general and cannot replace the full
class matrix.

## Exact checker

Run the built-in regressions with

```bash
.venv/bin/python plane-jc/cas/boundary_lattice_prefilter.py
```

They cover two boundary blowups of `P^2`, the marked-point matrices in
dimensions two and three, an incorrect doubled divisor class, and the
standard completion of `G_m x A^1`.

The same program accepts a JSON file.  A direct matrix input is

```json
{
  "name": "candidate boundary",
  "expected_unit_rank": 0,
  "boundary_matrix": [[1, 0], [-1, 1]]
}
```

For a boundary-blowup audit beginning with the line `L` in `P^2`, the program
can construct the class matrix from proximity centers:

```json
{
  "name": "two boundary blowups",
  "expected_unit_rank": 0,
  "boundary_blowups": [["L"], ["L", "E1"]]
}
```

Acceptance means only that the localization invariants match the declared
chart.  It does not prove that the indicated centers exist with the required
Newton initial forms, that the map extends with the claimed multiplicities,
or that the Keller coefficient equations are soluble.

For an `A2` chart, the full intersection matrix now has a strictly stronger
second gate: adjunction reconstructs the canonical class and requires
`K_X^2+rho(X)=10`; a target pole vector then reconstructs the ordinary and
logarithmic ramification divisors and selects dicritical components without
an external marking.  See
[`INTRINSIC_A2_BOUNDARY_GATE.md`](INTRINSIC_A2_BOUNDARY_GATE.md) and its
exact checker.  This rejects some unimodular matrices which the localization
test necessarily accepts.

## Application boundary for `(72,108)`

The published reduction begins with

\[
 (8,28)\longrightarrow(11/4,7)
\]

The primitive normal obtained by equating the weights of these two endpoints
is `(4,-1)`.  This identifies one monomial valuation sector, but it does not
identify every boundary component, the shared prefixes of later valuation
resolutions, or a global proximity graph.  The primary proof later moves to
`k[x,x^{-1},y]` and exposes Laurent translations of orders `2,3,4`.  At the
transformed `x=0,1/y=0` crossing these give the positive local rays
`(2,1),(3,1),(4,1)`, as audited in
[FRONTIER_LOG_SCALE_AUDIT.md](FRONTIER_LOG_SCALE_AUDIT.md).  Those resolve
the smooth equality branches, not the rational transformations: their
adapted base ideals are `(t,x^4),(t,x^6),(t,x^8)`, and the compiler resolves
their isolated source graphs in `4,6,8` blowups.  Additive composition with
the common order-four step proves that all cases have the same eight-blowup
translation graph.  The alternative factor residue is separated exactly.
The final monomial involution itself is absorbed by the
Hirzebruch `F_4` transition, which swaps the two base divisors.  Filling
pre-transition `Xinf` (post-transition `X0`) after the common eight-blowup
graph leaves a unimodular `10 x 10` boundary matrix with
unit rank zero, trivial Picard group, and a passing adjunction/Noether audit.
The full target pole vector on that graph has no dicritical component, so the
graph is provably nonexhaustive and at least one additional global cluster is
required.  A smooth point of `E3` is the unique one-blowup zero-pole
extension, but exact two-step witnesses show that positive-pole preparatory
chains are not selected by the lattice data alone.  The coefficient geometry
does select them.  The first weighted Wronskian forces a base-order-ten
blowup at `E3 intersect E4` followed by ten simple children.  At the other
corner the top bracket forces `A=a*r^2,C=c*r^3` for a quartic `r`.  All five
root-partition fans compile and pass the lattice and divisor-class gates; the
primary split-factor formula further forces the quadruple-root member in
Case 1 and supplies its transverse chart.  Thus both source-selected
terminal cases have the same 23-component boundary package, remaining
self-intersection `29`, and one degree-twelve dicritical.  The complete
matrices and local different/conductor data are in
[FRONTIER_LOG_SCALE_AUDIT.md](FRONTIER_LOG_SCALE_AUDIT.md).
That audit also applies the necessary correction for the final transformed
bracket `[P,Q]=X^2`: `K+3H` is only a boundary-supported class
representative there, while the actual ramification is
`K+3H+div(X^2)`.  The corrected dicritical normal indices are three and
five in the two terminal cases.
The two corners
therefore do not determine a full merged cluster of infinitely near points,
and the Laurent stage must satisfy (3), not (2).  A valid application
therefore requires:

1. reconstruct the divisorial valuations before the Laurent localization;
2. resolve them into blowup clusters and merge their common prefixes;
3. restore the original line at infinity and every nondicritical ancestor;
4. compute `B` and its Smith form;
5. only then compare pullback multiplicities and dicritical data.

The repository's exact exclusion shows that the surviving obstruction for
this chain lies in coefficient compatibility over a degree-five field.  No
claim is made here that the lattice filter recovers that contradiction.  Its
immediate value is to reject incomplete future translations before any large
coefficient ideal is built.
