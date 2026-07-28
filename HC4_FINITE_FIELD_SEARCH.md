# Collision-first finite-field search for `HC_4`

## Status

This is a bounded computational experiment.  It is not a proof of `HC_4`,
not a counterexample, and not a characteristic-zero lifting theorem.

An exact search over \(\mathbf F_{11}\) and \(\mathbf F_{13}\) found no
four-variable potential in the collision-normalized one/two-direction family
defined below.  The search covered degree bounds \(5,6,7,8\) and tested
45,181,194 coefficient choices.  Every point-filter survivor failed the
full coefficient expansion of the Hessian determinant.

The main output is therefore the reusable collision-first search
infrastructure, not a mathematical exclusion of general degree-\(5\) through
degree-\(8\) potentials.

A second exact experiment sampled 96 denser supports with 6, 8, 10, or 12
collision-kernel directions.  For each support it formed the full coefficient
ideal of \(\det\operatorname{Hess}\Psi-1\), adjoined the finite-field
equations, and reduced over both \(\mathbf F_{11}\) and \(\mathbf F_{13}\).
All 192 ideals were unit ideals, with no timeout or solver error.

A structural follow-up first identified the 3, 5, 5, and 7 respective kernel
directions in degree bounds 5, 6, 7, and 8 that can alter both forced
determinant defects on the \(x_0\)-axis.  It sampled 32 supports containing as
many of these directions as their support size allowed.  All 64
support-prime ideals were again unit ideals.

A principal-Hessian follow-up sampled another 128 supports.  Their
degree-\(d\) terms omit either \(x_2\) or \(x_3\), so the top homogeneous
Hessian determinant vanishes identically.  Lower-degree monomials involving
the omitted variable provide the fourth-variable bridges.  All 256
support-prime full coefficient ideals were unit ideals.

A non-coordinate follow-up used
\(u=x_2+\lambda x_3\) for \(\lambda=-1,1,2\).  Its 144 oblique
cone/complementary-bridge families again have identically zero principal
Hessian determinant.  All 288 full coefficient ideals over the two primes
were unit ideals, with no timeout or solver error.

## 1. Collision normalization

Write \(x=(x_0,x_1,x_2,x_3)\), normalize the two proposed collision points
to \(0,e_0\), and begin with

\[
q=x_0x_1+x_2x_3,\qquad \det\operatorname{Hess}q=1.
\]

For a degree bound \(d\), set

\[
\Psi_{\mathrm{base}}
=q-x_1x_0^{d-1}.
\]

Then

\[
\nabla\Psi_{\mathrm{base}}(0)
=\nabla\Psi_{\mathrm{base}}(e_0).
\]

Among monomials of any fixed degree \(k\), only

\[
x_0^k,\quad x_1x_0^{k-1},\quad
x_2x_0^{k-1},\quad x_3x_0^{k-1}
\]

contribute to the gradient difference at these two points.  Consequently,
the script constructs an explicit basis of the full kernel of

\[
\bigoplus_{k=3}^d \operatorname{Sym}^k(k^4)
\longrightarrow k^4,\qquad
h\longmapsto\nabla h(e_0)-\nabla h(0).
\]

For a lower-degree visible monomial, its basis vector subtracts the matching
degree-\(d\) carrier.  In the pure \(x_0\) channel the correcting coefficient
is \(-k/d\).  Taking primes \(p>d\) makes this basis valid without derivative
or denominator degeneracy.

The kernel dimensions for \(d=5,6,7,8\) are respectively

\[
107,\quad191,\quad311,\quad476.
\]

## 2. Exhausted family

For every kernel basis vector \(v_i\), and every pair \(v_i,v_j\), the
search exhausts

\[
\Psi_{\mathrm{base}}+a v_i,\qquad
\Psi_{\mathrm{base}}+a v_i+b v_j,
\]

with \(a,b\ne0\) in the selected finite field.  Thus “support two” refers to
support in this collision-kernel basis.  A carrier direction can itself
contain two ordinary monomials.

The exact counts are:

| prime | degree bound | directions | potentials |
|---:|---:|---:|---:|
| 11 | 5 | 107 | 568,170 |
| 11 | 6 | 191 | 1,816,410 |
| 11 | 7 | 311 | 4,823,610 |
| 11 | 8 | 476 | 11,309,760 |
| 13 | 5 | 107 | 817,908 |
| 13 | 6 | 191 | 2,615,172 |
| 13 | 7 | 311 | 6,945,252 |
| 13 | 8 | 476 | 16,284,912 |

At six deterministic points the script evaluates

\[
\det\operatorname{Hess}\Psi-1.
\]

A nonzero value is an exact rejection, not a probabilistic inference.  Every
choice surviving all six evaluations is then expanded as a sparse polynomial
over \(\mathbf F_p[x_0,x_1,x_2,x_3]\).  Acceptance requires literal equality

\[
\det\operatorname{Hess}\Psi=1.
\]

No exact candidate survived for any row of the table.

## 3. Sampled dense-support coefficient solving

The follow-up search uses the same collision normalization and kernel basis,
but replaces coefficient enumeration by exact Gröbner solving.  It samples
supports of sizes

\[
6,\quad8,\quad10,\quad12
\]

using three deterministic strategies:

- `uniform` samples from the full collision-kernel basis;
- `homogeneous` samples only single-monomial degree-\(d\) directions;
- `mixed` forces approximately one third of the support to be two-monomial
  carrier directions and fills the rest with single-monomial directions.

For every degree bound, support size, strategy, and two deterministic trials,
the script expands

\[
\det\operatorname{Hess}
\left(\Psi_{\mathrm{base}}+\sum_{i=1}^{k}a_iv_i\right)-1
\]

over \(\mathbb Q[a_1,\ldots,a_k][x_0,x_1,x_2,x_3]\).  After reduction modulo
\(p\), every coefficient in the \(x\)-variables is put into an ideal together
with

\[
a_i^p-a_i,\qquad 1\leq i\leq k.
\]

Thus the resulting ideal is nonunit exactly when that selected support has a
coefficient solution over \(\mathbf F_p\).  This is a full polynomial
identity test; it does not use point sampling.

There are 96 distinct supports and two primes, hence 192 ideals.  Singular
reduced every one to the unit ideal.  The largest determinant expansion had
2,400 terms, and the largest characteristic-zero coefficient system had
2,135 equations.

### Axis-obstruction support selection

For the base potential, restriction to the \(x_0\)-axis has forced
determinant-defect terms in degrees \(d-2\) and \(2d-4\).  Most random kernel
directions cannot change either coefficient.  The `axis` strategy computes
the exceptional directions that can change both:

| degree bound | eligible directions |
|---:|---:|
| 5 | 3 |
| 6 | 5 |
| 7 | 5 |
| 8 | 7 |

It then includes as many eligible directions as possible and fills the
remaining positions deterministically.  Two trials for each degree and each
support size \(6,8,10,12\) give 32 further supports.  Their 64 full
coefficient ideals over the two primes are all unit ideals, again with no
timeout or solver error.

### Principal-Hessian cone supports

The `cone2` and `cone3` strategies pass the highest-degree determinant gate
by construction.  For `cone2`, every selected degree-\(d\) monomial omits
\(x_2\); for `cone3`, every such monomial omits \(x_3\).  The base
degree-\(d\) collision carrier depends only on \(x_0,x_1\).  Hence the full
top homogeneous potential is a polynomial in at most three variables and

\[
\det\operatorname{Hess}\Psi_d=0
\]

identically in all selected top coefficients.

One half of the selected directions lie in this top cone.  The remaining
directions are lower-degree, collision-invisible monomials involving the
omitted fourth variable.  They are the bridges that can contribute to the
lower determinant layers without disturbing the principal identity.

Four deterministic trials for each degree bound, support size
\(6,8,10,12\), and omitted coordinate give 128 supports.  These are disjoint
from the preceding 128 supports.  Their 256 exact ideals over
\(\mathbf F_{11}\) and \(\mathbf F_{13}\) all reduce to the unit ideal, with
no timeout or solver error.  Thus these families fail below the principal
homogeneous Hessian layer.

### Non-coordinate principal cones

The oblique search uses the determinant-one linear coordinates

\[
u=x_2+\lambda x_3,\qquad v=x_3,
\qquad \lambda\in\{-1,1,2\}.
\]

In the adapted coordinates \((x_0,x_1,u,v)\), the quadratic base is

\[
x_0x_1+uv-\lambda v^2,
\]

which still has Hessian determinant one.  The degree-\(d\) correction lies
in \(k[x_0,x_1,u]\), so its four-variable Hessian determinant vanishes
identically.  Half of each selected family consists of these top cone
directions.  The other half consists of lower-degree, collision-invisible
monomials involving \(v\).

Returning to the original coordinates expands the top and bridge parameters
into tied collections of \(x_2,x_3\) monomials.  Thus this is not one of the
coordinate-cone supports above, although working in the adapted chart keeps
the exact determinant expansion sparse.

Three trials for each slope, degree bound, and support size \(6,8,10,12\)
give 144 families.  All 288 coefficient ideals over
\(\mathbf F_{11},\mathbf F_{13}\) are unit ideals, with no timeout or solver
error.  The largest determinant expansion has 2,268 terms and the largest
characteristic-zero coefficient system has 1,892 equations.

## 4. Scope and next search layer

The exhaustive result excludes every affine perturbation supported on at
most two vectors of the chosen collision-kernel basis.  The dense runs add
400 individual support-family exclusions, not complete support-size classes.
The combined experiment does not exclude:

- any unrecorded support of three or more kernel directions;
- a different quadratic Hessian or collision normalization over the finite
  field;
- potentials whose collision points are not rational over the searched
  prime field;
- solutions present only in other characteristics;
- rational or number-field lifts.

The sampled dense-support results add exact exclusions for 400 recorded
families, but they do not exhaust all supports of sizes 6 through 12.  The
useful next layer is a moving cone direction: make the omitted linear form
itself a coefficient parameter in the full projective three-space, or allow
successive homogeneous layers to use different cone directions.  Any family
recurring as nonunit over several primes can then be passed to rational
reconstruction and an exact characteristic-zero check.

## Reproduction

Run:

```bash
.venv/bin/python scripts/search_hc4_finite_field_potentials.py \
  --degrees 5 6 7 8 \
  --primes 11 13 \
  --support-bound 2 \
  --points 6 \
  --output artifacts/generated-results/hc4_finite_field_sparse_search.json
```

The generated JSON records every search count, the rejection points, all
exact modular candidates (none in this run), and supports occurring at more
than one prime.

Replay the sampled dense-support coefficient ideals with:

```bash
.venv/bin/python scripts/search_hc4_finite_field_dense_supports.py \
  --degrees 5 6 7 8 \
  --primes 11 13 \
  --support-sizes 6 8 10 12 \
  --strategies uniform homogeneous mixed \
  --trials 2 \
  --timeout 30 \
  --output artifacts/generated-results/hc4_finite_field_dense_support_search.json
```

This command requires Singular.  Its JSON ledger records every support,
coefficient-system size, modular Gröbner outcome, timeout status, and a
deterministic content hash that excludes runtime measurements.

Replay the axis-obstruction-guided supports with:

```bash
.venv/bin/python scripts/search_hc4_finite_field_dense_supports.py \
  --degrees 5 6 7 8 \
  --primes 11 13 \
  --support-sizes 6 8 10 12 \
  --strategies axis \
  --trials 2 \
  --timeout 30 \
  --output artifacts/generated-results/hc4_finite_field_axis_support_search.json
```

Replay the principal-cone/fourth-variable-bridge supports with:

```bash
.venv/bin/python scripts/search_hc4_finite_field_dense_supports.py \
  --degrees 5 6 7 8 \
  --primes 11 13 \
  --support-sizes 6 8 10 12 \
  --strategies cone2 cone3 \
  --trials 4 \
  --timeout 30 \
  --output artifacts/generated-results/hc4_finite_field_cone_bridge_search.json
```

Replay the non-coordinate cone/bridge families with:

```bash
.venv/bin/python scripts/search_hc4_oblique_cone_bridges.py \
  --degrees 5 6 7 8 \
  --primes 11 13 \
  --support-sizes 6 8 10 12 \
  --slopes -1 1 2 \
  --trials 3 \
  --timeout 30 \
  --output artifacts/generated-results/hc4_finite_field_oblique_cone_bridge_search.json
```
