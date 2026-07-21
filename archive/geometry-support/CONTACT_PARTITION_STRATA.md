# Contact-partition strata

This note generalizes the named dual-curve incidences to an arbitrary contact
partition

\[
\lambda=(\lambda_1,\ldots,\lambda_k),\qquad \lambda_i\ge2.
\]

The implementation is `contact_partition_incidence` in
`jcsearch.discriminant_geometry`.

## Uniform parameterization

For full contact, put

\[
M_\lambda(W)=\prod_{i=1}^k(W-r_i)^{\lambda_i}
\]

and define

\[
\Phi_\lambda=M_\lambda(1)-M_\lambda(0)-M_\lambda'(0),
\qquad
D_\lambda=M_\lambda'(1)-M_\lambda'(0).
\]

On `Phi_lambda=0` and `D_lambda!=0`, the normalization `c=1` gives

\[
a=-{1\over D_\lambda},\qquad
s={M_\lambda'(0)\over D_\lambda},\qquad
t=-{M_\lambda(0)\over D_\lambda},
\]

and

\[
H(W)={-M_\lambda(W)+M_\lambda'(0)W+M_\lambda(0)
       \over D_\lambda}.
\]

Thus `H(0)=H'(0)=H(1)=0` and `H'(1)=-1`. Scaling all four quantities by a
nonzero `c` restores arbitrary normalization.

For partial contact, the API multiplies the marked contact product by a monic
residual polynomial `Q`. Resultants `Res(Q,Q_mu)` ensure that `Q` does not
meet a marked root, so every requested multiplicity is exact.

## Equal-part quotient and saturation

Roots carrying the same multiplicity are quotiented before elimination. If a
multiplicity `mu` occurs `j` times, their contribution is

\[
Q_\mu(W)^\mu,
\]

where the coefficients of the monic degree-`j` polynomial `Q_mu` are the
elementary symmetric root coordinates. The quotient has generic degree
`product_mu (j!)`. Distinctness is imposed by:

- `Disc(Q_mu)` within each equal-part block;
- `Res(Q_mu,Q_nu)` between different blocks;
- `Res(Q_mu,Q)` against a residual factor.

The weighted admissibility factor also contains

\[
D_\lambda\bigl(M_\lambda''(1)-2D_\lambda\bigr).
\]

For normalized seed coordinates `h_3,...,h_(n-1)`, put

\[
h_n={-1-\sum_{k=3}^{n-1}(k-2)h_k\over n-2}.
\]

The coefficient-space elimination ideal has the exact Rabinowitsch
presentation

\[
\left\langle
\Phi_\lambda,
D_\lambda h_k+[W^k]M_\lambda\ (3\le k<n),
1-zA_\lambda
\right\rangle
\cap \mathbb Q[h_3,\ldots,h_{n-1}],
\]

where `A_lambda` is the full distinctness, residual, scale, and weighted-open
factor. The API returns this presentation without forcing an impractical
expanded Groebner basis; `incidence_elimination_generators` computes the
elimination when desired.

For full contact the root-quotient hypersurface has `k` coordinates and one
equation.  The weighted-Vandermonde theorem proves—not merely predicts—that
its seed image has dimension `k-1`; see
`UNIFORM_EXCEPTIONAL_SEEDS.md`.  Uniform nonemptiness of the exact admissible
open, finiteness of the top-coefficient map, and the exact closure converse
are proved together in
[COINCIDENT_ROOT_REBUILD.md](../../experimental/geometry/COINCIDENT_ROOT_REBUILD.md).

## Simultaneous omissions

`multiple_omission_incidence(degree, partitions)` identifies the normalized
seed coefficients for two or more full-contact factorizations.  For every
pair it returns the cleared numerators of `s_i-s_j` and `t_i-t_j`:

- their joint vanishing is the common-value collision locus;
- a two-generator Rabinowitsch gate forces the omitted values to be genuinely
  distinct.

This avoids conflating a coarser root partition of one inverse polynomial
with a seed having two different omitted targets.

The unique omitted-value theorem proves that the distinct-target
Rabinowitsch locus is empty in every degree.  The simultaneous API remains a
useful exact presentation of that emptiness and of common-value collision
boundaries.  Consequently the exact `E_lambda` form a disjoint stratification
of the nonsurjective seed locus.

For two partitions, `two_omission_incidence` additionally returns the uniform
affine-difference presentation

\[
aM_\lambda-bM_\mu=\alpha W+\beta.
\]

Polynomial Mason--Stothers proves that its off-collision saturation is empty
whenever `lambda!=mu` are maximal partitions using only twos and threes.  The
proof and collision partial order are stated in
`UNIFORM_EXCEPTIONAL_SEEDS.md`.

For a maximal partition, `maximal_two_three_phi(a,b)` returns the compact
quotient model `M=Q^2R^3`.  Its `Phi` hypersurface is irreducible for every
`a,b`: the stable proof is linear in `R'(0)` for `b>=3`, quadratic with an odd
valuation for `a>=3`, and has seven explicit endpoint-rank cases.  Hence the
maximal partition types index the irreducible components of the closure of
the exceptional seed locus.

The admissible quotient hypersurface remains smooth when collision diagonals
are retained.  Its finite map to the corresponding seed component is the
normalization, and `component_decomposition_count` counts its set-theoretic
fibers over collision partitions.

## Degree-five regression

For `lambda=(3,2)`, write

\[
M=(W-a)^3(W-b)^2.
\]

The API produces the root curve `Phi_(3,2)=0`. Its map to the existing
two-extra-root coordinates is

\[
R=3a+2b-1,
\qquad
P=[W^3]M-R.
\]

Eliminating `a,b` gives

\[
{5\over64}F(R,P),
\]

where `F` is exactly the polynomial in `OMITTED_VALUE_CLASSIFICATION.md`.
This identifies the old exceptional equation with the uniform full-contact
incidence, including its normalization and saturation data.

## Degree-six geometry

The normalized seed space has coordinates `(h_3,h_4,h_5)` and dimension
three. The exact stratum dimensions are:

| Partition | Dimension |
|---|---:|
| `(2,2,2)` | 2 |
| `(3,3)` | 1 |
| `(4,2)` | 1 |
| `(6)` | 0 |

For the main stratum, write

\[
Q(W)=W^3-e_1W^2+e_2W-e_3,qquad M=Q^2.
\]

Then

\[
\Phi_{2,2,2}
=(1-e_1+e_2-e_3)^2-e_3^2+2e_2e_3.
\]

It is irreducible and linear in `e_3`, with coefficient `2(e_1-1)`. The
coefficient map is birational: `e_1,e_2,e_3` are recovered successively from
the top three monic coefficients of `M`. Consequently:

- the `(2,2,2)` closure is irreducible;
- it is a rational surface;
- its normalized coefficient equation is irreducible of total degree four.

The full quartic is generated and checked by
`scripts/verify_contact_partition_strata.py`; its compact definition is the
numerator obtained after substituting the recovered `e_i` into
`Phi_(2,2,2)`.

The boundary relations are:

- `(4,2)` lies in the `(2,2,2)` closure through
  `Q=(W-u)^2(W-v)`;
- it is not contained in the singular locus of the quartic surface—the audit
  supplies an exact smooth point;
- generic `(3,3)` is not in the main closure—the admissible seed
  `(h_3,h_4,h_5)=(-35/27,-10/9,4/9)` gives a nonzero quartic value;
- `Phi_6` is a squarefree quartic, so the `(6)` stratum consists of four
  admissible points over the algebraic closure;
- every `(6)` point is a common collision boundary of `(4,2)`, `(3,3)`, and
  `(2,2,2)`;
- all four `(6)` points are smooth on the main quartic.
- the complete `(3,3)` pullback of the main quartic is the square of its
  two-root discriminant, up to invertible factors.  Hence the exact
  `(3,3)` and `(2,2,2)` strata are disjoint and their closures meet only on
  `(6)`; there is no seed with two distinct omitted values of these types.

Thus the lower curves describe different phenomena: `(4,2)` is a boundary of
the main square locus, whereas `(3,3)` is a separate exceptional curve sharing
the maximally ramified boundary points.

## Degree-seven geometry

The normalized seed space has dimension four. The exact dimensions are:

| Partition | Dimension | Codimension |
|---|---:|---:|
| `(3,2,2)` | 2 | 2 |
| `(4,3)` | 1 | 3 |
| `(5,2)` | 1 | 3 |
| `(7)` | 0 | 4 |

Every displayed root-space `Phi_lambda` is irreducible and coprime to its
weighted admissibility factor. The uniform weighted-Vandermonde certificate
proves that the leading `(3,2,2)` coefficient image has dimension two; the
older exact rank point at quotient coordinates `(-6,-3,-4)` remains as a
low-degree regression check.

Therefore the degree-seven nonsurjective locus has codimension two. In
particular, no single nonzero coefficient equation detects nonsurjectivity:
the first possible exceptional set is genuinely higher-codimensional.

Run:

```bash
.venv/bin/python scripts/verify_contact_partition_strata.py
.venv/bin/python scripts/verify_uniform_exceptional_seed_theorem.py
```
