# Fixed-rank symmetric/Hessian frontier

## Improved repository rank endpoint

The rank-directed circuit search now gives a homogeneous quartic
Hessian-nilpotent Vanishing counterexample in 44 variables with

\[
 \boxed{\operatorname{rank}_{\mathbb C(u,v)}\operatorname{Hess}P=37}.
\]

The separate 42-variable witness remains smaller in ambient dimension and
has exact Hessian rank 38.  Rank and dimension are therefore attained by
different witnesses:

| objective | quartic dimension | exact Hessian rank |
|---|---:|---:|
| smallest repository-certified homogeneous HN-VC witness | 42 | 38 |
| smallest repository-certified Hessian rank | 44 | 37 |

These are incumbents among the exact artifacts currently stored here.  The
table is not a literature-wide best-known claim.  Nor is the rank-37 artifact
a new logical disproof: the unrestricted quartic HN Vanishing Conjecture is
already falsified by any one valid witness, including the 42-variable
construction.  Its contribution is a different quantitative realization
with lower certified Hessian rank.

For the new source, the cubic circuit exposes the `qb` and `x2s` atoms,
performs twenty shared-factor cleanups, rank-compresses a cubic-output space
of dimension eight, and removes constant kernels of dimensions five and one.
The resulting cubic-homogeneous collision has dimension 22, exact generic
Jacobian rank 18, and exact polynomial nilpotency index 18.

Its cotangent Hessian has block ranks

\[
 \operatorname{rank}JH=18,\qquad
 \operatorname{rank}\sum_i y_i\operatorname{Hess}H_i=16,\qquad
 \operatorname{rank}\operatorname{Hess}(y^TH)=37.             \tag{0}
\]

Twelve exact polynomial syzygy generators contain seven generically
independent kernel columns, giving rank at most 37.  A rational
specialization has rank 37, giving the opposite inequality.  The generator
and an independent final-artifact replay are
[`verify_hessian_rank_reduced_bcw_22_route.py`](../scripts/verify_hessian_rank_reduced_bcw_22_route.py)
and
[`audit_hessian_rank_reduced_bcw_22_independent.py`](../scripts/audit_hessian_rank_reduced_bcw_22_independent.py).

## The 42-variable dimension baseline

For comparison, the homogeneous quartic Hessian-nilpotent witness in 42
variables has exact rank 38 and generic corank four.

It is best computed before expanding the orthogonal change.  For the
21-variable cubic-homogeneous collision `V=x+H`, put

\[
 p(x,y)=y^T H(x),\qquad
 \operatorname{Hess}p=
 \begin{pmatrix}
   \sum_i y_i\operatorname{Hess}H_i&(JH)^T\\
   JH&0
 \end{pmatrix}.                                      \tag{1}
\]

The change from `(x,y)` to `(u,v)` is an invertible congruence, so it
preserves the generic rank of (1), although it does not preserve ordinary
matrix powers.  Exact coefficient construction gives 309 nonzero entries.
At three deterministic good-prime specializations,

\[
 \operatorname{rank}JH=18,\qquad
 \operatorname{rank}\sum_i y_i\operatorname{Hess}H_i=15,
 \qquad \operatorname{rank}\operatorname{Hess}p=38. \tag{2}
\]

For the upper bound, a Gröbner syzygy computation over
`Q[x_0,...,x_20,y_0,...,y_20]` produces four polynomially independent kernel
columns of (1).  Their product with (1) is identically zero and their rank is
four at the same rational specialization.  Thus (2) is the exact generic
rank, not only a sampled lower bound.

The rank inflation along the present reduction is:

| stage | map dimension | generic correction-Jacobian rank | cotangent-Hessian rank |
|---|---:|---:|---:|
| foundational normalized map | 3 | 3 | 6 |
| terminal quadratic--cubic trace | 17 | 12 | 24 |
| rank-compressed homogeneous map | 24 | 19 | 38 |
| constant-kernel quotient | 21 | 18 | 38 |

The first row is not an HN quartic: the original correction is
nonhomogeneous and its cotangent Hessian is not nilpotent.  It nevertheless
shows that almost all rank inflation occurs in degree reduction and
homogenization, not in the final symmetric lift.

## Redundant and isotropic directions

The generic Hessian kernel has dimension four, but its constant kernel has
dimension only one.  In Witt `(x,y)` coordinates that constant direction is
`partial/partial y_20`, caused by the identity homogenizing output `H_20=0`.
It is isotropic for the quadratic form `x^T y`.  Deleting it alone makes that
form degenerate, so it does **not** directly produce a 41-variable witness for
the ordinary nondegenerate Laplacian.

The earlier 24-to-22-to-21 constant-Jacobian-kernel quotients do remove three
ambient variables from the cubic map, but all three stages have the same
cotangent-Hessian rank 38.  This illustrates an important distinction:
quotienting redundant variables can reduce essential dimension or corank,
but it cannot by itself improve the rank objective unless it changes the
nondegenerate part of the potential.

At a generic finite-field specialization after the orthogonal change, the
HN Hessian has power ranks

\[
38,36,34,32,30,29,28,\ldots,2,1,0,
\]

and nilpotency index 35.  This is only a specialization diagnostic; the
repository's existing determinant/homogeneity argument is the all-points HN
certificate.  It also shows that the current witness is not close to a
low-nilpotency-index frontier.

## Rank-directed search

[`search_fixed_rank_hessian.py`](../scripts/search_fixed_rank_hessian.py)
uses the same exact shared-factor BCW transition family as the existing
essential-dimension search, but profiles the cotangent Hessian of every
terminal rank-compressed homogenization and its constant-kernel quotient.
The bounded width-64 run examines 115 terminal traces.  Its minimum is 38;
the terminal histogram is 50 traces of rank 38 and 65 of rank 40, and no
constant-kernel quotient lowers either rank.

This is a negative search result only for the enumerated monomial
shared-factor traces.  It is not a lower bound for arbitrary symmetric
reductions.

The circuit-level successor
[`search_restricted_bcw_circuits.py`](../scripts/search_restricted_bcw_circuits.py)
now implements polynomial DAG gates, multi-term target shears, partial
rank and index winners lower the cubic Jacobian rank to 17 and index to 18.
A wider width-64 beam evaluated 396 terminal traces and found the separate
22-variable source certified above,
lowering the quartic Hessian rank from 38 to 37.  Its 44-variable HN lift is
larger than the 42-variable dimension incumbent.  The earlier index winner
has sampled Hessian index 34 rather than 35, but this remains a
specialization diagnostic.
See the [restricted-minima frontier](RESTRICTED_MINIMA_FRONTIER.md).

The next rank-directed search should therefore optimize the operations that
increase rank:

1. track the exact block rank (1) after every candidate stable circuit;
2. allow polynomial DAG gates and multi-term target shears, rather than only
   monomial factor exposures;
3. penalize the induced pairing on `ker JH`, since
   `rank Hess(p)` is controlled by both `rank JH` and that residual pairing;
4. test direct two-parameter quadratic--cubic homogenizations before adding
   the six cubic-output variables;
5. freeze and independently replay only candidates with rank below 37.

## Literature boundary

The classical vanishing formulation and its equivalence with quartic HN
polynomials are due to [Zhao](https://arxiv.org/abs/math/0409534).  The
symmetric reduction is due to
[de Bondt and van den Essen](https://doi.org/10.1090/S0002-9939-05-07570-2).
Known positive islands include homogeneous symmetric nilpotent-Jacobian maps
through dimension five
([de Bondt--van den Essen](https://doi.org/10.1016/j.jpaa.2004.08.030)) and
geometric regular-intersection cases
([van den Essen--Zhao](https://arxiv.org/abs/0704.1690)).  More directly for a
rank parameter, cubic-homogeneous Keller maps with Jacobian rank at most two
are invertible
([de Bondt--Sun](https://arxiv.org/abs/1803.05551)); this includes symmetric
cubic corrections whose quartic potential has Hessian rank at most two.
These results do not put rank 38 near a known sharp boundary.  The computation
here therefore identifies a baseline and a search objective, not yet a
counterexample to a specific published fixed-rank theorem.

## Reproduction

```bash
.venv/bin/python scripts/audit_fixed_rank_hessian_witness.py
.venv/bin/python scripts/verify_hessian_rank_reduced_bcw_22_route.py
python3 scripts/audit_hessian_rank_reduced_bcw_22_independent.py
.venv/bin/python scripts/search_fixed_rank_hessian.py --width 64 --max-steps 24
```

The first three commands give exact certificates for the rank-38 dimension
baseline and the rank-37 circuit witness.  The final command is an
exploratory search and takes several minutes.
