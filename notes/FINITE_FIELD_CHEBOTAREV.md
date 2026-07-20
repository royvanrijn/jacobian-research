# Finite-field Chebotarev law for weighted seeds

This note derives the finite-field statistic common to every weighted seed
from the universal `S_n` monodromy theorem. It is an asymptotic good-reduction
result. Exact formulas in exceptional or small characteristics remain separate
computations, as in `FINITE_FIELD_VALUE_DISTRIBUTION.md`.

## Pencil Chebotarev theorem

Let `H` be a characteristic-zero primitive of degree `n>=3`, defined over a
number field, and consider

\[
E_{s,t}(W)=H(W)-sW+t.
\]

After choosing an integral model, exclude the finite set of primes where a
coefficient denominator, the degree, separability, or the `S_n` cover has bad
reduction. For every remaining residue field `F_q`, let `U` be the complement
of the discriminant curve in the `(s,t)`-plane.

The geometric monodromy of the degree-`n` root cover over `U` is `S_n`.
Effective Chebotarev for finite morphisms of varieties over finite fields,
using Lang--Weil on the corresponding twists, gives for every conjugacy class
`K` in `S_n`

\[
\#\{(s,t)\in U(\mathbb F_q):\operatorname{Frob}_{s,t}\in K\}
={|K|\over n!}q^2+O_H(q^{3/2}).
\]

The cycle lengths of Frobenius are the degrees of the irreducible factors of
`E_{s,t}`. In particular, its fixed points are exactly the simple rational
roots.

This use of Chebotarev for finite quasi-projective morphisms follows the
finite-field formulation proved via twisted varieties and Lang--Weil in
[Steve Meagher, *A simple proof of Chebotarev's density theorem over finite
fields*](https://doi.org/10.1017/S0004972718000448). The specialization from
factorization types to symmetric-group cycle types is also developed for
families `f(W)+sW+t` by [Pär Kurlberg and Lior Rosenzweig, *The Chebotarev
density theorem for function fields -- incomplete
intervals*](https://arxiv.org/abs/1901.06751).

## Fixed-point distribution

Let `D_m` denote the number of derangements in `S_m`. The number of
permutations in `S_n` with exactly `j` fixed points is

\[
{n\choose j}D_{n-j}.
\]

Consequently, if `M_j(q)` counts all pencil parameters `(s,t)` with exactly
`j` simple rational roots, including the lower-order discriminant locus, then

\[
M_j(q)=p_{n,j}q^2+O_H(q^{3/2}),
\qquad
p_{n,j}={{n\choose j}D_{n-j}\over n!}.
\]

The discriminant curve has only `O_H(q)` rational points, so adding its
ramified specializations does not change the error term.

The first laws are:

| `n` | Nonzero probabilities `p_{n,j}` |
|---:|---|
| 3 | `p_0=1/3`, `p_1=1/2`, `p_3=1/6` |
| 4 | `p_0=3/8`, `p_1=1/3`, `p_2=1/4`, `p_4=1/24` |
| 5 | `p_0=11/30`, `p_1=3/8`, `p_2=1/6`, `p_3=1/12`, `p_5=1/120` |
| 6 | `p_0=53/144`, `p_1=11/30`, `p_2=3/16`, `p_3=1/18`, `p_4=1/48`, `p_6=1/720` |

There is never a permutation with exactly `n-1` fixed points. For every
`1<=k<=n`, the falling-factorial moment satisfies

\[
\mathbb E[(\#\operatorname{Fix})_k]=1.
\]

These identities are the group-theoretic limits of the ordered distinct
fiber-product counts.

## Transfer to the three-dimensional weighted map

For fixed `C!=0`, the change of target coordinates

\[
(A,B)\longmapsto(s,t)=(BC,cAC^2)
\]

is a bijection over `F_q`. The affine fiber of the weighted map has one point
for every simple rational root of `E_{s,t}`. Hence the count on `C!=0` is
exactly `(q-1)M_j(q)`.

The plane `C=0` contains only `q^2` targets. Its direct boundary fibers, all
ramified specializations, and every omitted-value lift therefore contribute
below the Chebotarev error after the extra `C` parameter is included. If
`N_j(q)` is the full number of targets in `F_q^3` with `j` rational source
points, then

\[
N_j(q)=p_{n,j}q^3+O_H(q^{5/2}).
\]

Thus the limiting rational image density is

\[
\lim_{q\to\infty}{|G_H(\mathbb F_q^3)|\over q^3}
=1-p_{n,0}
=1-{D_n\over n!}.
\]

As `n` grows, this tends to `1-e^{-1}`. This does not conflict with geometric
surjectivity over the algebraic closure: a target can have preimages without
having a rational preimage.

The limiting mean fiber size is one, and in fact the mean over all targets is
exactly one for every finite field because the source and target both have
`q^3` points.

## Exact incidence check and diagnostics

There is one useful identity before taking a limit. For each candidate root
`r` and each slope `s!=H'(r)`, exactly one intercept `t` makes `r` a simple
root. Therefore

\[
\sum_{s,t}\#\{\text{simple roots of }E_{s,t}\}=q(q-1).
\]

The module `jcsearch.chebotarev` implements derangements, exact `S_n`
fixed-point probabilities, and prime-field pencil histograms. Run:

```bash
.venv/bin/python scripts/verify_weighted_chebotarev.py
```

The verifier proves the finite group identities through `S_8`, checks the
exact simple-root incidence, and reports finite-sample discrepancies for
inverse degrees three, four, and five. Those samples are diagnostics; the
asymptotic theorem comes from monodromy plus finite-field Chebotarev.

## Scope

The theorem applies after excluding finitely many bad primes and then letting
`q` grow through good residue fields. It does not replace:

- exact formulas for a specified small characteristic;
- analysis where the inverse degree drops modulo the characteristic;
- inseparable reductions or reductions where the geometric group changes;
- the seed-specific `C=0` fiber tables.

These distinctions explain why the original cubic has a universal `S_3`
limit while characteristic two has a different exact distribution.
