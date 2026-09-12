# Low-complexity profile of the essential cubic Keller collision

## Result

The repository's 21-variable cubic-homogeneous collision is already the
lowest-Jacobian-rank terminal found in the current shared-factor BCW beam, but
it is not close to the low-nilpotency boundary.  Exact computation gives

\[
 F=X+H,\qquad \operatorname{rank}_{\mathbb Q(X)}JH=18,
 \qquad (JH)^{19}=0,\qquad (JH)^{18}\ne0.
\]

Five elementary linear conjugations reduce the expanded correction from 60
component-monomial terms to 53.  They also reduce the number of distinct
scalar cubic monomials from 54 to 45 and give a constructive vector-Waring
decomposition with 164 linear forms, down from the direct 194-form
polarization of the original coordinates.

The generator and exact checker is
[`generate_low_complexity_bcw_21.py`](../scripts/generate_low_complexity_bcw_21.py).
It writes the transported sparse collision to
[`low_complexity_bcw_21_counterexample.json`](../artifacts/generated-results/low_complexity_bcw_21_counterexample.json).

## 1. Exact structural profile

For the sparse conjugate, the verified measures are:

| measure | value |
|---|---:|
| variables | 21 |
| nonzero components of `H` | 20 |
| component-monomial terms | 53 |
| distinct scalar monomials | 45 |
| component coefficient rank | 20 |
| generic rank of `JH` | 18 |
| nilpotency index of `JH` | 19 |
| derivative-flattening rank | 42 |
| vector-Waring rank | between 42 and 164 |

The exact rational power ranks at the deterministic point
`x_i=i^2+3i+5` are

\[
 18,17,16,\ldots,2,1,0.
\]

Thus the specialized Jordan type is `(19,1,1)`.  The rank over
`QQ(x_0,...,x_20)` is computed exactly by Singular, not inferred from the
specialization.  Nilpotence is inherited from the certified
cubic-homogeneous Keller map.  A nilpotent 21-by-21 matrix of rank 18 has at
least three Jordan blocks, so its largest block has size at most 19; the
nonzero eighteenth power at the displayed point makes the generic
nilpotency index exactly 19.

This rules out the hoped-for immediate index-four witness in the existing
coordinates.  Linear conjugacy cannot improve either rank or nilpotency
index, and the previously certified invariant-module audit rules out a
further collision-preserving linear quotient.

## 2. Sparse conjugation

For `S=I+aE_ij`, replace `F` by `S^-1 F S`.  The frozen trace is

| `i` | `j` | `a` | terms after the move |
|---:|---:|---:|---:|
| 4 | 6 | `-3/2` | 58 |
| 5 | 8 | `-3` | 56 |
| 10 | 8 | `1` | 55 |
| 4 | 3 | `-3` | 55 |
| 13 | 7 | `6/7` | 53 |

The fourth move is deliberately term-neutral.  It exposes the two
cancellations in the fifth move; a strictly greedy descent stops at 55.
Every collision point and the common image are transported by the same
inverse linear maps and are checked over `QQ` in the generated artifact.

There is a tradeoff: the number of nonzero coordinate entries of `JH`
increases from 92 to 96.  The optimization here is for expanded cubic terms
and scalar-monomial support, not entrywise Jacobian sparsity.

## 3. Tensor and linear-form bounds

Use "vector-Waring rank" for the least `r` in a decomposition

\[
 H(x)=\sum_{\nu=1}^r v_\nu\ell_\nu(x)^3,
 \qquad v_\nu\in\mathbb Q^{21}.
\]

If such a decomposition has `r` terms, then

\[
 JH=3\sum_\nu (v_\nu\otimes\ell_\nu)\ell_\nu^2.
\]

Flattening with rows indexed by `(output,input)` and columns by quadratic
monomials therefore has rank at most `r`.  Its exact rank is 42, proving
`r>=42`.

For the upper bound, the checker simultaneously polarizes the coefficient
vector of every supported monomial using

\[
 a^2b=\frac{(a+b)^3-(a-b)^3-2b^3}{6}
\]

and

\[
 abc=\frac{(a+b+c)^3-(a+b-c)^3-(a-b+c)^3-(-a+b+c)^3}{24}.
\]

After identical forms are merged and zero coefficient vectors removed, 164
distinct forms remain, and the script re-expands them to the exact 53-term
map.  Hence

\[
 \boxed{42\le \operatorname{rank}_{\mathrm{vW}}(H)\le164}.
\]

This is not yet a square Druzkowski representation.  It is the relevant
linear-form count for a vector-valued power-linear expansion and can feed a
power-linear stabilization, but the extra pairing and square-map identities
must be constructed before calling 164 a Druzkowski dimension.

## 4. Bounded low-rank search

[`search_low_rank_bcw.py`](../scripts/search_low_rank_bcw.py) changes the
terminal diagnostic from essential dimension alone to the post-quotient
power-rank profile.  With width 24 and 17 steps it profiles 39 terminal
traces at one deterministic point modulo 1,000,003:

| terminals | quotient dimension | sampled rank | sampled index | terms |
|---:|---:|---:|---:|---:|
| 24 | 21 | 18 | 19 | 60 |
| 2 | 22 | 18 | 19 | 62 |
| 13 | 23 | 19 | 19 | 64 |

The modular values are search diagnostics, not characteristic-zero proofs.
The winning 21-variable trace is separately frozen and audited exactly as in
Section 1.  No terminal in this beam lowers the rank below 18 or the index
below 19.

Run the bounded search with

```bash
.venv/bin/python scripts/search_low_rank_bcw.py --width 24 --max-steps 17
```

## 5. Circuit-level successor

The current architecture carries a length-19 generic Jordan chain through
all stored reductions.  The dimensions 21, 22, and 24 differ mainly by zero
directions and triangular quotient data; optimizing the final basis cannot
shorten that chain.

That search has now been implemented.  Exposing the two polynomial gates
`A=w+xy`, `Q=s+xy^2` and applying the multi-term target shear
`F_1 <- F_1-12 Q(2+3A)` before monomial cleanup produces a different
22-variable cubic-homogeneous collision.  Exact computation gives

\[
 \operatorname{rank}JH=18,\qquad
 (JH)^{17}\ne0,\qquad (JH)^{18}=0.
\]

Thus the circuit change lowers the certified upper bound for the minimum
nilpotency index from 19 to 18.  It does not lower the rank or the best
ambient dimension.  A second two-atom circuit gives a 24-variable collision
with exact generic rank 17 and the same exact index 18, lowering the rank
upper bound as well.  The constructions, exact polynomial-power
certificates, and search architecture are in the
[restricted-minima frontier](RESTRICTED_MINIMA_FRONTIER.md).

For context, the rank-two cubic-homogeneous Keller case is known invertible
([de Bondt--Sun](https://arxiv.org/abs/1803.05551)), while an inverse formula
is known for homogeneous power-linear maps with `(JH)^3=0`
([de Bondt--Yan](https://arxiv.org/abs/1302.5864)).  These results concern
different restricted classes.  The exact frontiers are now
`3<=r_cub<=17` and `3<=nu_cub<=18`; the full cubic-homogeneous index-three
invertibility case is not settled by the power-linear or symmetric-Jacobian
theorems.  The stronger
uniform degree-nine bound is false already in dimension five: van den Essen's
generic-rank-three weak-index-three automorphism has inverse degree thirteen.
