# Fixed-row Newton ramification

This note extracts the strongest Galois-theoretic consequence that currently
follows uniformly from the fixed-derivative-order argument of
[Borisov--Filaseta--Lam--Trifonov](https://doi.org/10.4064/aa-90-2-121-153).
It is stronger than density-one irreducibility, but it is not yet a
density-one alternating-group theorem.

Fix `r>=1`, put

\[
 h=(m+1)r,\qquad d=h-r=mr,
\]

and write

\[
 B_{d,r}(x)=\sum_{j=0}^{d}\binom{r+j}{r}x^j
 =\frac1{r!}G_h^{(r)}(x),
 \qquad G_h(x)=1+x+\cdots+x^h.                         \tag{1}
\]

The splitting fields of `B_(d,r)` and the cancellation polynomial
`M_(m,r)` are identical after the affine-reciprocal change recorded in
[ARITHMETIC.md](ARITHMETIC.md).

## 1. The extracted density-one theorem

### Theorem

For every fixed `r>=1`, all but

\[
 O_r\!\left(\frac{X\log\log X}{\log X}\right)          \tag{2}
\]

integers `m<=X` have all of the following properties:

1. `B_(mr,r)` is irreducible over `Q`;
2. its discriminant is not a rational square; and
3. there is a prime

   \[
   (\log X)^{10}<p\le mr-3                            \tag{3}
   \]

   that divides the order of its Galois group.

Consequently the transitive Galois group contains an element of order `p`.
In its degree-`d=mr` root action this element has cycle type

\[
 p^t1^f,\qquad t\ge1,\qquad f\ge p-1.                 \tag{4}
\]

The exponent `10` is only the convenient cutoff used in the cited proof; no
optimization is intended.

### Proof

The first assertion is Theorem 2 of Borisov--Filaseta--Lam--Trifonov,
restricted from all geometric degrees `h<=r(X+1)` to the progression
`h=r(m+1)`.  Their proof gives the quantitative exceptional count (2).
The nonsquare assertion follows after removing the
`O_r(sqrt(X))` square-discriminant parameters proved in
[the discriminant note](../archive/cancellation-components/PARAMETER_DISCRIMINANT.md).

It remains to extract the ramification statement from the same proof.  Form
the monic reciprocal numerator

\[
 W_{h,r}(x)=(x-1)^{r+1}x^dB_{d,r}(x^{-1}).             \tag{5}
\]

In the notation of the cited paper, take its integer near the derivative
degree to be

\[
 a=h-r+1=d+1.                                         \tag{6}
\]

This choice occurs in every one of their cases, including the separately
treated derivative orders `2` and `4`.  If `p\parallel a`, `p>r`, and
`a=pa_0`, then

\[
 W_{h,r}(x)\equiv (x^a-1)x^r\pmod p.                  \tag{7}
\]

Outside their bad-pair set, for every nontrivial `a_0`-th root of unity
`zeta`, Proposition 1 gives

\[
 v_p(W_{h,r}(\zeta))=1
\]

and the Newton polygon of `W_(h,r)(x+zeta)` begins with

\[
 (0,1)\longrightarrow(p,0).                           \tag{8}
\]

Thus a root has valuation `1/p` over the unramified field
`Q_p(zeta)`.  The ramification index of the local splitting field is
therefore divisible by `p`.  Unramified base change does not alter inertia,
so `p` divides the order of the global Galois group.  The repeated rational
root `x=1` in (5) contributes no ramification; the order-`p` action is
nontrivial on the roots of `B_(d,r)`.

For fixed `r>=2`, Lemmas 16 and 18, the bad-pair set `T`, and the
large-squarefree-kernel set `T_(00)` in the proof of their Theorem 2 show
that, apart from (2) exceptions, the
integer `a=d+1` has a prime `p\parallel a` above
`z=(log X)^10` for which the exact valuation above holds.  The cases where
`a` itself is prime contribute only `O_r(X/log X)` further exceptions by the
upper-bound sieve.  Otherwise `a_0>1`, and a nontrivial root `zeta` exists.

For `r=1` the same conclusion is elementary from their displayed

\[
 W_{h,1}(x)=x^{h+1}-(h+1)x+h.
\]

If `p\parallel h`, `h=pa_0`, and `zeta!=1` is an `a_0`-th root of unity,
then

\[
 W_{h,1}(\zeta)=h(1-\zeta),
\]

so its valuation is exactly one.  The same large-squarefree-kernel estimate
removes `O(X log(log X)/log X)` integers without a suitable `p`; primes and
powerful integers are smaller exceptional sets.

Finally, when `a=d+1` is composite, `p<=a/2`, hence `p<=d-3` once `d>=7`.
Cauchy's theorem supplies an order-`p` element.  Every order-`p` permutation
has type `p^t1^f`, and

\[
 f\equiv d\equiv-1\pmod p.
\]

Therefore `f>=p-1`, proving (3)--(4).  The finitely many smaller degrees are
absorbed in (2).  QED

## 2. What this does and does not prove

The new conclusion is a growing wild prime divisor of the Galois group, not
yet a Jordan prime cycle.  Edge (8) occurs around every relevant residue
root of unity.  An inertia element of order `p` can therefore act as a
product of several disjoint `p`-cycles.  The Newton polygon alone does not
show that `t=1` in (4).

This is the precise support gap in the tempting argument

\[
 \text{Newton edge of length }p
 \quad\Longrightarrow\quad
 \text{a }p\text{-cycle}.
\]

The implication is valid if one proves that a single ramified cluster can be
separated from the others in the local splitting field, or if an independent
Frobenius factorization isolates one `p`-cycle.  Neither separation is
currently uniform on the diagonal.

Even a pure prime cycle would leave a second issue: Jordan's theorem requires
primitivity.  Transitivity follows from irreducibility, but transitivity does
not imply primitivity in composite degree.  The two genuinely missing steps
are therefore:

1. **local support:** isolate one `p`-cycle rather than only an order-`p`
   product; and
2. **blocks:** rule out every nontrivial block system in degree `mr`.

If both are proved for a density-one set, Jordan's theorem gives containment
of `A_(mr)`.  The square-discriminant estimate then upgrades the result to
`S_(mr)` for density-one many `m`.

There are two useful reductions once local support is available.  Since the
Newton prime satisfies `p | d+1`, one has `gcd(p,d)=1`.  A pure `p`-cycle in
a transitive imprimitive action must lie inside one block: if it cycled `p`
blocks of size `c`, its support would have size at least `pc`.  Hence every
nontrivial block would have size

\[
 c>p,\qquad c\mid d.                                  \tag{9}
\]

Thus primitivity only has to eliminate the divisors in (9), not every
abstract block size.  In particular, if `d+1=2p`, any order-`p` permutation
in degree `d=2p-1` is automatically one `p`-cycle and no divisor in (9)
exists.  The support and block gaps are both closed on this sparse subfamily
whenever the Newton edge, irreducibility, and `p<=d-3` hypotheses hold.

The detailed block lemma, the exact local-separation target, and a modular
factorization experiment are developed in
[FIXED_R_GALOIS_UPGRADE.md](FIXED_R_GALOIS_UPGRADE.md).  That note also
records why finding a favorable auxiliary factorization experimentally is
not yet a density-one proof: Chebotarev cannot be used to manufacture the
desired cycle type before its presence in the Galois group is known.

## 3. Exact regression

The companion script

```bash
.venv/bin/python scripts/verify_fixed_r_newton_ramification.py
```

checks (1), (5), and (7) on a bounded grid.  For derivative orders
`1<=r<=8` it also checks representative cyclotomic residue fields
coefficientwise: the Taylor coefficients through degree `p-1` are divisible
by `p`, the constant coefficient has exact valuation one, and the degree-`p`
coefficient is a unit.  These are exact certificates for the edge (8).  The
script verifies the algebraic and local identities, not the cited analytic
exceptional-set estimate.
