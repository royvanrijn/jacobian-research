# Two-tail minimal-excess closure in the scalar `HC4` packet

## Status

This continues `HC4RSD48--49`.  Write

\[
 f=A^2B,\qquad o=\deg B,\qquad k=\deg A=\frac{d-o}{2},
\]

with `h=0`.  The same-weight tail degrees are

\[
 \deg R_j=o-(j-2)k.
\]

> **Theorem HC4RSD50 — complete `h=0` closure for `d>=5o/3`.**
> Every synchronized scalar minimal-excess packet satisfying
> \[
> d\ge\frac{5o}{3}
> \]
> has a fixed ruling and reduces to `HC2` or the exact `JC2` cotangent
> endpoint.  There is no bound on total degree.

The region `d>=2o` is HC4RSD49.  It remains to treat

\[
 \frac{5o}{3}\le d<2o.
\]

## 1. Strict two-tail strip

Suppose

\[
 \frac{5o}{3}<d<2o.
\]

Then

\[
 0<r:=o-2k<k
\]

and the only possible higher tails are

\[
 z^3R_3+z^4R_4,
 \qquad
 \deg R_4=r,
 \quad
 \deg R_3=r+k.
\]

HC4RSD48 first gives, after a constant linear change,

\[
 R_4=x^r.
\]

The highest mixed coefficient then forces

\[
 (R_3)_{yy}=0,
\]

so

\[
 R_3=x^{r+k-1}(\alpha x+\beta y).                 \tag{1.1}
\]

Let `q=\kappa a^2B`; factorial normalizations are irrelevant to the zero
sets below.  The next bordered coefficient is

\[
 4\kappa a^2r(r+4)x^r B_{yy}
 +D_{k,r}\beta^2x^{2r+2k-2}=0,                    \tag{1.2}
\]

where

\[
 D_{k,r}=16k^2-32k-3r^2-12r+16.                  \tag{1.3}
\]

If `\kappa\ne0`, equation (1.2) gives

\[
 B_{yy}=\gamma x^{o-2}
\]

for some scalar `\gamma`; hence

\[
 B=x^{o-2}Q_2(x,y).
\]

But this strip has `k>r>=1`, hence `o=r+2k>=5`; this contradicts
squarefreeness of `B`.

If `\kappa=0`, equation (1.2) becomes `D_{k,r}\beta^2=0`.  Since `k>=r+1`,

\[
 D_{k,r}
 =r(13r-12)+16(k-r-1)(k+r-1)>0,                  \tag{1.4}
\]

so `\beta=0`.  Thus the two upper tails depend only on `x,z`.  The highest
coefficient linear in the next nonzero lower term `P=AB` is

\[
 -4r(r+4)x^{3r-2}P_{yy},                           \tag{1.5}
\]

up to a nonzero scalar.  Therefore `P_{yy}=0`; but `P` is divisible by the
squarefree degree-`o>=5` polynomial `B`, impossible.  This closes the strict
strip.

## 2. Boundary `d=5o/3`

Here `o=3k`, `d=5k`, and admissibility forces `k>=2`.  The complete higher
tails are

\[
 z^3R_{3,2k}+z^4R_{4,k}+t z^5.
\]

If `t=0`, the strict two-tail argument applies.  Suppose `t\ne0`.  The top
coefficient is a nonzero scalar multiple of

\[
 t^2\det H_{R_4},
\]

so normalize

\[
 R_4=x^k.
\]

The next coefficient forces `(R_3)_{yy}=0`, hence

\[
 R_3=x^{2k-1}(\alpha x+\beta y).
\]

The following coefficient is

\[
 k(k-1)\kappa a^2x^kB_{yy}
 -(2k-1)^2\beta^2x^{4k-2}=0                       \tag{2.1}
\]

up to a common nonzero scalar.  For `\kappa\ne0`, this again gives
`B=x^{3k-2}Q_2`, impossible because `o=3k>=6`.  For `\kappa=0`, equation
(2.1) forces `\beta=0`; the next coefficient linear in `P=AB` is a nonzero
multiple of

\[
 k(k-1)t^2x^{k-2}P_{yy},
\]

again impossible because `B|P` is squarefree of degree at least six.

This proves HC4RSD50.

## 3. Emerging general recurrence

The proofs now cover all `h=0` packets through two nontrivial higher tails.
The pattern is stable:

1. the highest tail has zero binary Hessian and becomes a pure linear power;
2. the next tail has transverse degree at most one;
3. the next equation either forces `B` to have transverse degree at most two,
   contradicting squarefreeness, or (when `\kappa=0`) aligns the next tail and
   forces `(AB)_{yy}=0`.

For a longer chain, the natural conjectural invariant is

\[
 \deg_y R_{J-\ell}\le\ell.
\]

If this reaches `q\propto B`, then

\[
 \deg_y B\le J-2,
 \qquad
 J=2+\left\lfloor\frac{o}{k}\right\rfloor.
\]

Since `k>=2`, one has `J-2<o-1` for every `o>=3`, contradicting
squarefreeness.  Proving this recurrence would close **all** `h=0` packets at
once.

## 4. Verification

Run

```bash
.venv/bin/python scripts/verify_hc4_minimal_excess_two_tail.py
```

The checker verifies the strict-strip and boundary identities symbolically;
no total degree or root partition is instantiated.
