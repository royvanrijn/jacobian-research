# Stage D: Newton-directed edge bands

## What is encoded

For the translated chart

\[
x=uv,\qquad y=u-(uv)^{-1},
\]

the implementation now translates every proposed ordinary monomial exactly:

\[
x^a y^b=\sum_{k=0}^b(-1)^k\binom bk
u^{a+b-2k}v^{a-k}.
\]

The union of these Laurent terms is only the initial ansatz.  Exact common-
denominator division independently recomputes the full polynomialization
kernel for `P` and `Q`.  The two coordinates may therefore have different
Newton supports and different kernel dimensions.  No rectangular completion
is inserted.

All systems impose `[A,B]=-u`, equivalent to a polynomial Keller map after
composition with the chart, and the normalization

\[
F(0,0)=F(1,0)=0,\qquad DF(0,0)=I.
\]

The inverse-chart pair `(A,B)=(uv,u-(uv)^{-1})` is present in every kernel and
is checked as a Keller control.  Hence emptiness below is specifically an
incompatibility with the normalized collision, not a failure to represent any
constant-Jacobian map.

## Executed families

| family | Laurent terms P/Q | kernel dimensions P/Q | equations/variables | result |
|---|---:|---:|---:|---|
| `(3,2)` bands, weighted levels 12/18 | 25/43 | 13/21 | 132/34 | exact `Q` basis `[1]` |
| `(2,3)` bands, weighted levels 12/18 | 20/34 | 13/21 | 94/34 | exact `Q` basis `[1]` |
| primitive `(2,7)` square/cube | 24/31 | 11/13 | 234/24 | exact `Q` basis `[1]` |
| reversed `(7,2)` square/cube | 14/16 | 11/13 | 84/24 | exact `Q` basis `[1]` |
| total-degree edges 6/9 | 31/61 | 18/34 | 158/52 | first modular F4 run exceeded 30 s |

For each of the four thin families, the exact characteristic-zero result agrees
with runs at `1000003`, `1000033`, and `1000037`.  The total-degree control is
not declared empty: it is retained as an F4 scaling timeout, and the harder
fields were skipped after the first modular timeout.

Every family passes the coarse Newton test that exponent `(1,0)` can occur in
the outer bracket and every pole-cancellation kernel is nonempty.  Thus none of
these failures is an "earliest uncancellable pole" or a missing Newton edge.
The obstruction is the complete collision-normalized nonlinear ideal.

## Scope and next refinement

These are deliberately small prototypes, not a transcription of the
admissible `(72,108)` chain.  In particular, placing the square and cube of the
primitive `(2,7)` monomial tests an orientation and exponent ratio; it does not
encode all corners, valuations, dicritical branches, or the initial corner
`A0=(8,28)`.

The result supports two concrete decisions:

1. keep thin valuation bands—dense edge completion already crosses the present
   F4 scaling boundary at only 52 parameters;
2. transcribe the excluded `(9,27)` chain completely and make its published
   terminal contradiction a regression test before generating the unresolved
   `(8,28)` chain.

Reproduction:

```bash
.venv/bin/python scripts/stage_d_newton_bands.py --timeout 30
.venv/bin/python scripts/verify_newton_translation.py
.venv/bin/python scripts/verify_certificates.py
```
