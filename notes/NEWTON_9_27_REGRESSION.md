# Published `(9,27)` exclusion as a regression oracle

Primary source: Guccione, Guccione, Horruitiner, and Valqui,
[Increasing the degree of a possible counterexample to the Jacobian Conjecture
from 100 to 108](https://arxiv.org/abs/2204.14178), especially Proposition 4.1,
Theorem 5.1, and equations (5.9)--(5.11).

## Reduction supplied by the paper

The paper shows that the `(9,27)` corner case would yield polynomials with
Newton polygons

\[
N(P)=\{(0,0),(1,1),(6,16),(6,18),(0,18)\},
\]

\[
N(Q)=\{(0,0),(1,0),(9,24),(9,27),(0,27)\}
\]

and bracket `[P,Q]=x`.  Section 5 introduces a Laurent series, changes to
polynomial coefficients `d_i`, and selects nine coefficient equations: six
from the square at exponents `-1,-2,-3,-4,-5,-7`, and three from the cube at
`-1,-2,-4`.

## Reproduced elimination

Eliminating

`d_-10,d_-8,d_-7,d_-6,d_-5,d_-4,d_-3,d_-2`

from those nine equations over `Q` produces a one-element reduced elimination
basis:

\[
18C_3^{23}d_1d_{-1}^{6}F_{-4}
+8C_3^{69}F_{-4}^{3}
+27d_0d_{-1}^{9}=0.
\]

This is exactly equation (5.9), up to term order.  The same elimination ran at
three large primes, but the characteristic-zero output is the proof artifact;
finite fields are only a computation check.

The script also verifies that the displayed quartic

\[
g(y)=-\frac{35-42y+54y^2-81y^3+243y^4}{910}
\]

is separable and is divisible by neither `y` nor `y+1`, and checks the
differential equation from which it arises.

After the paper's valuation substitutions, equation (5.9) becomes (5.11).  If
`dtilde_-1=(y+1)^k`, the two exhaustive cases are machine checked:

- `k>=8`: the three orders at `y=-1` are at least `70`, exactly `66`, and at
  least `72`, so the unique order-66 term cannot cancel;
- `k<=7`: the three degrees are at most `76`, exactly `78`, and at most `75`,
  so the unique degree-78 term cannot cancel.

## Scope

This is a machine reproduction of the explicit Section 5 terminal argument.
It does not independently prove all results cited by the paper to reduce an
arbitrary hypothetical Keller pair to Proposition 4.1.  That distinction is
recorded in the JSON result and must remain visible when this is used as a
validation rung for `(72,108)`.

Run:

```bash
.venv/bin/python scripts/newton_9_27_regression.py
.venv/bin/python scripts/verify_certificates.py
```
