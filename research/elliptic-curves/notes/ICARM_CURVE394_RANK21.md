# ICARM curve 394: rank at least 21 at `log(N)=166.252...`

## Result

ICARM curve 394 is the compact Elkies rank-17 K3 fibre at `t=3/8`.  The
repository-local checker
[`verify_icarm_curve394_rank21.py`](../cas/verify_icarm_curve394_rank21.py)
specializes the pinned compact family and its seventeen exact sections,
globally minimizes the fibre, and obtains the public model coefficient for
coefficient:

```text
y^2 + x*y = x^3
              - 354803089674674467206754048738095*x
              + 2558194545892203175112161719607326368645810580537.
```

All 21 public points are replayed exactly.  For an independent rank proof, the
checker uses the seventeen specialized generic sections and public points
1--4.  Their images in the good-reduction quotients `E(F_p)/2E(F_p)` have
combined binary rank 21 at

```text
11,17,31,37,41,43,53,61,67,73,83,97,101,107,109,131,137,149,157.
```

Good reduction at 19 has group order 25, excluding rational 2-torsion.  The
usual infinite-descent argument therefore proves unconditionally

```text
rank E(Q) >= 21.
```

This is a lower bound only; no Selmer upper bound or exact-rank claim is made.

## Exact conductor replay

The minimal discriminant factors as

```text
2^9 * 3^8 * 5^6 * 7^4 * 13^5 * 23^3 * 29^4 * 89^2
* 43207 * 226549
* 22823593909227592035983291
* 44013936637595415741483513793.
```

PARI proves every displayed factor prime, returns the identity local minimal
change at all twelve bad primes, and reconstructs the exact conductor

```text
N = 1593562111507190066539814084004447718921281851572777685020200143306222910
log(N) = 166.2520985277272016652232895273070674463...
```

All reductions are multiplicative:

| `p` | `v_p(Delta)` | Kodaira | `f_p` | `c_p` | local sign |
| ---: | ---: | --- | ---: | ---: | ---: |
| 2 | 9 | `I9` | 1 | 9 | -1 |
| 3 | 8 | `I8` | 1 | 8 | -1 |
| 5 | 6 | `I6` | 1 | 6 | -1 |
| 7 | 4 | `I4` | 1 | 4 | -1 |
| 13 | 5 | `I5` | 1 | 5 | -1 |
| 23 | 3 | `I3` | 1 | 3 | -1 |
| 29 | 4 | `I4` | 1 | 4 | -1 |
| 89 | 2 | `I2` | 1 | 2 | +1 |
| 43207 | 1 | `I1` | 1 | 1 | -1 |
| 226549 | 1 | `I1` | 1 | 1 | -1 |
| 22823593909227592035983291 | 1 | `I1` | 1 | 1 | +1 |
| 44013936637595415741483513793 | 1 | `I1` | 1 | 1 | -1 |

The Tamagawa product is `207360` and the global root number is `-1`.
Multiplying `p^f_p` gives the displayed conductor exactly.  Independently of
the decimal approximation, `N<10^73` and `log(10)<231/100` give
`log(N)<168.63<173.25`.

Curve 394 therefore improves the repository-local rank-at-least-21 conductor
anchor from curve 285's `173.2515...` to `166.2521...`.  It remains above
ICARM 245's rank-at-least-20 conductor `150.6689...`; finding one more
independent point below that line is still the stronger target.

## Replay

```sh
.venv/bin/python elliptic-curves/cas/verify_icarm_curve394_rank21.py --check
```

The generated certificate is
[`icarm_curve394_rank21_v1.json`](../../artifacts/generated-results/elliptic-curves/icarm_curve394_rank21_v1.json).
