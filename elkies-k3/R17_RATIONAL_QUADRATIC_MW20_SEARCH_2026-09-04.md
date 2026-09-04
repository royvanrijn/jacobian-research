# Rational-quadratic MW20 search on R17 (2026-09-04)

Operational continuation is paused and recorded in
[`R17_RATIONAL_QUADRATIC_MW20_HANDOFF_2026-09-04.md`](R17_RATIONAL_QUADRATIC_MW20_HANDOFF_2026-09-04.md).

## Status

The requested construction remains **UNKNOWN**.  No rational quadratic
character has been proved to give twist rank at least three, so no
MW20-over-`QQ(t)` surface and no tail-survival claim is promoted here.

The search did close one previously untested same-trace collision gate, added
a native-`074d9` trace-zero subchart, and made the five-control transport test
exactly replayable.  Every negative computation below is kept within its
stated chart.

## Native `074d9` rational-cover screens

The five requested controls are fibres of the `norm12-orbit-074d9` lineage,
not rational fibres of the compact published equation.  The heuristic search
therefore uses the exact eight rational-`PGL2` coordinates in the `074d9`
class.

Three finite-prime screens were run over six disjoint eight-prime blocks from
211 through 491:

| family | covers scored, with coordinate multiplicity | best weakest block | best mean block |
|---|---:|---:|---:|
| compact published coordinate, `q=t^2+b*t+c`, `|b|,|c|<=100` | 40,380 | 0.251 | 0.357 |
| eight exact `074d9` lineage coordinates, same box | 323,040 | 0.256 | 0.370 |
| each of the five controls used as a rational cover-point anchor, eight coordinates, `|b|,|c|<=50` | 407,440 | 0.436 | 0.499 |

The best anchored item occurs in chart `norm12-orbit-08aaa`, above curve 356,
and has

```text
q(z) = 1 + 46*(z-z_356) + 40*(z-z_356)^2.
```

Its exact cover preflight has rational points above 356 and no rational point
above 351, 376, 377, or 385.  This is not an MW20 candidate: no new section is
known, and its score is not a rank bound.  In particular, the heuristic pass
did not authorize a characteristic-zero section solve.

The exact diagnostic preflight is stored as
`artifacts/generated-results/elkies-k3-r17-mw20-control-transport-anchored-heuristic-top-v1.json`.

The artifacts are:

- `artifacts/generated-results/elkies-k3-r17-rational-quadratic-twist-nagao-h100-v1.json`;
- `artifacts/generated-results/elkies-k3-r17-rational-quadratic-twist-pgl8-nagao-h100-v1.json`;
- `artifacts/generated-results/elkies-k3-r17-rational-quadratic-twist-control-anchored-pgl8-nagao-h50-v1.json`.

All three artifacts state explicitly that their scores are heuristic and do
not imply a twist-rank lower bound.

## Exact low-genus collision gate on `11952`

The existing genus-two and genus-three collision ledgers compared distinct
trace masks but did not test repeated covers inside one trace mask.  The new
checker groups the complete smooth atlas and every affine genus-two and
genus-three survivor by the full quadratic cover key, including the scalar
squareclass.

At primes `17,23,29,31,37`, the number of trace masks still capable of a
same-mask triple falls

```text
65, 4, 0, 0, 0.
```

The distinct-mask genus-two gate was already empty.  The four distinct-mask
genus-three pairs surviving the original five primes are all absent in the
targeted `43,47,53` screens.  Consequently no quadratic cover occurs three
times in the processed smooth/genus-two/genus-three affine layers on
`norm12-orbit-11952` with the stated simultaneous-good-reduction hypotheses.

This is an exact low-genus chart result, not a twist-rank upper bound and not a
result on `074d9`.  Bad-reduction denominators, parameter-at-infinity charts,
and higher arithmetic genus remain open.  The certificate is

`artifacts/generated-results/elkies-k3-r17-norm12-11952-low-genus-rank3-cover-collision-v1.json`.

## Native `074d9` trace-zero boundary

The rootless lattice has no norm-two row, so after arithmetic genus three the
next trace-minimum layer is the unique norm-zero, arithmetic-genus-five row.
For the polynomial-x trace-zero subchart

```text
x(u)=x0+x1*u+...+x4*u^4,
x^3+A*x+B=S(u)^2*Q(u),
deg(S)=5, deg(Q)=2,
```

the complete `p^5` census is inexpensive.  At the good prime `p=19`, all
`19^5=2,476,099` polynomials were checked.  Exactly three survive, their three
scalar-sensitive quadratic cover keys are distinct, and hence there is no
pair or triple collision in this subchart.

This does **not** exhaust trace-zero `P.O=0` sections.  A general polynomial
section on the short twist corresponds to `x=X/Q`, not necessarily polynomial
`x`; the missing denominator form is the next exact gate.  The relevant
certificate is

`artifacts/generated-results/elkies-k3-r17-074d9-tracezero-genus5-polynomial-x-p19-v1.json`.

## Five-control transport gate

For a future character `q(z)=a*z^2+b*z+c` in any of the eight exact lineage
coordinates, a control parameter `z_i` lifts over `QQ` exactly when `q(z_i)`
is a rational square.  The script

`elkies-k3/scripts/check_r17_mw20_control_transport.py`

performs those five exact square tests from the certified native parameters.
If twist rank three is first certified, the absolute fibre lower bounds and
the numerical budget beyond an assumed exact generic rank 20 will be

| curve | old jump beyond 17 | rank budget beyond 20 |
|---:|---:|---:|
| 351 | 8 | 5 |
| 356 | 12 | 9 |
| 376 | 5 | 2 |
| 377 | 6 | 3 |
| 385 | 12 | 9 |

The square test is only the transport preflight.  To prove that the displayed
tail itself survives, the three new generic sections must be specialized at
each rational lift and compared exactly with the certified quotient bases
`P18,...`.  If the certified generic rank is larger than 20, the displayed
budget must first be recomputed.

The alternate-Q80 controls (curve 12/395 in class `11952` and curves
363/364/378 in class `08f72`) are not base parameters of the `074d9`
fibration.  Their transport is therefore undefined without an additional
explicit birational common-K3 transport; they are not silently substituted
for the five requested controls.

## Next theorem-directed step

The remaining native target is the full trace-zero denominator chart, most
conveniently written on the short twist as

```text
Y^2 = X^3 + q^2*A*X + q^3*B,
deg(q)=2, deg(X)<=6, deg(Y)<=9.
```

It must be searched for repeated scalar-sensitive `q`, followed by exact
characteristic-zero lifting and a height-pairing rank-three certificate.  A
bounded parameter search or another small-prime count cannot replace those
gates.
