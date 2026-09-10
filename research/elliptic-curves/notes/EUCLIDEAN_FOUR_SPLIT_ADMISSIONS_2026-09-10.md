# The four certified-M17 Euclidean split UNKNOWNs are inherited

All four split points in the corrected frozen roster are
`INHERITED_RATIONAL_SPAN`. Independent exact replay verifies **integer
relations** in their original seventeen specialized generic sections.
Three chains close after two successful halves, and one after three,
within the predeclared eight-layer cap. No candidate enters seeded V3.

This applies the existing
[halving-or-cycle theorem](DET1092_SPLIT_SEED_HALVING_DICHOTOMY_2026-09-08.md).
It is not a specialized rank upper bound. The
[status authority](../../MATH_STATUS.json) and
[four-case summary](../../artifacts/generated-results/elliptic-curves/euclidean_split_admission_v2/summary.json)
bind the result.

## Exact receipts

Write `P1,...,P17` for each fibre's original ordered, signed generic points
on the archived short model. Write `Q` for the positive-square-root conic
branch, transported by `(x,y) -> (d^4 x,d^6 y)`, where `d` is the denominator
of the reduced parameter. The checker reconstructs these conventions from
the frozen generic parent and conic maps.

| Parameter | Orbit mask | Sealed places / quotient rows | Successful halves | Exact integer relation |
|---|---:|---:|---:|---|
| `2` | 82931 | 160 / 133 | 2 | `Q = P2-P10-P11+P13-P14+P16` |
| `-4/3` | 65035 | 160 / 124 | 2 | `Q = P1+P5+P6-P9-P10+P12-P14-P17` |
| `-5/7` | 82931 | 157 / 123 | 2 | `Q = -P3-P6-P12+P13+P16+P17` |
| `2/7` | 30223 | 159 / 132 | 3 | `Q = -P1+P2+2P3-P4+P5+P6-2P7+P8+P9-P11+P12-P14-P15+P16-P17` |

The initial uniquely compatible binary words are:

```text
2:     (0,1,0,0,0,0,0,0,0,1,1,0,1,1,0,1,0)
-4/3:  (1,0,0,0,1,1,0,0,1,1,0,1,0,1,0,0,1)
-5/7:  (0,0,1,0,0,1,0,0,0,0,0,1,1,0,0,1,1)
2/7:   (1,1,0,1,1,1,0,1,1,0,1,1,0,1,1,1,1)
```

Every parity word, exact subtraction target, duplication quartic and rational
half is retained in the receipts:
[2](../../artifacts/generated-results/elliptic-curves/euclidean_split_admission_v2/cases/4/receipt.json),
[-4/3](../../artifacts/generated-results/elliptic-curves/euclidean_split_admission_v2/cases/17/receipt.json),
[-5/7](../../artifacts/generated-results/elliptic-curves/euclidean_split_admission_v2/cases/43/receipt.json),
[2/7](../../artifacts/generated-results/elliptic-curves/euclidean_split_admission_v2/cases/46/receipt.json).
Each has a neighbouring `frame.json` and `independent-replay.json`.

For `Q_0=Q`, the recurrence is `2Q_(n+1)=Q_n-T_n` and
`W_n=sum_(j<n) 2^j T_j`. The checker verifies
`Q=2^n Q_n+W_n` at every step. The first three cases have `Q_1=Q_2`;
the fourth has `Q_2=Q_3`. Consequently
`(2^(b-a)-1)Q=2^(b-a)W_a-W_b` has multiplier **one** in all four cases,
giving exactly the displayed relations. The companion branch has inherited
trace minus `Q`, also checked by exact group law.

## Bounds and independent replay

The [immutable protocol](../../artifacts/generated-results/elliptic-curves/euclidean_split_admission_v2/protocol.json)
copies the original plan, generic-only parent, input manifest, selected conic
receipts, incidences, admission requests, M17 packets and supervisor seals.
At freezing, every archived rank-17 packet equals its original generic packet.
No later points enter the reference or candidate.

Every finite place and quotient row comes from those packets. The adapter
neither adds primes nor shortens the footprint. Each complete footprint has
column rank 17 and an odd-order reduction excluding rational 2-torsion.
Thus each binary parity is unique. **An in-span finite column alone
contributes no dependence conclusion.**

Bounds: eight successful halvings per case, 120 seconds per classifier,
180 seconds per independent replay, 1.5 GiB RSS, at most two workers.
All classifiers completed in under one second each; the independent
replays completed in under four seconds each on Sage 10.9. Execution receipts
retain measured times and memory use.

The producer uses the existing exact quartic-halving engine. The independent
checker imports neither it nor the prospective wrapper. Sage finite elliptic
groups replace the producer's integer finite-group kernel; manual Fraction
arithmetic replaces Sage's rational elliptic group law. The checker rebuilds
the generic sections and trace words, verifies the complete Euclidean lift,
rational square, original branch and scaling, checks all finite quotients and
nine exact halves, and derives and verifies the final relations.
Nonhalving replay also supports complete exact rational quartic
factorizations; none is needed for these inherited conclusions. Polynomial
arithmetic uses Sage in both implementations.

## Reuse and reproduction

The [prospective routine](../cas/prospective_split_admission.py) accepts an
explicit sealed frame and rational point. It returns exactly the three
requested classifications and retains the pending unique subtraction at the
step cap. Invalid or deficient frames raise an error, never an independence
obstruction. The caller must impose execution limits, as the frozen runner
does.

```python
from prospective_split_admission import frame_from_packet, ProspectiveAdmission
from replay_split_admission import ReplayFrame

frame = frame_from_packet(generic_only_packet)
receipt = ProspectiveAdmission(frame).consider(split_point, max_steps=8)
independent = ReplayFrame(frame).replay(receipt)
```

The committed inputs replay without the ignored original pilot directory:

```sh
for index in 4 17 43 46; do
  timeout 120s sage -python research/elliptic-curves/cas/run_frozen_euclidean_admissions.py case --index "$index"
  timeout 180s sage -python research/elliptic-curves/cas/verify_frozen_euclidean_admissions.sage --index "$index"
done
timeout 25s sage -python -m unittest discover -s research/elliptic-curves/tests -p test_prospective_split_admission.py -v
timeout 25s sage -python -m unittest discover -s research/elliptic-curves/tests -p test_split_seed_descent.py -v
```

Ten new admission/replay tests and eight existing descent tests pass. They
cover finite escape, global nonhalving despite an in-span column, independence
hidden by two halvings, an odd-index cycle, integral relations, zero, UNKNOWN
bounds, deficient frames, preserved places, tampering and no artifact reads
in the arithmetic APIs.

The five generic-gate failures at `1,1/2,-3/5,4,-7/3` remain unchanged.
Earlier search reporting described the four finite-column misses too strongly
as certified failures to add a direction. Their original receipts correctly
retained UNKNOWN; the relations above now supply the missing proof.
This removes all four certified-M17 candidates from this frozen seed intake
and gives no reason to expand this lane on their account. It does not exclude
other points on these fibres or productive splits elsewhere in the Euclidean
atlas. No parameter/orbit expansion or seeded V3 run was performed.
