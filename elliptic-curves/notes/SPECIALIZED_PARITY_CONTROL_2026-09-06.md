# Specialized parity sampling recovers the known 29th direction at lower search height

A fixed sample of parity classes in the already recovered 28-point subgroup
of the native `11952` control recovers a 29th independent point at height
125,000. Recovery occurs on chart 18 of 49; all 18 boxes complete, and the
remaining 31 are unsearched after the declared target stop. Exact admission
replay and complete-cloud proofs modulo 2, 3 and 5 confirm lower bound 29.
The union with the original 49-chart control contains 73 points up to sign.

This is a previously known curve. The new centre selection and execution use
only the 28 points already recovered by the earlier generic-17 search. They
read no unrecovered public point, oracle word, or oracle visibility minimum.
The choice of this control and the research hypothesis use its earlier
history, so this is not independent population validation or a new rank record.

## The coverage gap and the fixed policy

The earlier 301-chart follow-ups paired a high-scoring generic parity label
with one cyclically assigned quotient word. They did not explore all parity
classes in the specialized subgroup. Their completed boxes cannot exclude
other rational points or establish specialized covering optimality.

Here the first 2,048 distinct 28-bit SHA256-derived masks with a nonzero
quotient above the generic seventeen form the fixed candidate set. The hash
domain and order are frozen before geometry. Numerical canonical heights at
384 bits are rounded at scale 10^6. A unimodular LLL basis change supports
numerical CVP; every representative is then transported exactly back to the
original point subgroup, with its parity and rounded quadratic norm checked.
The 49 largest computed norms, with mask ties, determine the chart roster.
All rational maps are frozen before searching at height 125,000 and ten
seconds per chart.

The metric and CVP calculations choose search representatives. Exact parity,
unimodular transport, norm recomputation and rational chart identities do not
turn numerical canonical heights or CVP into proofs of global optimality.
Only separate exact point-independence certificates establish rank bounds.

The earlier generic-centre control recovered 28 at the same height. Its
retrospectively chosen chart at height 1,000,000 and sixty seconds recovered
29. The new result therefore shows that changing the centre policy can expose
the missing direction at the lower height on this control. It does not prove
a universal height reduction or that every new curve benefits.

## Prospective test

A separate frozen protocol selects all six rank-27 curves in inventory V12,
sorted by stable ID: 40, 41, 48, 71, 72 and 90. It uses only each curve's own
27 certified points, whose prefix is its exact generic 16- or 17-point basis.
For each curve it samples 2,048 parity classes with a nonzero quotient,
selects 49 by the same numerical policy, and freezes every map file before
any prospective point search. The maximum is 294 boxes at height 125,000 and
ten seconds per chart, with one worker and a per-curve target stop at 28.
The separate protocol introduces no new parameter, trace, public point,
score tuning or result-dependent refill.

The [aggregate certificate](../../artifacts/generated-results/elliptic-curves/specialized_parity_experiment_v1.json)
records the terminal prospective outcomes and links the complete point clouds.
A bounded miss cannot show that a curve has exact rank 27. The control's
success is evidence of one repaired exposure gap, not a guarantee of a new
near-record curve.

All six trials complete all 49 boxes and exact history replay, for 294
completed boxes total. Each trial's seed-plus-returned-point cloud still
certifies 27 modulo 2, 3 and 5. These clouds concern this policy's inputs
and outputs; they are not a union of every historical search on each curve.

| Inventory ID | Family | Parameter | Points in this trial cloud | Lower bound |
| --- | --- | --- | ---: | ---: |
| 40 | `074d9` | `2818/1535` | 125 | 27 |
| 41 | `11952` | `-2448/11` | 115 | 27 |
| 48 | `11952` | `2828/2015` | 39 | 27 |
| 71 | `103b2` | `3726/881` | 97 | 27 |
| 72 | `11952` | `2012/211` | 93 | 27 |
| 90 | `a1-fibration-01` | `-1867/270` | 49 | 27 |

No inventory entry changes. Preparation of all six map files, all point
workers and exact history replays took 325.151 seconds. Full-cloud checks
took 29.878 seconds and exact geometry replay 2.215 seconds. A further
prospective campaign requires a separate mathematical gate; no larger sample
or higher search box is automatically inferred from these misses.

The [portable evidence manifest](../../artifacts/generated-results/elliptic-curves/specialized_parity_evidence_v1.json)
names twenty base archives and the supplement.
`verify_specialized_parity_bundle.py` checks the seven histories and clouds,
the exact specialized parity transports and maps, and the aggregate in an
extracted workspace without new point searching.
