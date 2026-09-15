# Four exact new-generator accessibility witnesses

**Four unused rebuilt charts make a withheld direction visible through an
anchor that requires the newly admitted generator.** The target is outside
every original-bank box tested on the same finite representative dictionary.
Deleting the new generator from the successful anchor also loses visibility.
These are exact retrospective mechanism witnesses, not recovered search gains
or a performance improvement. The preceding
[24-pair comparison](BASIS_AWARE_AMPLIFICATION_2026-09-14.md) remains failed.

The [summary and explicit point/anchor words](../../artifacts/generated-results/elliptic-curves/cancellation_basis_accessibility_v3/summary.json)
cover all18 successful rank18-to19 rebuilds from that comparison. Independent
checks validate4164 maps and60606 coordinate-square evaluations. They certify
44 endpoint targets **individually above their respective rank19 subgroups**;
this does not assert44 mutually independent new directions. Three additional
endpoint classifications remain UNKNOWN. No point searches were added.

Let `M` be the original rank18 subgroup, `G` its newly admitted generator and
`M' = M + ZG` the independently certified rank19 subgroup. For each retained
endpoint target `P` independently certified above `M'`, the finite dictionary
is

\[
\mathcal R_P=\{\pm(P+kG):k=-1,0,1\}.
\]

All original exported anchors, all rebuilt exported anchors with a nonzero
`G` coefficient, and their coefficient-deletion counterfactuals are prepared
before endpoint evaluation in each case. Each uses the same factor-free map
and at most two retained prime neighbours. There are859 original,271 new and
271 counterfactual anchor entries across the18 cases. This is not a full
subgroup, coset, chart or neighbour enumeration.

For an anchor `Q` and horizontal matrix `[[a,b],[c,d]]`, the exact coordinate
of a representative `R` is the primitive pair

\[
[m:n]=[d\,t_Q(R)-b:-c\,t_Q(R)+a],
\qquad t_Q(R)=\frac{y_R+y_Q}{x_R-x_Q},
\]

with the tangent endpoint treated separately. The box height is
`max(|m|,|n|)`, with parameter infinity included. Each witness has an exact
quartic square, a verified inverse point map and native Sage independence.

| Case / family | Minimum in original bank | New anchor, factor-free height | Best new-chart height | Minimum after deleting `G` |
|---|---:|---:|---:|---:|
| `9a7287102d5dc305a8bc` /07ca9 | 193300803059524 | 95474 | 32435 | 146085642906 |
| `6f5afe31891f7dc1938b` /08234 | 597032192702056 | 4457 | 1542 | 406840411 |
| `1f344a0d9056734f634c` /08234 | 440492405376957 | 11800 | 4287 | 418728887891 |
| `259d0615516598ba9cdc` /103b2 | 2954213075 | 93915 | 93915 | 97249789033 |

Every entry compares the same six representatives at bound125000. The new
factor-free column evaluates `P` itself, which is also the winning
representative in all four cases. The best coordinates are respectively
`[-32435:21013]`, `[169:1542]`, `[-3226:4287]`, `[-16837:93915]`.
The first three minima use neighbours at primes5,7,3 respectively; the fourth
uses the factor-free map. Neighbours are not necessary for these four positive
visibility witnesses because their factor-free boxes also contain `P`.

The successful anchors have exact relations `Q=Q0-G`, `Q=Q0-G`, `Q=Q0+G`,
`Q=Q0+G` respectively, with `Q0 in M`. The full integral words and coordinates
are retained. Rank19 independence proves these anchors are unavailable in
`M`. Each native rank20 certificate proves the target still adds a direction
after `G` was admitted. Counterfactual maps use `Q0`, prepared without looking
at the target, and their displayed minima remain above125000. This identifies
the generator, target and exact representation change.

These four charts were exported but **none was searched** in the original
comparison. They occur5th,7th,1st and14th among the respective new-generator
anchors. The audit neither repeats those searches nor credits hypothetical
recoveries. Among the44 target rows,17 are visible in an original-bank chart,
four in a new-generator chart, and those four new-visible rows are
disjoint from the17 old-visible rows. The remaining finite misses do not
exclude another representative of any target coset.

## Relation to the existing height diagnostics

The audit reuses [`rank_growth.py`](../../elkies-k3/scripts/rank_growth.py),
with PARI numerical canonical-height matrices at384-bit precision. In the
table order, the numerical Schur residual heights decrease as follows:

| Case | Above `M` | Above `M'` |
|---|---:|---:|
| `9a7287102d5dc305a8bc` | 21.042189 | 16.925963 |
| `6f5afe31891f7dc1938b` | 21.119721 | 14.693181 |
| `1f344a0d9056734f634c` | 16.731504 | 15.222452 |
| `259d0615516598ba9cdc` | 17.931570 | 14.706684 |

These values describe real-span projection, not integral-coset coordinate
minima. The [existing midpoint theorem](../../elkies-k3/SPECIALIZATION_QUOTIENT_AND_RANK_JUMP_THEOREMS.md#4-midpoints-and-orthogonal-decomposition)
separates the quotient-height contribution from misalignment with `Q/2`.
The exact coordinate witnesses above supply the missing finite accessibility
check; the numerical height reductions do not imply it.

`last_increment_corr=1` is automatic for a nonzero projection onto this
one-dimensional newly added block. It is not evidence of unusually strong
cascade coupling. No target-aware numerical diagnostic enters a selector.

## Reproduction, failures and scope

The first sealed audit stopped before target evaluation: the existing `short`
normalizer also removed an integral scaling from an endpoint already on the
search model. Its map preparation and3.328140 CPU-second failed worker remain.
Version2 requires the supplied endpoint equation to equal the search model;
all18 inputs meet this condition. Its first scan completed, but the checker
passed native Sage rationals into Python Fraction, causing a type failure.
That4.903003 CPU-second worker and its complete scan remain unchanged.

Version3 preserves the first scan verbatim, fixes only the checker's rational
serialization and executes the17 still-unstarted cases. All18 independent
replays pass. Their isolated workers cost123.202766 CPU seconds, including
imports and verification. Predecessor failures and supplemental internal
checker/regression costs are recorded separately in the
[completion receipt](../../artifacts/generated-results/elliptic-curves/cancellation_basis_accessibility_v3/completion.json).
Manual smoke/import and packaging costs are outside these diagnostic component
meters. This is not a complete development-programme timing comparison.

Two [negative replays](../../artifacts/generated-results/elliptic-curves/cancellation_basis_accessibility_v3/regressions.json)
reject removal of the required generator and alteration of a claimed minimum
coordinate, even after updating the temporary outer map hash. All original
protocols and failed prefixes remain in the
[portable replay bundle](../../artifacts/generated-results/elliptic-curves/cancellation_basis_accessibility_v3/portable-replay.tar.gz).
The first byte-checked archive omitted one imported supplemental module; it is
retained. The portable successor adds that dependency and checks all source
hashes, imports and18 case input bindings in an empty root, without repeating
the arithmetic. Its [manifest](../../artifacts/generated-results/elliptic-curves/cancellation_basis_accessibility_v3/portable-manifest.json)
records the precise addition.
After extracting it into an empty directory, a standalone arithmetic replay is:

```sh
sage -python elliptic-curves/cas/retain_cancellation_basis_accessibility.py \
  replay --output /tmp/new-basis-accessibility-replay.json
```

The next performance gate must expose anchors using newly admitted generators
early, with the same initial subgroup, full-cloud certification, box vocabulary
and complete CPU accounting in both arms. It must use previously untested
whole-j controls and evaluate later gains separately from first recovery.
The six-representative oracle dictionary and numerical target diagnostics are
retrospective evidence only. Do not refit residues or suppress unsearched
neighbour coordinates because old radius-model support is zero. A repeatable
later-gain advantage, prediction without missing points, and comparison against
an established external search implementation remain open.
