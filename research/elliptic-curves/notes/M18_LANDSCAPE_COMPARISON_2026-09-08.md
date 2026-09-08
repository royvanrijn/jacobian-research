# What distinguishes the initial M18 landscapes?

**Curve 302 has a conspicuously deep extension tail, but the retained features
provide no large-margin separator for amplification across fibres.** In
particular, the successful `1926/2699` state overlaps the bounded null controls.
This is a retrospective diagnostic, not a learned model or a selection rule.

The [full comparison table](../../artifacts/generated-results/elliptic-curves/m18_landscape_comparison_v1/comparison.md),
[103-feature CSV](../../artifacts/generated-results/elliptic-curves/m18_landscape_comparison_v1/comparison.csv)
and [self-contained replay package](../../artifacts/generated-results/elliptic-curves/m18_landscape_comparison_v1/manifest.json)
cover 15 initial states on 12 distinct rational `j`-invariants. All four repeated
states have the identical curve-302 equation, so they are one fibre, not four
independent examples. Four states have completed certified gains, ten have
completed bounded no-gain runs, and one has no sealed outcome at the snapshot.
The gains represent only two distinct fibres.

## Findings

Norms below are divided by the median diagonal entry of the initial generic
M17 height matrix, keeping raw norms in the feature files. The comparison uses
all 64 refined parity cosets, including those whose charts were never executed.

| Initial state | Verified outcome | Normalized CVP maximum | Masks with norm ≤2 | Median quartic coefficient bits |
|---|---:|---:|---:|---:|
| 302, recovered-strict-02 | 18→31 | 4.438 | 13 | 594 |
| 302, recovered-strict-03 | 18→31 | 7.720 | 13 | 594 |
| 302, residual-strict-03 | 18→31 | 7.835 | 13 | 594 |
| 1926/2699 | 18→21 | 2.7022 | 16 | 415 |
| 2980/1967 | 18→18 | 2.6986 | 13 | 459 |
| 1117/2193 | 18→18 | 2.566 | 15 | 429 |
| 5193/35630 | 18→18 | 2.450 | 30 | 483 |
| 2953/1671 | 18→18 | 2.546 | 22 | 460 |
| Six orbit8044 factory seeds | 18→18 each | 2.551–2.697 | 9–19 | 855–1009 |

The curve-302 maxima come from the seed-containing extension (`extension=1`).
Its generic-only extension has the same maximum, 2.369, in all four states.
The three completed 302 states have maximum/minimum ratios 3.075–5.429;
all ten null states have ratios 1.421–1.740. However, the rank-21 success has
ratio 1.568, inside that null range. “More cheap masks” also fails: it has
16, whereas the null small-conic seed has 30.

The seed's orthogonal squared-height component, computed as the Schur
complement of the rounded generic M17 Gram block, gives a similar diagnostic:
2.254–5.709 for the three completed 302 states, but 0.3899 for `1926/2699`
and 0.3900 for the null `2980/1967`. Thus a large orthogonal seed component
accompanies the deep 302 tail; it does not distinguish the two latter fibres.
Raw seed height is additionally sensitive to the chosen representative.

Every one of the 960 cosets has exactly two signed minimizing vectors, hence
one minimizing pair. Shell-8 and shell-10 each contribute 32 pairs per state.
Minimum multiplicity therefore gives no discrimination here. Competitive
anchors and parity concentration overlap too: `1926/2699` has 31 masks on
23 anchors within 1.25 times its minimum; a null factory seed has 40 on 27.
Its smaller quartics may make searching cheaper, but coefficient complexity
is presentation dependent and is not an arithmetic success certificate.

Scanning all frozen scalar features for strict separation of **any gain**
from **bounded no gain** produces only three correlated maxima:
overall CVP maximum, shell-10 maximum, and seed-containing maximum. The overall
boundary is 2.7022 versus 2.6986, only about 0.13% of the adjacent null value.
This is a narrow post-hoc gap, not a useful large-margin result. No ML model,
p-value, tuned cutoff or prospective validation is claimed.

## Feature boundary and exposure

The roster was frozen from available M18 packets and epoch-zero landscapes in
six named campaigns, without filtering on outcome. This adds `2953/1671`, the
newly completed original-funnel `2980/1967`, and two available universality-panel
states to the requested controls. The sealed panel outcome for
`residual-strict-03` is 18→31 in 385 charts; `recovered-strict-01` remains
`UNKNOWN_NO_SEALED_REPLAY` in this snapshot. Its features are retained but it
is excluded from outcome comparisons.

The feature extractor reads exactly the 18-point packet, initial landscape
and its 32 full-score NPZ files. It never reads chart transcripts, later
points, later epochs or outcome summaries. All vectors were hashed into the
feature seal before the separate outcome join. The analyst already knew the
user's outcome descriptions, and curve-302 seed construction was retrospective;
this is an enforced data boundary, not a claim of analyst blinding.

At M18 unchanged V3 has two extension masks at each of 32 selected anchors:
64 exact CVPs, not 200 or 128. The 128th order statistic is explicitly null.
Coverage is complete for those selected anchors, not the full generic shell
or all M18 parities. Shell names refer to generic norms. Quartic statistics
use all 64 refined parity profiles, not just successful or executed charts.
The canonical lane is excluded from quartic statistics consistently.

All campaigns use the same initial anchor, metric, height-125000 and
10-second chart policies. The 302 runs stop at target 31, versus target 32
for the other fibres. All ten nulls exhaust their 114-centre initial schedule;
rank-21 exhausts a later schedule after 568 total charts. These are bounded
outcomes and subgroup lower bounds, not exact ranks or proofs of absent points.

## Reproduction and next use

The package contains all input snapshots and original verified outcome
summaries in `frozen-inputs.tar.gz`, individual feature vectors, exact integer
norm rows, policy, source hashes, feature seal and label bindings. No point
search was launched. All 960 CVPs replayed by exhaustive rational LDL
enumeration, bounded at two million nodes per coset. The exactness is for the
retained integer Gram rounded from numerical heights; this is not an exact
canonical-height theorem. Quartic profiles are retained initial V3 metadata;
this audit does not independently regenerate quartic maps or height matrices.

Run from the repository root:

```sh
python3 research/elliptic-curves/cas/check_m18_landscape_comparison.py
python3 research/elliptic-curves/cas/check_m18_landscape_comparison.py --replay-cvp
python3 -m unittest discover -s research/elliptic-curves/tests -p test_m18_landscape_comparison.py
```

Six focused tests cover rejection of M19 input and mismatched M18 bases,
feature reads restricted to the frozen NPZ inputs, Gram tampering, label
access before sealing, feature tampering, and quantile/entropy conventions.
The lightweight checker verifies archive hashes and original label bindings.

For a next prospective panel, retain the deep-tail and orthogonal-seed-height
features as hypotheses, while admitting controls with ordinary tails as well.
The `1926/2699` versus `2980/1967` pair is the immediate counterexample any
proposed universal explosion score must explain. No production ranking or V3
policy has been changed on the strength of this small comparison.
