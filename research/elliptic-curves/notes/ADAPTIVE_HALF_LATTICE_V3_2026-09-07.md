# Calibration-only V3: stratified terminal shortlists

V2 is immutable. V3 changes the two shortlist gates, not the retained anchors,
canonical lane, 125000 search box, ten-second chart budget, incremental exact
rank audits, or stop-on-first-no-gain-epoch rule. This is a calibrated experiment,
not an unconditioned discovery claim or a rank upper bound.

## Frozen rule

Protocol SHA256:
`eed77c4ac724c0132bad7cd39348da75e6b7ebace2a9986f6715b51ff18d5123`.

For every retained anchor, score every extension coset. Select top four and
bottom four Babai norms, two nearest each of sixteen score-quantile centres,
eight deterministic semantic-coset coverage samples, and eight coverage samples
from the lower-score half. All metric boundary ties survive. Exact-CVP every
selected coset, retaining the unchanged V2 minimum-centre convention. Construct
factor-free maps to measure coefficient bit sum and maximum and record exact
minimum multiplicity.

Final selection unions shallow/deep exact-norm Pareto frontiers against small
coefficient bit profiles and high multiplicity, eight exact-norm quantile
centres, eight coefficient-bit quantile centres, and four coverage samples.
All ties survive; cardinalities are not fixed top-k cuts. Shell diversity is
provided by the unchanged sixteen anchors from each of shells eight and ten.
Chart execution ordering and actual-centre deduplication remain V2's.

Coverage uses SHA256 of injective finite-reduction coset fingerprints, never
the coordinate mask integer. Reversing enumeration leaves selection unchanged.
This does **not** prove invariance of LLL/Babai scoring or chart conventions
under arbitrary subgroup rebasing.

First execution starts from a sealed copy of V2's certified M30. An independent
M17 replay is authorized only if that execution certifies rank 31. Execution
denies reads of retrospective artifacts and unredacted exceptional-point data.

The combined protocol hashes both seed files during validation. The M17 branch
initializes only the generic seventeen points; it does not parse the M30 seed
or read the first trial's charts/certificates. It reads that trial's terminal
summary only to enforce the required success gate. The finalizer checks this
cross-trial read boundary. This is a source/dataflow boundary, not a claim that
the M30 seed bytes were unavailable to the process for hash validation.

## Pre-freeze calibration

The terminal diagnostic independently replayed 2353 exact witnesses and verified
4585 sealed V2 files. The retained finite height-one example has Babai rank
285/8192 within its anchor (9602/262144 terminal scores), and would rank 16/17 by
exact centre norm when appended to V2's sixteen refined candidates. It was
excluded at both gates. Its availability does not establish general
anti-correlation between the score and visibility.

In the **oracle-enriched**, nearest-target diagnostic population of 511 distinct
cosets, Spearman correlations are:

| Metrics | Correlation |
| --- | ---: |
| Babai / exact centre norm | 0.242209 |
| Babai / quartic coefficient bits | 0.076543 |
| Babai / retrospective witness height | -0.039311 |
| Exact centre norm / witness height | -0.051319 |
| Quartic coefficient bits / witness height | 0.015240 |

The height-one example ranks 449th descending (63rd ascending) in exact norm
within that population. Its bit sum is 596, ascending competition rank 290
with 100 ties. Multiplicity is identically one, hence has no ranking power.
These are conditional diagnostic ranks, **not full-parity-space ranks**.

In V2's prospectively selected but top-16-truncated population of 512 terminal
cosets, Babai / exact norm correlation is 0.063627. Earlier winning-chart ranks
are recorded individually; their selection bias prevents using them to assess
the success frequency of excluded mediocre-score charts. No population-wide
anti-correlation claim is made. Full-space exact-norm and quartic ranks remain
unmeasured: computing them would defeat the cheap-pass/shortlist experiment.

## Reproduction and evidence

- Runner: [`adaptive_visibility_cascade_v3.sage`](../cas/adaptive_visibility_cascade_v3.sage).
- Pure selector: [`visibility_selection_v3.py`](../cas/visibility_selection_v3.py).
- Independent replay: [`check_visibility_cascade_v3.sage`](../cas/check_visibility_cascade_v3.sage).
- Pre-freeze diagnostics: [`quantify_visibility_shortlists_v3.py`](../cas/quantify_visibility_shortlists_v3.py).
- Terminal decomposition: [`v2_terminal_failure_decomposition_v1.json`](../../artifacts/generated-results/elliptic-curves/v2_terminal_failure_decomposition_v1.json).
- Local frozen protocol, inputs and execution checkpoints:
  `research/artifacts/local/elliptic-curves/adaptive-visibility-cascade-v3/`.
- Local pre-freeze metric table:
  `research/artifacts/local/elliptic-curves/v2-terminal-retrospective-v1/v3-prefreeze-metrics.json`.

Run with Sage Python: `adaptive_visibility_cascade_v3.sage run --start 30`.
Only after certified success: `adaptive_visibility_cascade_v3.sage run --start 17`.
Check each with `check_visibility_cascade_v3.sage --start 30` (or `17`).

## Fixed-M30 outcome

The frozen fixed-state execution certified **30→31 on chart 101**, then
cancelled 1394 stale charts. It scored 262144 cosets, refined 1667 exactly,
and scheduled 1495 charts (the all-ties requirement retains many coincident
coefficient-bit profiles). All 101 executed searches completed. The gain has
independent finite mod-2, mod-3 and mod-5 rank-31 certificates.

The winning finite coordinate is `-42239/7324`, within the unchanged 125000 box.
It is **not** extension 1832. Its displayed extension is 4033, norm-10 anchor
126935, and Babai rank **2286/8192** within that anchor. It entered solely through
the deterministic semantic-coset coverage arm, not a top/bottom or quantile
arm. After exact CVP its norm 181063985 ranks fifth among that anchor's 53
refinements. These labels were read only retrospectively after the frozen
execution found and certified the point.

If inserted into V2's old sixteen exact refinements, this particular winner
would rank third of seventeen and pass its final top-eight cut. Thus this
successful witness directly identifies the **first** shortlist gate as its
obstruction; it does not establish that changing both gates was necessary.
The distinct height-one counterexample above failed both. No gate-ablation
experiment is claimed.

The independent full fixed-M30 schedule replay passed, including every exact
minimum, tie convention, chart transcript and rank certificate. The conditional
M17 execution subsequently completed: **17→18→…→31 in 1169 charts**.
Both runs have independent replay and certificate-to-transcript bindings for
every executed chart in the [final packaged result](../../artifacts/generated-results/elliptic-curves/adaptive_visibility_cascade_v3.json).
Thus the frozen, retrospectively calibrated policy autonomously recovers a
rank-31 subgroup from generic MW17, without exceptional points as execution
inputs. This is not an unconditioned discovery, exact-rank proof, or transfer
guarantee. V1 and V2 remain unchanged.

The subsequent [two-seed amplifier and transfer outcomes](CURVE302_SEEDED_V3_RESULTS_2026-09-08.md)
separate successful known-seed amplification from unsuccessful finite
bootstrap/control exposures.
