# Determinant1092: parameter intake, one seed, unchanged V3

**Production intake and full replay complete; seed confirmation running.**
The [first results](DET1092_FUNNEL_FIRST_SEEDS_2026-09-08.md) certify one M18
from the original90 and two M18 fibres from a separately frozen conic-splitting
follow-up. One of the latter has a standalone verified rank21 subgroup from
an active V3 cascade. A separate
[small-conic example](DET1092_SMALL_CONIC_SEED_2026-09-08.md) now supplies an
independently proved M18 fibre; its114-chart V3 run has no further certified gain.
The objective is to find a fibre with certified rank at least32 through

\[
10^7\text{ parameters}\longrightarrow
10^4\text{ arithmetic candidates}\longrightarrow
10^2\text{ bounded seed attempts}\longrightarrow
\text{certified }M_{18}\longrightarrow\text{V3}.
\]

Ten seeded fibres is a desired yield, not a quota, prediction, stopping
theorem, or permission to refill the population. The intermediate counts
are budget ceilings. This protocol does not extend any existing pilot.
The [machine-readable specification](../protocols/det1092_search_funnel_v1.json)
records these boundaries. Live progress and terminal outcomes belong to the
production controller, not the launch snapshot on this page.

The completed intake contains exactly10million primitive addresses. Its full
draw, equation, fingerprint, deduplication and selection replay passes. The
frozen height quotas retain8,999 distinct arithmetic candidates and90 seed
inputs across bands7 through16; band17 is empty. Intake and replay cost751.96
and789.96 supervised seconds. None of the retained candidates splits the
frozen orbit8044 conic; their bounded quartic seed attempts continue unchanged.
The full data contain three conic-splitting fibres outside the original
shortlist. Exact follow-up certifies two extra directions; both conic points
on the third fibre have verified words in its inherited17-point subgroup.
This follow-up preserves the original selection and its yield accounting.

The [post-search results adapter](../cas/det1092_funnel_results.py) independently
rechecks any completed V3 epoch with rank above22 and compares its equation
with frozen catalogue/inventory snapshots. Its conductor stage uses exact
local Tate and independent PARI exponents through prime10,000, retaining a
proved conductor divisor, upper bound and unresolved cofactor. Large upper
bounds alone never exclude an improvement. A nonminimal conductor37 control
and three comparison-bound tests pass. The first failed Sage uniformizer
attempt is retained; requesting principal prime generators over `Q` fixes
the generic local-reduction branch without global factorization.

## Implemented run and validation

The [controller](../cas/run_det1092_funnel.py) supervises one worker, preserves
every attempt, and shares finite job/campaign allocations across resumes.
The production configuration has ten million accepted primitive parameters,
at most50million draws, and parameter-height shells with upper endpoints
`2^12,2^16,2^20,2^24,2^28,2^32`. Arithmetic-height bands7 through17 each retain
up to810 score-selected and90 SHA-selected controls. Each band supplies at
most eight score-selected seed attempts and one control: ceilings9,900 and99,
respectively. Empty strata are not refilled.

The current arithmetic interface computes smoothness, exact `j` heights,
fixed-prime traces/cubic splitting at primes5 through97, and splitting of
the frozen orbit8044 conic after exact parameter transport. Strict, Selmer,
extra-direction Kummer and `xi` fields are explicitly `UNKNOWN`; only the
frozen local trace sum schedules fibres. This does not claim that the missing
discriminator has been solved or that this score predicts seed yield.

Seed confirmation first specializes the conic, then uses the retained
2,048-parity/49-centre policy with the384-bit rounded metric, height125,000
and ten seconds per chart. Its geometry freezes before quartic search.
The first certified extra point ends admission and cancels later charts.
Independent replay recomputes the conic transport, geometry, actual returned
point witnesses and finite-group proof. V3 receives exactly18 points and
keeps the original numerical source files and settings. Queue labels are
implemented; no additional deep/aggressive allocation is silently attached
to them in this first run.

Declared limits are3GiB RSS,8GiB intake database storage, four hours each for
intake and full replay,1,200seconds each for a seed attempt and its replay,
and14,400seconds combined for a V3 attempt and replay. The controller has a
48-hour total supervised allocation. Interrupted reservations conservatively
consume their allocation when final supervision data are unavailable.

The [validation package](../../artifacts/generated-results/elliptic-curves/det1092_funnel_validation_v1/manifest.json)
retains the final source snapshot and certificates. Twenty-one focused tests
pass, including interruptions on both sides of a transaction commit,
duplicate and singular intake, tampered checkpoints, exact height/parameter
transport, nonminimal local residues, unresolved-feature/queue handling, rejection of an
off-curve point, and stopping before consuming another candidate after the
first certified extra point.
The final256-parameter smoke run passes all nine supervised stages: independent
PARI checks on1,078 projective table entries, full intake replay, formula-derived
conic rank18 certification and repeated exact V3-map preflight, and two bounded
seed misses with independent replay. The smoke point budget is one chart at
height100/two seconds with64 sampled parities. It does not run a full V3
cascade or contribute a new seed to prospective yield.

The initial smoke's JSON tuple/list replay failure and all outputs remain
under their original local run path. The final smoke uses a fresh namespace;
the earlier bound and equations were not silently changed.

From the repository root:

```sh
python3 research/elliptic-curves/cas/run_det1092_funnel.py status \
  --directory research/artifacts/local/elliptic-curves/det1092-funnel-production-v1
python3 research/elliptic-curves/cas/run_det1092_funnel.py resume \
  --directory research/artifacts/local/elliptic-curves/det1092-funnel-production-v1
```

`resume` reuses committed intake and certified stages within remaining frozen
limits. It rejects changed source/input/runtime bindings. A fresh smoke replay
uses `freeze --profile smoke --directory NEW_LOCAL_PATH`, then `run` against
that path. Its output must remain separate from production statistics.

## Mathematical entry and exit gates

The generic parent has full geometric and arithmetic MW rank17. Use the
[verified reduced parameter chart](DET1092_REDUCED_PARAMETER_CHART_2026-09-07.md)
for affordable evaluation, retaining its exact map from the scan coordinate
`s=a/b` to the original parent coordinate `t`. Bind every fingerprint,
section, point and certificate to that coordinate map and exact equation.

Write `H_s=sp_s(M17)`. Seed confirmation succeeds only when an explicit
rational point `P` and an independently replayable certificate prove

\[
\operatorname{rank}H_s=17,
\qquad \operatorname{rank}\langle H_s,P\rangle=18.
\]

Nonmembership in a displayed integral subgroup is insufficient: `P` might
only saturate it. Strictness at a prescribed collection of local places is
not required either. The
[historical first-unlock certificate](DET1092_FIRST_UNLOCK_JACOBIAN_CLASS_2026-09-08.md)
records that the first302 gain is not strict for the frozen bad-place test.
The success criterion is independence in the rational Mordell--Weil span.

If specialization loses inherited rank, retain `INHERITED_RANK_UNRESOLVED`
or `INHERITED_SPAN_BELOW_17`, as the evidence warrants. Such a fibre cannot
enter this single-extra-point route until its rank17 input is established;
it is not thereby excluded from high rank or other search routes.

## Stage1: equation-only mass intake

Generate at most ten million distinct primitive scan addresses, using a
versioned deterministic stream across logarithmic parameter-height shells.
Use disjoint shells, a fixed domain/seed, positive denominators, exact gcd
reduction, and a finite draw cap. Persist rejected draws and their reasons
through reproducible counter intervals and block digests. Parameter height
controls sampling coverage; it is not the final arithmetic comparison.

For each fibre compute the smoothness test and reduced rational `j=N/D`.
Stratify by `max(bit_length(abs(N)),bit_length(D))`, in frozen64-bit bands.
Also record numerator and denominator heights separately and the integral
model coefficient heights. This avoids treating a small rational coordinate
or a small `j` numerator as an intrinsically cheap fibre. Equal `j` values
form comparison buckets, not proofs of rational isomorphism: retain twists
unless an exact rational model transport identifies them. Keep deduplication
and any source symmetries explicit.

The per-fibre record has the following fields:

| Data | Interpretation and cost boundary |
| --- | --- |
| Exact model, address/transport, discriminant nonzero | Mandatory arithmetic identity and smoothness; no full discriminant factorization or global minimal model required. |
| Fixed-prime cubic splitting and trace data | Cache parameter-residue tables where valid. Ambiguous/nonminimal residues require exact local scaling or an explicit unresolved marker. |
| Strict, Selmer and Kummer fingerprints | Versioned mathematical definition, local support, inherited subgroup reference, completeness status and certificate provenance. Unimplemented or incomplete fields are `UNKNOWN`. |
| Low-degree multisection splitting | Evaluate an already frozen atlas; record exact rational square/splitting conditions. Splitting alone is not independence. No new per-fibre multisection enumeration. |
| Condition characterizing `xi` | Reserved until an applicable, proved condition is available. A retrospective class or a supplied-point recognizer is not a prospective selector. |
| Scheduling score and cost estimate | Heuristics, with all components and missing-data semantics recorded. Neither a rank lower bound nor an upper bound. |

Local calculations must have declared small bounds. Full number-field
class groups, full Selmer groups, unbounded bad-prime discovery and quartic
point searches are outside this mass stage. A costly promising fingerprint
belongs on the retained arithmetic shortlist with a separate finite budget.

Existing work imposes useful negative controls. Every smooth member of the
historical RR net is already rationally pointed, and the
[generic-point RR controls](DET1092_RR_GENERIC_POINT_SPECIFICITY_2026-09-08.md)
produce extra Jacobian classes from generic elliptic points. Consequently RR
solubility, Selmer-set nonemptiness, or a Jacobian rank increment cannot be
used as an established extra-elliptic-point criterion. The
[trace/twist transfer](DET1092_TRACE_TWIST_KUMMER_OBSTRUCTION_2026-09-08.md)
also returns an inherited Kummer class on the old comparison panel.

Retain at most10,000 arithmetic candidates across the frozen height strata.
The feature definitions, weights, per-band quotas and deterministic tie
rules must freeze before production intake. A small benchmark can measure
cost and test exact table/scalar agreement before that freeze. If a later
`xi` theorem changes the selector, version a new selector and replay the
equation-only data; do not silently edit the existing ranking.

## Stage2: bounded confirmation of exactly one extra direction

Freeze at most100 distinct seed-search inputs from the arithmetic shortlist.
Thus the arithmetic shortlist is at most `10^-3` of the mass population,
and the point-search set at most `10^-5`. Keep a small, predeclared
height-matched SHA-selected control allocation within those100 to measure
whether the fingerprints enrich seed yield. Controls use the same budgets.
No point outcomes may refill the selection or tune its weights.

Specialize and certify the inherited17 points before point construction.
Try exact specialization of already available low-degree multisections
first. Otherwise run a bounded extra-point constructor with a frozen centre
policy, chart count, height, per-chart timeout and total per-fibre cap.
The implemented numerical choices and their bounds are frozen above; there
is no implicit unlimited search.

On the first replay-certified extra point, seal exactly `M17 + P18`, cancel
remaining seed charts, and enqueue amplification. If a batch returns several
points, preserve the complete raw batch but select the first independent
point in a frozen order; later returned points do not enter the initial
amplifier basis. Stopping occurs at the first certification checkpoint.

The seed packet includes the ordered17-point specialization, `P18`, exact
equation, model transports, point identities, rank certificate and independent
replay, plus all input/source hashes and search bounds. Existing finite-group
certificates are sufficient when they close independence; failure of a chosen
mod-2 certificate is unresolved certification, not a proof of dependence.

Use distinct outcomes:

- `CERTIFIED_M18`: eligible for V3 after replay.
- `BOUNDED_NO_CERTIFIED_SEED`: completed finite attempt, no rank upper bound.
- `INHERITED_RANK_UNRESOLVED` / `INHERITED_SPAN_BELOW_17`: unavailable rank17 input.
- `CENSORED` / `ERROR`: unfinished computation, never a mathematical miss.

The [proved conic progression](DET1092_CONIC_SEED_PROGRESSION_2026-09-08.md)
already supplies infinitely many certified rank18 specializations. Use it
as a separately labelled constructive lane and a positive handoff test.
It has a different parameter distribution and height profile, so its seeds
must not be counted as successes of the broad population ranking. Its
sufficient progression criterion is not a characterization of all seeds.

## Stage3: amplification from the sealed M18

Feed only the exact certified18-point packet into the unchanged V3 numerical
landscape, chart selection, bounded backend and adaptive update policy.
Bind the engine source hashes and the calibrated numerical settings before
execution. Use the normal source-bound checkpoints and independent replay
of actual point clouds. Do not carry a retrospectively known later point
into the amplifier seed.

Queue classification uses independently certified lower bounds:

| Certified lower bound | Queue |
| --- | --- |
| 18--19 | Retain bounded outcome. |
| 20--27 | Interesting; preserve the20--24 milestone explicitly. |
| 28--29 | Deep queue. |
| 30--31 | Aggressive queue. |
| At least32 | Target lower bound achieved; export and independently replay the certificate. |

Queues schedule resources; they do not change V3 numerical policy. Deeper
or aggressive treatment needs its own declared finite allocation and
checkpoint continuation rule. Do not restart an exhausted run under the
same protocol, relabel old boxes as fresh work, or imply that a queue entry
authorizes unbounded computation. Rank at least32 needs no completed descent;
exact rank32 would require a matching unconditional upper bound.

## Checkpoints, audit and implementation boundary

The intake stores append-only compressed equation blocks and checkpoints in
SQLite. A single WAL transaction commits the block, monotone draw cursor,
accepted-address count, complete deduplication index, height histograms and
retained heaps. Digests bind every block and its running chain. A full replay
reconstructs all draws, rejections, records, selections and the complete
deduplication index in a separate bounded disk-backed database. The crash
tests verify equality with uninterrupted execution. Changing a frozen
protocol or source stops resume rather than silently changing selection.

Before launch, freeze the listed unresolved configuration, all finite wall,
memory and worker limits, software versions, source hashes, exact inputs,
and the replay plan. Test resume, duplicate-address handling, singular
fibres, nonminimal local residues, inherited specialization loss, missing
fingerprints, a false/non-independent seed, and first-seed cancellation.
Run a small, separately labelled end-to-end smoke test, including one known
constructive seed and one negative/control intake. These are implementation
checks, not prospective discoveries or calibration data for selection.

Existing components to reuse, without changing their frozen campaigns:

- [`det1092_record_scale_selection.py`](../cas/det1092_record_scale_selection.py): homogeneous equations, exact `j` heights and corrected residue-table scoring. Its fixed population, paths and non-resumable intake are not yet this funnel.
- [`det1092_conic_seed_v2.sage`](../cas/det1092_conic_seed_v2.sage): exact constructive positive seed lane; its original-coordinate output needs exact transport before reduced-coordinate use.
- [`det1092_v3_contract.py`](../cas/det1092_v3_contract.py): exact specialization and source/model bindings.
- [`det1092_v3_worker.py`](../cas/det1092_v3_worker.py): checkpointed amplifier loop.
- [`run_curve302_seeded_v3_amplifier.py`](../cas/run_curve302_seeded_v3_amplifier.py): prior M18 handoff mechanics, with a retrospective oracle seed constructor that must not be reused for prospective selection.

The new [intake](../cas/det1092_funnel.py),
[first-seed/V3 adapter](../cas/det1092_funnel_worker.py) and
[prospective M18 replay](../cas/det1092_funnel_replay.py) implement these stages
in a separate namespace. The old numerical engine, controllers, protocols
and point outcomes are unchanged. The mathematical `xi` selector remains
open, and complete production cascade performance is not established by the
smoke checks. No mathematical-status entry changes: no new theorem, seed,
or rank is claimed.
