# Wide-search equation-only arithmetic profile

Status: **six-worker BASE+LOCAL census launched after Sage commissioning**. The attached
patch was applied as commit `604528fa`. The source is the completed
`artifacts/local/elliptic-curves/broad-rank-v1` campaign, not the earlier live
triangle snapshot. No full arithmetic-population result is claimed yet.

The bound census has 2,080 distinct equations, final tail `23=25,24=4,25=1`,
and 157 follow-up improvements. Every initial rank is bound to its actual
batch-000 packet and replay receipt. Final endpoints are checked against the
completion receipt. This is a provenance replay, not a new point-independence
calculation. The frozen [plan](../../artifacts/local/elliptic-curves/wide-arithmetic-profile-v1/plan.json)
records all 9,185 referenced input/receipt hashes.

The [three-fibre smoke replay](../../artifacts/local/elliptic-curves/wide-arithmetic-profile-v1/SMOKE_VALIDATION.json)
passes all three BASE transports/cubic identities. Two LOCAL rows complete
with Brumer–Kramer local terms 8 and 6, and replay through the integral monic
2-division polynomial, certified PARI maximal order, and local reductions.
The third LOCAL row timed out at 180 seconds and remains UNKNOWN. These are
not three successful local computations or a representative tail result.

The [detached launch receipt](../../artifacts/local/elliptic-curves/wide-arithmetic-profile-v1/launch.json)
records six workers, BASE 60 seconds and LOCAL 180 seconds per curve. The
controller reuses all six smoke checkpoints, including the timeout. CLASS
is disabled, and no point searches are launched. After the full census the
detached wrapper runs `check`; the checker was already tested to reject the
partial smoke. Live progress is in `/tmp/wide-arithmetic-profile.log` and
the stage directories. `SUMMARY.md` is still the partial smoke summary until
the full controller finishes or `summarize` is explicitly invoked.

While workers are active, use the race-free observer from the repository root:

```sh
python3 research/elliptic-curves/cas/status_wide_arithmetic_profile.py
```

The original `status` glob can see a worker's transient `*.tmp.json` file and
race its removal. This observer reads only the frozen population's final
checkpoint names. It leaves the active arithmetic sources and their hashes
unchanged; the post-completion checker has no concurrent worker files.

This controller turns the completed wide R17 specialization search into a frozen
arithmetic census.  It does **not** launch rational-point searches and it does
not interpret the final certified lower bounds as exact Mordell--Weil ranks.

## Purpose

The point-search campaign produced a large, adaptively followed population.
The arithmetic profile asks a different question: do equation-derived
2-descent / cubic-field features differ systematically across the certified
lower-bound strata?

The primary exact local quantity is the Brumer--Kramer term

\[
  \dim_{\mathbf F_2}\operatorname{Sel}_2(E/\mathbf Q)
  \le g(E)+u(E)+n(E),
\]

for curves with no rational 2-torsion and irreducible 2-division cubic, where
`g(E)` is the cubic field class-group 2-rank,

* `u(E)=1` for negative minimal discriminant and `2` for positive;
* `Phi_m` is the set of multiplicative primes with even minimal discriminant
  valuation;
* `Phi_a` is the additive set;
* `n_p` is the number of primes of the cubic field above additive `p`;
* `n(E)=#Phi_m + sum_{p in Phi_a}(n_p-1)`.

The default all-population run **does not compute class groups**.  It records
maximal-order / local data and the equation-derived `u+n`.  A separate bounded
`class-probe` command runs provisional `proof=False` computations only on a
frozen high-rank-lower-bound panel plus deterministic controls.  Such outputs
are explicitly conditional/provisional and are not unconditional rank bounds.

## Commands

First point the normalizer at the frozen wide-search campaign folder:

```bash
sage -python elliptic-curves/cas/run_wide_arithmetic_profile.py probe \
  --source /path/to/frozen-wide-search
```

The probe only reports what it can bind.  Preparation is fail-closed and by
default requires the equation/final-tail fingerprint: exactly 2080 distinct
equations and final lower-bound tail `23=25,24=4,25=1`:

```bash
sage -python elliptic-curves/cas/run_wide_arithmetic_profile.py prepare \
  --source /path/to/frozen-wide-search --expected-improved 157
```

Before the long run, use the deterministic first three rows as a Sage/API smoke:

```bash
sage -python elliptic-curves/cas/run_wide_arithmetic_profile.py run --jobs 1 --limit 3
```

Then launch all six workers into the same folder; completed smoke checkpoints are
reused:

```bash
sage -python elliptic-curves/cas/run_wide_arithmetic_profile.py run --jobs 6
```

Default per-curve time budgets are 60 seconds for BASE and 180 seconds for
LOCAL.  A hard address-space memory cap is disabled by default because Sage's
virtual-memory mappings are host-dependent; `--*-memory-gb` enables one when
desired.  Timeouts and arithmetic failures remain `UNKNOWN` and are retained in
the stage status counts.  Final `check` refuses a partial smoke: BASE and LOCAL
must each have a checkpoint for all 2080 frozen rows.

Optional provisional class-group work is separate:

```bash
sage -python elliptic-curves/cas/run_wide_arithmetic_profile.py class-probe \
  --jobs 2 --class-timeout 120 --class-memory-gb 2
```

Finally:

```bash
sage -python elliptic-curves/cas/run_wide_arithmetic_profile.py check
```

`check` verifies source/population hashes and deterministically rebuilds the
aggregate CSV/JSON/Markdown from the immutable per-curve checkpoints.

## Outputs

`artifacts/local/elliptic-curves/wide-arithmetic-profile-v1/` contains:

* `plan.json`, `population.json` — frozen input/provenance;
* `base/*.json` — exact equation and 2-division data;
* `local/*.json` — exact cubic-field/local Brumer--Kramer terms where completed;
* `class/*.json` — optional provisional class-group probes;
* `profiles.csv`, `profiles.json`, `summary.json`, `SUMMARY.md`, `REPORT.json`.

The summary is intentionally descriptive.  It reports medians by certified
lower-bound stratum and Spearman diagnostics, but no p-values.  The search was
adaptive and the ranks are lower bounds, so ordinary iid significance language
would be misleading.

## Mathematical boundaries

1. `forced_class_2rank_lower_from_known_rank = max(0, rank_LB-(u+n))` is a
   consequence of the known lower bound plus Brumer--Kramer.  It is not an
   equation-only predictor of rank.
2. A provisional `proof=False` class group may support a GRH-conditional
   diagnostic only.  This controller never upgrades it to an unconditional
   class-group or rank certificate.
3. Rational 2-torsion / reducible 2-division cases are profiled but the stated
   Brumer--Kramer implementation is marked not applicable.
4. Missing/timeout local arithmetic is UNKNOWN, never zero.
5. Follow-up improvement is a searchability phenotype.  It is kept separate
   from final certified lower-bound strata.

Primary reference for the local bound: Klagsbrun--Sherman--Weigandt,
*The Elkies Curve has Rank 28 Subject only to GRH*, Section 3.1,
restating Brumer--Kramer Proposition 7.1.

`probe` separately reports how many curves have a bound pre-follow-up stage and
how many are observed to improve.  If that count is the expected **157**, rerun
`prepare --expected-improved 157` to freeze the history phenotype.  Generic
starting rank is deliberately not treated as the initial-search result.
`--expected-tail ''` disables the tail gate only when intentionally profiling a
different frozen population.

## Integration checks and retained commissioning failure

The original eleven ordinary regressions passed. The expanded suite has eighteen
passing tests, covering worker invocation/resume, timeout process-group cleanup,
incorrect worker output, endpoint bindings, altered packets/replays, and missing
completion evidence. Integration corrected
an undefined `memory` variable in the worker launcher. It also removed the
unlabelled-minimum-observation fallback: a collection of ranks without an
identified initial endpoint does not certify a pre-follow-up phenotype.

The first smoke is retained at
`artifacts/local/elliptic-curves/wide-arithmetic-profile-precommissioning-v1`.
Its BASE worker timed out because Sage's inherited `global_minimal_model()`
uses the general number-field path and factors the entire discriminant.
The controller and its active child group were stopped. Source snapshots,
the timeout checkpoint and a stop receipt are preserved there.
The corrected rational-curve `minimal_model()` path completed the same
calculation in about 0.007 seconds excluding process startup. A fresh plan
was frozen before repeating the three-fibre smoke; no failed checkpoint was
silently converted to PASS.

Runtime limits and Sage version are frozen separately per stage; changing the
worker count from one to six does not change per-curve budgets. Workers have
isolated process groups, final checkpoints are written atomically, and altered
worker inputs fail the plan check. The final `check` recomputes JSON, CSV and
Markdown in a temporary directory so a mismatch cannot overwrite its evidence.
The human-readable tables include completed-case denominators, signature/root
counts, conductor sizes, and the number of paired observations per correlation.

Absence of a descriptive local-term association would not by itself locate a
missing signal in g: censoring, adaptive exposure, parent composition and
incomplete local arithmetic remain alternative explanations. Conversely, a
local-term association is not a measurement of the full Selmer dimension.
