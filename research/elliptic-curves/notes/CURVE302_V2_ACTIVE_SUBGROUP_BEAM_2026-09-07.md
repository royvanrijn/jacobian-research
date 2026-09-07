# Curve 302 — target-free V2 active-subgroup beam policy

## Scope

This is a frozen **geometry-only policy design**, not a point-search result.
It uses V2's generic shell/orbit landscape, a certified active subgroup, and
already discovered independent points in a current wave.  It excludes known
exceptional directions, residual labels, catalogue ranks, point-search
transcripts, and all point-search execution from selection.  The sealed run
records `point_searches_run = 0`.

The complete artefacts are

- [selector](../../artifacts/generated-results/elliptic-curves/curve302_v2_active_subgroup_beam_v1.json),
- [post-sealing diagnostic evaluation](../../artifacts/generated-results/elliptic-curves/curve302_v2_active_subgroup_beam_retrospective_v1.json), and
- [the prior fixed `2^14` diagnostic](../../artifacts/generated-results/elliptic-curves/curve302_exceptional_subgroup_landscape_v1.json).

The latter is only an evaluator.  It was not an input to the selector.

## Frozen prospective rule

For a wave of independent points `Q_1,…,Q_k`, form every nonempty intermediate
active subgroup

```text
M + <Q_i : i in S>,     empty != S subseteq {1,…,k}.
```

For each candidate, score a fixed generic anchor catalogue: the four least
finite-reduction-fingerprint actual parity classes from each of V2's shell-8
and shell-10 catalogues.  For every anchor, enumerate every binary extension
parity of that candidate subgroup.  An exact rational LDL ellipsoid CVP
certificate gives the minimum and multiplicity in each parity coset; every
shortest actual centre is then mapped to obtain its reduced quartic
coefficient-bit complexity.

The selection key is lexicographic:

1. maximize the fraction at or below the active baseline's lower-quartile normalized CVP minimum;
2. maximize the fraction at or below its median;
3. minimize mean determinant-normalized CVP minimum;
4. maximize mean `log2` multiplicity of shortest representatives;
5. minimize median, then 90th-percentile, reduced-chart bits;
6. use the finite-reduction fingerprint tuple only to break an exact tie.

The baseline thresholds are fixed before comparing branches.  The beam retains
the best candidate at every represented nonzero increment rank, then fills by
the same key to a hard maximum of four states.  Thus it does not automatically
absorb a multi-point wave and it never gives a coordinate-prefix privilege to
the first newly listed point.

For a finite implementation, entrywise rounding of a decimal height Gram
matrix is *not* itself covariant.  The policy therefore rounds once per
certified active subgroup; a rebase by `U` uses exactly `U G U^T`, rather than
re-evaluating and rounding decimal pairings.  This is part of the frozen
definition, not a post-hoc adjustment.

## Offline `k=3` calibration

V2 certified three consecutive gains after a rank-24 state.  To exercise the
rule without a new search, those three already-certified directions were
coalesced as one hypothetical `k=3` wave.  All seven nonzero intermediate
subgroups were scored.  The frozen beam has masks `4`, `6`, `7`, and `2`: one
rank-25 branch, the best retained rank-26 branch, the rank-27 full branch,
and a second rank-25 state used to fill the width-four beam.  It evaluates
`16`, `32`, or `64` generic parity cosets according to the increment rank and
maps 82 distinct exact-CVP centres.

This synthetic calibration does **not** claim that V2 historically discovered
a simultaneous three-point wave or that these are known-good point-search
branches.  It only seals what the next implementation must do if a real wave
has `k > 1` certified independent directions.

## Basis-invariance check

The selector performed 15 deterministic random unimodular rebases: all seven
candidates once and the four retained beam states twice more.  Every check
passed exactly: the parity-coset CVP norms, multiplicities, actual centre
identities, derived chart-bit values, and whole score signature all agree
after transport of the frozen integer metric.

## Strictly retrospective diagnostic

After sealing the selector, a different program transported the exact rounded
metric already fixed by the 16,384-state exceptional-direction diagnostic to
the seven off-grid V2 subgroups.  This evaluates the documented local-CVP
upper bound and reports its location in the frozen diagnostic distributions.
It does not change the selector or schedule a search.

At the V2 rank-24 baseline, 11 fixed diagnostic directions remain outside the
active integral subgroup and the largest stage-local numerator is
`553,595,235`.  The full coalesced branch (`mask 7`, retained in the beam)
reduces that maximum to `343,550,813`, improves nine fixed directions, and
has a retrospective `log(1+C)` unlock sum of `21.2315733232`.  Other branches
have different retrospective returns—for example mask `3` has a similarly
large `20.4521125421`—which is useful negative evidence: the generic policy
is not tuned on, nor asserted optimal for, the diagnostic labels.

These are finite, nearest-plane-plus-coordinate-descent upper bounds in a
fixed rounded metric.  They are neither exact CVP minima, pointed-quartic
coordinate minima, a predicted search cost, nor a rank statement.

## Replay

```bash
sage -python elliptic-curves/cas/design_curve302_v2_active_subgroup_beam.sage --check
sage -python elliptic-curves/cas/evaluate_curve302_v2_active_subgroup_beam.sage --check
python3 -m unittest elliptic-curves.tests.test_curve302_v2_active_subgroup_beam
```
