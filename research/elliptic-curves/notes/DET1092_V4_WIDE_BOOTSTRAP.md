# V4 wide-bootstrap follow-up on the eight determinant-1092 fibres

## Question

The completed V3 same-parent pilot searched 82 M17 charts on each of eight
prospective determinant-1092 fibres. All eight independently replayed at rank
17 with no gain and no censorship. V4 asks one narrower question:

> Was the first extra direction simply outside V3's small M17 exact-CVP bootstrap schedule?

This is deliberately parallel to the separate rank-jump/incidence investigation.
It does not attempt to explain why a jump exists.

## Repaired quotient construction

The first V4 launch stopped during preflight before any point-search chart. The
failure was useful: `curve302_parent_degree2_multisection_orbits_v1.tsv` is **not**
a serialization of all `2^17` quotient classes. The complete quotient was
enumerated internally when that artifact was built, but the TSV intentionally
exports only the rational/genus-one survivor classes. In the certified parent
basis the TSV therefore contains exactly:

```text
norm 8:  63,922
norm 10: 40,917
norm 12:    139
----------------
total:   104,978
```

The omitted classes include the unique zero class and the norm-4/norm-6
nonsurvivors. Also, `orbit_mask` names the class in the enumeration's reduced
basis, whereas `parent_MW17_w` is transported to the certified parent marking.
Treating the TSV as a complete parent-basis mask table was the preflight bug.
No curve or search result was affected because zero charts had run.

The failed namespace is preserved unchanged. The repaired experiment writes to:

- `artifacts/local/elliptic-curves/det1092-v4-wide-bootstrap-v2/`
- `artifacts/local/elliptic-curves/det1092-v4-wide-controller-v2/`

## Frozen repaired experiment

For each of the same eight fibres, V4 starts from the exact ordered generic MW17
seed and excludes every parent-basis mod-2 parity touched by that fibre's
completed V3 82-chart schedule. V3 charts and parity classes are counted
separately: two centres differing by `2M` may share parity.

The complete quotient is now constructed directly as the integer masks
`0..2^17-1`; it does not depend on a TSV row existing. In the specialized
rounded canonical-height lattice V4 then:

1. constructs the integral LLL transport;
2. Babai-scores **all 131,072 parity classes**;
3. removes zero and every parity already touched by V3;
4. forms an outcome-blind diversified candidate union from:
   - 64 deepest Babai classes;
   - 64 shallowest Babai classes;
   - 128 global Babai-norm quantile representatives;
   - 32 quantile representatives from each exported survivor shell 8, 10, 12;
   - 32 quantile representatives from the complement of the exported survivor set;
   - 96 SHA coverage representatives;
   - deterministic SHA fill to 512 classes;
5. runs the same exact rational-LDL CVP solver used by V3 on those classes.

The TSV is used only as an independent cross-check that its 104,978 exported
parent-basis parities are unique and have the exact 63,922/40,917/139 survivor
counts. It does not create or delimit the complete quotient.

Exact CVP remains fail-closed. If a candidate exhausts the frozen V3 node bound,
it is recorded as a CVP-preparation censor and replaced by the next deterministic
SHA-ordered fresh parity. No point or rank outcome participates in this
replacement. Exactly 512 classes with completed exact-CVP certificates are
mapped to pointed quartics and searched.

Search order prefers smaller exact quartic coefficient profiles, then deeper
exact metric norm and a fixed SHA tie-break. No rank, exceptional point,
catalogue label, or previous search success enters the selection.

Thus the clean null budget is exactly **4,096 fresh exact-CVP charts** across the
eight fibres, in addition to the already completed 656 V3 charts. Chart height,
per-chart timeout, GP binary, generic seed, exact-CVP node bound and point-search
backend are inherited from the frozen V3 jobs.

V4 stops a fibre immediately when an exact finite-reduction certificate proves
rank at least 18. It does **not** implement a new cascade. An independently
replayed gain is exported as `v4-discovery.json`, explicitly marked as eligible
for the already-frozen V3 cascade.

## Run

From the repository root:

```sh
git pull --ff-only
cd research
python3 -m unittest discover -s elliptic-curves/tests -p 'test_det1092_v4_bootstrap.py' -q
python3 elliptic-curves/cas/run_det1092_v4_bootstrap.py launch --hours 24
```

The repaired controller uses its new `-v2` evidence namespace automatically;
do not delete or edit the old failed preflight directory. The launcher returns
after starting one detached controller. A 24-hour value is a resource lease,
not an estimate. Resume without deleting checkpoints:

```sh
python3 elliptic-curves/cas/run_det1092_v4_bootstrap.py status
python3 elliptic-curves/cas/run_det1092_v4_bootstrap.py diagnose
python3 elliptic-curves/cas/run_det1092_v4_bootstrap.py resume --hours 24
python3 elliptic-curves/cas/run_det1092_v4_bootstrap.py stop
```

## Evidence and replay

Search checkpoints every executed chart. It does not run a full finite-rank
audit after a chart that returned no unseen rational coordinate. Whenever a new
coordinate appears, it runs the exact mod-2 audit immediately; this allows an
immediate certified stop on a genuine new direction. A canonical final cloud is
always audited, followed by independent mod-3 and mod-5 checks.

Replay independently reconstructs the direct `2^17` mask atlas, rechecks the
exported survivor TSV as a subset, recomputes the specialized Babai landscape,
the exact-CVP selection and all mapping profiles, traverses chart indices
numerically, reconstructs every exact map and point witness, reconstructs the
returned point cloud, reruns the finite-rank certificate checker and odd-prime
checks, and only then writes `v4-verified.json`.

## Interpretation

Completed result: **eight independently verified 17→17 outcomes**, 512 fresh
charts per fibre (4096 total), and zero certified gains. No automatic V3
cascade or new-parameter expansion was released. The retained summary is
`artifacts/local/elliptic-curves/det1092-v4-wide-bootstrap-v2/summary.json`.

A positive case means that widening exact-CVP coverage at M17 exposes a first
independent direction that V3's 82-chart bootstrap missed. The resulting
certified subgroup is a natural seed for the unchanged V3 cascade.

These eight V4 nulls mean another 4,096 fresh, parity-distinct, specialized
exact-CVP M17 charts failed to bootstrap the eight fibres. Widening this
particular schedule did not repair bootstrap sensitivity. This still does
**not** prove that any fibre has rank exactly 17, that no
rank jump exists, or that all possible M17 pointed charts have been searched:
512 fresh parity classes per fibre are exposed to the finite point box, not the
entire quotient, and alternative `2M` translations remain another visibility
dimension.

Resource stops and timeouts remain censored, not negative arithmetic evidence.
No automatic expansion to the forty reserve fibres occurs.
