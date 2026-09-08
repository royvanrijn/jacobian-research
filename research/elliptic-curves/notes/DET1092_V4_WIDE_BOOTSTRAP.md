# V4 wide-bootstrap follow-up on the eight determinant-1092 fibres

## Question

The completed V3 same-parent pilot searched 82 M17 charts on each of eight
prospective determinant-1092 fibres. All eight independently replayed at rank
17 with no gain and no censorship. V4 asks one narrower question:

> Was the first extra direction simply outside V3's small M17 bootstrap schedule?

This is deliberately parallel to the separate rank-jump/incidence investigation.
It does not attempt to explain why a jump exists.

## Frozen experiment

For each of the same eight fibres, V4 starts from the exact ordered generic MW17
seed and excludes every mod-2 parity class touched by that fibre's completed V3
82-chart schedule. V3 charts and parity classes are counted separately: two
centres differing by `2M` may share parity.

From the complete nonzero degree-two quotient it admits generic minimum shells
4, 6, 8, 10 and 12. In the specialized rounded height metric it forms an
outcome-blind union of:

- 64 deepest classes;
- 64 shallowest classes;
- 128 global metric-quantile representatives;
- 32 metric-quantile representatives from each of the five generic shells;
- 96 SHA coverage representatives;
- deterministic SHA fill after deduplication until exactly 512 classes remain.

Only those 512 are mapped to pointed quartics. Search order prefers smaller
quartic coefficient profiles, then a fixed SHA tie-break. No rank, point,
exceptional-direction, catalogue, or prior-search success information enters
that ordering.

Thus the clean null budget is exactly **4,096 fresh charts** across the eight
fibres, in addition to the already completed 656 V3 charts. Chart height,
per-chart timeout, GP binary, generic seed and exact point-search backend are
inherited from each frozen V3 job.

V4 stops a fibre immediately when an exact finite-reduction certificate proves
rank at least 18. It does **not** implement a new cascade. An independently
replayed gain is exported as `v4-discovery.json`, explicitly marked as eligible
for the already-frozen V3 cascade. This separation prevents bootstrap changes
from being confused with later-cascade changes.

## Run

From the repository root:

```sh
git pull --ff-only
cd research
python3 -m unittest discover -s elliptic-curves/tests -p 'test_det1092_v4_bootstrap.py' -q
python3 elliptic-curves/cas/run_det1092_v4_bootstrap.py launch --hours 24
```

The launcher returns after starting one detached controller. A 24-hour value is
a resource lease, not an estimate; the job can be resumed without deleting
checkpoints:

```sh
python3 elliptic-curves/cas/run_det1092_v4_bootstrap.py status
python3 elliptic-curves/cas/run_det1092_v4_bootstrap.py diagnose
python3 elliptic-curves/cas/run_det1092_v4_bootstrap.py resume --hours 24
python3 elliptic-curves/cas/run_det1092_v4_bootstrap.py stop
```

Local evidence is written under:

- `artifacts/local/elliptic-curves/det1092-v4-wide-bootstrap/`
- `artifacts/local/elliptic-curves/det1092-v4-wide-controller/`

No V3, 302, 11952, selection, or generated proof artifact is edited.

## Evidence and replay

Search checkpoints every executed chart. It does not run a full finite-rank
audit after charts that returned no previously unseen rational point. Whenever
a new rational coordinate appears, it runs the exact mod-2 audit immediately;
this permits an immediate certified stop on a genuine new direction. A canonical
final cloud is always audited, followed by independent mod-3 and mod-5 checks.

Replay independently recomputes the complete V4 selector and all 512 mapping
profiles, traverses chart indices numerically, reconstructs every exact map and
point witness, reconstructs the returned point cloud, reruns the finite-rank
certificate checker and the odd-prime checks, and only then writes
`v4-verified.json`.

## Interpretation

A positive case means V3's original M17 schedule was too narrow on that fibre.
The resulting certified subgroup is a natural seed for the unchanged V3
cascade.

Eight clean V4 nulls mean that another 4,096 fresh, parity-distinct M17 charts
failed to bootstrap the eight fibres. This would substantially weaken the
"shortlist depth alone" explanation. It would still **not** prove that any fibre
has rank exactly 17, that no rank jump exists, or that all possible M17 pointed
charts have been searched. The complete quotient is used for candidate
selection, but only 512 new classes per fibre are exposed to the finite point
box.

Resource stops and timeouts remain censored, not negative arithmetic evidence.
No automatic expansion to the forty reserve fibres occurs.
