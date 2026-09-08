# V4 wide-bootstrap follow-up on the eight determinant-1092 fibres

## Question

The completed V3 same-parent pilot searched 82 M17 charts on each of eight
prospective determinant-1092 fibres. All eight independently replayed at rank
17 with no gain and no censorship. V4 asks one narrower question:

> Was the first extra direction simply outside V3's small M17 exact-CVP bootstrap schedule?

This is deliberately parallel to the separate rank-jump/incidence investigation.
It does not attempt to explain why a jump exists.

## Frozen experiment

For each of the same eight fibres, V4 starts from the exact ordered generic MW17
seed and excludes every mod-2 parity class touched by that fibre's completed V3
82-chart schedule. V3 charts and parity classes are counted separately: two
centres differing by `2M` may share parity.

The complete generic degree-two quotient contains all `2^17 = 131072` parity
classes with the certified minimum-shell histogram

```text
0:      1
4:   1218
6:  24875
8:  63922
10: 40917
12:   139
```

V4 first constructs the specialized rounded canonical-height Gram matrix and
its integral LLL transport. It then **Babai-scores all 131,072 parity classes**
in that specialized lattice. After removing zero and every parity touched by
V3, it forms an outcome-blind union of:

- 64 deepest Babai classes;
- 64 shallowest Babai classes;
- 128 global Babai-norm quantile representatives;
- 32 Babai-norm quantile representatives from each of the five nonzero generic shells;
- 96 SHA coverage representatives;
- deterministic SHA fill after deduplication until exactly 512 classes remain.

Every one of those 512 classes is then run through the **same exact rounded-
metric CVP solver used by V3**, with V3's frozen node limit. The exact minimum is
transported back to the original MW17 marking, parity is rechecked integrally,
and the same semantic minimum-centre tie convention is used. Only then is the
pointed quartic constructed.

Search order prefers smaller exact quartic coefficient profiles, then deeper
exact metric norm and a fixed SHA tie-break. No rank, exceptional point,
catalogue label, or prior-search success enters selection or ordering.

Thus the clean null budget is exactly **4,096 fresh exact-CVP charts** across the
eight fibres, in addition to the already completed 656 V3 charts. Chart height,
per-chart timeout, GP binary, generic seed, exact-CVP node bound and point-search
backend are inherited from the frozen V3 jobs.

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

Replay independently recomputes the full 131,072-class Babai landscape, the 512
selected exact CVPs, semantic minimum centres and all 512 mapping profiles. It
then traverses chart indices numerically, reconstructs every exact map and point
witness, reconstructs the returned point cloud, reruns the finite-rank
certificate checker and odd-prime checks, and only then writes `v4-verified.json`.

## Interpretation

A positive case means that widening exact-CVP coverage at M17 exposes a first
independent direction that V3's 82-chart bootstrap missed. The resulting
certified subgroup is a natural seed for the unchanged V3 cascade.

Eight clean V4 nulls mean that another 4,096 fresh, parity-distinct, specialized
exact-CVP M17 charts failed to bootstrap the eight fibres. This would
substantially weaken the explanation that V3 merely exact-CVP'd too few M17
classes. It would still **not** prove that any fibre has rank exactly 17, that no
rank jump exists, or that all possible M17 pointed charts have been searched:
512 fresh parity classes per fibre are exposed to the finite point box, not the
entire quotient, and alternative `2M` translations of those parities remain an
additional visibility dimension.

Resource stops and timeouts remain censored, not negative arithmetic evidence.
No automatic expansion to the forty reserve fibres occurs.
