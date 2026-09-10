# Generic-rank ceiling fibration pivot

The A1 parent foundry is stopped after its current bounded jobs. Its new
fibrations are genuine, but their first panels have not supplied evidence for
a specialization tail that justifies expanding that MW16 lane.

The next production question is therefore which **ceiling-attaining
fibrations** exist on each active K3 surface. For a Jacobian elliptic K3,
Shioda--Tate gives

```text
rank MW = rho(Xbar) - 2 - rank(reducible-fibre root lattice).
```

Both active surfaces have exact geometric Picard rank 19. Thus their generic
Mordell--Weil ceiling is 17, and every fibration attaining it is rootless.
Neither surface can carry an MW18 fibration. An MW18 construction needs a
different Picard-rank-at-least-20 K3, with a separate arithmetic marking and
equation-realization proof.

## Current ceiling atlas

| Surface | Ceiling | Attainment | Lattice-type census | Production state |
|---|---:|---|---|---|
| `X948` | MW17 | rootless | complete at `O(NS)/J2`: two rank-17 rootless frame types, the published R17 and alternate-Q80 types | equation realization may proceed only through marked-`U` and exact equation gates |
| `X1092` | MW17 | recovered determinant-1092 rootless curve302 parent | complete: 19 rootless J2 types, including the known frame as class 6 | class 1 selected for marked-U realization; equations pending |

The complete `X948` census is an enumeration of frame/isometry types, not an
exact classification of fibration orbits under surface automorphisms: the
known finite J1 count is between two and eight. It nevertheless supplies the
right level for a fibration construction queue. The recovered `X1092` parent
proves that its ceiling is attained. The now-complete
[X1092 census](DET1092_ROOTLESS_J2_CENSUS_2026-09-10.md) independently enumerates
208 primitive embeddings in the Weyl cover and 19 complement types. The 18
additional types require rational marked-U and equation realization.

[`build_fibration_ceiling_atlas.py`](../cas/build_fibration_ceiling_atlas.py)
binds the source certificates and creates a fail-closed preflight artifact.
The historical v1 packet records `X1092_J2_PENDING`; v2 binds the completed
census. New types may enter marked-U construction, while specialization stays
blocked pending an exact rational equation and generic MW17 certificate.

## Completed X1092 proof computation

The exact rank-seven Nishiyama auxiliary with opposite discriminant form is
now retained in
[`det1092_nishiyama_auxiliary_v1.json`](../../artifacts/generated-results/elliptic-curves/det1092_nishiyama_auxiliary_v1.json).
Its complete primitive embedding cover through all Niemeier cases is now
recorded in the linked census, with Weyl accounting, representative Grams,
minimum data, automorphism orders and integral isometries. A clean replay
reproduces all sixteen anchor packets. The equation layer may now attempt
a marked nef `U` for each new type over `Q`.

The existing generic `Genus(...).representatives()` interface is not an
acceptable substitute here: in this environment it requires Magma, which is
unavailable. The census must remain Sage/PARI/Niemeier-replayable, as the
determinant-948 proof is.

## Sources

- [Shioda--Tate and frame classification](../../elkies-k3/RANK_MUTATION_AND_LIFT_THEOREMS.md)
- [complete determinant-948 rootless J2 certificate](../../artifacts/generated-results/elkies-k3-rootless-j2-niemeier-first.json)
- [finite determinant-948 J1 bound](../../artifacts/generated-results/elkies-k3-rootless-j1-uniform-bound-v1.json)
- [exact geometric Picard19 proof for the determinant-1092 surface](CURVE302_PARENT_SEARCH_AND_GEOMETRY_2026-09-07.md)
- [determinant-1092 Nishiyama auxiliary](../../artifacts/generated-results/elliptic-curves/det1092_nishiyama_auxiliary_v1.json)
