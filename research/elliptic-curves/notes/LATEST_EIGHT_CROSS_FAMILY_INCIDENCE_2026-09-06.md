# Exact cross-family incidence for the latest eight curves

The eight additions in the [47-curve inventory](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v4.json)
have now been checked against all twelve recorded presentations: six compact
MW17, five compact MW16, and the published R17 model. These are presentations,
not twelve asserted independent generic structures.

The [96-pair certificate](../../artifacts/generated-results/elliptic-curves/latest8_cross_family_j_incidence_v1.json)
and [Sage-free replay](../../artifacts/generated-results/elliptic-curves/latest8_cross_family_j_incidence_replay_v1.json)
prove 88 rational-preimage exclusions. The remaining eight cases have complete
rational-root factorizations with only the known parameter in the original
family, and no preimage at infinity. In particular, neither new rank-27 curve
acquires another recorded generic point structure by this route.

For a modular exclusion, the reduced rational j-map is checked to remain a
morphism of the same projective degree; its image then excludes the target.
For each surviving equation, the rational factors and all residual factors
are multiplied back exactly. A finite prime without projective roots excludes
a rational root of each residual factor. Equal j alone is never promoted to
rational isomorphism or point independence.

The bounded audit took 3.551 seconds and its exact replay 0.386 seconds. The
original 32-curve proof, seven-curve extension and this eight-curve extension
cover **384 + 84 + 96 = 564 = 47 × 12** pairs. This closes the missing coverage
for the current inventory only. Other families, future specializations,
nongeneric rational points and rank upper bounds remain outside the claim.

```sh
python3 elliptic-curves/cas/replay_latest8_cross_family_incidence.py \
  --input artifacts/generated-results/elliptic-curves/latest8_cross_family_j_incidence_v1.json \
  --output /path/to/new-replay.json
```

The [small evidence supplement](../../artifacts/generated-results/elliptic-curves/latest8_cross_family_evidence_v1.json)
is self-contained for the incidence check. It retains the required family
models and inventory, the two new sources, raw certificate, replay, protocol
and supervised logs. Its isolated replay supplies no new point search or
rank claim.
