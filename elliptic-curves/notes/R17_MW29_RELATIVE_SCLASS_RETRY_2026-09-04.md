# MW29-relative S-class retry for record fibres 356 and 385

Date: 2026-09-04  
Status: bounded exact relation experiment; no S-class, Selmer, or rank bound

## Target

The calculation was performed only in the quotient by the 29 certified point
half-ideals, with every prime ideal above the declared bad rational primes
inverted.  No full BNF was requested.  The bounded factor-base presentations
have dimensions 38 for curve 356 and 35 for curve 385 after inserting the
canonical rational principal rows and S-columns.  These dimensions are not
global upper bounds because factor-base generation has not been certified.

## Exact additions

`run_r17_kummer_quotient_sclass_collector.sage` now checks for equality of
reduced-ideal HNFs before attempting norm factorization.  Equality, together
with the two retained reduction multipliers, is an exact principal relation;
only such proof-bearing collisions are stored.

`run_fermigier_rank20_minkowski_specialq.py` now permits a deterministic cap
on the expensive batch-GCD stage.  It selects the lowest-bit-length unresolved
cofactors first, records the selected and omitted counts, and retains the full
unresolved generator/cofactor cache.  The ledger also records all collection
settings needed by subsequent runs.

## Bounded results

The dense signed half-ideal experiment used factor-base bound 240, target
widths one through three, 12--20 signed point companions, up to 20 S-ideal
companions, and a strict 90-second budget per curve.

| curve | attempts | cached reduced ideals | exact collisions | quotient-rank gain |
| ---: | ---: | ---: | ---: | ---: |
| 356 | 14,341 | 14,341 | 0 | 0 |
| 385 | 13,213 | 13,213 | 0 | 0 |

The full three-dimensional short-vector pass then enumerated every selected
degree-one special ideal through rational prime 239, combination bound four,
and the fixed determinant-one shape twists.  Trial division stopped at
20,000; unresolved cofactors were retained rather than generically factored.

| curve | exact candidates | batch-GCD tranche | GCD split operations | affected records | quotient-rank gain |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 356 | 88,648 | 30,000 | 4,142 | 3,847 | 0 |
| 385 | 94,134 | 10,000 | 698 | 662 | 0 |

The exact local artifacts and SHA-256 hashes are:

- `relative-collision-v1/curve356-b240-dense90.json`:
  `257559288c972d731a99fa47db0fe35bd52d1cd8e52d37106180bc802ca195c6`;
- `relative-collision-v1/curve385-b240-dense90.json`:
  `2eaa6833fd691580057983d8e895dedea2b6a79f9d66b99fd73fc72b49616a4f`;
- `relative-minkowski-v1/curve356-b240-allq-c4-batchgcd30k.json`:
  `d528c0f757fd7ea57ae763f7b8ccb80794b7d8b845bec417cd07d002ede7e50c`;
- `relative-minkowski-v1/curve385-b240-allq-c4-batchgcd10k.json`:
  `6093ed5a18f4fd2c31b67a4e980391d1d61390911b3dfeaec6162868f6046136`.

All paths are below
`artifacts/local/elliptic-curves/r17-kummer-quotient-sclass/`.

## Interpretation and next gate

Dense coset randomization, single reduced-ideal representatives, and generic
smoothness should not receive a larger budget: three independent exact
channels produced zero quotient-rank gain.  Batch GCD does expose many exact
shared rational factors, but insisting that each remaining cofactor be fully
factored discards all of them.

## Residual-ideal and targeted closure pass

That next implementation was completed.  For every affected record,
`refine_r17_unresolved_ideal_vertices.sage` reconstructs the principal
generator, removes all exact factor-base valuations and every proved shared
prime-ideal valuation, verifies the product, and retains the reduced
unfactored remainder as an exact ideal-HNF vertex.  It does not call a generic
factorization on that remainder.

| curve | proved-prime edges | principal residual tails | opaque residual tails | incidence dependencies | quotient-rank gain |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 356 | 3,607 | 540 | 3,067 | 0 | 0 |
| 385 | 636 | 89 | 547 | 0 | 0 |

Every opaque reduced tail was distinct in each ledger.  In particular, the
shared rational factors did not turn into repeated residual ideal classes at
this sample size.

The second-stage collector
`close_r17_residual_ideal_vertices.sage` then attacked the leaves rather than
resampling blind special ideals.  For a residual ideal `I`, it enumerates
projectively distinct short `beta in I`, factors the much smaller quotient
`J=(beta)/I`, and verifies the exact identity and every prime-ideal valuation.
For a certified point half-ideal the source class `I` is killed, so the same
operation is an actual MW29-relative relation.  Projective normalization is
essential: without it, rational multiples of the original generator create
tautological zero rows.

A declared tranche of 200 opaque ideals per curve used one new projective
representative per ideal and a 230-bit exact norm-factor limit.  It produced
29 direct known-half-ideal edges on each curve and 193/200 new tail edges.
The cache-only canonical postprocessor then inserted the exact principal
relations `(p)=product_{P|p}P^e` for all 3,660/1,454 outside rational primes
actually present in the respective graphs.

| curve | source edges | MW29 target edges | new tail edges | outside `(p)` rows | final edges / vertices | dependencies | quotient-rank gain |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 356 | 3,607 | 29 | 193 | 3,660 | 7,489 / 12,253 | 0 | 0 |
| 385 | 636 | 29 | 200 | 1,454 | 2,319 / 4,175 | 0 | 0 |

The additional exact artifacts and SHA-256 hashes are:

- `relative-ideal-vertices-v1/curve356-b240.json`:
  `0378eab099604d5a2b3d90595c13ba77896967cf7cfba365d1924e1ff056cffa`;
- `relative-ideal-vertices-v1/curve385-b240.json`:
  `d7052a6bfad12e8f597b9407a5163c40214416a3cea7ba3c114323f7ba164a10`;
- `relative-targeted-closure-v1/curve356-lll-tails200.json`:
  `36b28a632039abff0e0053494955f4677d387df53159d188d269bd37df4d1219`;
- `relative-targeted-closure-v1/curve385-lll-tails200.json`:
  `51503a57a3e256ef2eeb0e3022c02d42e85aef7f90d6bcce339e569fb7df9a7a`;
- `relative-targeted-closure-v1/curve356-lll-tails200-canonical.json`:
  `f05fb21f2cc4d6049469a174a45cdc4b8d72082aa2533ba54f85ad4ed26b7ccb`;
- `relative-targeted-closure-v1/curve385-lll-tails200-canonical.json`:
  `26929c0b463a6d7f1c2c28b25498be64288d403dd9aa73e30bc8f0d7c9e72c43`.

All paths are below the same local artifact root.  These exact graph sizes
also give a stopping rule: the new relations introduce roughly three fresh
outside prime-ideal vertices each, and even the free outside `(p)` rows leave
the incidence matrices full row rank.  Scaling this particular family keeps
the graph underdense; it is not a rational next expenditure.

The remaining decisive gate is therefore genuinely global: a certified
F2-only generation/upper bound for

```text
Cl(K) / (2*Cl(K) + <bad-prime ideals> + <known MW29 half-ideals>),
```

or an equivalent ray-class/complete relative-descent provider.  The complete
MW29-relative 2-Selmer quotients remain `UNKNOWN`; none of the bounded
dimensions `38` and `35` is an upper bound.

## Replay

For each curve, refine the retained cofactor tranche without a new search:

```bash
sage -python elliptic-curves/cas/refine_r17_unresolved_ideal_vertices.sage \
  --input artifacts/local/elliptic-curves/r17-kummer-quotient-sclass/relative-minkowski-v1/curve356-b240-allq-c4-batchgcd30k.json \
  --output artifacts/local/elliptic-curves/r17-kummer-quotient-sclass/relative-ideal-vertices-v1/curve356-b240.json \
  --max-inputs 30000 --overwrite

sage -python elliptic-curves/cas/refine_r17_unresolved_ideal_vertices.sage \
  --input artifacts/local/elliptic-curves/r17-kummer-quotient-sclass/relative-minkowski-v1/curve385-b240-allq-c4-batchgcd10k.json \
  --output artifacts/local/elliptic-curves/r17-kummer-quotient-sclass/relative-ideal-vertices-v1/curve385-b240.json \
  --max-inputs 10000 --overwrite
```

Run the declared targeted tranche by replacing `CURVE` with `356` and `385`:

```bash
sage -python elliptic-curves/cas/close_r17_residual_ideal_vertices.sage \
  --input artifacts/local/elliptic-curves/r17-kummer-quotient-sclass/relative-ideal-vertices-v1/curveCURVE-b240.json \
  --output artifacts/local/elliptic-curves/r17-kummer-quotient-sclass/relative-targeted-closure-v1/curveCURVE-lll-tails200.json \
  --max-tail-ideals 200 --max-samples-per-tail 1 \
  --max-tail-factor-bits 230 --max-samples-per-known-ideal 1 --overwrite
```

The current closer includes outside `(p)` rows directly.  To reproduce the
cache-only augmentation of the retained pre-augmentation target files, run:

```bash
sage -python elliptic-curves/cas/augment_r17_targeted_closure_canonical.sage \
  --input artifacts/local/elliptic-curves/r17-kummer-quotient-sclass/relative-targeted-closure-v1/curveCURVE-lll-tails200.json \
  --output artifacts/local/elliptic-curves/r17-kummer-quotient-sclass/relative-targeted-closure-v1/curveCURVE-lll-tails200-canonical.json \
  --overwrite
```
