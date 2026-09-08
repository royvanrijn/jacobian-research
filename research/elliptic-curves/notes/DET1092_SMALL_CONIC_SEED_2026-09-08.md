# A smaller equation-derived determinant1092 seed

**Verified example.** The reduced-parent parameter

\[
\boxed{s=5193/35630}
\]

has an explicitly certified rank18 subgroup consisting of its17 inherited
sections and one point from the orbit8044 conic. Its original parent address is

\[
t=-1663732982925420458917240939000/
41904862529800516500295571942909.
\]

The reduced rational `j` has551 numerator bits and534 denominator bits. The
previous equation-derived `u=0` progression control had6,176 numerator bits.
This is a different, substantially smaller fibre; no model change of the old
control could change its `j` height.

The [certificate package](../../artifacts/generated-results/elliptic-curves/det1092_small_conic_seed_v1/manifest.json)
contains the equation, ordered18 points, exact parameter/model transports,
conic solution and independent finite-group proof. The fibre is unmatched by
exact rational isomorphism in the pinned630-entry public catalogue and201-entry
local inventory. This is a snapshot comparison, not literature-wide novelty.

## Construction

Use the existing
[rank18 conic base change](DET1092_RANK18_BASE_CHANGE_AND_INITIAL_UNLOCK_2026-09-08.md)
and the verified reduced parameter chart. Substitution of the chart matrix
into the original quadratic, followed by division by its exact square content,
gives

\[
z^2=4940136120823281-21931265870123820s
       +51437496474840100s^2.
\]

This primitive conic has coefficients of at most56 bits. One bounded PARI
`qfsolve` call on its ternary quadratic form supplies an isotropic integer
vector. Exact substitution verifies the rational point and the address above.
The frozen residual quadratic and line equations lift both rational branches
to the elliptic fibre. The first branch that passes independent rank
certification is retained; the two branches are not counted as two directions.

No exceptional curve302 point, public rank label, quartic point search or
catalogue enters construction. The input conic was already constructed from
the generic parent. This new specialization is a separate constructive lane
and does not count as a success of the
[ten-million-parameter selector](DET1092_SEARCH_FUNNEL_2026-09-08.md).

## Independent proof

The [standalone checker](../cas/verify_det1092_funnel_small_conic_seed.sage)
imports no constructor, repository rank backend or cubic Kummer implementation.
It checks the reduced/original parameter relation, all eighteen elliptic point
identities, equality of the first seventeen points with the specialized generic
sections, and the last point's residual-conic/line equations.

At the recorded good primes it enumerates each complete finite elliptic group
and its doubled subgroup using elementary integer addition. The combined
quotient matrix has column rank18; the inherited submatrix has rank17. An
odd-order good reduction excludes rational2-torsion. Every integral relation
is therefore divisible by2, and infinite descent proves independence.
The complete finite groups, cosets and point reductions are retained in the
[standalone replay](../../artifacts/generated-results/elliptic-curves/det1092_small_conic_seed_v1/standalone-replay.json).

Construction took8.67 supervised seconds, a repeated constructor/preflight
took1.52 seconds, and the independent proof took0.75 seconds. The respective
caps were90,90 and60 seconds; all inputs and supervision records are retained.

## Bounded amplification outcome

The sealed M18 entered unchanged V3. All114 selected charts completed, with
no certified extra direction. The independent map, point-cloud, landscape and
mod-2/3/5 replays pass. The terminal reason is `COMPLETE_FINITE_NO_GAIN`.
Thus the certified lower bound remains18; no exact-rank upper bound, conductor
result or rank above18 is asserted.

This supplies a practical seeded control at a height comparable to the lower
production bands. It does not prove that a small seeded fibre amplifies, or
that other fibres in the conic family lack additional points.

## Reproduction

The [constructor](../cas/construct_det1092_funnel_small_conic_seed.sage) freezes
a separate protocol before solving the conic. Fresh construction requires a
fresh local output directory; the broad population and its selected inputs
remain unchanged. The [amplifier wrapper](../cas/run_det1092_constructive_amplifier.py)
permits one finite attempt and retains the original V3 policy.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_funnel_small_conic_seed.sage \
  --run research/artifacts/local/elliptic-curves/det1092-small-conic-seed-v1 \
  --output research/artifacts/local/elliptic-curves/det1092-small-conic-seed-v1/standalone-replay.json
```

The result is an independently proved rational specialization of an existing
rank18 base-change construction. The search for higher-rank related curves
and improved conductors remains open.
