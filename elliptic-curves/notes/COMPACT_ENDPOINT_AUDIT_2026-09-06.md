# Twenty-one nonsingular endpoint curves omitted from the compact parameter boxes

The six compact R17 and five compact MW16 parameter scans explicitly exclude
zero and infinity. The [exact endpoint audit](../../artifacts/generated-results/elliptic-curves/compact_atlas_endpoints_v2.json)
now evaluates all 22 cases. Twenty-one are nonsingular and mutually
nonisomorphic over Q; none matches the pinned 593-equation catalogue or the
528 previously measured address-equations. The `a1-fibration-03` fibre at
infinity is singular in the endpoint chart and is not counted as an elliptic
curve.

Only existing sections have been evaluated. No endpoint rational-point search
has run, and no curve is added to the high-rank inventory. The
[summary certificate](../../artifacts/generated-results/elliptic-curves/compact_endpoint_summary_v1.json)
binds all transport and finite-independence checks.

At zero, the short coefficients are `A(0),B(0)`. At infinity put `t=1/v`,
`X=v^4 x`, `Y=v^6 y`, giving endpoint coefficients `A_8,B_12`. Rational
section coordinates are evaluated by exact Laurent valuations and leading
coefficients. Compatible poles map to the elliptic group identity; finite
values must satisfy the endpoint equation exactly.

The first audit incorrectly limited the accepted section format to constant
denominators. It failed before producing a rank certificate. Its frozen
protocol, source and failure log remain retained. Version 2 handles rational
section denominators, with explicit valuation checks, and completes all
22 cases. No failed output is overwritten.

The retained finite sections prove lower bounds from 11 to 17. Every
nonsingular endpoint's full list of transported sections is independently
checked modulo 3 and 5; both give the same bound as modulo 2. This does not
prove that the remaining sections are rationally dependent or establish an
exact specialization rank. In particular, specialization of a generically
independent set must not be treated automatically as an independent basis.

The following remain open: point searches on these 21 omitted curves,
additional independent directions, whole-curve exact ranks, and any
near-record result. A point campaign must use each endpoint's separately
certified subgroup and numerical geometry, rather than assume the full
generic 16/17-section Gram matrix stays independent after specialization.
