# One planner-ready foundry route: exact H3 A1 to published R17

Status: **exact certified positive control; one route only**.

This note records one end-to-end route and does not bulk-fill the 936 foundry
route rows.  The selected target is the determinant-948 published rootless
frame `NS0001-F001`.  Its source is the exact characteristic-zero H3 final
reverse model with fibre type `A1` and Mordell--Weil rank 16 (assuming
`rho=19`).  The physical final edge has old-fibre degree two, so this is the
lowest-cost exact equation source currently available for the experiment.

The six hand-curated inputs are in
`data/lattice-foundry/planner-ready-h3-a1-r17-v1.json`:

1. a common marked `NS=U+W(-1)` convention, the source `U`, and the exact
   `QQ` Weierstrass/pointing artifacts;
2. an explicit rank-15 core plus determinant-47 binary bridge decomposition;
3. the order-47 graph-glue generators and multiplier;
4. the single good-prime plan `p=5`;
5. root-death/rootless-birth data plus the marked parent--target shell
   intersection counts at norms 4, 6, and 8;
6. the verification-only exact relative-`U` transport for H3 edge 13.

## Blind planner protocol and result

The planner input contains the parent core Gram matrix, binary bridge, graph
state, prime, empty root-survival set, desired rootless signature, and shell
overlap counts.  It does **not** contain the isotropic neighbour line, adjusted
lift, child Gram matrix, target frame Gram matrix or ID, target `U`, or the
elliptic-neighbour path.

With maximum parameter support 3 and 20,000 deterministic dense probes, the
planner proposed 7,954 isotropic lines.  It rejected 1,552 by the prescribed
root incidence and 6,401 by the shell fingerprint, then materialized one
neighbour.  That neighbour is rootless.  Only after this hit does the checker
load the endpoint data; exact PARI integral isometry identifies the completed
frame with `NS0001-F001` and rejects the alternate rootless frame
`NS0001-F002`.  The completed frame has determinant 948 and 2,622 signed
norm-four vectors.

The post-hit transport check then reads the independent relative-`U`
certificate.  H3 edge 13 has `q=6`, orbit 2247, old-fibre degree 2, saturation
index 1, and root-rank transfer `1 -> 0`.  Thus the certified chain is

```text
exact QQ A1/MW16 equation
  -> blind inverse-ADE core planner
  -> rootless NS0001-F001 and its U'
  -> exact degree-2 fibration hop (q6 orbit2247)
  -> MW17.
```

## Reproduction

```bash
sage -python elkies-k3/scripts/certify_single_planner_ready_foundry_route.sage
sage -python elkies-k3/scripts/certify_single_planner_ready_foundry_route.sage --check
```

The generated certificate is
`artifacts/generated-results/elkies-k3-single-planner-ready-foundry-route-v1.json`.
It serializes the actual planner input, records the endpoint fields that were
withheld, and records the discovered line only in the post-hit result.
The readiness audit now reports
`PASS_EXACTLY_ONE_CURATED_PLANNER_READY_ROUTE`: one curated ready route, zero
ready bulk routes, and all 936 original bulk rows preserved.

## Boundary

This is a positive control for the complete workflow, not a theorem that any
of the existing 936 bulk pairs has acquired an equation source.  `NS0001` had
no rootful companion in the original Niemeier-shell route ledger, so this exact
H3 source is represented as one separate curated route.  The bulk rows remain
unchanged and unready.  The computation certifies the lattice discovery and
uses the already-certified exact relative-`U`/elliptic-neighbour transport; it
does not replace those equation-level certificates.
