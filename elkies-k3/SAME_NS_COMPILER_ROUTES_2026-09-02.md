# Same-NS compiler routes for Golay-720 and NS0031 — 2026-09-02

<!-- status-consumer: EC-K3-GOLAY-DET720-PHYSICAL-CORRIDOR 0868df67fe8c37ad -->
<!-- status-consumer: EC-K3-NS0031-F017-PHYSICAL-CORRIDOR 67f5a71fa0733cdb -->
<!-- status-consumer: EC-K3-NS0031-MARKED-FORMAL-BRANCH b31e99bce4edac0a -->

## Outcome

The source-to-MW17 middle layer is now explicit at the marking/lattice level
for both promoted source classes.  Routes were selected by the compiler key

```text
(maximum horizontal P.O., maximum q, maximum old-fibre degree,
 edge count, physical Weyl repairs),
```

not by graph distance alone.

| NS/source | selected MW17 target | compiler cost | route length |
|---|---|---:|---:|
| Golay-720 `G720-S0128` (`3A5/MW2`) | `G720-F001` | `(4,4,2,6,2)` | 6 |
| NS0031 `(prescribed-root-sources-all-ns-3e8-all-a-v1.json, NS0031-S001)` (`A1+2A7/MW2`) | trisection-first `NS0031-F017` | `(6,8,2,5,1)` | 5 |

Every displayed edge has a primitive isotropic fibre, passes the complete
finite/affine component and all-section nef gates, has no negative
finite-degree horizontal wall, and carries a determinant-one integral
Neron--Severi transport.  The terminal frames are integrally isometric to the
named rootless MW17 frames.

These are exact selected-route certificates, not graph-optimality theorems.
The searches use bounded forward beams.  The NS0031 terminal meeting search
does completely enumerate the 63,902 minimum-norm-eight degree-two classes on
`F017`, finding six matches to the retained low-root forward frontier.

## Golay-720

The selected chain is

```text
3A5 -> 4A2+A5 -> 3A1+2A2+A3 -> 4A1+A2
    -> 3A1 -> 2A1 -> rootless G720-F001.
```

All six arrows have old-fibre degree two and `q=4`.  Their horizontal pole
orders are `(0,0,2,2,2,4)`.  Only the final arrow needs a chamber correction:
two component-Weyl reflections.  This gives compiler cost `(4,4,2,6,2)`.

The stronger `P.O.<=2` end-to-end corridor remains open.  A longer
zero-repair route with terminal cost `(6,8,2,13,0)` is lexicographically worse
and is therefore not the selected compiler route.

Replay the exact certificate with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_golay_det720_physical_corridor.sage --check
```

The compact certificate is
[`../artifacts/generated-results/elkies-k3-golay-det720-3a5-to-mw17-physical-corridor-v1.json`](../artifacts/generated-results/elkies-k3-golay-det720-3a5-to-mw17-physical-corridor-v1.json).

## NS0031

The selected chain is

```text
A1+2A7 -> A1+2A3+A5 -> 5A1+A3 -> 4A1 -> A1
       -> rootless NS0031-F017.
```

The edge `(q,P.O.)` pairs are

```text
(6,0), (8,2), (4,2), (8,6), (8,6).
```

All old-fibre degrees are two.  The first four presentations need no Weyl
repair.  The final `A1` bridge is obtained by matching the complete
minimum-norm-eight `F017` reverse shell to the low-root forward frontier; one
finite-simple-root reflection makes it physical.  The resulting compiler
cost is `(6,8,2,5,1)`.

`F017` is retained for target interest as well as accessibility.  Its exact
degree-two counts are 41,885 rational and 63,902 genus-one bisection
translation orbits.  In the deterministic 256-coset pilot it has 47 rational
and 80 genus-one trisection candidates, the largest rational-trisection count
among the thirteen catalogued rootless NS0031 frames.  The sampled trisection
ranking is a heuristic; the route certificate does not promote it to
arithmetic effectivity or rank gain.

Replay the exact certificate with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_ns0031_f017_physical_corridor.sage --check
```

The compact certificate is
[`../artifacts/generated-results/elkies-k3-ns0031-a1-2a7-to-f017-physical-corridor-v1.json`](../artifacts/generated-results/elkies-k3-ns0031-a1-2a7-to-f017-physical-corridor-v1.json).

## Proof boundary

Both corridors begin at exact prescribed-root lattice markings.  Neither is
yet a compiled algebraic neighbour sequence from a rational source equation.

- The simple rational Golay `3I6` specialization saturates to determinant 20,
  not 720, so it cannot be the start of the certified Golay route.
- NS0031 has exact marked finite-field evidence, a one-parameter formally
  smooth `ZZ_7` marked branch, and the finite lift as a regression, but no
  algebraized rational characteristic-zero source has yet been proved.

Thus the same-NS fibration-graph middle is now explicit and auditable.  The
remaining work is equation-level realization of the marked source and then
edge-by-edge pencil compilation.
