# Lattice foundry report — 2026-09-01

## Outcome

The first deterministic foundry shell found 136 rootless rank-17 frame
classes on determinant-changing Picard-19 lattice classes after exact
saturation and ternary-realizability gates. A subsequent exact source hunt
and neighbour replay promotes `NS0024` (discriminant 950) to the preferred
lattice-level outcome:

```text
5A1+A2+A5 / MW5
    -- eleven primitive nef degree-two edges, q only 4 or 6 -->
new rootless frame / MW17.
```

The source is a saturated complement in `N(A11+D7+E6)`. Every route edge has
complete component, all-section, and finite horizontal-wall gates; all eleven
physical Weyl reflection counts are zero; and the composed 19-dimensional
integral marking is unimodular. The terminal positive frame is not integrally
isometric to any of the three rootless `NS0024` frames in the initial shell,
so the route itself discovers a fourth exact class.

This is success level A at the lattice/arithmetic/route layer. It is not yet
an equation: no effective equation-side zero, resolved Riemann--Roch pencil,
or characteristic-zero Weierstrass model is claimed.

## Exact arithmetic class

The new surface class has

```text
disc(NS) = -950,
T = [[0,0,1], [0,950,0], [1,0,0]],
signature(T) = (2,1),
disc(T) = -950.
```

The finite quadratic module is cyclic of order 950, with the stored normal
form `diag(3/2, 2/25, 2/19)` in `Q/2Z`. The exact opposite-discriminant-form
certificate glues `NS` and `T` to the even unimodular lattice of signature
`(3,19)`, hence to the K3 lattice. This is the same certified `NS0024/T`
class for source and target; it is not the determinant-948 H3 surface.

The equation-friendly source has root rank 12, MW rank 5, trivial torsion,
and determinant-predicted MW regulator `475/288`. Its root lattice is
primitive. It has 46 norm-two vectors, 2,712 norm-four vectors, and
automorphism group order 276,480. Exact cyclic discriminant gluing embeds the
rank-seven auxiliary primitively into `N(A11+D7+E6)` and recovers the source
as its saturated orthogonal complement.

## Certified route

| edge | q | orbit | source | target | old-fibre degree | P.O | Weyl reflections |
| ---: | ---: | ---: | --- | --- | ---: | ---: | ---: |
| 1 | 4 | 1 | `5A1+A2+A5/MW5` | `2A1+A2+D6/MW7` | 2 | 0 | 0 |
| 2 | 4 | 1 | `2A1+A2+D6/MW7` | `A3+A6/MW8` | 2 | 0 | 0 |
| 3 | 4 | 2 | `A3+A6/MW8` | `A1+A7/MW9` | 2 | 0 | 0 |
| 4 | 4 | 121 | `A1+A7/MW9` | `A3+D4/MW10` | 2 | 0 | 0 |
| 5 | 4 | 78 | `A3+D4/MW10` | `A6/MW11` | 2 | 0 | 0 |
| 6 | 6 | 2,933 | `A6/MW11` | `A5/MW12` | 2 | 1 | 0 |
| 7 | 6 | 13,293 | `A5/MW12` | `A4/MW13` | 2 | 1 | 0 |
| 8 | 6 | 33,561 | `A4/MW13` | `A3/MW14` | 2 | 1 | 0 |
| 9 | 6 | 131,287 | `A3/MW14` | `A2/MW15` | 2 | 1 | 0 |
| 10 | 6 | 274,686 | `A2/MW15` | `A1/MW16` | 2 | 1 | 0 |
| 11 | 6 | 493,546 | `A1/MW16` | `rootless/MW17` | 2 | 1 | 0 |

The route cost vector is `(maximum degree, maximum q, reflections,
root-rank area, edges, sum q) = (2,6,0,67,11,56)`. No zero-changing loop is
used. The planned resolved-RR dimension is two at each step, but that number
is deliberately stored as an equation-planning heuristic rather than a
lattice theorem. Global optimality among every possible neighbour route is
not claimed.

At each late node the complete cheaper `q=4` shell was exhausted before
using `q=6`. The selected first-growth searches are deterministic but do not
make a completeness claim for unvisited higher-cost shells after their first
hit.

## New rootless target score

The routed terminal has determinant 950, minimum squared norm 4, 2,644
norm-four vectors (1,322 unoriented pairs), automorphism group order 4, and
1,322 mod-2 short cosets hit through norm 4. Its initial theta series is
`1 + 2644 q^4 + ...`, and its Hermite invariant is approximately
`2.6723896387`. The 1,322 unoriented minimal representatives have absolute
pairing histogram

```text
|<v,w>| = 0: 344815 pairs
|<v,w>| = 1: 433452 pairs
|<v,w>| = 2:  94914 pairs.
```

There are 2,644 oriented degree-two, `q=2` isotropic candidates obtained from
these norm-four vectors before automorphism/chamber quotienting. These
statistics are discovery scores, not evidence for exceptional specialization
rank.

## Equation promotion

The source root type admits the semistable planning profile
`I6 + I3 + 5I2 + 5I1`, whose Euler numbers total 24. Root type alone does not
prove those Kodaira symbols, so additive alternatives remain in the ansatz
ledger. Since `rho=19`, the K3 moduli dimension is one. A rational or
low-genus parameter is plausible enough for modular reconnaissance, but no
such parametrization has yet been proved.

The next authorized computation is therefore a small-prime source-equation
ansatz, followed by zero/component marking on edge 1. Expensive
characteristic-zero reconstruction remains closed until those modular and
resolved-RR gates pass. Direct construction from the rootless equation is
not the preferred entry point.

## Positive controls and search boundary

The initial run contains 48 auxiliary classes, 509 frame classes, and 138
rootless frames including the two locked H3 controls. The complete H3 branch
is replayed independently: 3,220 sixth-vector candidates reduce to 167
seventh-vector labels, twelve primitive rootless embeddings, and exactly two
frame classes—published R17 and alternate Q80—with no third rootless H3 `J2`
frame. The known `N(2A7+2D5)` controls are checked before the determinant is
varied.

All 768 signed ambient-root mutations of the two stored control embeddings in
the declared mutable rows are tested, saturated, and filtered exactly. This
does not enumerate every rank-four-to-seven root subsystem, every short-glue
extension, or every primitive rank-seven Niemeier auxiliary below determinant
5,000. The auxiliary generator is a deterministic heuristic shell, not a
bounded-classification theorem. The exact H3 result remains complete only
for the pinned determinant-948 auxiliary.

## Reproduce

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_lattice_foundry.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/hunt_lattice_foundry_rootful_source.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_lattice_foundry_route.sage \
  --manifest elkies-k3/data/lattice-foundry/ns0024-nef-route-v1.json
```

Each command also supports `--check` after its output has been generated.
The eleven neighbour-search commands and their exact exhausted-shell bounds
are retained in the route manifest's input artifacts.

Configuration:
[`data/lattice-foundry/one-root-control-shell-v1.json`](data/lattice-foundry/one-root-control-shell-v1.json)

Route specification:
[`data/lattice-foundry/ns0024-nef-route-v1.json`](data/lattice-foundry/ns0024-nef-route-v1.json)

Database:
[`artifacts/generated-results/elkies-k3-lattice-foundry-v1.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-v1.json)

Source certificate:
[`artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-source-hunt.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-source-hunt.json)

Route certificate:
[`artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-nef-route.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-nef-route.json)
