# NS0024 new rootless frame: equation/source route

Date: 2026-09-03.

## Outcome

The new determinant-950 completion is an exact catalog-external rootless frame
candidate for `NS0024`: it has rank 17, minimum 4, no norm-two
vectors, 2,634 norm-four vectors, and automorphism-group order 2.  Its
norm-four count differs from the three rootless frames in the foundry catalog:

| frame | norm-four vectors |
|---|---:|
| `NS0024-F001` | 2,640 |
| `NS0024-F002` | 2,630 |
| `NS0024-F005` | 2,632 |
| new completion | 2,634 |

The exact Gram matrix and replay are in
[`../artifacts/generated-results/elkies-k3-ns0024-new-rootless-source-route-v1.json`](../artifacts/generated-results/elkies-k3-ns0024-new-rootless-source-route-v1.json).
This identifies a new integral frame class by an exact isometry invariant; it
does not yet produce a Weierstrass equation for the corresponding fibration.

There is, however, a clean equation **source**.  The canonical foundry source
`NS0024-S001` is

```text
2E8 / MW1,  height Gram [950],  torsion 1.
```

It is the frame of an Inose fibration attached to a non-isomorphic pair of
elliptic curves joined by a cyclic isogeny of degree 475.  This gives a
symbolic Weierstrass source over the level-475 isogeny locus.  The missing
part is an equation-level fibration change from that Inose model to the new
rootless frame.

## Exact completed-core route

Apply the same order-191 bridge class 4 to the canonical masked-core seed and
to the three stored Kneser neighbours.  The integral glue multipliers change
with the core basis, but the completed frames have the following exact
sequence:

```text
D5+E8 / MW4
   -- core 17-neighbour --> 3A1+A2 / MW12
   -- core 13-neighbour --> 3A1+A2 / MW12
   -- core  7-neighbour --> rootless / MW17 (new frame).
```

The selected multipliers are respectively `59, 50, 76, 83`; changing sign in
the cyclic group gives `132, 141, 115, 108`.  Each completed frame has
determinant 950.  This is an exact path in the genus of positive-definite
frames.  It is **not** an elliptic-neighbour corridor between marked primitive
copies of `U`, and the primes `17,13,7` must not be read as pencil degrees.

Its practical value is source discovery: `D5+E8/MW4` is the nearest
equation-facing source currently attached to the new frame.  A useful next
compiler target is therefore an explicit `II*+I1*` model for that fibration,
followed by a marked elliptic-neighbour search toward the two
`3A1+A2/MW12` frames and the rootless endpoint.

## Inose equation source

Write

```text
E1: y1^2 = x1^3 + a2*x1^2 + a4*x1 + a6,
E2: y2^2 = x2^3 + a2'*x2^2 + a4'*x2 + a6',
```

and let `Delta1, Delta2` be their discriminants.  Put

```text
A = (a2^2 - 3*a4) * ((a2')^2 - 3*a4'),
B = (32/27)
    * (2*a2^3 - 9*a2*a4 + 27*a6)
    * (2*(a2')^3 - 9*a2'*a4' + 27*a6').
```

Then the Inose surface is

```text
Y^2 = X^3 - (A/3)*X + (Delta1*s + B + Delta2/s)/64.       (1)
```

It has two `II*` fibres at `s=0,infinity`.  For non-isomorphic `E1,E2`, its
geometric Mordell--Weil lattice is `Hom(E1,E2)<2>`.  Thus a primitive
degree-475 isogeny supplies a section of height `2*475=950`, exactly the
`NS0024-S001` source lattice.  Equation (1), the Mordell--Weil identification,
and a general procedure for constructing the section attached to an isogeny
are given by Kazuki Utsumi, [*The Mordell-Weil lattice of an Inose surface
arising from isogenous elliptic curves*](https://arxiv.org/abs/2209.02463),
equation (2.7), Proposition 3.1, and Section 5.

This is a genuine equation route at the source-family level, not a numerical
ansatz.  To obtain one explicit specialization, one still has to provide an
explicit non-CM point of the desired field on `X0(475)`, construct the two
475-isogenous curves, and substitute their coefficients in (1).  The present
calculation does not claim such a point over `QQ`.

## Why the canonical source is not a cheap compiler start

Exact Sage modular-symbol calculations give

```text
X0(475): index 600, genus 45, 12 cusps, no elliptic points;
trace(W19,W25,W475) on differentials = (-3,1,-7);
genus X0(475)/<W475> = 19;
genus X0(475)/<W19,W25> = 9.
```

Consequently no rational global one-parameter chart for the level-475 source
should be assumed.  This does not preclude isolated rational or number-field
specializations.

There is also a large section-complexity barrier.  The exact affine-CVP audit
for `NS0024-S001` gives minimum nonzero section frame norm 950 and pole order

```text
(950 - 4)/2 = 473.
```

The root-adapted low-`q` search was extended with an explicit
`--include-zero-mw` option, because the default search deliberately omits
root-supported divisors.  Exhaustive Weyl-orbit searches then found:

| old-fibre degree | tested `q` | dominant orbits | primitive neighbours | rank-growing |
|---:|---|---:|---:|---:|
| 2 | every even `q=4,...,40` | 2,574 | 2,474 | 0 |
| 3 | every multiple `q=3,...,30` | 708 | 699 | 0 |

Below norm 950 the Mordell--Weil projection is necessarily zero.  Every
primitive neighbour in these bounded searches retains root rank 16.  This is
an exact bounded negative result, not a proof that a higher-`q` route does not
exist.  It does show that the canonical `2E8/MW1` fibration is a poor start for
the usual low-pole neighbour compiler.

## Recommended continuation

The most credible route is now:

1. Use (1) to anchor the surface and source identity on the degree-475
   isogeny locus.
2. Locate the `D5+E8/MW4` primitive `U` embedding on that Inose surface and
   compile its `II*+I1*` equation.  This is the first missing equation step.
3. Search actual degree-two/three elliptic neighbours from that marked model,
   using the two `3A1+A2/MW12` completions and the 2,634-vector Gram matrix as
   exact landing tests.
4. Once the rootless equation is reached, recover a 17-section basis and
   verify its height Gram against the stored new frame.

The alternative is to construct the degree-475 Inose section directly by
Utsumi's general-degree method and start a high-pole fibration change from
`2E8/MW1`; the pole-473 certificate makes that route substantially less
attractive.

## Replay

```bash
sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-ns0024-2e8-source-root-adapted.txt \
  --root-rank 16 --degree 2 --include-zero-mw \
  $(for q in $(seq 4 2 40); do printf -- '--q %s ' "$q"; done) \
  --output artifacts/generated-results/elkies-k3-ns0024-2e8-zero-mw-degree2-q40-v1.json

sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-ns0024-2e8-source-root-adapted.txt \
  --root-rank 16 --degree 3 --include-zero-mw \
  $(for q in $(seq 3 3 30); do printf -- '--q %s ' "$q"; done) \
  --output artifacts/generated-results/elkies-k3-ns0024-2e8-zero-mw-degree3-q30-v1.json

sage -python elkies-k3/scripts/certify_ns0024_new_rootless_source_route.sage --check
```

The generated certificate records SHA-256 hashes of every lattice, source,
pole, and bounded-neighbour input.  Its `proof_boundary` field is normative.
