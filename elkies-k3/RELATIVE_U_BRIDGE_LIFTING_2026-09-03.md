# Cross-Gram reconstruction of relative `U`-embeddings and the first NS0024 search

<!-- status-consumer: EC-K3-RELATIVE-U-BRIDGE-LIFTING 800e22abf69b91aa -->
<!-- status-consumer: EC-K3-NS0024-RELATIVE-U-FIRST-EDGE-OBSTRUCTION d57544697149506f -->

Date: 2026-09-03

## Outcome

The proposed matrix identity is correct and gives a complete finite
parameterization in every declared box of the four marked intersections.  The
canonical statement and proof are Lemma H-1 and Corollary H-1a in
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md).
The necessary qualification is that `F.O'`, `O.F'`, and `O.O'` have their
physical compiler interpretation only after both derived `(-2)` classes have
been certified as the actual effective zero curves.  A splitting mate that has
not passed the chamber audit is only a lattice pseudo-zero.

The identity was replayed on all 42 recorded fibration hops.
Every stored old-fibre degree is recovered from `A_11`, and the saturation of
the two projected vectors recovers the stored bridge in both orientations.
All 84 oriented projected pairs happen to have saturation index one.  Thus the
regression verifies the square-index formula but does not supply an example
with a nontrivial correction.  The theorem note includes a separate exact
binary example with index four.  The exact regression record is
[`elkies-k3-relative-u-bridge-lifting-regression-v1.json`](../artifacts/generated-results/elkies-k3-relative-u-bridge-lifting-regression-v1.json).

## The finite theorem

For a fixed ordered source splitting `NS=U+W(-1)`, an ordered target `U'` with
cross-pairing matrix `A` projects to an ordered pair in `W` with Gram

```text
G_A = A^t*J*A-J.
```

Conversely, every ordered integral representation of `G_A` in `W` reconstructs
the literal target basis by `u'_j=(J*A)_(bullet j)+w_j`.  Its Gram is `J`, so
the resulting `U'` is unimodular, primitive, and an integral direct summand.
If `B=<w_1,w_2>` and `C=saturation_W(B)`, then

```text
det(G_A) = [C:B]^2 det(C).
```

Positive definiteness of `W` makes every fixed norm shell finite.  Therefore a
finite box of `(F.F',F.O',O.F',O.O')` gives a terminating exhaustive search of
literal ordered `U'` markings, not merely of abstract bridge classes.  A bound
on a downstream notion called “physical cost” has this conclusion only if it
provably bounds those four integers; that coercivity hypothesis cannot be
omitted.

The literature search found the ingredients separately rather than this
relative `2 by 2` packaging.  Primitive `U` embeddings are the standard
lattice encoding of Jacobian elliptic K3 fibrations; Brandhorst--Elkies,
Lemmas 2.5--2.6, connect their relative position to neighbouring frames; and
saturation/overlattice corrections belong to Nikulin's discriminant-form
formalism.  No source checked states the displayed cross-pairing identity as a
named theorem, but the proof is elementary integral lattice algebra and the
identity is treated here as a tailored lemma, not a novelty claim.  The new
results are the 84-presentation replay and the bounded NS0024 obstruction.

## NS0024 prospective result

The first requested edge starts from the completed `D5+E8/MW4` frame and seeks
the completed `3A1+A2/MW12` frame.  Before searching a second basis vector, it
is stronger and much cheaper to enumerate the possible primitive fibres
`F'`.  In the source basis `(F,O+F)`, write

```text
F' = (d+t)F + d(O+F) + w,
q = d(d+t),
```

where `d=F.F'` and `t=O.F'`.  The root-adapted Weyl enumerator exhausts these
classes modulo the source root Weyl group.  If no child has root rank five,
then no choice of `O'`, hence no representation of a full `G_A`, can produce
the desired MW12 target.

The exact results are:

| degree `d` | audited `t` | audited `q=d(d+t)` | dominant orbits | primitive fibres | largest child MW rank |
|---:|---:|---:|---:|---:|---:|
| 2 | 0 through 18 | even 4 through 40 | 431174 | 429877 | 5 |
| 3 | 0 through 4 | 9, 12, 15, 18, 21 | 13711 | 13704 | 7 |
| 4 | 0 through 4 | 16, 20, 24, 28, 32 | 79701 | 79375 | 8 |

The required child MW rank is 12.  Hence the first edge, and therefore the
three-edge chain, has no relative-`U` lift in these declared `(d,t)` ranges.
This is an exact bounded nonexistence result, not evidence against a higher
degree or larger-`t` lift.  Since the fibre obstruction already applies, no
prospective primitive `U'` marking exists to print in the audited range; a
positive-definite Kneser path has not been substituted for one.

The completed frames were also compared by exact integral isometry with every
stage of the previously known 13-edge NS0024 route.  None of the four completed
frames matches a route frame, including the stages with the same ADE/MW label.
Thus an old marking cannot simply be transplanted.  The comparison and compact
search summaries are:

- [`elkies-k3-ns0024-relative-u-first-edge-obstruction-v1.json`](../artifacts/generated-results/elkies-k3-ns0024-relative-u-first-edge-obstruction-v1.json), the compact combined certificate;
- [`elkies-k3-ns0024-completed-frame-comparison-v1.json`](../artifacts/generated-results/elkies-k3-ns0024-completed-frame-comparison-v1.json);
- [`elkies-k3-ns0024-relative-u-degree2-fibre-summary-v1.json`](../artifacts/generated-results/elkies-k3-ns0024-relative-u-degree2-fibre-summary-v1.json);
- [`elkies-k3-ns0024-relative-u-degree3-fibre-summary-v1.json`](../artifacts/generated-results/elkies-k3-ns0024-relative-u-degree3-fibre-summary-v1.json);
- [`elkies-k3-ns0024-relative-u-degree4-fibre-summary-v1.json`](../artifacts/generated-results/elkies-k3-ns0024-relative-u-degree4-fibre-summary-v1.json).

The `--summary-only` switch performs the same exact enumeration and root
classification while omitting hundreds of megabytes of individual negative
witnesses.  On a hit, summary-only mode must not be used: the full fibre
witness should be retained and passed to
[`search_ns0024_relative_u_bridge_lifts.sage`](scripts/search_ns0024_relative_u_bridge_lifts.sage)
for the second-vector representation, literal primitive-`U` construction,
target isometry, and marking export.

## Reproduction

First rebuild and compare the completed frames while exporting the root-adapted
source frame:

```bash
sage -python elkies-k3/scripts/search_ns0024_relative_u_bridge_lifts.sage \
  --compare-known-only \
  --export-adapted-source artifacts/generated-results/elkies-k3-ns0024-completed-d5e8-root-adapted-frame-v1.txt \
  --output artifacts/generated-results/elkies-k3-ns0024-completed-frame-comparison-v1.json
```

Then run
[`search_root_adapted_weyl_neighbors.sage`](scripts/search_root_adapted_weyl_neighbors.sage)
with `--root-rank 13 --adapt-mw-at-least 12 --rank-growth-only
--include-zero-mw --summary-only`.  Use:

```text
degree 2: q = 4,6,8,...,40
degree 3: q = 9,12,15,18,21
degree 4: q = 16,20,24,28,32
```

The input frame is
[`elkies-k3-ns0024-completed-d5e8-root-adapted-frame-v1.txt`](../artifacts/generated-results/elkies-k3-ns0024-completed-d5e8-root-adapted-frame-v1.txt).

## Proof boundary

Proved here: the relative-`U` parameterization, automatic primitivity, the
square-index bridge correction, bounded finiteness, the 42-hop regression,
and the stated exact first-edge nonexistence boxes.

Not proved here: a global degree bound, a bound on `t`, existence outside the
audited boxes, nefness or zero effectivity for an un-audited marking, the
horizontal-wall certificate, an equation, or a rational map for the new
NS0024 completion.
