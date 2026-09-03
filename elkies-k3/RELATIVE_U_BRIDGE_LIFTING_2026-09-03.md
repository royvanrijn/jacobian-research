# Cross-Gram reconstruction of relative `U`-embeddings and the first NS0024 search

<!-- status-consumer: EC-K3-RELATIVE-U-BRIDGE-LIFTING 800e22abf69b91aa -->
<!-- status-consumer: EC-K3-LOCAL-BRIDGE-MUTATION-H1C 7421328afadcf61f -->
<!-- status-consumer: EC-K3-PRIME-LOCAL-BRIDGE-MUTATION-H1D d4c6c84967a8fbc5 -->
<!-- status-consumer: EC-K3-NS0024-RELATIVE-U-FIRST-EDGE-OBSTRUCTION d57544697149506f -->
<!-- status-consumer: EC-K3-R17-NONCYCLIC-4A1-DIRECT-EQUATION f657620e07f8f3f0 -->

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

The two-sided strengthening is Theorem H-1c.  The saturation index is the
same in both directions, the saturated bridges have isomorphic discriminant
forms, and their graph-glue orders agree.  If their common determinant is
`c`, the common glue order is `h`, and `D=|det NS|`, then

```text
c/h divides gcd(c,D).
```

Thus 35 of the 42 stored maximal glues are forced by coprimality.  Seven
shared-prime cases require local graph data beyond the support theorem.  They
are now completely classified by Theorem H-1d.  In all seven, the raw bridges
are already saturated, the common-core determinant forces maximal local glue,
and every marked maximal graph has the required local NS discriminant form.
The finite quadratic form therefore fixes the order but not the graph label.
Exact norm-two coset enumeration distinguishes the alternatives: conditional
on the stored good-prime glue and the declared old/new ADE systems, the
historical transported graph is unique on every edge.

The parity conclusion also needs an odd-primary qualification.  Even relative
determinant kills the 2-primary bridge discriminant, but does not force
cyclicity at odd primes.  A saturated degree-two lattice counterexample has
`A=[[2,3],[6,8]]` and bridge Smith form `diag(3,21)`.  Conversely, odd
relative determinant and odd saturation index force two 2-primary generators.
The exact R17 control realizes the latter branch geometrically: it gives a nef
`4A1/MW13` fibration with maximal non-cyclic `ZZ/4+ZZ/8` bridge glue.  Its
frame is not isometric to either stored H3 `4A1` frame, so it is a new local
frame rather than a shortcut between the historical nodes.

The equation-facing follow-up is now complete: the same marking gives an
explicit `4I2+16I1` model over `QQ`, a saturated arithmetic MW13 basis, and a
target-free reverse hop to the literal published R17 equation.  See
[`R17_NONCYCLIC_4A1_DIRECT_FIBRATION_2026-09-04.md`](R17_NONCYCLIC_4A1_DIRECT_FIBRATION_2026-09-04.md).

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

Theorem H-1c is likewise recorded as a tailored consequence of Nikulin graph
glue and unimodular splitting.  Its R17 instance, the 35/42 forced-maximal
count, and the comparison with the two stored H3 `4A1` frames are new exact
computations.  Theorem H-1d and its seven-edge/R17 enumeration are the next
prime-local layer.

## Prime-local normal form and seven-edge census

Theorem H-1d separates three finite operations at each prime:

1. index-`ell^v_ell(m)` isotropic subgroups of the raw bridge discriminant
   form classify its possible saturation;
2. isotropic graph subgroups of `A_(K,ell)+A_(C,ell)` are filtered by
   `H^perp/H isomorphic to -q_(NS,ell)`;
3. the rank-four bridge isometry transports the old graph to the new graph,
   so the two sides are coupled rather than chosen independently.

The first two operations classify unmarked finite modules from the proposed
input list.  Root births require one extra marked layer: the saturation
embeddings/transport `tau_ell` and the coset theta minima of the metric
blocks.  Local Smith invariants by themselves cannot determine ADE births;
the historical ambient markings supply the missing decoration in the census.

For fixed `K,C` the graph order is forced by

```text
2 v_ell(|H|)=v_ell(det K)+v_ell(det C)-v_ell(|det NS|).
```

Away from `det NS`, the quotient form is zero and the graph is the full graph
of an anti-isometry.  This is one unmarked local orbit.  At the shared bad
primes the marked graph set can be larger, and the root effect is not encoded
by `q_NS`: it is decided by the theta/coset minima after the primary graphs
are assembled.

The seven exact marked graph counts are

```text
NS0024-8: 10     Golay720-3: 18     Golay720-4: 18
H3-6:       6     Q80-1:       4     Q80-5:       6     Q80-7: 4.
```

Their distinct transported ADE-pair counts are respectively
`10,12,16,6,3,5,4`.  In every case the actual historical ADE pair occurs
once.  Conversely, all marked graphs pass the local `q_NS` filter, so root
decoration is essential.

The R17 control has `A_(K,2)=(4,4,8)`, `A_(C,2)=(4,8)` and 32 marked maximal
graphs.  They yield eight ADE-pair types; the actual rootless-to-`4A1` pair
has multiplicity four.  It therefore disproves both cyclic one-generator
glue and uniqueness from the ADE pair.

There is one correction to the proposed global cutoff.  All seven historical
shared-prime cross matrices have `det(A)=0`, while their raw bridge Grams are
positive definite.  Hence `ell | 2D det(A)` is not a finite condition.  The
safe exact normal form is

```text
global bridge/glue completion
  = product_(ell | 2 D det(A^t J A-J)) X_ell,
```

with singleton unmarked factors away from `D`; after the saturated bridges
are fixed, possible graph-glue defect is supported on
`ell | gcd(det C,D)`.

Generate and byte-check the census with

```bash
sage -python elkies-k3/scripts/certify_prime_local_bridge_mutation.sage
sage -python elkies-k3/scripts/certify_prime_local_bridge_mutation.sage --check
```

The output is
[`elkies-k3-prime-local-bridge-mutation-v1.json`](../artifacts/generated-results/elkies-k3-prime-local-bridge-mutation-v1.json),
with SHA-256
`c158be14bbb32c8426ff90486194aff0e3fa8feafa66662f9e385d7c48db6909`.
It records the raw Grams, saturation and valuation ledgers, both local bridge
forms, local NS forms, every marked graph count, every exact ADE transition,
and the R17 noncyclic control.

## R17 local mutation replay

Generate and byte-check the exact certificate with

```bash
sage -python elkies-k3/scripts/certify_r17_local_bridge_mutation.sage
sage -python elkies-k3/scripts/certify_r17_local_bridge_mutation.sage --check
```

The output is
[`elkies-k3-r17-local-bridge-mutation-v1.json`](../artifacts/generated-results/elkies-k3-r17-local-bridge-mutation-v1.json).
Its pinned SHA-256 is
`71a49d17e95ac5861822530b3476276ab3fee4844f2727cb5af51f45f6bcad96`;
the record names SageMath 10.9 and hashes every input.
It checks the symbolic raw-Gram and parity identities, the 42-edge support
count, the primitive R17 bridge, both saturations and graph-glue orders, the
nef/physical-zero gate, root span and MW rank, and exact frame invariants.  It
does not construct a Weierstrass equation or decide a `J1` orbit.

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
square-index bridge correction, its two-sided saturation/glue strengthening,
the corrected 2-primary parity law, the prime-local saturation/graph normal
form, good-prime unmarked rigidity, bounded finiteness, the 42-hop regression,
the exact seven-edge bad-prime census, the exact R17 maximal non-cyclic local
mutation and negative-control census, and the stated exact first-edge
nonexistence boxes.

Not proved here: a global degree bound, a bound on `t`, existence outside the
audited boxes, nefness or zero effectivity for an un-audited marking, the
horizontal-wall certificate, an equation or `J1` identification for the new
R17 fibration, a rational map for the new NS0024 completion, variation of the
stored good-prime graph decorations in the seven-edge census, or automatic
lifting of every finite discriminant-form automorphism to an integral
automorphism of the metric blocks.
