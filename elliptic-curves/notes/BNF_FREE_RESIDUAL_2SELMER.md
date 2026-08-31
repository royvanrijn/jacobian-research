# BNF-free residual 2-Selmer bookkeeping

Status: implementation infrastructure.  This note proves no class-group,
2-Selmer, Cassels--Tate, or exact-rank result.

The large cubic fields arising from the fixed high-rank fibres need only the
modulo-squares part of an (S)-class computation.  A generic BNF calculation
is therefore the wrong intermediate target.  The relation collectors for
curve 273 use products of degree-one prime ideals, full three-dimensional
Minkowski lattices, several determinant-one archimedean twists, and ordinary
plus large-prime support.  See
[`search_curve273_ideal_lattice_relations.sage`](../cas/search_curve273_ideal_lattice_relations.sage)
and its exact replay
[`verify_curve273_full_ideal_descent_chain.sage`](../cas/verify_curve273_full_ideal_descent_chain.sage).

## Exact relation records

[`residual_selmer_quotient.py`](../cas/residual_selmer_quotient.py) supplies
the CAS-independent GF(2) layer.  Every relation has both its ordered
prime-ideal valuation-parity row and a reproducible principal generator.  Its
sparse elimination returns a dependency as the actual list of principal
generators whose product is a unit modulo squares.  The curve-273 pool audit
now prints this product for each independent LP-free relation; it no longer
reports an anonymous incidence-vector dependency.

The implementation does not factor a class group.  The field-specific
collector is responsible for exact ideal factorization and generator
coordinates; this layer preserves those results and cannot convert a bounded
relation plateau into a theorem.

For curve 273, persist this bridge from the collector with:

```bash
sage -python elliptic-curves/cas/analyze_curve273_relation_pool.py \
  --include-full-ideal-chain --include-crt-cycle-logs \
  --include-ideal-lattice-logs \
  --write-principal-relations \
  artifacts/local/elliptic-curves/curve273_principal_relations.json
```

The resulting local ledger contains ordered factor-base and large-prime ideal
columns, parity rows, sources, and the exact power-basis generator for every
row.  Its schema explicitly calls it an exact relation checkpoint, not a
class-group-completeness certificate.

The Fermigier quadratic special-q collector uses the same principle. Its
`--relation-ledger` output retains every sampled `a+b*theta+c*theta^2`
generator and, for each single-/double-large-prime closure, the mod-two list
of generators whose product supplies that relation. It is deliberately a
per-run exact ledger: its checkpoint remains a compact parity-row cache and
is not retroactively misrepresented as retaining generators from older runs.
The original one-dimensional special-q collector now writes the same ledger.
Its reproducible `5689:5096` calibration closes 335 relations in one bounded
30,000-candidate run; each retained double-large-prime cycle has its exact
one- or two-generator product. This is an exact relation replay, not a
class-group-completeness result.

[`run_fermigier_rank20_minkowski_specialq.py`](../cas/run_fermigier_rank20_minkowski_specialq.py)
now supplies the intended full three-dimensional alternative: it takes an
actual special-q ideal basis, applies determinant-one archimedean shape
twists, LLL-reduces each Minkowski lattice, and enumerates short combinations.
Every sampled algebraic integer and every merged large-prime cycle is written
to the same audited ledger schema. It is now parameterized by a monic cubic
and a declared Selmer-prime set; for signature `(1,1)` it uses one real
coordinate plus the real and imaginary coordinates of a complex embedding.
The `--curve273` preset loads the pinned ICARM cubic and all 14 declared
Selmer rational primes, so the same full ideal-lattice collector can begin a
curve-273 ledger directly:

```bash
sage -python elliptic-curves/cas/run_fermigier_rank20_minkowski_specialq.py \
  --curve273 --factor-base-bound 5000 \
  --special-q-min 5003 --special-q-max 15000 --max-special-q 20 \
  --lattice-combination-bound 2 \
  --shape-twists=-30:0,-20:0,-10:0,-5:0,0:0,5:0,10:0,20:0,30:0 \
  --trial-prime-bound 5000 \
  --relation-ledger artifacts/local/elliptic-curves/curve273_minkowski_relations.json
```

The `--elkies-rank28` preset loads the published `t=-9529/5471` cubic and all
twelve rational bad primes from the proved bad-place ledger. The same factor
hints are consumed by `NumberField` construction, canonical-row augmentation,
and the final audit, avoiding a repeated factorization of the 168-digit
polynomial discriminant. The pinned factor-base-1000 paired-special-ideal
pilot samples 10,288 algebraic integers and closes no noncanonical relation.
After adding 172 exact canonical `(p)` generators, its 327-column model has 26
`S` columns, relation rank 172, and displayed quotient dimension 141. This is
`UNCERTIFIED_FACTOR_BASE`: 1,000 is below its Bach/ERH generation bound
1,202,640, so the displayed dimension is not an `S`-class or Selmer upper
bound. The exact command chain is in [`../REPRODUCE.md`](../REPRODUCE.md).

This bounded full-S tranche has 1,225 factor-base ideals, 30 S-columns, and
21,916 exact sampled generators. It closes 60 cycles of relation rank one.
Augmenting its ledger with all 671 free canonical `(p)` rows raises the exact
relation rank to 671 and lowers the displayed S-quotient model from 1,195 to
538 dimensions. This is still **uncertified**, because its ERH generation
threshold is 1,231,857. The current non-S projection kernel has 74 explicit
products; its 59-dimensional square-norm kernel is entirely explicit global
squares. Thus that bounded norm-compatible span is trivial. It is therefore a
useful portability and relation-replay check, not new Selmer information.
A wider declared run (120 special ideals, coefficient bound four, eleven
shape twists) samples 955,510 exact generators and closes 840 cycles. It has
the same zero gain beyond the canonical quotient baseline; its 854
non-S-kernel products have an 839-dimensional norm-square kernel, all explicit
global squares. This is a bounded negative result for that exact family, not
a completeness claim. The collector now removes the associate pair `alpha`
and `-alpha` before merging—otherwise their product gives the tautological
global square `-alpha^2`—and rejects an element divisible by every prime above
its selected rational special prime, since that is merely a canonical `(p)`
multiple. It can sample degree-two special ideals or products `Q_i Q_(i+1)`
around a special-ideal graph cycle. The latter modes are the next
quotient-novel relation families; their smoke runs make no Selmer claim.

The bridge from closed relations to explicit squareclasses is:

```bash
sage -python elliptic-curves/cas/extract_bnf_free_squareclasses.py \
  --relation-ledger artifacts/local/elliptic-curves/r20_5689_relations.json \
  --output artifacts/local/elliptic-curves/r20_5689_ks2_candidates.json
```

It performs sparse GF(2) elimination on the complement of `S`, retaining a
basis of combinations whose remaining factor-base parity is supported at `S`.
Their stored generator products therefore have even valuations outside `S`.
Use `--individual-s-supported-only` only to reproduce the older, weaker
single-row diagnostic.
On the reproducible `5689:5096` calibration, the non-S projection kernel has
dimension 305 (from 335 closed rows). Its 296-dimensional square-norm kernel
has an explicitly verified global-square basis, so the full norm-compatible
span of this larger candidate space has no nontrivial class. The older
individual-row diagnostic had only 15 candidates and raw residual signature
rank seven modulo the known rank-20 image; that raw signature was not a
descent result. This remains a bounded relation-collection result, not a
complete `K(S,2)` or a rank conclusion.

Before evaluating local signatures, enforce the exact cubic norm condition:

```bash
sage -python elliptic-curves/cas/filter_bnf_free_norm_condition.py \
  --candidates artifacts/local/elliptic-curves/r20_ks2_candidates.json \
  --generate-norm-kernel \
  --relation-ledger artifacts/local/elliptic-curves/r20_relations.json \
  --output artifacts/local/elliptic-curves/r20_norm_square_candidates.json
```

For a monic short model, a Kummer representative `x-theta` has norm `y^2`.
The kernel mode uses only valuations at the declared rational Selmer primes
plus the norm sign to find all products with rational-square norm in the span
of the supplied candidates. It then discards exact global squares, recording a
verified square root for each discarded kernel basis class. On the bounded
complete-base ten-special calibration, the non-S projection kernel has 464
explicit products and its 455-dimensional square-norm kernel consists
entirely of global squares. The entire norm-compatible span is therefore
trivial. The previously reported one-dimensional signature remainder was a
fixed-precision real-sign cancellation artifact in two individual squares;
real signs are now evaluated as exact algebraic-real comparisons. This remains
only a bounded collection result, not a complete `K(S,2)` or a rank
conclusion.

The next exact hand-off is to materialize the corresponding 2-covers, before
any point search or Cassels--Tate pairing:

```bash
sage -python elliptic-curves/cas/build_bnf_free_two_covers.py \
  --candidates artifacts/local/elliptic-curves/r20_norm_square_candidates.json \
  --output artifacts/local/elliptic-curves/r20_two_covers.json
```

For each `alpha`, the output records the two homogeneous quadrics in
`[u:v:w:z]` obtained from `x-theta = alpha*(u+v*theta+w*theta^2)^2`, together
with `x` as the constant coefficient divided by `z^2`.  The output is an
explicit arithmetic target only: its local solubility must be checked at every
relevant place before it can enter a residual Selmer basis or a
Cassels--Tate calculation.

At selected manageable finite primes, the reduction audit gives a genuine
one-sided certificate: a smooth projective `F_p` point proves a `Q_p` point by
Hensel lifting, while an empty projective reduction proves no `Q_p` point.
For singular residue points it also tests a valuation-Hensel minor and, if
needed, exhaustively lifts normalized projective classes modulo increasing
powers of `p`; an empty lift tree is again a certified obstruction. The
two-adic trivial cover, for example, is certified by this singular-lift path.
Primes above the enumeration limit and exhausted/state-capped lift trees are
retained as inconclusive, never silently accepted:

```bash
sage -python elliptic-curves/cas/audit_bnf_free_two_cover_reduction.py \
  --covers artifacts/local/elliptic-curves/r20_two_covers.json \
  --primes 2,3,5,7,13 --output artifacts/local/elliptic-curves/r20_cover_local.json
```

It is a finite-place filter, not a full local-Selmer certificate: all relevant
finite places and the real place still need complete local analysis.

An exact finite obstruction can already be fed into the early quotient layer.
This strict mode keeps a candidate only when every finite place listed in the
cover audit is positively certified; it separates inconclusive covers rather
than admitting them:

```bash
python3 elliptic-curves/cas/filter_bnf_free_local_selmer.py \
  --images artifacts/local/elliptic-curves/r20_candidate_images.json \
  --cover-local-audit artifacts/local/elliptic-curves/r20_cover_local.json \
  --output artifacts/local/elliptic-curves/r20_finite_locally_filtered.json
```

## Quotient early

The JSON interface accepts exact candidate global-squareclass images in two
separate target spaces:

1. local squareclasses at the selected (S)-places; and
2. a fixed auxiliary-prime fingerprint.

It concatenates these spaces and reduces every candidate modulo the supplied
Kummer images of the known Mordell--Weil subgroup.  A surviving coordinate is
therefore an explicit target for a cover or Cassels--Tate calculation; a zero
coordinate says only that the selected targets do not distinguish it from the
known image.  It is not a global-square or Selmer conclusion.

For a manifest with field-specific values, run:

```bash
python3 elliptic-curves/cas/residual_selmer_quotient.py \
  --input artifacts/local/elliptic-curves/curve273_squareclass_images.json \
  --output artifacts/local/elliptic-curves/curve273_squareclass_quotient.json
```

The input holds `local_dimension`, `fingerprint_dimension`,
`known_mw_images`, and `candidate_images`.  Each image retains its exact
`generator`, `local`, and `fingerprint` bit masks.  It also contains the
separately scoped `class_quotient_certification` record.

The Fermigier rank-20 calibration can now generate its known-image side of
this manifest directly, rather than retaining only rank counts:

```bash
sage -python elliptic-curves/cas/run_fermigier_rank20_auxiliary_fingerprints.py \
  --prime-bound 59 \
  --output artifacts/local/elliptic-curves/fermigier_rank20_signature_map.json

python3 elliptic-curves/cas/residual_selmer_quotient.py \
  --input artifacts/local/elliptic-curves/fermigier_rank20_signature_map.json \
  --output artifacts/local/elliptic-curves/fermigier_rank20_signature_audit.json
```

This produces 51 Selmer-relevant local coordinates and 24 auxiliary witness
coordinates.  Its known-MW target has rank 20 using exactly
`11,19,23,29,59`; the witness coordinates remain non-Selmer coordinates.

Curve 273 has the same concrete handoff:

```bash
sage -python elliptic-curves/cas/analyze_curve273_kummer_fingerprint.py \
  --prime-bound 5000 \
  --output artifacts/local/elliptic-curves/curve273_signature_map.json

python3 elliptic-curves/cas/residual_selmer_quotient.py \
  --input artifacts/local/elliptic-curves/curve273_signature_map.json \
  --output artifacts/local/elliptic-curves/curve273_signature_audit.json
```

Its reproducible known-image map has 59 local and 54 witness coordinates and
faithfully separates the certified rank-30 Kummer image. It prepares the
quotient target for a candidate (K(S,2)) representative but creates no
residual Selmer class by itself.

To evaluate actual generators from the relation ledger, provide a JSON list
whose entries have `label` and three ascending
`generator_coefficients` (for `1,theta,theta^2`), then run:

```bash
sage -python elliptic-curves/cas/evaluate_bnf_free_signature_map.py \
  --signature-map artifacts/local/elliptic-curves/curve273_signature_map.json \
  --candidates artifacts/local/elliptic-curves/curve273_candidate_generators.json \
  --output artifacts/local/elliptic-curves/curve273_candidate_images.json

python3 elliptic-curves/cas/residual_selmer_quotient.py \
  --input artifacts/local/elliptic-curves/curve273_candidate_images.json \
  --output artifacts/local/elliptic-curves/curve273_candidate_quotient.json
```

The evaluator preserves every newly encountered 2-adic direction by extending
the local target (with zero coordinates for the old known image). Its result
is still only a signature image: separate local-solubility and global
class-quotient arguments are required before calling it a Selmer class.

Local images must be filtered against an independently computed basis of the
full local Kummer image, not merely the images of known rational points.
[`filter_bnf_free_local_selmer.py`](../cas/filter_bnf_free_local_selmer.py)
accepts a versioned `elliptic-curves.bnf-free-local-kummer-map.v1` JSON map
with `local_dimension`, `method`, and packed `allowed_local_images`. It
checks that every known MW image lies in this supplied local space, rejects
locally obstructed candidates, and only then recomputes the fingerprint/MW
quotient:

```bash
python3 elliptic-curves/cas/filter_bnf_free_local_selmer.py \
  --images artifacts/local/elliptic-curves/r20_candidate_images.json \
  --local-kummer-map artifacts/local/elliptic-curves/r20_local_kummer_map.json \
  --output artifacts/local/elliptic-curves/r20_locally_filtered.json
```

This remains a local filter: its survivors still need the global S-class
bound and the descent norm condition before they form a residual 2-Selmer
basis.

`audit_bnf_free_local_kummer_coverage.py` provides a separate certification
check for the odd places already represented in a signature map. At an odd
prime it bounds `dim E(Q_p)/2E(Q_p)` by the 2-part of the Tamagawa number times
the nonsingular special-fibre group. If known Kummer images attain that bound,
they are certified to span the full local image at that place; otherwise the
place stays unresolved. It treats the real component similarly and explicitly
does **not** claim two-adic coverage:

```bash
sage -python elliptic-curves/cas/audit_bnf_free_local_kummer_coverage.py \
  --signature-map artifacts/local/elliptic-curves/fermigier_rank20_signature_map.json \
  --output artifacts/local/elliptic-curves/r20_local_coverage.json
```

For Fermigier this certifies full known-point coverage at five of eleven odd
bad primes and at the real place; primes `3,5,7,31,79,1049` and the two-adic
place remain unresolved. Thus it narrows the local-descent work without
silently replacing it by the known Mordell--Weil image.

The published Elkies rank-28 fibre is now a stronger positive control for
this distinction. The exact builder
[`build_elkies_2026_rank28_local_coverage.py`](../cas/build_elkies_2026_rank28_local_coverage.py)
recomputes 53 bad-place coordinates for both the generic seventeen and the
certified eleven-point public complement. The generic local-signature rank is
15, and the rank after adjoining all eleven globally independent exceptional
directions is still 15. Every individual bad-place block has incremental rank
zero. Thus a signature map can completely miss genuine Mordell--Weil quotient
gain; it cannot be used as either a lower or an upper bound for that quotient.

The corresponding coverage audit proves equality with the full local Kummer
image at odd primes `3`, `19`, `20650099`, and
`315574902691581877528345013999136728634663121`, and at the real place. It
leaves seven odd primes and the two-adic place unresolved. These exact local
equalities reduce later membership tests, but do not address the global
`S`-class quotient.

For explicit norm-one candidates, use
[`run_bnf_free_two_cover_local_supervisor.py`](../cas/run_bnf_free_two_cover_local_supervisor.py)
instead of a single monolithic cover audit. It gives each `(cover,p)` pair its
own process group, wall/RSS limits, and metadata-bound cache block, so a hard
singular Hensel tree cannot erase completed local certificates. Its pinned
rank-28 pilot tests 12 of 49 bounded norm-one covers at seven odd primes. It
certifies 60 local points; 19 singular lift trees hit the state cap and five
workers time out. No local obstruction is found. All 24 incomplete
mathematical cases remain inconclusive, and untested finite places, the real
place, and global completeness remain open.

There is now a genuine positive-control suite alongside that synthetic pilot.
[`build_elkies_2026_rank28_public_selmer_controls.py`](../cas/build_elkies_2026_rank28_public_selmer_controls.py)
takes the eleven public complement points already certified independent modulo
the generic seventeen. For each point `Q`, it records the exact Kummer class
`alpha=X(Q)-theta`, verifies `Norm(alpha)=Z(Q)^2`, and supplies the associated
cover point `[1:0:0:1]`. The generic cover builder and local audit independently
recheck all eleven rational witnesses. Consequently these are genuine Selmer
classes and prove residual 2-Selmer dimension at least 11 on the published
rank-28 fibre. This is the expected positive-control floor, not the missing
ambient upper bound: rank 32 still requires residual dimension at least 15.

Those certified place factors can already eliminate a candidate before the
unresolved local work. The coverage mode uses only these exact local spans and
does not treat a survivor as locally soluble everywhere:

```bash
python3 elliptic-curves/cas/filter_bnf_free_local_selmer.py \
  --images artifacts/local/elliptic-curves/r20_candidate_images.json \
  --local-coverage-audit artifacts/local/elliptic-curves/r20_local_coverage.json \
  --output artifacts/local/elliptic-curves/r20_covered_places_filtered.json
```

## Cassels--Tate handoff

Once those arguments provide a certified residual Selmer basis, record the
explicit pairing matrix and audit it before searching any cover:

```bash
python3 elliptic-curves/cas/audit_residual_cassels_tate.py \
  --input artifacts/local/elliptic-curves/curve273_residual_pairing.json \
  --output artifacts/local/elliptic-curves/curve273_residual_pairing_audit.json
```

The input names `known_mw_rank`, the `residual_basis`, an alternating
`cassels_tate_matrix`, the pairing algorithm, and whether the basis has been
certified by the preceding global/local computation. The audit derives the
radical and hence the post-pairing rank upper bound. It calls a nondegenerate
residual pairing an exact-known-rank result only with that basis certification;
otherwise it remains an unchecked pairing audit. A nonempty `cover_searches`
list is rejected unless it is explicitly recorded as post-pairing.

## Completeness gate

The audit refuses to label relation-rank stabilization as a certificate.  A
non-null remaining mod-2 (S)-class bound must name a valid source, such as a
relation-lattice index bound, an enumeration of quadratic extensions
unramified outside (S), an analytic class-number bound, or an explicitly
labelled GRH 2-class computation.  The output records a hypothesis when one
is used.

The curve-273 pool also implements the conservative relation-lattice route.
With `--complete-factor-base`, it materializes every prime ideal over rational
primes through the requested bound, computes the cubic Minkowski bound, and
computes the quotient dimension of the collected pure factor-base relations
after killing (S)-columns. It labels that number
`CERTIFIED_MINKOWSKI_FACTOR_BASE` only when the materialized factor base reaches
the Minkowski bound; otherwise it explicitly reports an uncertified
factor-base model, not a bound.
This is a safe one-sided upper bound when certified, since adding omitted
principal relations can only shrink the quotient. For curve 273 the raw
Minkowski bound is currently far beyond the practical relation factor base, so
the current run correctly remains uncertified rather than making a false
completeness claim.

[`audit_bnf_free_s_class_quotient.py`](../cas/audit_bnf_free_s_class_quotient.py)
adds the practical, explicitly conditional alternative.  It reconstructs and
checks every closed relation from the retained principal generators, verifies
that the ledger contains every prime ideal over rational primes through its
declared factor-base bound, and computes the quotient of the resulting
GF(2) prime-ideal space by the verified relation rows and `S` columns.  It
then permits an upper bound only through one of two generation gates:

1. the unconditional cubic Minkowski bound; or
2. Bach's `12(log |Delta_K|)^2` prime-ideal bound after an explicit
   `--assume-erh` opt-in.

For the Fermigier cubic field the latter threshold is `262523`, versus an
impractical Minkowski threshold of about `2.9e31`.  Materializing the ERH
factor base gives 42,251 prime-ideal columns and requires no BNF, regulator,
or odd class-group calculation.  The resulting bound is conditional on
ERH/GRH and is labelled as such in JSON; it becomes useful only after enough
exact relations lower the quotient dimension.

Before the lattice collector contributes a single relation, every rational
prime in that base has a free exact principal relation `(p)`. The augmentation
tool stores these 23,034 generators and factorization rows explicitly. On the
ERH-complete Fermigier factor base, they lower the certified mod-two S-class
quotient upper bound from 42,207 to 19,204 (and the `K(S,2)` envelope from
42,235 to 19,232). This is a substantial certified reduction, but it is still
only an ERH-conditional ambient bound, not a Selmer bound.

A bounded complete-base Fermigier calibration used ten special primes just
above `262523`, `b`-bound `5000`, and the hybrid norm path with trial bound
`10000`. It sampled 300,000 principal generators and closed 483 exact
single-/double-large-prime relations. The verified relation lattice has rank
27 and lowers the **ERH-conditional** S-class quotient upper-bound model from
42,226 to 42,207 before the canonical `(p)` rows. Adding the 23,034 exact
canonical rows reduces that certified model to 19,204. Eliminating the
non-S columns of the collector rows yields 464 explicit
squareclass products; the 455-dimensional rational-square-norm kernel has an
explicit global-square basis. The older individual-row report had 251
products and raw signature rank 11 modulo the known rank-20
Mordell--Weil image. Neither the large class bound nor a raw signature rank is
a residual Selmer dimension.

The source-separated audit exposes the next bottleneck: those 483 collected
noncanonical Fermigier rows yield **zero** further mod-two quotient reduction
after the canonical `(p)` and `S` baseline is inserted. With the canonical rows
included, the non-S projection has 495 products and its 482-dimensional
square-norm kernel is again entirely global squares. The next collector
tranche must therefore be judged by rank gain after that baseline, not by raw
relation rank or cycle count.

The curve-273 pool now emits the same standard generator/closure fields after
its large-prime elimination and reports progress after the canonical `(p)` and
`S` quotient baseline. Its 20-special bounded run has zero gain beyond that
baseline, correctly identifying an unproductive relation family. It remains
`UNCERTIFIED_FACTOR_BASE`: the ERH generation threshold is `1231857`, far
above its 5,000-prime base. Thus the same gate can be applied to a future
complete curve-273 relation collection without changing the global/local
quotient machinery.

The complementary arbitrary-ideal pool has an exact sparse-incidence handoff.
It ranks unresolved large-prime ideals by elimination-pivot status, incidence,
and rational-prime size, then writes only the prime and degree-one residue
needed by the full three-dimensional lattice collector. The pinned checkpoint
[`curve273_large_prime_target_plan_20260822.json`](../../artifacts/generated-results/elliptic-curves/curve273_large_prime_target_plan_20260822.json)
records 444 exact rows on 962 large-prime columns, of rank 444 and nullity
zero; its twelve targets are therefore **search targets**, not relation
dependencies or Selmer classes. It was generated with:

```bash
sage -python elliptic-curves/cas/analyze_curve273_relation_pool.py \
  --include-full-ideal-chain --include-crt-cycle-logs \
  --include-ideal-lattice-logs \
  --write-large-prime-target-plan \
  artifacts/generated-results/elliptic-curves/curve273_large_prime_target_plan_20260822.json \
  --target-plan-count 12
```

This was run with SageMath 10.9; the resulting JSON has SHA-256
`0109839ec26c89840c7a2bd6da821c98c8e39383c66330854a82ab86279c8690`.
Pass one selected target directly into the lattice collector, for example,
with `--target-plan ... --target-plan-rank 4`; plan entries are intentionally
run independently, rather than multiplied into one unrelated target ideal.

For example, the smallest high-incidence planned ideal
`1858485403:649378816` had a bounded radius-24, five-shape pass with no exact
large-prime closure. Its only reported smooth norm was a rational-prime
multiple, which has no quotient value after the canonical `(p)` baseline.
That negative result is deliberately not promoted to any class-group,
`K(S,2)`, or Selmer assertion.

For a relation ledger whose factor base reaches that threshold, run:

```bash
sage -python elliptic-curves/cas/run_fermigier_rank20_fixedfb_specialq.py \
  --factor-base-bound 262523 \
  --norm-factor-mode hybrid --trial-prime-bound 10000 \
  --relation-ledger artifacts/local/elliptic-curves/r20_grh_relations.json

sage -python elliptic-curves/cas/augment_bnf_free_canonical_principal_relations.py \
  --relation-ledger artifacts/local/elliptic-curves/r20_grh_relations.json \
  --output artifacts/local/elliptic-curves/r20_grh_relations_canonical.json

sage -python elliptic-curves/cas/audit_bnf_free_s_class_quotient.py \
  --assume-erh \
  --relation-ledger artifacts/local/elliptic-curves/r20_grh_relations_canonical.json \
  --output artifacts/local/elliptic-curves/r20_grh_s_class_audit.json
```

The evaluator can carry this verified class-quotient record into the early
local/fingerprint quotient, while refusing a mismatched cubic field:

```bash
sage -python elliptic-curves/cas/evaluate_bnf_free_signature_map.py \
  --signature-map artifacts/local/elliptic-curves/fermigier_rank20_signature_map.json \
  --candidates artifacts/local/elliptic-curves/r20_ks2_candidates.json \
  --class-quotient-audit artifacts/local/elliptic-curves/r20_grh_s_class_audit.json \
  --output artifacts/local/elliptic-curves/r20_candidate_images.json
```

This certifies only the `K(S,2)` envelope.  The norm equation and all local
descent conditions are still required before a residual 2-Selmer basis or
rank upper bound can be claimed.

Even a certified (S)-class bound must still be joined to the local descent
conditions before it gives a 2-Selmer bound.  If a residual Selmer basis is
obtained, the next arithmetic step is the Cassels--Tate pairing before any
cover-point search; neither operation is implemented or claimed by this
bookkeeping module.

Run its regression checks with:

```bash
PYTHONPATH=elliptic-curves:elliptic-curves/cas \
  python3 -m unittest discover -s elliptic-curves/tests \
  -p 'test_residual_selmer_quotient.py' -v
```
