# Recovering the other rank-17 fibration

<!-- status-consumer: EC-K3-H3-OTHER-R17-J2-CANDIDATE f1884d1f6168a934 -->
<!-- status-consumer: EC-K3-ELKIES-2026-BISECTION-VISIBILITY-RECORD-CURVES 1c39220ee5fedc77 -->
<!-- status-consumer: EC-K3-ELKIES-2026-R17-SMALL-ISOGENY-EXCLUSIONS fc2c4caaa79fb36c -->

## Outcome

The immediate structural result has two exact parts.

First, the published `R17` fibration is excluded for every current high-rank
target in the requested test set.  The primitive equations

```text
j_R17(t) = j(E)
```

are irreducible of degree 24 for the 2024 Elkies--Klagsbrun rank-29 curve and
ICARM curves 273, 302, and 398.  The newly ingested ICARM curves 399 and 400
have the same outcome.  The rank-28 control instead has the exact factor
`5471*t+9529`, recovering its published parameter.  Thus the negative result
is calibrated and survives quadratic twisting.

Second, the repository already contained an equation-open lattice candidate
for the missing fibration.  Replaying and composing its transports produces a
primitive rootless `U` embedding in the pinned H3/R17 Neron--Severi lattice.
Its positive frame has rank 17 and determinant 948, has no norm-two roots, and
is not integrally isometric to the published `rank17_gram.txt`.  Therefore the
pinned lattice has at least two distinct rootless `J2` frame classes.

This is not yet the rank-29 identification.  The second frame still lacks a
generic characteristic-zero Weierstrass equation and hence has no `j`-map to
test against the rank-29 curve.

## Gate A: intrinsic fingerprint and rank-29 embedding search

The alternate frame now has an exact intrinsic fingerprint.  With
`theta_L(q)=sum_v q^(v.v/2)`, its exactly enumerated initial theta expansion is

```text
1 + 2626*q^2 + 53290*q^3 + 460360*q^4 + O(q^5).
```

The `O(q^5)` is a declared computation boundary, not a claim that the full
infinite theta series is determined.  There are 1,313 unoriented minimal
norm-four lines.  Joining two lines when their absolute pairing is two gives
a connected graph with 92,676 edges and the complete degree distribution
stored in the certificate.  The lattice automorphism group has order four;
its action on the cyclic discriminant group `Z/948Z` has image
`{1,473,475,947}`.  The minimal-line count already differs from the published
frame's 1,311, in addition to the exact `qfisom` non-isometry.

Specialization distortion was calibrated in the correct basis: the seventeen
known generic sections on each published rank-25--28 fibre.  After fitting the
scalar in `H_t approximately h*G_published`, the observed errors are:

| control | fitted `h` | relative Frobenius error | maximum `|residual|/h` |
| --- | ---: | ---: | ---: |
| rank 25 | 9.790535 | 0.174808 | 0.829330 |
| rank 26 | 10.335812 | 0.158529 | 0.788986 |
| rank 27 | 10.547799 | 0.170800 | 0.704497 |
| rank 28 | 13.624052 | 0.135631 | 0.732249 |

The replacement rank-29 test enumerates all 11,692 unoriented integral vectors
of displayed height at most 70 in the full 29-dimensional public-point
lattice.  It searches integral `29 by 17` matrices, rather than subsets of the
29 generators.  Across 1,000 deterministic seeded construction/refinement
trials, the best matrix has exact rational column rank 17 and fitted scale
`h=15.394989`, but relative Frobenius error `0.210330` and maximum
`|residual|/h=0.882957`.  It therefore misses both positive-control envelopes.

This is a calibrated bounded negative experiment.  It closes the old subset
comparison as evidence, but it neither proves non-embedding nor excludes a
generic subgroup requiring a vector above height 70 or a better discrete
search path.  In particular it does not demote the alternate equation as the
main reconstruction target.

The exact invariant artifact is
[`elkies-k3-other-rank17-invariants.json`](../artifacts/generated-results/elkies-k3-other-rank17-invariants.json),
and the calibrated bounded-search ledger is
[`elkies-k3-other-rank17-rank29-gate-a.json`](../artifacts/generated-results/elkies-k3-other-rank17-rank29-gate-a.json).

## Exact published-fibration exclusions

| target | rank lower bound | degree | irreducible-mod-prime witness |
| --- | ---: | ---: | ---: |
| Elkies--Klagsbrun 2024 / ICARM 12 | 29 | 24 | 461 |
| ICARM 273 | 30 | 24 | 367 |
| ICARM 302 | 31 | 24 | 397 |
| ICARM 398 | 30 | 24 | 1009 |
| ICARM 399 | 29 | 24 | 83 |
| ICARM 400 | 28 | 24 | 157 |

The rank-28 control factors with degrees `1+23`; its degree-23 cofactor is
irreducible modulo 197.  The exact certificate and the primitive integer
coefficient vectors are in
[`elkies_2026_bisection_visibility_record_curves_v1.json`](../artifacts/generated-results/elliptic-curves/elkies_2026_bisection_visibility_record_curves_v1.json).
The public target equations and source hashes are pinned in
[`elkies_2026_r17_j_recognition_targets.json`](../elliptic-curves/data/elkies_2026_r17_j_recognition_targets.json).

The separate small-prime isogeny gate also excludes rational published-R17
parameters satisfying

```text
Phi_ell(j_R17(t),j(E))=0,  ell=3,5,7,11,
```

for all six targets.  Every test has a clean projective no-root witness; the
witness table and complete modular factor degrees are in
[`PUBLISHED_R17_SMALL_ISOGENY_EXCLUSIONS_2026-09-01.md`](PUBLISHED_R17_SMALL_ISOGENY_EXCLUSIONS_2026-09-01.md).
This does not test the alternate Q80-derived fibration or other isogeny
degrees.

The 2024 announcement is the decisive provenance clue: it says the rank-29
curve was found by specializing a rank-17 fibration of the same K3 surface and
then searching outside the generic rank-17 subgroup.  It does not publish the
other equation or its marking.  The local statement “same pinned H3 surface”
therefore remains a reconstruction hypothesis until a candidate equation
specializes to that curve.

## The recovered second frame

The existing alternate q80 ending is

```text
alternate A1/MW16 --q6--> rootless/MW17.
```

The final raw fibre has `(a,b)=(2,3)`.  Reflection in the old zero and the
single `A1` component gives a degree-two nef representative with

```text
D.F=2,  D.O=1,  MW norm=23/2.
```

The exact horizontal section satisfies

```text
D = O + S - F,  S.O=4,
```

so the generic-fibre pencil is generated by

```text
1, (y+y(S))/(x-x(S)).
```

The complete wall audit gives minimum section pairing one and excludes a
negative bisection by parity.  A nef ambient line bundle has `h0=5`; three
disjoint elementary transforms should cut it to the required two-plane.  This
is materially cheaper than generic unrestricted Riemann--Roch reconstruction.

The new canonicalization certificate composes the alternate q80 basis with
the exact q80-to-pinned-R17 transport.  It stores the fibre, isotropic mate,
zero, and the full determinant-one 19-by-19 change of marking in pinned NS
coordinates:

[`elkies-k3-other-rank17-candidate.json`](../artifacts/generated-results/elkies-k3-other-rank17-candidate.json).

The distinctness proved here is `J2` distinctness of frame lattices.  The
complete `J1` classification up to automorphisms of the K3 surface remains
open, but Corollary H2a of
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md)
now proves that the rootless `J1` count is finite and lies between two and
eight.

## Candidate ranking by actual equation cost

| candidate | structural value | old-fibre degree | RR/local cost | equation status | priority |
| --- | --- | ---: | --- | --- | ---: |
| published R17 | exact control and exclusion map | 2 on final construction edge | completed | exact over `QQ` | control |
| alternate q80 q4--q6 frame | unique non-published rootless `J2` class; strongest rank-29 provenance candidate currently present | 2 | `h0=5`, three local transforms, final kernel 2 | generic characteristic-zero parent/child open; exact `GF(73)` specialization shadow | 1 |

The apparent simplicity of the final q6 is conditional on reaching its own
generic parent equation.  Transporting the q80 corridor directly from the
current H3 `D12` equation has enormous degree and is not a viable shortcut.
Construction should resume from the q80 characteristic-zero checkpoints, not
by pulling this fibre through the published R17 equation.

## Classification programme

<!-- status-consumer: EC-K3-H3-ROOTLESS-J2-COMPLETE c6f054948b04b507 -->

The correct completeness target is a Kneser--Nishiyama `J2` classification,
not another unrestricted neighbor walk.  Nishiyama's method classifies frame
lattices through primitive embeddings of a discriminant-compatible auxiliary
lattice into the 24 Niemeier lattices; rootless complements select rank-17
fibrations immediately.  See [Nishiyama's original
paper](https://doi.org/10.4099/math1924.22.293) and a modern summary of the
[Kneser--Nishiyama construction](https://doi.org/10.1002/mana.12018).

For this Picard-19 lattice the implementation gates are:

1. pin one explicit positive-definite rank-seven auxiliary lattice whose
   discriminant form is opposite to the positive frame form (equivalently,
   compatible with the sign convention for `NS(H3)`);
2. enumerate primitive embeddings into every Niemeier root lattice modulo its
   Weyl group and glue automorphisms;
3. compute each orthogonal complement including primitive closure and glue;
4. retain complements with no roots and determinant 948;
5. canonicalize by integral frame isometry, then retain the full pinned `U`
   marking because frame isometry alone is insufficient for equation work;
6. refine from `J2` to `J1` only after the relevant surface automorphism action
   is known.

The two exact frames were used as mandatory positive controls.  The complete
classifier returns exactly those two classes.

The exact first-pass audit is now recorded in
[`ROOTLESS_J2_COMPLETENESS_TRACK_2026-08-31.md`](ROOTLESS_J2_COMPLETENESS_TRACK_2026-08-31.md).
The target genus has mass
`77731517730627488307787/925557271717478400`, forcing at least 167,967
isometry classes in the full genus.  Thus an unfiltered genus traversal is not
the cheap route.  The same audit rejects all 65 old files under
`seeds/target-genus-rootless-pneighbor`: they occupy a different 2-adic and
79-adic genus and contribute no new `J2` control.  Niemeier filtering should
therefore precede, not follow, exhaustive positive-definite enumeration.

The auxiliary and both Niemeier positive controls are exact.  The pinned
rank-seven lattice has cyclic quadratic module `[1267/948]` and root system
`D5`.  All eight discriminant anti-isometries for each control glue to
`N(2A7+2D5)` and recover the requested saturated complement.  Modulo
`Aut(K) x Aut(R)`, the published control has two primitive embedding orbits
and the alternate has one.

The all-24-Niemeier enumeration is also complete.  The Leech lattice and ten
rooted classes fail the exact `D5` gate; thirteen rooted classes give sixteen
full-automorphism `D5` anchor orbits.  Residual-Weyl enumeration produces
3,220 primitive sixth-vector representatives, of which only 21 can support a
rootless dominant seventh vector.  The 167 positive-label cases reduce to 25
exact fixed-space ellipsoid solutions and twelve integral primitive
embeddings.  Their saturated complements form exactly two integral isometry
classes: the published frame (eight representatives in the enumerated cover)
and the alternate frame (four), all in `N(2A7+2D5)`.  The cover counts are not
embedding-orbit counts because anchor-stabilizer duplicates are retained.
Thus there is no third rootless `J2` equation candidate; the alternate frame
is the unique non-published construction target.

## Finite J1 boundary

<!-- status-consumer: EC-K3-H3-ROOTLESS-J1-UNIFORM-BOUND b71330a75ad2c9ad -->

Picard rank 19 forces the rank-three transcendental Hodge-isometry group to be
exactly `{+id,-id}`.  The cyclic determinant-948 quadratic form has eight
isometries,

```text
{1,157,317,473,475,631,791,947} modulo 948.
```

Braun--Kimura--Watari Proposition C' therefore bounds the `J1` multiplicity
of each `J2` class by four.  Since the Niemeier classification gives exactly
two rootless `J2` classes, the rootless `J1` count lies in `[2,8]`.  This is a
usable finite search boundary, not the exact classification.  The smaller
quotients suggested by the frame automorphism images are not promoted because
they still require ample-cone stabilizer control.

## Rank-29-first construction gate

For each constructed candidate equation, compute its reduced `j`-map and form
the primitive recognition polynomial for ICARM 12 first.  A rational linear
factor is only a candidate parameter; the specialized minimal model must then
be checked for equal `j`, twist class, and finally a `QQ`-isomorphism to the
rank-29 model.  If the rank-29 curve is found, specialize all seventeen generic
sections and certify that they span the expected subgroup before replaying the
outside-subgroup search.

Curves 398--400 should then be tested against the same candidate `j`-map.
Their common authorship and discovery dates are useful prioritization evidence,
not proof of a common fibration.

## Replay

```bash
.venv/bin/python \
  elliptic-curves/scripts/analyze_elkies_bisection_visibility_and_record_curves.py \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python -c \
  'import sys; from sage.all import *; sys.argv=["verify_q80_alternate_fifth_q6_rootless.sage","--write-artifact"]; globals()["__file__"]="/home/royvanrijn/src/jacobian-research/elkies-k3/scripts/verify_q80_alternate_fifth_q6_rootless.sage"; load(__file__)'

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python -c \
  'from sage.all import *; globals()["__file__"]="/home/royvanrijn/src/jacobian-research/elkies-k3/scripts/verify_q80_alternate_final_q6_nef.sage"; load(__file__)'

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python -c \
  'from sage.all import *; globals()["__file__"]="/home/royvanrijn/src/jacobian-research/elkies-k3/scripts/canonicalize_other_rank17_candidate.sage"; load(__file__)'

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_rootless_j2_niemeier_controls.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/classify_rootless_j2_niemeier_first.sage --check

python3 elkies-k3/scripts/certify_rootless_j1_uniform_bound.py --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_elkies_2026_r17_isogeny_exclusions.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compute_other_rank17_invariants.sage

python3 elkies-k3/scripts/gate_a_alternate_rank17_rank29.py \
  --short-bound 70 --trials 1000 --seed 291317 \
  --top-choices 8 --refinement-passes 2
```

The first command proves the `j`-exclusions and control factor.  The next two
replay rootlessness, non-isometry, nefness, and the cost model.  The following
command pins the alternate `U` embedding in H3/R17 NS coordinates, and the
first Niemeier replay certifies the auxiliary and both positive-control
gluings.  The complete Niemeier replay proves that those controls are exactly
the two rootless `J2` frame classes.  The dependency-free `J1` replay checks
the four-coset uniform multiplicity bound and the resulting interval `[2,8]`.
The final two commands compute the exact bounded intrinsic invariants and
replay the calibrated full-lattice Gate A search.  None constructs the missing
characteristic-zero equation or completes the exact `J1`
surface-automorphism classification.
