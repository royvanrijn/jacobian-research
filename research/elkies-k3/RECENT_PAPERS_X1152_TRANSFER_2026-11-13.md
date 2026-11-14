# Recent papers and exact transfer to the 1152 parents

This is the canonical source for the November 13–14 adoption work. Mathematical
status remains in [MATH_STATUS.json](../MATH_STATUS.json). This work reproduces
external results and makes scoped computations; it claims no new rank record.

## Sources and distinct objects

- C. R. Hill, [Two Rational Families of Elliptic Curves with Trivial Torsion and
  Generic Rank at Least Eighteen](https://hillaxiom.com/papers/two-rational-rank-18-families.html),
  November 11, 2026.
- Hill, [Nineteen Independent Sections from Two Quadratic Covers of an Elliptic
  K3 Surface](https://hillaxiom.com/papers/k3-19-sections.html), October 24, 2026.
- [Rank Ledger](https://github.com/crhillresearch/rank-ledger), immutable root tree
  1c63174ed3b85002b229c8b5a08b53466077886a. The imported files are pinned by
  Git blob identity and SHA-256 in the [source manifest](data/external/hill2026/SOURCES.json).
- Elkies–Klagsbrun–Li, [ICARM4239](https://elliptic-rank.icarm.cloud/curve/4239),
  submitted October 15, 2026. The [raw snapshot](../elliptic-curves/data/external/icarm_curve4239_20261113.json)
  carries its equation, 32 witness points and published parent equation.

Hill's rootless parent has 17 displayed independent sections with Gram
determinant 1152. ICARM4239 is reported to come from a rank-17 fibration of the
K3 surface of discriminant 1152. Equality of these numbers does not identify
the fibrations or the underlying surfaces. Their geometric elliptic fibrations are now [proved distinct](X1152_JMAP_FIBRATION_DISTINCTION_2026-11-14.md) by intrinsic j-critical values. Their underlying K3 surface relationship remains UNKNOWN.
Neither is identified with our determinant948, determinant1092 or
determinant1020 parent by this import.

The generic [ICARM4239 equation export](data/fibrations/icarm4239_x1152_generic_parent.json)
contains no exceptional parameter or witness points. Its full generic rank17
basis has not yet been independently imported. Four rational sections are
verified, with a dependence among three; no four-point independence is claimed.

## Published replay and record control

The byte-preserved upstream [rank18 checker](data/external/hill2026/verify_rank18.py)
and [rank19 checker](data/external/hill2026/verify_rank19.py) both pass under the
local Sage10.9 environment:

- 17 exact parent sections; positive height Gram of determinant1152;
- both published rational covers: 18 exact sections, determinant452984832,
  inherited block twice the parent Gram, extra Schur complement3;
- the paired genus-one base has an explicit nontorsion rational point;
- the two independent quadratic characters give 19 independent sections
  over that paired base, and hence infinitely many rank-at-least19 fibres;
- the published finite-field torsion arguments replay.

This is a rerun of the authors' checkers, not an independent full implementation.
The upstream SHA256SUMS file has a stale README entry at the fetched tree;
all mathematical data and checker bytes match both their Git blob identities
and the mathematical checksums. No checksum was rewritten.

Our [finite-quotient checker](../elliptic-curves/cas/verify_icarm_curve4239_rank32.py)
separately certifies the public rank-at-least32 curve. All32 points lie on the
equation and its integral short model exactly. A33-by32 binary matrix from24
good primes has rank32; the2-division cubic is irreducible modulo7. Infinite
descent proves independence unconditionally. The certificate was replayed with
an initially empty arithmetic cache using standard Python. This proves neither
exact rank32 nor a new programme discovery.

The [parent checker](scripts/verify_icarm4239_x1152_parent.py) verifies the public
specialization at t=-2868768/362209 by an exact Q-isomorphism and transports all32
points back. It is a KNOWN_RECORD_CONTROL, excluded from prospective selection.

## Frozen generic-only transfer pilot

The [producer](scripts/pilot_hill1152_transfer.py) freezes all choices before
running. It selects one sign-preferred representative per parity among generic
basis words of support at most2 and norm8 or10. No exceptional fibre, target
point or rank label selects a word. Limits are120 CPU seconds,180 external wall
seconds,4GiB address space; the completed v2 run used less than1 CPU second.

The [frozen packet](../artifacts/generated-results/hill1152-transfer-20261113-v2/input.json)
contains35 norm8 words and13 norm10 words. A common Mobius chart moves trace
poles away from infinity. Exact trace/chord identities are retained for every
row. The earlier v1 stopped on a degenerate residual; its inputs, partial
checkpoints and source are preserved locally. The frozen v2 successor retains
all nine repeated-root quadratic residuals as DEGENERATE_RESIDUAL, rather than
silently dropping them. The completed bank contains no linear residual.

### Complete pencils: a finite negative result

For each of35 norm8 traces, the regular residual chord slopes form a complete
one-parameter pencil, and its branch quartic is represented by an invertible
5-by5 rational matrix applied to [1,lambda,...,lambda^4]. This includes the
projective lambda=infinity endpoint.

All595 cross-pairs have disjoint rational projective branch images:

| Good prime | Newly excluded pairs | Remaining |
|---|---:|---:|
|101|559|36|
|103|27|9|
|107|9|0|

The [independent standard-library checker](scripts/verify_hill1152_transfer.py)
reconstructs the labelled traces from the imported rational basis, checks every
chord coefficient identity and recomputes the finite projective images. It
does not import Sage or the producer.

Disjoint projective branch images exclude equality of the quadratic extensions
represented by these smooth binary quartics. Projective equality alone would
still require the rational constant squareclass to match. This is only the
frozen35-pencil bank, not an exclusion of all bisections on this parent.
Singular images and larger trace words remain outside the claim.

### Four certified rational rank18 covers

Four of the13 norm10 traces give squarefree quadratic covers with rational
conic points; the other nine residuals are degenerate. The four accepted indices
are0,1,2,4. Their branch polynomials are pairwise distinct and disjoint from each
other and from the two published branch quadratics.

The [family loader and checker](scripts/verify_hill1152_four_covers.py) returns
the rational parameter map, polynomial Weierstrass model and18 explicit
rational-function points for each accepted cover. Each degree48 discriminant is
squarefree, so all singular fibres are I1. Exact18-by18 height matrices are
positive definite, have determinant452984832, restrict to twice the parent
Gram and have Schur complement3. Thus each family has generic rank at least18.

These are additional reproducible presentations from this selected bank; no
priority, exact-rank or saturation claim is made. Generic independence must be
checked again at a chosen specialization. The certificates are in
[four-generic-covers.json](../artifacts/generated-results/hill1152-transfer-20261113-v2/four-generic-covers.json).

Pairing the four new branch quadratics with each other and the two published
ones gives15 disjoint biquadratic genus-one bases. Necessary primitive local
square conditions modulo2^3,3^2,5^2,7^2,11^2 reject none. This initial finite gate supplied no global conclusion. The later [paired-Jacobian and complete first-pair descent](X1152_PAIRED_JACOBIANS_AND_FIRST_DESCENT_2026-11-14.md) proves infinite rational points for pair[0,1], using a complete CAS descent and independently checked all-place local witnesses. The other14 bases retain UNKNOWN global solubility.

## Eight-fibre commissioning, November14

The four new generic loaders were specialized at the fixed parameters0 and1,
giving eight predeclared fibres. The first frozen controller asked for two lazy
centres with a2000-node CVP allowance. Seven worker outputs reached the CVP
budget with zero searched charts; the eighth worker stopped at12 seconds.
A reporting bug then read MWState's nested reductions as a top-level basis.
Those original receipts remain unchanged and are not eight point-search
negatives.

A separately frozen successor used the identical eight fibres and two direct
centres: the first inherited section and the new cover section. Before search,
it exactly undid the redundant polynomial pullback scale: if u=N/D and the old
parameter is c+1/u, the specialized cover coordinates differ from the original
parent by x_cover=N^4*x_parent, y_cover=N^6*y_parent. Removing this square/cube
scale needs no factorization and preserves the labelled subgroup.

Every starting subgroup certified rank18. All16 declared pointed-quartic
searches completed their height10000 boxes at the0.5-second-per-chart allowance.
Every map, rational hit and point-admission transition replayed; no additional
independent direction was certified. Charged wall time was70.639 seconds.
This is a commissioning result, not a family-productivity comparison, exact
rank, new record or justification to enlarge the search.

The [portable bundle](../artifacts/generated-results/elliptic-curves/hill1152_specializations_20261114_v1.zip)
retains both frozen input protocols, all eight exact source and transported
requests, the16 search witnesses and replay receipts.
[verify_hill1152_specialization_panel.py](scripts/verify_hill1152_specialization_panel.py)
replays it without a point search, CVP enumeration or Sage. The
[compact replay](../artifacts/generated-results/elliptic-curves/hill1152_specializations_replay_20261114_v3.json)
records all eight lower bounds and zero added directions. Initial and final finite-reduction certificates were rebuilt in newly empty per-fibre caches, rather than trusted from the discovery snapshot.

The later [j-map distinction](X1152_JMAP_FIBRATION_DISTINCTION_2026-11-14.md) resolves fibration identity, and the [first complete descent](X1152_PAIRED_JACOBIANS_AND_FIRST_DESCENT_2026-11-14.md) advances global cover arithmetic. Further work concerns underlying surfaces, the missing record-parent basis and explicit paired-base points, rather
than repeating this completed allowance. The four generic loaders and the
eight certified fibre subgroups remain reusable inputs.

## Mestre obstruction transfers to both 1152 equations

For E:y^2=x^3+A(t)x+B(t), the retained Mestre identity at auxiliary u=2 is

    D=-5*A*B*(9261*B^2+400*A^3),
    P=(-21*B/(5*A),sqrt(D)/(25*A^2)),
    Q=(-21*B/(20*A),sqrt(D)/(200*A^2)).

The exact polynomial and sum identities are the same as in our
[existing proof](R17_MESTRE_CORRELATED_SECTIONS_AND_GENUS_GATE_2026-09-13.md).
The following new witnesses transfer its full root-field valuation argument:

| Parent | A irreducible mod | A simple local root | B irreducible mod | B simple local root | Smooth disjoint branch witness |
|---|---:|---|---:|---|---:|
|Hill1152|11|4 mod83|37|84 mod103|179|
|ICARM4239 parent, short model|41|20 mod71|107|62 mod151|197|

The A local primes are11 mod12, excluding i and cube roots of unity from its
degree8 root field; the B local primes are3 mod4, excluding i from its degree12
root field. Both polynomials are irreducible by their full-degree irreducible
reductions. For every rational auxiliary u(t), the valuation table in the
referenced proof forces all20 geometric roots of A and B to stay branched.
Consequently this identity has covering genus at least9 on either parent.
At u=2 the degree44 branch divisor is squarefree, disjoint from the squarefree
degree24 parent discriminant, and the genus is21. The anti-invariant section
Gram is diag(24,24). Thus the two directions exist but this route cannot provide
an infinite rational base. The obstruction is to this identity, not all covers.

## Rational-point methods and next gates

See the [separate Xns121 method transfer](../elliptic-curves/notes/XNS121_METHOD_TRANSFER_2026-11-13.md).
The next useful tasks are:

1. retain the proved geometric fibration distinction; determine the underlying
   K3 surface relation and import a full basis for the record parent;
2. reuse the completed eight-fibre commissioning and its rank18 certificates;
   require a new mathematical selection or arithmetic gate before another panel;
3. use the15 explicit paired Jacobians and the completed first-pair descent;
   extract an explicit point on pair[0,1] before a further family commissioning;
4. apply quotient and Selmer-direction arguments to third-cover incidence
   curves or marking curves of genus>1, where rational-point finiteness is
   relevant. A positive-rank genus-one base is meant to have infinitely many
   points and is not itself a finite-point target.

No old stalled-fibre campaign is restarted by this note. The stronger goal of
two independent new directions on one quadratic extension of an actual MW17
parent, with infinitely many rational base points, remains open.

## Replay

From the repository root:

    sage -python research/elkies-k3/data/external/hill2026/verify_rank18.py
    sage -python research/elkies-k3/data/external/hill2026/verify_rank19.py
    python3 research/elkies-k3/scripts/verify_hill1152_transfer.py \
      --input-dir research/artifacts/generated-results/hill1152-transfer-20261113-v2
    sage -python research/elkies-k3/scripts/verify_hill1152_four_covers.py \
      --output /tmp/hill1152-four-covers-fresh.json
    python3 research/elliptic-curves/cas/verify_icarm_curve4239_rank32.py --check

The checkers' output options use exclusive creation. Choose a fresh output path.
Missing historical artifacts elsewhere in the checkout still block global
navigation audits; no absent evidence was regenerated or reclassified.
