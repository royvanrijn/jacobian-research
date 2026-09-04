# Determinant-aware re-ranking of the rank-seven foundry

<!-- status-consumer: EC-K3-NS0031-MARKED-FORMAL-BRANCH b31e99bce4edac0a -->
<!-- status-consumer: EC-K3-NS0031-QQ-MARKING-OBSTRUCTION 8e2dc35cdf9b6bc3 -->
<!-- status-consumer: EC-K3-GOLAY-DET720-QQ-MARKING-OBSTRUCTION 972f591d2885f9ba -->
<!-- status-consumer: EC-K3-DIFFERENT-NS-ARITHMETIC-GATE-RERANK 252991e141c42e55 -->
<!-- status-consumer: EC-K3-RANK19-ARITHMETIC-MARKING-CLASSIFIER 93e6c5626d369572 -->
<!-- status-consumer: EC-K3-ARITHMETIC-FIRST-MARKED-T-FOUNDRY 6b9d34ae8d722280 -->

## Result

The rank-seven surface catalogue now has a fail-closed pre-solver ranking:

[`../artifacts/generated-results/elkies-k3-rank7-determinant-aware-ranking-v1.json`](../artifacts/generated-results/elkies-k3-rank7-determinant-aware-ranking-v1.json).

It changes the old policy in two ways.

1. A surface is not scored as a rootless MW17 target until it passes four
   exact or necessary lattice gates: an explicit rootless rank-17 frame, the
   Hermite--Blichfeldt determinant bound, discriminant length at most three,
   and an explicit even rank-17 Gram witness with the catalogue discriminant
   form.
2. Determinant is no longer a scalar objective to minimize.  The surviving
   surfaces are placed in determinant regimes, then ordered by typed
   arithmetic/source/corridor readiness and exact short-vector quality.
   Rank-jump mechanism evidence is a separate coordinate.  Raw bisection,
   trisection, and quadrisection counts do not enter the priority key.

The current exact accounting is:

```text
827 imported surface classes
761 deferred before scoring
 66 explicit even rootless MW17 witnesses
 11 in the observed dense determinant band 500--948
  3 above-band short-vector-quality exceptions
 52 above-band sparse-pressure rows
  1 fully ready for expensive equation scoring
```

Here “deferred” is deliberate.  The catalogue is not complete in a
determinant band, so absence of a rootless frame is not a nonexistence
theorem.

## The lower boundary

For a positive-definite rank-17 lattice with Gram determinant `D` and minimum
at least four, the Hermite invariant satisfies

```text
gamma(L) >= 4 / D^(1/17).
```

Blichfeldt's bound

```text
gamma_17 <= (2/pi) Gamma(2+17/2)^(2/17)
```

therefore gives

```text
D >= (4/gamma_17)^17 >= 28.865863123698173...
```

and hence the necessary integer boundary `D >= 29`.  This is only a
necessary packing bound.  The new scorer does not promote `D >= 29` to the
existence of an even lattice, the realization of a prescribed discriminant
form, or a K3 embedding.

The discriminant-length gate is independent.  A Picard-rank-19 K3 has a
rank-three transcendental complement, so its discriminant group has length
at most three.  The positive-definite rank-17 existence gate is stronger than
this numerical length test: it requires an actual even rootless Gram witness
in the imported catalogue.  The relation with primitive embeddings and
discriminant forms is the one developed by Nikulin in
[*Integral symmetric bilinear forms and some of their applications*](https://doi.org/10.1070/IM1980v014n01ABEH001060).

## Determinant regimes

The present imported catalogue first realizes a rootless MW17 frame at
determinant 500.  The published R17 reference has determinant 948.  The
ranking therefore records three logically different regions:

```text
29--499   low-determinant design frontier, if an explicit form is realized
500--948  observed dense rootless sweet spot
>948      above-band sparse pressure, with exact short-vector exceptions
```

The interval `29--499` is not asserted empty.  It is precisely the range in
which a new construction can improve the observed boundary.  Within a
regime, raw determinant is not minimized.  The exact number of unoriented
norm-four pairs is used as a lattice-density coordinate, with the published
R17 value 1,311 as a named benchmark; it is explicitly not a rank-jump
predictor.

The three imported determinant-384 surfaces currently have maximum
catalogued MW ranks `15,13,13`, with root systems `A2,4A1,4A1`.  They fail the
pre-solver gate because no rootless rank-17 Gram has yet been attached.  This
does not exclude another fibration or another lattice in the required genus.
Any determinant-384 construction that supplies the missing even rootless
rank-17 witness will automatically enter the low-determinant design-frontier
tier on the next replay.

## Genus-theoretic replacement for the witness gate

<!-- status-consumer: EC-K3-ROOTLESS-GENUS-MASS e7589727ca8f7e50 -->

The present ranking is witness-first because the imported catalogue supplies
Gram matrices rather than complete local-genus certificates.  The next
version should insert a genus gate before explicit class generation.  For
each discriminant form, compute the rootless mass `mu_0(G)` from weighted
local ADE representation averages and assign one of three exact states:

```text
mu_0(G)=0       genus excluded before class generation,
mu_0(G)>0       genus proved rootless-capable and prioritized,
not computed    fail-closed UNKNOWN.
```

The degree-one shortcut `[q]Theta_G<2` is a cheap sufficient certificate;
selected higher ADE averages give intermediate linear-programming bounds,
and the full triangular root-system inversion decides the sign of `mu_0`.
This changes the logical input from “an explicit rootless class was already
found” to “the local genus is proved capable of containing one.”  Explicit
Gram witnesses remain necessary later for construction and equation work.

The first exact replay does not yet change the 827-row ranking.  Its weighted
root means for the determinant-78, 948, and 950 controls are approximately
`49.13`, `14.15`, and `14.08`, all above two.  Thus these controls require
higher ADE data; determinant 78 remains the mandatory zero-mass calibration.
The rank-at-most-four LP is now also known to be insufficient on determinant
78: even after restricting to its 621 census-realized ADE types, it permits a
rootless fraction `0.00275682...` although the true value is zero.

Large compatible Kneser neighbours turn an exact rootless mass fraction into
an asymptotic line-hit probability.  The complete determinant-948 rootless
classification gives exact fraction `8.930328...e-6`, or an unguided baseline
of about `111,978` compatible lines per hit.  Four distinct determinant-950
classes give only the lower probability bound `1.020293...e-5`, hence an
upper baseline of about `98,011` lines per hit.  These scores apply to the
unmarked one-spinor-genus controls and are asymptotic; they are not measured
small-prime frequencies or speedup claims.  The theorem, exact values, and
implementation boundary are in
[`ROOTLESS_GENUS_THEORY_2026-09-03.md`](ROOTLESS_GENUS_THEORY_2026-09-03.md).

## Readiness gates and current priority

After the lattice gates, the scorer attaches three fail-closed readiness
coordinates:

1. arithmetic field-of-definition evidence;
2. a source-equation precursor;
3. an exact compiler corridor.

The strict expensive-scoring queue requires an exact rational
field-of-definition record, a positive source precursor, and at least a
marking-level exact corridor.  On present evidence it contains only the
determinant-948 H3/R17 surface.  This is expected: it has exact equations and
markings over `QQ`, a complete equation-level corridor, and the separate
Pasten--Salgado non-thin rank-jump theorem recorded in
[`PASTEN_SALGADO_NONTHIN_RANK_JUMPS_2026-08-31.md`](PASTEN_SALGADO_NONTHIN_RANK_JUMPS_2026-08-31.md).

The former first unresolved row, the determinant-720 Golay surface, is now
arithmetically excluded. Its literal stable marking curve is `X_0(60)`, whose
rational points are cusps. Its one-parameter formal `Z_7` precursor and exact
marking-level corridor remain geometric controls only.

The arithmetic pre-screen now has an explicit rejection ledger. Of 827
catalogue surfaces, 66 pass the exact lattice-theoretic MW17 filters;
determinant `950`/`NS0024` and determinant `1184`/`NS0031` are then removed by
their rational-marking obstruction certificates, together with determinant
`720`, leaving 63 arithmetic candidates. The determinant-948 control remains
the only fully ready row and is outside the different-NS milestone. There is
therefore no unresolved rootless-frame row authorized for equation work.

The upstream arithmetic-marking classifier types these as one
`ARITHMETICALLY_POSSIBLE`, three `ARITHMETICALLY_EXCLUDED`, and 62 `UNKNOWN`.
Only the first type can enter expensive equation scoring; the determinant-948
positive control is already realized, so the new different-NS equation-agent
handoff is empty.

The operational priority now comes from the global `T`-first planner. It
orders all 827 transcendental rows before NS/rootless inspection and has an
823-row arithmetic research queue. Its first exact-coarse diagnostics are
determinants `378`, `256`, and `512`; their `X_0(7)`, `X_0(2)`, and `X_0(4)`
labels do not replace the still-missing literal stable-kernel calculations.

The determinant-1184 NS0031 surface has a one-parameter formally smooth
`Z_7` marked branch and marking-level corridor evidence, but it is now
arithmetically excluded over `QQ`: the split-Clifford modular curve is
`X_ns(4) x X_0(37)` over the `j`-line, and neither noncuspidal rational point
of `X_0(37)` passes the mod-4 Cartan Frobenius test. It is retained as a
geometric/local control, not an equation candidate. Future rankings must put
this rational-marking screen ahead of equation cost.

## Rank-jump coordinate

The rank-jump coordinate has its own ordering and is excluded from the main
priority tuple.  At present only H3/R17 has theorem-level evidence.  The
other rows remain `UNKNOWN`; low-degree multisection abundance is not used as
a substitute.  This follows the calibration in
[`R17_MULTISECTION_DIVERSITY_2026-09-02.md`](R17_MULTISECTION_DIVERSITY_2026-09-02.md),
where global cover counts and fibre-specific exceptional directions are
shown to be different data.

## Reproduction

Run:

```bash
sage -python elkies-k3/scripts/build_rank19_arithmetic_marking_classifier.sage
python3 elkies-k3/scripts/build_rank7_determinant_aware_ranking.py
python3 elkies-k3/scripts/build_arithmetic_first_marked_t_foundry.py
sage -python elkies-k3/scripts/build_rank19_arithmetic_marking_classifier.sage --check
python3 elkies-k3/scripts/build_rank7_determinant_aware_ranking.py --check
python3 elkies-k3/scripts/build_arithmetic_first_marked_t_foundry.py --check
```

The typed evidence overlay is
[`data/lattice-foundry/determinant-aware-ranking-evidence-v1.json`](data/lattice-foundry/determinant-aware-ranking-evidence-v1.json).
It names every non-default arithmetic, equation, corridor, and rank-jump
claim and links it to an existing exact artifact or canonical note.  Missing
evidence remains `UNKNOWN` and cannot authorize the expensive queue.

The Blichfeldt bound is a classical geometry-of-numbers filter.  For the
elliptic-K3 determinant and Mordell--Weil conventions used here, see
Elkies--Kumar,
[*K3 surfaces and equations for Hilbert modular surfaces*](https://msp.org/ant/2014/8-10/ant-v8-n10-p01-s.pdf),
especially the Shioda--Tate and discriminant formulas summarized in its
elliptic-fibration preliminaries.
