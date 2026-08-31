# Elkies K3 / neighbor programme

This directory contains the K3-surface and elliptic-neighbor work used to search for high-rank elliptic fibrations and explicit equation routes.

For the chronological account of rank exchange, cancellations, corrections,
equation lifts, and reusable methods, start with
[`ELKIES_K3_PROCESS_ATLAS.md`](ELKIES_K3_PROCESS_ATLAS.md). For the
general rank, determinant, fibration, specialization, and lift statements,
see [`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md).
Before extending the certified equation route or using its endpoint, read
[`PRIOR_WORK_SHORTCUT_AND_ENDPOINT_CERTIFICATION_2026-08-25.md`](PRIOR_WORK_SHORTCUT_AND_ENDPOINT_CERTIFICATION_2026-08-25.md).
It records what the historical construction supplies, the minimum
intermediate construction record, and the now-completed direct R17 endpoint
gate. Optional intermediate completeness must not delay arithmetic use of the
certified endpoint.
For the
human-readable map of how the reverse searches, H3 source construction, H2
comparison, Q80 route, and CM24 scaffold fit together, see
[`CONSTRUCTION_ROUTES.md`](CONSTRUCTION_ROUTES.md). The exact lossless reverse
transport from the pinned recovered 17-by-17 lattice to H3, including the q8
marking distinction, is documented in
[`RANK17_TO_H3_REVERSE_TRANSPORT_2026-08-23.md`](RANK17_TO_H3_REVERSE_TRANSPORT_2026-08-23.md).
For executable entry points and the research history behind
failed/superseded scripts, see [`scripts/README.md`](scripts/README.md) and
[`SCRIPT_ROUTE_AND_FAILURE_LEDGER.md`](SCRIPT_ROUTE_AND_FAILURE_LEDGER.md).

## One-minute status

- Think of one K3 surface with several coordinate systems. A neighbour changes
  the elliptic fibration, moving rank between fibre roots and sections.
- The generic budget is always `root rank + MW rank = 17`.
- The physical equation route continues through q8/orbit376 and the preferred
  q12/orbit5867 edge to the pinned rootless `MW17` endpoint; q12/orbit4484 is
  a certified but unnecessary fallback.
- The q4/orbit230 equation and q4 return are exact, but the formerly promoted
  q6/orbit1315 continuation is withdrawn as an equation-cost target: its cheap
  score used a chamber pseudo-zero rather than the effective P230 section.
- The older q6/orbit1307 horizontal survives physical Weyl reduction, but its
  advertised component-10 zero does not: after correction component 10 has
  degree zero, while components 3, 5, and 9 have degree one.  Its 10,334-point
  continuation is therefore also withdrawn pending a physical-zero rerank;
  see
  [`../artifacts/local/elkies-k3/q24-2a5-q6o1307-physical-nef-audit.json`](../artifacts/local/elkies-k3/q24-2a5-q6o1307-physical-nef-audit.json).
- Physical q4/orbit208 is complete: `P.O=0`, exact RR dimensions
  `4 -> 2 -> 2`, quartic degree 4, and a minimal `3I4+12I1` Jacobian.  The
  q4/orbit1584 edge then gives `D4+A3+3A1/MW7`, and q4/orbit164 gives
  `2A3+2A1/MW9`, both with exact characteristic-zero equations and the
  effective markings needed by the next edge.  The exact q8/orbit376
  horizontal and resolved Riemann--Roch compiler are complete.  The preferred
  q12/orbit5867 edge now gives an exact rootless Jacobian with `24I1` and
  seventeen exact independent sections whose determinant-948 height Gram is
  integrally pinned to R17. An exact point at the old `v=0` I2 support proves
  source identity; good-reduction counts at 131 and 137 prove geometric
  Picard rank 19; and the unique index-two candidate has odd norm 73, proving
  that the full geometric MW lattice is saturated R17 of exact rank 17.
- The separately certified first-q8 and D13 lattice detours remain available,
  but combined totals using the withdrawn 4,199 suffix are not valid compiler
  comparisons.  See
  [`A11_EQUATION_COST_ROUTE_HANDOFF_2026-08-24.md`](A11_EQUATION_COST_ROUTE_HANDOFF_2026-08-24.md).

<!-- status-consumer: EC-K3-H3-A11-Q8-QQ-2A5 b7aaf4bf483eac68 -->
<!-- status-consumer: EC-K3-H3-A11-R17-PHYSICAL-Q4O208-PROMOTED-ROUTE 86e7a3affbf35a9e -->
<!-- status-consumer: EC-K3-H3-A11-R17-PHYSICAL-Q10-PROMOTED-ROUTE e8dbe599e076f13d -->
<!-- status-consumer: EC-K3-H3-D13-R17-Q4O11-PROMOTED-LATTICE-ROUTE b648ea30fd562496 -->
<!-- status-consumer: EC-K3-H3-Q4O208-Q4O1584-QQ 8463d62d0e9f2b83 -->
<!-- status-consumer: EC-K3-H3-Q4O1584-Q4O164-QQ bafe854a24d6762b -->
<!-- status-consumer: EC-K3-H3-Q4O164-Q8O376-MOD131-HORIZONTAL 2249c509c1217d7c -->
<!-- status-consumer: EC-K3-H3-Q4O164-Q8O376-QQ-HORIZONTAL 688f0a5f6d989e9c -->
<!-- status-consumer: EC-K3-H3-Q4O164-Q8O376-QQ-SMOOTH-RR c0bc752848961743 -->
<!-- status-consumer: EC-K3-H3-Q12O5867-DEGREE1-SECTION-QQ 28056772646e9fc7 -->
<!-- status-consumer: EC-K3-H3-Q12O5867-QQ-ROOTLESS 6a3dbee5942ddd0a -->
<!-- status-consumer: EC-K3-H3-Q12O5867-QQ-R17-BASIS a2097150acf00645 -->
<!-- status-consumer: EC-K3-H3-Q12O5867-ENDPOINT-QQ a83b08acd921c32b -->
<!-- status-consumer: EC-K3-H3-Q4O208-R17-CURRENT-MARKED-ROUTE 432b34c44c78bcb9 -->

Use `R17` for the recovered rootless endpoint and `H3 source` for the level-474
Kumar polarization. Low-q, E6, H2 and Q80/CM24 are comparison or regression
routes, not substitutes for the selected H3 marking. See
[`CONSTRUCTION_ROUTES.md`](CONSTRUCTION_ROUTES.md).

## Current priority

The **H3 source family** is the primary characteristic-zero starting point.
The currently selected **H3 degree-two corridor** is one certified route from
that source to the pinned recovered rootless/MW17 frame; the endpoint
identification and complete inverse NS transport are exact.  The route is not
proved shortest, globally optimal, or cheapest to compile.  The equation
route now reaches the rootless endpoint:

```text
H3 E7+E8/MW2
 --q6--> E8+E6/MW3
 --q8--> D13/MW4
 --q24 orbit85--> D12/MW5
 --q6 orbit42--> A11/MW6
 --q8 orbit12--> 2A5/MW7
 --q4 orbit208--> 3A3/MW8
 --q4 orbit1584--> D4+A3+3A1/MW7
 --q4 orbit164--> 2A3+2A1/MW9
 --q8 orbit376--> 4A1/MW13
 --q12 orbit5867--> rootless/MW17.
```

The q8 repair found two independent bugs in the previous child-side compiler: a binary-quartic 2-cover multiplier was applied twice, and the q-frame CRT normalizer omitted the `Dx` factor of `x(S)=Nx/Dx`. After both corrections the q8 Riemann--Roch problem collapses to an exact `13 -> 2` characteristic-zero intersection and its Jacobian has one `I9*` fibre plus nine `I1`, hence `D13/MW4` exactly.

The selected `D13/MW4 --q24 orbit85--> D12/MW5` equation arrow is now closed. Its exact q24 horizontal section is exported as `artifacts/local/elkies-k3/q8-q24-horizontal-section-qq.json` with status `PASS_EXACT_Q24_HORIZONTAL_SECTION`, and the resolved component-valuation RR lift is exported as `artifacts/local/elkies-k3/q24-d13-to-d12-component-valuation-qq.json` with status `PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR`. The RR profile is `ambient=56`, `collision=48`, `post=8`, `resolved=6`, `kernel=2`; the compiled quartic has degree four and its Jacobian has root rank `12`, root determinant `4`, Euler number `24`, and MW rank `5`.

At lattice/planning level, the candidate D13 detour is

```text
D13 --q4 orbit11--> D5+D9 --q4--> D13(new zero) --q24--> current D12.
```

The landing is the exact stored D12 full basis and the composite endpoint is
the pinned R17 full basis. Its historical planning score is 25,323 against
the calibrated direct score 27,885. Because the q4/orbit230 return exposed a
physical-chamber pseudo-zero, this changed-zero score is not treated as an
equation target until its own return equation identifies an effective zero.
A complete scan of all 1,119 declared-nef A11 q8
candidates, including all 24 torsion-root presentations via saturated adapted
frames, found no cheaper inherited-explicit operational-cost zero loop.

For a replay beginning one stage earlier, the corresponding exact lattice
detour is

```text
E8+E6 --q4 orbit11, zero=old_E8E6_component_1--> A2+D5+E7
      --q4--> E8+E6(new zero) --q4--> equation D13.
```

Its historical planning score is 3,961 versus 5,802 for the direct q8. The
D13 landing and pinned R17 endpoint
are identified by full integral unimodular bases, not ADE/MW labels. An exact
second-loop q4/q6/q8 scan finds no further improvement from the returned
E8+E6 marking. Widened q10 degree-two, q6/q9/q12 degree-three, and
q8/q12/q16 degree-four first-edge scans also have no winner. The machine certificate is
[`../artifacts/generated-results/elkies-k3-h3-first-q8-q4o11-promoted-route-certificate.json`](../artifacts/generated-results/elkies-k3-h3-first-q8-q4o11-promoted-route-certificate.json).
The changed-zero D13 landing was also searched directly toward current D12;
its degree-two q60 target and all 73 exact q4/q6/q8 crossover presentations
are operationally more expensive than the canonical D13 continuation. As
with the D13 detour, its changed-zero compiler score awaits a physical return
chamber audit.

The following formerly promoted lattice splice is no longer an equation-cost
target. It was

```text
2A5 --q4 orbit230--> A1+A4+A5 --q4--> 2A5(new zero)
     --q6 orbit1315--> 3A2+A3 --q4--> 2A5(second new zero) --q4--> 3A3,
```

using `old_A11_component_10` and then `old_A5A5_component_1` as explicit
zeros, then resume the existing suffix. Its marked fibres and determinant-one
transports to stored `3A3` and pinned R17 remain exact. However, the 4,199
score is withdrawn: the stored return zero has `P.O=26` and differs from the
effective P230 section by vertical roots; q6/orbit1315 meets effective P230 in
degree 54. Thus this is an exact lattice path, not a compiler-cost promotion.

For orbit42, use `mw=(-1,0,-1,-1,0)`, height `7`, correction `3`,
`P.O=3`, fibre twist `0`, and `D42=O+P+V`. The physical I8* marking and
complete twenty-point zero-pole shell are exact. Four exact section candidates
are also known. The resolved-RR artifact certifies the A11 child, and the
shell-degree marking binds it to equation-side orbit64/mapping7 in the chosen
C10 orientation.  The next q8 lift is also exact: its resolved RR profile is
`ambient=14`, `collision rank=12`, `h0=2`; its quartic Jacobian has
`2I6+12I1`, and old A11 component 9 is the exact child zero. See
[`ORBIT42_EQUATION_LIFT.md`](ORBIT42_EQUATION_LIFT.md)
and the version-locked [`scripts/success-path/`](scripts/success-path/) ledger.

Useful negative lessons: the fast q6 transport has degrees `435/703`; the
rational-halving candidates have no A11 chord; and the old correction-one,
`P.O=2`, easy zero-pole and archived Hensel paths are not the selected lift.
These close shortcuts, not every possible construction.

Marking warning: the pinned dominant D13 representative used by the lattice corridor and the component-nef D13 representative used by the q8 equation compiler are not the same stored 17-by-17 frame. Their exact full-NS bridge changes the embedded `U`; use the reverse-transport ledger rather than substituting one frame for the other by ADE/MW label.

The **Q80 Low-q Compiler Route** is now a completed secondary compiler route through its terminal characteristic-zero specialization shadow. Its generic lattice corridor reaches rootless/MW17, while the exact CM24 characteristic-zero terminal equation is `4A2+A3+A5/MW2`. The former final-q6 MW-marking problem is closed; see [`Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md`](Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md).

## H3 q8 exact certificate

Reproduce the corrected marking with

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_marking.sage
```

Expected geometry:

```text
MW coordinate = (-2,-2,0)
height = 24
O-intersection = 10
smooth collision degree = 10
II*, IV* = identity component
```

Then reproduce the complete second neighbour with

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_physical_root_target.sage \
  --representative component-nef \
  --output artifacts/local/elkies-k3/q8-target-component-nef-audit.json
cmp artifacts/local/elkies-k3/q8-target-component-nef-audit.json \
  elkies-k3/data/fibrations/h3_q8_component_nef_physical_root_target.json
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_corrected2cover_qq.sage
```

Expected endpoint:

```text
Q8QQRR|ambient=13|rows=11|rank=11|kernel=2|...
Q8QQQUARTIC|degree=4
Q8QQCHILD|finite=[(1,[2,3,15],'I9*'),(9,[0,0,1],'I1')]|infinity=((0,0,0),'smooth')|root_rank=13|root_euler=24|root_det=4|MW_rank=4|status=PASS_EXACT_CORRECTED_Q8_D13_CHILD
```

The detailed diagnosis and supersession boundary is in [`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md).

## Q80 status: generic route and exact CM24 characteristic-zero shadow

The generic Q80 alternate route is certified to a new rootless `MW17` frame:

```text
E6+D5+A3/MW3
 --q4--> D9+A4/MW4
 --q4--> D7+D5/MW5
 --q6--> D7+D4/MW6
 --q4--> A6+A4/MW7
 --q4--> A6+A3/MW8
 --q6--> A4+A2+A1/MW10
 --q4--> A3+A2/MW12
 --q4--> 4A1/MW13
 --q4--> A1/MW16
 --q6--> rootless/MW17.
```

Every retained new divisor from `D7+D5/MW5` onward has chamber-reduced old-fibre degree `2`.

The old suffix labels `7774`, `1938`, `6855`, and `candidate1` are search-result identifiers rather than intrinsic names. New notes should use the canonical stage names in [`CONSTRUCTION_ROUTES.md`](CONSTRUCTION_ROUTES.md), while [`Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md`](Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md) records the exact historical search lineage and the successful terminal-marking reconstruction.

The former characteristic-zero obstruction at the `A1/MW16` parent is closed. The selected generic terminal divisor is replayed with exact unimodular neighbour transports, the exact final horizontal is reconstructed as a difference of easier high-incidence MW sections, and the historical CM24 `P2-P3` point is used only to identify the correct difference modulo `73`.

The exact final resolved Riemann--Roch pencil over `QQ(sqrt(-3))` has

```text
ambient = 4
A4 quotient rank = 1
A5 quotient rank = 1
condition rank = 2
kernel = 2
h0(D) = 2.
```

Its exact degree-four binary quartic compiles to a characteristic-zero CM24 child with

```text
finite fibres = 4 I3 + I4 + I6 + 2 I1
infinity = smooth
root lattice = 4A2+A3+A5
root rank = 16
root determinant = 1944
MW rank = 2
Euler number = 24.
```

The late specialized sequence is therefore exact in characteristic zero:

```text
2A6+3A1/MW3
 --q6--> A5+2A4+2A1/MW3
 --q4--> 2A4+2A3+A1/MW3
 --q4--> A1+2A3+2D4/MW3
 --q4--> A1+A2+A3+A4+A5/MW3
 --q6--> 4A2+A3+A5/MW2.
```

This does **not** identify the specialized child with the generic endpoint. The generic terminal child remains the independently certified rootless/MW17 frame; CM24 has higher Picard rank and acquires extra algebraic classes.

Reproduce the terminal reconstruction with

```bash
sage elkies-k3/scripts/trace_q80_candidate1_marked_transport.sage
sage elkies-k3/scripts/recover_q80_final_q6_via_basis_sections.sage
sage elkies-k3/scripts/certify_q80_final_q6_char0_rr_from_basis.sage
sage elkies-k3/scripts/compile_q80_final_q6_char0_child.sage
```

The final pinned certificate/model are

```text
data/fibrations/q80-final-q6-char0/Q80_CHAR0_FINAL_Q6_CERTIFICATE.md
data/fibrations/q80-final-q6-char0/q80_char0_final_q6_child.sage
```

The historical GF73 final certificate remains useful as a regression:

```text
data/fibrations/kumar_q80_final_q6_cm24_equation_gf73.txt
```

and the machine-readable stage ledger is

```text
data/fibrations/kumar_q80_cm24_equation_progress.tsv.
```

## Reusable equation/compiler results

The combined H3/Q80 work now gives these durable rules:

1. **Binary-quartic covariant maps are 2-covering maps.** Differences of transported quartic points must not be interpreted as primitive MW differences without the factor-of-two check.
2. **Clear the full rational expression before CRT normalization.** For the corrected H3 q8 frame, the `Nx` residue is `Ny*Dx/(h*Dy)`; omitting `Dx` leaves a hidden vertical pole.
3. **Specialize the actual divisor before searching equations.** Generic `P.O`, MW height, twist, and vertical support can change dramatically at CM24.
4. **Connected ADE corrections compile as resolved quotient-line/module conditions**, not one independent row per listed exceptional component.
5. `root_component_data()` may return arbitrary integral root-lattice bases; use discriminant groups rather than assuming Cartan coordinates.
6. For an A3 correction `(-2,-1,-1)`, the exact local module is the middle-component double-vanishing condition.
7. For a D4 correction `(-1,0,-1,-1)`, the ramified-chart outer-complement condition is the deterministic quotient residue `c=0`.
8. In the final A5 correction `(-1,0,-1,-1,0)`, the exact quotient line is the `+/-4` residue pair for the two horizontal signs.
9. **A singular modular marking need not be lifted directly.** Reconstruct easier exact sections first and identify the required section by exact group-law difference plus modular regression.
10. **Keep generic and specialized endpoints separate.** A characteristic-zero equation on CM24 certifies the specialization shadow; it does not replace the generic rootless lattice certificate.
11. **Optimize equation paths separately from lattice paths.** The first rank-growing or smallest-q neighbour need not minimize pole order, coefficient growth, RR dimension, or marking complexity.

The field-generic exact module compatibility layer remains in

```text
scripts/elliptic_neighbor_compiler_field_generic.sage
scripts/verify_elliptic_neighbor_compiler_field_generic.sage.
```

## Start here

- [`ELKIES_K3_PROCESS_ATLAS.md`](ELKIES_K3_PROCESS_ATLAS.md) — short chronology, rank flow, cancellations, and reusable lessons.
- [`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md) — proved deductions, conditional lift theorem, and open navigation conjectures.
- [`AGENTS.md`](AGENTS.md) — compact handoff rules, hints, and known dead ends.
- [`CONSTRUCTION_ROUTES.md`](CONSTRUCTION_ROUTES.md) — human-readable geometric route map and naming scheme.
- [`scripts/README.md`](scripts/README.md) — current proof/compiler/search entry points and archive policy.
- [`SCRIPT_ROUTE_AND_FAILURE_LEDGER.md`](SCRIPT_ROUTE_AND_FAILURE_LEDGER.md) — exact path-selection reasons, failed lanes, supersessions, and remaining gaps.
- [`Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md`](Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md) — completed Q80 terminal characteristic-zero reconstruction and supersession boundary.
- [`Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md`](Q80_LOWQ_SUFFIX_PROVENANCE_2026-08-23.md) — exact Q80 suffix search lineage and successful marking recovery.
- [`Q80_CM24_EQUATION_LEDGER_2026-08-22.md`](Q80_CM24_EQUATION_LEDGER_2026-08-22.md) — Q80 CM24 equation ledger, updated with the characteristic-zero closeout.
- [`H3_Q8_CURRENT_FRONTIER.md`](H3_Q8_CURRENT_FRONTIER.md) — concise H3 q8 exact frontier.
- [`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md) — two-bug diagnosis and exact H3 q8 repair.
- [`RESEARCH_UPDATE_2026-08-22.md`](RESEARCH_UPDATE_2026-08-22.md) — historical 2026-08-22 repository-wide checkpoint; later Q80 status is superseded by the closeout note above.
- [`KUMAR_E7E8_BACKTRACK.md`](KUMAR_E7E8_BACKTRACK.md) — source recovery, H2/H3 distinction, and exact H3 lattice corridor.
- [`data/fibrations/kumar_q80_cm24_equation_progress.tsv`](data/fibrations/kumar_q80_cm24_equation_progress.tsv) — machine-readable Q80 CM24 stage summary.
- [`Q80_TO_ROOTLESS_PATH_2026-08-21.md`](Q80_TO_ROOTLESS_PATH_2026-08-21.md) — Q80 rootless lattice certificate and historical detailed ledger.

## Next strategic gate

For **H3**, the active lifting target is the physical degree-two
`q8/orbit376` pencil from the exact C8-pointed `2A3+2A1/MW9` equation to
`4A1/MW13`.  The q4/orbit1584 predecessor has an exact QQ Riemann--Roch plane,
quartic, `I4+3I2+8I1+I0*` minimal Jacobian, and second-I6-affine pointing.
The q4/orbit164 edge has an exact `4 -> 2 -> 2` plane, quartic,
`I4+2I2+12I1` finite fibres plus `I4` at infinity, and an exact C8 pointing.
The inherited-`P1` construction now reaches a complete modular Abel trace on
this model: the three pointed q4 maps have degrees `3 -> 6 -> 7`, and at
`p=131` the fibrewise `7 x 8` kernels interpolate an exact section with
`x=(32,28)` and `y=(47,42)`.  Its difference from the marked q8 horizontal
has residual MW tail `(0,0,0,0,-1,-1,2,0,0)`, so the immediate gate is the
last local marking ambiguity in the already exact equation-side rank-eight
section subgroup.  Fourfold pole growth gives its corrected height Gram with
determinant `459/8`.  Testing all eight component-compatible embeddings
against the full modular trace leaves one section with the certified
`P.O=4` profile:
`H=T-C8opp-B0+2B1+B2-3B3-B4-2B5+B7`, with compact degrees
`x=(12,8)`, `y=(18,12)`.  Direct CRT of this smaller compact section over 22
good primes gives a 566-bit modulus; simultaneous projective reconstruction
recovers exact `QQ(t)` coordinates with primitive-vector maxima 363 and 526
bits.  Exact substitution and reduction to every input prime pass.  Thus the
horizontal is proved in characteristic zero without lifting the larger Abel
trace; exact fourfold pole growth also gives its marked canonical height 11.
The resolved q8 calculation is now exact: the chord ambient has dimension 12,
smooth saturation leaves dimension 4, and the two connected old-`I4` quotient
rows leave `h0=2`. The resulting quartic has a minimal Jacobian with
`4I2+16I1`, Euler number 24, and `4A1/MW13`. The selected marked embedding and
an exact `B0` tangent orient the finite `T=0` I4 chain, identifying the first
pointed origin as old component 6. Directly pointing instead at `P1229`, the
nonidentity component of the finite I2 at `T=25281/168246841`, selects quartic
sign `-1` and reproduces the same child invariants exactly. Thus the q8 edge is
equation-marked with the preferred `P1229` zero.
The marked q8/orbit376 and q12/orbit5867 lattice edges are primitive, nef,
old-fibre-degree two, unimodular in both directions, and end at pinned R17.
The q12/orbit5867 nominal four-section lattice word has four physical
`P.O=0` branches, parent degrees `(3,2,1,2)`, and parent
`a-b=(2,2,1,1)`, improving on the q12/orbit4484 fallback as a planning proxy.
It is the preferred optional final edge. The stable certificate label is
`q12o5867_after_q8o376` with exact fibre fingerprint `d676cab5...`; the
expanded-frontier index 30357 is sample-local. Complete polynomial shells at
four good primes contain no realization of its nominal Q1 branch, so the
equation compiler instead uses classes `499+500+69+511` with correction
`-489+933-913`; its actual parent degrees are `(4,2,1,5)` and `a-b` values
`(4,2,1,4)`. All sections, the q12 horizontal, the exact `22 -> 2` resolved
plane, quartic, and rootless `24I1` Jacobian are complete over QQ. A complete
mod-131 polynomial shell then selects a regular 17-section basis; all seventeen
lift exactly over QQ, and exact intersections reproduce the pinned
determinant-948 height Gram. Thus rank at least 17 and trivial torsion are
unconditional. Switch to the Picard/discriminant and source-identity gates in
[`PRIOR_WORK_SHORTCUT_AND_ENDPOINT_CERTIFICATION_2026-08-25.md`](PRIOR_WORK_SHORTCUT_AND_ENDPOINT_CERTIFICATION_2026-08-25.md).

The former 4,199, 10,334, and 13,518 targets remain withdrawn: respectively
they used a pseudo-zero, a non-section component-10 zero, and a non-nef q104
representative.  The q10 certificates remain exact provenance, but q10 is no
longer the lifting target.

For **Q80**, the terminal marking/RR/equation gate is closed. Treat the direct two-parameter final-section resultants, digit-by-digit `73`-adic lift, and local-73 singularity probes as diagnostics only. Any new Q80 work should start from the exact certificates above rather than reopening the final q6 reconstruction.
