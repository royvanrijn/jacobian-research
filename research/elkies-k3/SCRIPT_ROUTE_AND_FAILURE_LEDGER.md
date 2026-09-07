# Elkies K3 script route and failure ledger

Status: **research-history and claim-boundary index**, updated through
2026-09-04.  Mathematical status remains exclusively in `../MATH_STATUS.json`.

This note records why the script tree contains several apparently competing routes. It
separates four questions that were repeatedly conflated during the reconstruction:

1. Which K3/source polarization is the correct one?
2. Which lattice path proves that the source and rootless `MW17` fibrations lie on the
   same surface?
3. Which path is easiest to execute with exact characteristic-zero equations?
4. Which specialization/compiler experiments are useful without being the generic
   construction?

The answers are not the same:

- **H3 is the recovered source family.**
- **The preferred physical H3 equation route is complete through q8/orbit376
  and q12/orbit5867 to the rootless `24I1/MW17` endpoint.**
- **The endpoint is source-identified, has geometric Picard rank 19, and has
  full saturated determinant-948 Mordell--Weil lattice R17.**
- **The older native D12 -> A1 suffix remains a lattice/root-Weyl and reverse-lift
  reference route; it is not the current arithmetic path.**
- **The selected H3-to-MW17 chain has not been proved shortest, globally optimal, or
  cheapest to compile.**
- **Q80 is an independent comparison/compiler route, not the source family.**
- **q32 is retained as an alternate/regression route; it is no longer preferred over
  the viable q24-native H3 path.**

The later direct norm12/orbit11952 hop changes the operational conclusion
without changing that history: the alternate-Q80 rootless class now has an
exact `QQ` equation, `24 I1`, seventeen saturated rational sections, and
arithmetic generic rank 17.  The degree-11511 Q80 transport and third-q12
reconstruction remain provenance.  The determinant-948 noncyclic
`4A1/MW13` bridge is also equation-explicit in both directions.  The active
foundry milestone has therefore moved to an arithmetic MW17 fibration on a
different Neron--Severi lattice.

The executable index is [`scripts/README.md`](scripts/README.md); archived experiments
are classified in [`scripts/archive/README.md`](scripts/archive/README.md). The current
route priority and arithmetic handoff are in [`README.md`](README.md); the older
[`H3_PREFERRED_PATH.md`](H3_PREFERRED_PATH.md) is retained as a dated lattice-route
decision.

## 1. What we started with

The initial hard datum was the determinant-948 rootless rank-17 Mordell--Weil lattice.
Directly solving for seventeen sections produced systems that were too large and too
poorly marked to be a sensible first reconstruction. The first strategic change was to
search the same Neron--Severi lattice for primitive isotropic classes and hence alternate
elliptic fibrations with fewer MW generators and more reducible fibres.

This reverse search produced valid lattice transports, but not automatically the
historical/source construction.

### Low-q MW2 Backtrack

```text
rootless/MW17
 --q25--> A3+7A1/MW7
 --q4 --> D4+A3+2A2+2A1/MW4
 --q4 --> A5+D4+2A2+A1/MW3
 --q4 --> E6+D4+2A2+A1/MW2.
```

Reason selected: it rapidly reduced the number of free sections and supplied exact
inverse divisor classes. Result: a strong proof that the rootless fibration can be
rewritten as a low-MW fibration. Limitation: the endpoint is not the Kumar source
polarization identified by the primary construction.

### E6 Backtrack

```text
rootless/MW17
 --q90--> MW7
 --q4 --> MW4
 --q4 --> E6+2A3+2A1/MW3.
```

Reason selected: an `E6/MW3` frame looked small enough for explicit section equations.
Result: the integral NS transport is valid. Limitation: the split Weierstrass chart later
attacked was inferred from fibre type rather than obtained by geometrically executing
this chain.

Authoritative lattice replay: `scripts/verify_e6_neighbor_chain.sage`.

## 2. Recovering the actual beginning

The source audit restored the construction order from the quaternionic moduli problem:

```text
principally polarized QM abelian surface
  -> Dolgachev--Kumar K3
  -> E7+E8 fibration with MW rank 2
  -> elliptic neighbours
  -> rootless MW17 fibration.
```

Exact height-form and discriminant-glue classification left three compatible Kumar
frames:

```text
H1 = [ 5/2   1 ]    H2 = [ 4       0 ]    H3 = [21/2   3]
     [  1  190 ]         [ 0   237/2]         [   3   46].
```

H2 was initially attractive because its extra involution matches `w2=w237`, and its
height directions split cleanly. The exact pullback of `H21` to the `H92` chart at CM24,
however, identifies the level-474 genus-two component with **H3**, on
`H21 intersect H92` (also lying on `H101`). The source family is therefore

```text
y^2 = -27*x^6 + 198*x^4 - 171*x^2 + 576
```

on the normalized base, with an exact `E7+E8/MW2` K3 equation and marked section
descent.

Current source replays include:

```text
classify_kumar_e7e8_anchors.sage
deconstruct_x0679_quotients.sage
verify_h21_h92_level474_branch.sage
reconstruct_h21_h92_level474_qq.sage
export_h3_level474_source_family.sage
verify_h3_noncm_q6_source_anchor.sage
verify_h21_q6_section_descent.sage
lift_h21_p1_modular.sage
verify_h92_section_descent.sage
```

This establishes **where construction begins**. It does not choose the globally easiest
route from that source to the rootless fibration.

## 3. The selected H3 degree-two corridor

The selected lattice/chamber chain is

```text
H3 E7+E8/MW2
 --q6 --> E8+E6/MW3
 --q8 --> D13/MW4
 --q24--> D12/MW5
 --q6 --> A11/MW6
 --q8 --> 2A5/MW7
 --q4 --> 3A3/MW8
 --q4 --> A3+2A2/MW10
 --q4 --> 5A1/MW12
 --q4 --> 4A1/MW13
 --q4 --> 3A1/MW14
 --q4 --> 2A1/MW15
 --q4 --> A1/MW16
 --q6 --> rootless/MW17.
```

### Why this path was selected

The discovery objective was approximately:

```text
small q
+ old-fibre degree two
+ nefness that can be proved exactly
+ immediate MW-rank growth
+ manageable root classification.
```

Those are excellent lattice-search criteria. They are not necessarily the right cost
function for symbolic equation compilation, where pole order, coefficient growth,
marking transport, local quotient complexity, and the availability of easy sections can
dominate.

### What is exhaustive and what is not

- The first H3 q6 shell was exhaustively enumerated under its stated marking and proper
  factor presentation. It gives the preferred `E8+E6/MW3` child.
- From that child, early q4/q6 second-step searches were bounded experiments, not global
  obstructions. The q8 orbit analysis then selected the exact `D13/MW4` child.
- From D13, the stated proper presentations through `q=20` were exhaustively checked for
  rank growth. q21 and q22 also have maximum MW rank four, and q23 has no proper factor
  presentation. The first rank growth is q24.
- The q24 shell has **three** primitive orbits with `D12/MW5` root data. Orbit 85 was
  preferred because of its root/MW projection and chamber behaviour; that does not prove
  it is the easiest possible equation child.
- The downstream preferred suffix is no longer justified merely by matching ADE labels
  to an older/historical frame. Commit `bfb96bf526b0b289e12d531e7aff697cae5f39c6`
  exports the equation-side q24 D12 frame as `D24eq` and replays the suffix directly from
  that frame.

Consequently the chain proves existence, integral transport, and nef degree-two pencils.
It does not rule out another route that is symbolically cheaper, but route selection is
no longer blocked by uncertainty over which D12 marking to continue from.

### q24-native D12 -> A1 replay

The latest native replay uses

```text
scripts/export_h92_q24_native_d12_frame.sage
scripts/run_h92_q24_native_suffix_to_a1.py
```

and starts from `D24eq`, explicitly described by the exporter as the `D12/MW5` frame
obtained directly from q24/D13 equation-side geometry, **not** from q32 and **not** from a
later historical D13/D12 frame.

The replay uses the suffix q-sequence

```text
6, 8, 4, 4, 4, 4, 4, 4, 4
```

and reaches

```text
A11/MW6
2A5/MW7
3A3/MW8
A3+2A2/MW10
5A1/MW12
4A1/MW13
3A1/MW14
2A1/MW15
A1/MW16.
```

The terminal q4 step reports root data `(1,2,2)`, MW rank 16, and selected orbit `17593`.
The runner's terminal status is

```text
PASS_Q24_NATIVE_D12_TO_A1
```

with local artifact

```text
artifacts/local/elkies-k3/q24-native-suffix/q24-native-d12-to-a1.json
```

**Claim boundary:** this is a pure lattice/root-Weyl replay downstream of D12. It pins
the preferred sequence and the native marking from which equation compilation should
continue; it does not claim that the nine downstream characteristic-zero equation
neighbors have already been executed.

## 4. Equation progress on H3

### H3 q6: source to `E8+E6/MW3`

The first arrow is exact at equation level. The q6 shell, marking, section descent,
chamber reduction, and child fibre data are all pinned.

### H3 q8: `E8+E6/MW3` to `D13/MW4`

The q8 equation route is exact after two independent compiler repairs:

1. **2-cover multiplier error.** The binary-quartic covariant map is a 2-covering map.
   The old calculation interpreted mapped point differences as primitive MW differences
   and effectively doubled the horizontal twice. The corrected section is
   `S=Pmap+Qmap`, with height 24 and `S.O=10`; the old height-96 `2S` marking is
   withdrawn.
2. **Missing denominator in CRT normalization.** For `x(S)=Nx/Dx`, the q-frame
   congruence must clear the complete rational expression. The corrected relation is
   `R*h*Dy == Ny*Dx (mod Nx)`. Omitting `Dx` left a hidden vertical pole and produced
   oversized false intersections.

After both fixes the exact RR problem is

```text
ambient dimension = 13
condition rank     = 11
kernel dimension   = 2
quartic degree     = 4
child              = D13/MW4.
```

Authoritative scripts:

```text
derive_h92_q6_child_q8_marking.sage
derive_h92_q6_child_q8_physical_root_target.sage
derive_h92_q6_child_q8_corrected2cover_qq.sage
certify_h92_q6_child_jacobian.sage
```

**Replay-packaging boundary (2026-09-04).**  The final corrected-q8 checker is
not currently standalone from a clean checkout: its q6 child-Jacobian,
transported-zero, and E7-infinity JSON inputs are ignored intermediates and are
absent.  Rebuild the complete `EC-K3-H3-Q6` chain before replaying it.  Do not
interpret a missing intermediate as a failed equation identity, and do not
advertise the one-line q8 command as portable until a compact sufficient input
projection is tracked.

### H3 q24: `D13/MW4` to native `D12/MW5`

The lattice/chamber arrow is exact and the preferred child is orbit 85.  The
equation arrow was the nearest characteristic-zero frontier at the time of the
calculations summarized below; it was subsequently closed as part of the
exact H3 route to q12/orbit5867.  The retained scripts remain valuable because
they record how the high-pole horizontal was recovered and how the physical
marking was attached.

`scripts/close_h92_q8_q24_by_q6_translation.sage` closes the exact marking/NS bridge.
`scripts/recover_h92_q24_exact_by_qq_trace_interpolation.sage` transports the selected
horizontal through exact specialization and interpolation. Its working construction uses

```text
W = Qmap - S3,
q24 horizontal = AJ_II*(W) + 2*G1,
```

with the high pole profile handled through trace interpolation rather than a direct giant
symbolic group-law expression. The script enforces denominator-square and RHS-square
structure and an independent finite-field check.

The new native suffix result changes the **route decision**, not this proof boundary:
once the q24 D12 equation is fully certified, continue from its `D24eq` marking along the
native q6/q8/q4 suffix rather than switching to q32 or identifying with a historical D12
frame.

Those formerly missing gates are supplied by the later exact endpoint chain
and its canonical certificates.  Do not restart this interpolation route
merely because the historical working description records its earlier proof
boundary.

## 5. Q80: why it exists and what it proves

H2 is a valid symmetry comparison polarization, although not the recovered H3 source.
A bounded stability search found its q80 neighbour

```text
H2 E7+E8/MW2 --q80--> Q80 E6+D5+A3/MW3.
```

At CM24 this chart retains all three generic MW directions and acquires only one extra
root, making it an excellent compiler laboratory.

The generic lattice route is

```text
Q80 E6+D5+A3/MW3
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

The generic endpoint is certified at lattice level. The exact characteristic-zero CM24
shadow ends instead at

```text
4A2+A3+A5/MW2,
```

because CM24 has larger Picard rank. The two endpoints must never be identified.

The final Q80 q6 route originally became stuck on a singular direct section lift. It was
closed by reconstructing easier `P.O=0` sections, identifying the required pair modulo
73, and taking their exact group-law difference. This is a durable compiler lesson:
preserve easy markings or reconstruct them; do not insist on lifting the hardest marked
section directly.

## 6. Failure and supersession ledger

| programme/assumption | what happened | evidence retained | current replacement |
|---|---|---|---|
| Directly solve the rootless MW17 equation from 17 sections | Systems remained too large even after motif compression. | Early section-system and lattice-search scripts. | Backtrack to simpler fibrations, then reconstruct forward. |
| Guessed split E6 chart represents the transported E6 backtrack | The chart had the right-looking fibre type but no geometric transport certificate. Exact glue added the missing off-diagonal height gate; all 69 old survivors from 406,655,040 tested cores fail it. | E6/MW3 root launchers and archived component/state scans. | H3 source family and exact q6/q8 transported equations. |
| All three `A2` fibres can be additive `IV` | The resulting family has constant `j=0`; order-three CM forces even geometric MW rank, while the target frame needs rank three. | Historical `E8+A2^3` scripts and correction note. | Mixed `II*+2IV+I3+3I1` model as a comparison; H3 is the primary route. |
| H2 is the source because of its clean involution | H2 correctly explains `w2=w237`, but the exact Humbert branch for level 474 is H3 on `H21 cap H92`. | H2/Q60/Q80 scripts and Kumar comparison notes. | H3 source family; retain H2/Q80 for comparison/compiler work. |
| Old q8 degree-46 marking | It carried an unnecessarily large section representation and obscured the true small RR intersection. | q8 audit/trace scripts in `archive/`. | Correct q8 2-cover marking and exact `13 -> 2` RR system. |
| `true1600` source-side q8 route | Large cap/jet systems attacked an inflated presentation. | `probe_h92_q8_true1600_*` archive cluster. | Correct child-side q-frame route. |
| `corrected1278` hand-built q8 route | Local conditions were useful diagnostics but still solved the wrong enlarged module problem. | `probe_h92_q8_corrected1278_*` archive cluster. | Field-generic resolved quotient modules and exact corrected q8 compiler. |
| Binary-quartic map returns primitive MW points | False: it is a 2-covering map; differences were doubled. | q8 marking-height audits. | Explicit factor-of-two check in the authoritative q8 scripts. |
| CRT normalize modulo `Nx` without `Dx` | Left a hidden vertical pole and fake degree growth. | old q-frame probes. | Clear the whole rational expression before normalization. |
| Direct singular Hensel lift of the hard q6/q8/Q80 section | Non-transverse coordinates caused unstable or very slow lifting. | q6 third-section and Q80 local-73 diagnostics. | Exact section descent, trace interpolation, or differences of easier exact sections. |
| CM24 endpoint can stand in for the generic Q80 endpoint | Extra CM algebraic classes change the fibre/root data. | exact CM24 equation ledger. | Keep generic rootless lattice certificate and CM24 equation shadow as separate claims. |
| First rank-growing neighbour is automatically best | Lattice search optimized immediate rank growth and degree, not equation complexity. | D13 q24 and later first-hit scripts. | Keep alternative searches available, but use the now-native q24 suffix as the default H3 route while it remains viable. |
| q32 D12 should become the preferred H3 continuation because its modular compiler is convenient | q32 gave useful modular quartic/spinor-marking diagnostics, but it introduces a different D12 route and does not improve the source provenance. The q24-native `D24eq` frame now replays exactly to A1. | q32 signature, spinor-quartic and pointed-section experiments. | Demote q32 to alternate/regression work; compile the native q24 `D24eq` suffix. |
| A historical D12 frame with the same ADE/MW label can stand in for the q24 equation-side D12 | Equal root data does not identify the marked elliptic fibration. | historical suffix searches and frame comparisons. | `export_h92_q24_native_d12_frame.sage` + `run_h92_q24_native_suffix_to_a1.py`. |
| The alternate rootless frame is only lattice-accessible | Minimum-incidence `U` search found norm12/orbit11952 directly on the compact published R17 equation. | Historical Q80 transport and third-q12 reconstruction remain as regression/provenance. | Use the direct degree-two compiler, exact alternate equation, and saturated rank-17 basis. |
| Published R17 controls automatically calibrate alternate Q80 | Exact alternate-j preimage polynomials for all four rank-25--28 controls have no rational roots. | The four misses are retained as negative controls. | Use the complete norm-twelve public-curve atlas; curve 12 and the `08f72`/`11952` fibres are native alternate controls. |
| Target-fitted visibility predicts unseen exceptional directions | Low-genus pencils can be forced through known quotient directions after the targets are supplied. | Exact fitted quartics and declared bounded splitting searches. | Separate rigid, predeclared, and post-hoc visibility; calibrate prospective predictors on held-out controls. |
| Incremental CRT height makes a finite beam a completeness sieve | Shortest rational-representative height is nonmonotone: an exact four-prime width-one run keeps a height-1409 endpoint while a pruned branch completes to `48/53`, of height 53. | The exact EC regression and bounded q12/orbit5867 beam outputs remain valid candidate generators. | Label beam output heuristic, retain omitted states as untested, and require a full Cartesian product or proved lossless reduction for exhaustion. |
| Maximum bridge minimum is a safe or useful hard sieve | The retrospective fixed-core rule retained only four of five rootless classes. On an untouched determinant-78 shell all 277 candidates received the same score, and computing it was slower than direct root classification. | The benchmark, held-out negative control, and fixed-core theta-convolution certificates remain exact. | Use bridge minimum only as descriptive data. Reject only with a proved monotone root-budget bound; use exact theta convolution only after the complete fixed-core bridge universe is known. |
| A geometric NS0024 MW17 route can close the arithmetic `QQ` milestone | A full rational NS0024 marking would force a forbidden rational non-CM degree-475 Fricke-quotient point. | Exact geometric frames, modular source work, and physical corridors remain valid over algebraic closures or larger fields. | Park arithmetic NS0024 over `QQ`; require the Galois/rational-marking gate before equation compilation. |
| A smooth modular point or long Hensel lift supplies a rational source | NS0031 model 157 has a unit Jacobian minor and hence a formal `ZZ_7` branch, but its split-Clifford curve has no noncuspidal rational point with the required marking. | The GF(7), finite-lift, formal-smoothness, and five-edge corridor certificates are retained as local/geometric controls. | Stop the NS0031 `QQ` lift campaign; run the rational-marking arithmetic screen before equation work on the next NS. |
| The cheapest native bisection prefix can be used as an exhaustive sieve | The 1,024 native alternate-Q80 covers were selected by cost and do not represent all native classes. | Exact priority keys and the full inherited/smooth-character domains. | Keep ranked prefixes and complete enumerations separately typed; never discard outside-prefix candidates mathematically. |
| Norm-eight inversion exhausts every zero-Tate-class carrier | Exact quotient analysis exposed 49 norm-twelve residual traces in addition to the 63,917 norm-eight minima. | The original inversion remains a valid first layer. | Enumerate the full minimum-norm quotient spectrum and close all 833 residual trace/target cases before claiming zero-class exclusion. |
| Successful chunks or aggregate counts imply complete shard coverage | Independent chunks can omit intervals, exceptional traces, malformed outputs, or an entire named ambient while preserving plausible totals. | All raw singular-search chunks, prescribed-root shards, and exact exceptions are retained. | Require half-open ranges or an exact expected key set, current source hashes, gap/overlap checks, fail-closed parsing, and one merged coverage certificate. |
| PARI `ellrank`'s third field is the full `Sha[2]` dimension | In the four-field convention `[r1,r2,s,L]`, `s` is the even Cassels-pairing quotient rank `dim(Sha[2]/2Sha[4])`; conflating it with all of `Sha[2]` obscures the distinction between Selmer, Mordell--Weil, and Tate--Shafarevich contributions. | The completed-backend formula `dim Sel_2(E)=r2+dim E(Q)[2]+s` remains valid. | Record the field as `pari_cassels_pairing_quotient_rank`, validate the interval and parity, and derive the total dimension through the tested fail-closed helper before applying the residual gate. |
| A URL plus a historical hash is a durable public-database replay | The ICARM database grew from 474 to 556 records, so the original `wgxli` checker could no longer retrieve the byte-identical September 1 snapshot even though its target records remained stable. | The original fingerprint and later exact lineage atlas remain valid. | Preserve a claim-sufficient exact local projection.  The repaired checker joins the committed ids-1--474 equation projection to all thirteen `wgxli` point records and reproduces the original artifact byte-for-byte offline; live retrieval is an optional drift audit only. |

## 7. What was deliberately kept

The audit retained unique failed scripts for three reasons:

1. They record exact bounded negative results that should not be recomputed or forgotten.
2. They explain compiler rules that otherwise look arbitrary.
3. They can serve as regression fixtures when a new implementation claims to solve the
   same local problem.

Only five byte-identical ` (1)` archive copies were removed. No unique historical file
was deleted, and no broad root-to-archive move was made because existing notes and shell
commands still refer to several historical root scripts.

The q32 experiments should likewise remain available: the spinor-point/pointed-quartic
work is a useful general marking-recovery technique even though q32 is no longer the
preferred H3 continuation.

## 8. What should be done next

### Primary arithmetic work

The determinant-948 equation routes are complete controls.  Use the compact
published chart and the native alternate-Q80 fibres for specialization, but
keep the proof and search-budget gates separate.  A certified actual quotient
`Sel_2(E_t)/<P1,...,P17>` of dimension below 15 rejects a rank-32 candidate.
A timeout, incomplete class ledger, or norm-one signature is scheduling
information only: it neither passes nor rejects the proof gate and does not
block a separately finite, checkpointed point search.  Fifteen exactly
certified quotient directions prove rank at least 32 without a Selmer upper
bound.

The first new construction milestone is a planner-found arithmetic MW17
fibration on a different NS.  Arithmetic NS0024 over `QQ` is closed
negatively, and the split-Clifford/`X_0(37)` obstruction now closes NS0031
negatively as well. Rerank the remaining frames through the rational-marking
arithmetic gate. Determinant `720` is the strongest lattice/corridor control,
but its known rational `3A5` point saturates to determinant `20` and is not a
valid source.

The alternate-Q80 product laboratory has closed only the zero Tate class and
the complete rational intersection-one deep-trace screen.  Nonzero quotient
classes, section existence, and product-twist ranks remain open.  Smooth and
singular arithmetic-genus-one character collisions are closed on the two
direct norm-twelve charts; the recorded genus-two miss is bounded.

### Parallel route optimization

Alternative routes remain scientifically useful, but they are secondary to
the different-NS arithmetic source gate, exact residual descent, and the
surviving quotient-class problems. Reopen them only when they supply a needed
source identity, specialization certificate, regression control, or reusable
compiler technique.

A route-cost model should still include:

```text
old-fibre degree
horizontal pole order and coefficient height
RR ambient dimension
resolved-component quotient complexity
field-extension degree
marking transport size
availability of easy exact sections
expected quartic/Jacobian complexity.
```

### Documentation discipline

For every new attack, record in the script header and in this ledger:

- exact source frame and marking;
- whether enumeration is exhaustive or bounded;
- selection cost/reason;
- exact PASS claim;
- exact failure or timeout boundary;
- replacement script/note when superseded.

That prevents a useful negative result from becoming an unexplained filename and prevents
a bounded search from being retold later as an obstruction theorem.
