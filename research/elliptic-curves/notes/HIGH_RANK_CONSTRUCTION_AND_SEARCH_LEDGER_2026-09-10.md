# High-rank construction, evidence and search ledger

Consolidated 2026-09-10; source review through `b883277c00c6bfe8d8193d9d5ada809010d91b35`.
This is a historical synthesis and experiment ledger, **not mathematical-status authority**.
`research/MATH_STATUS.json` and the linked canonical certificates retain that role.
The initial version at `17e8d3e` is preserved in Git history. This revision corrects
scope errors, distinguishes user-reported local work from inspected commits, and
replaces the proposed shallow broad run with an executable successor.

**Current operational handoff:** [broad search and launch instructions](BROAD_RANK_SEARCH_2026-09-10.md).
No search was launched by this documentation/implementation update. The user reports
that earlier experiments are stopped; historical notes saying “running” are dated
snapshots, not a live process check.

## 1. Executive conclusion

We did **not discover ICARM curve302**. It entered as a public curve with 31 points.
We independently certified its rank lower bound, reconstructed a different explicit
K3 parent whose generic rank is exactly 17, identified the primitive generic core,
and developed a calibrated search that recovers all fourteen displayed extra
directions. We still do not know the original discoverers' unpublished construction
or a prospective arithmetic condition forcing another comparable specialization.

The most useful operational evidence is not the newest parent by itself. The
[sixty-fibre X948 R17 panel](R17_SIXTY_SEED_COMPLEMENT_PANEL_2026-09-09.md) found 48
seeds, 38 further cascades, one rank-at-least-27 subgroup and three at least26.
The [subsequent foundry publication](FOUNDRY_CURVE_LEDGER_2026-09-09.md) added 30
curves at bounds22–26. Those are prospective discoveries relative to pinned
exclusion snapshots. By contrast, curve302's spectacular 17-to-31 recovery is
retrospectively calibrated, and the new X1092 class1's reported shallow panel has
no gain. These are different kinds of evidence, not competing estimates of one rate.

**Recommendation:** reuse the productive six-presentation R17 seed/amplification
pipeline for most of a wide, fixed campaign; give the two realized X1092 MW17
presentations a smaller comparison arm. Keep arithmetic-class research optional.
A new theorem, complete class group, exact conductor or unconditional rank upper
bound is not required to search for and certify additional rational points.

## 2. Curve302: provenance, rank and what “works” means

### External input, not our new record

The [rank31 note](ICARM_CURVE302_RANK31.md) records the ICARM source, its attribution,
the public equation and all31 points. It records an August23 pointer by David
Renshaw. The public conditional BSD+GRH exact-rank statement was not reproduced
by this repository; our result is **unconditional rank at least31**. Original
family/parameter/section provenance remains unknown. No outreach is authorized.

The exact checker verifies point membership, a rational short-model transport,
trivial torsion and a full-rank matrix in products of finite groups modulo2.
A separate quadratic-character implementation corroborates independence. It also
checks the minimal model, exact local reduction and conductor, with a later
primality-proof supplement replacing mere probable-prime checks.

For a hypothetical integer relation, the injective finite signatures force all
coefficients even. Absence of rational2-torsion permits repeated division of the
relation, proving independence. This argument does **not** supply a rank upper
bound. Numerical canonical heights schedule search; they are not the rank proof.

### The successful inverse reconstruction

[Construction investigation](ICARM_CURVE302_CONSTRUCTION_INVESTIGATION.md),
[historical recovery](CURVE302_CONSTRUCTION_RECOVERY.md),
[determinant1092 reconstruction](CURVE302_DET1092_RECONSTRUCTION_2026-09-07.md),
and [the completed parent](CURVE302_RECOVERED_MW17_PARENT_2026-09-07.md) document the route:

1. **Exclude inappropriate identifications.** The published X948 R17 j-map has no
   rational parameter for302; the degree24 inverse polynomial is irreducible
   modulo397. Broader finite atlases were also tested. Those exclusions concern
   specified fibrations, not all possible parents. The submitted first17 points
   were not automatically the latent generic subgroup.
2. **Construct a baseline parent.** Point-configuration methods produced a K3 of
   arithmetic generic rank9 specializing to302. Its own rational Picard ceiling
   prevents improving it to MW17 by merely changing that surface's fibration.
3. **Recognize a better lattice.** Of401 integral-point words in a candidate
   rank17 core, exactly one is the norm6 outlier; the other400 determine a
   minimum4 determinant1092 form. This is an inverse, target-conditioned step.
4. **Realize the lattice geometrically.** A degree4 polarization `O+P+Q`, of old
   fibre degree3, supplies15 predicted line sections. Linear incidence equations
   have a four-dimensional kernel: three concurrent-line directions and one
   useful configuration. All27 required incidences and78 nonincidences check.
   The lines and node determine a quartic `C L + W Q3`; the pencil `W=tL`
   leaves a cubic with a rational origin.
5. **Compile and complete sections.** Tangent/Cremona transformations give a
   Weierstrass equation. Fifteen lines plus the first two conics give17 displayed
   sections but only rank15. Two further conics raise the19-section collection
   to rank17; a unimodular subcollection gives the desired basis.
6. **Prove, rather than assume, the endpoint.** Exact height pairings, the
   discriminant-form saturation obstruction and Frobenius/Picard upper bounds
   close generic rank17. The fibre `t=0` is literally302, not just an equal-j twist.

The resulting equation has degrees8 and12, twenty-four geometric I1 fibres,
smooth infinity, trivial torsion and full saturated height determinant1092.
The [two-prime follow-up](CURVE302_PARENT_SEARCH_AND_GEOMETRY_2026-09-07.md)
closes geometric Picard19 and geometric generic MW17. X1092 is distinct from the
production determinant948 K3; recovering it is genuinely different-parent progress.

Write `D` for the **displayed** independent rank31 subgroup and `M17` for the
specialized generic group. Exact coordinate transport and Smith factors establish

`D / M17 = Z^14`, with `M17` primitive in `D`.

That explains the origin of seventeen directions, not the arithmetic cause of the
remaining fourteen. The [parity-domain proof](CURVE302_CONSTRUCTION_RECOVERY.md#the-parity-domain-has-no-missing-rational-halves)
closes2-saturation of the displayed group; further small-prime saturation work
still does not prove full saturation at every prime or exact rank31.

### Why the recovery engine succeeds on302

[V1](ADAPTIVE_HALF_LATTICE_VISIBILITY_2026-09-07.md) reached28;
[V2](ADAPTIVE_HALF_LATTICE_V2_2026-09-07.md) reached30 in906 charts;
[V3](ADAPTIVE_HALF_LATTICE_V3_2026-09-07.md) reached31 from the generic17 in1169 charts.
V3's M30-to31 winner entered via deterministic coverage, not the extreme-score
shortlist: its Babai rank was2286/8192 within its anchor. This directly motivates
coverage reserves rather than repeatedly tightening a top-only chart gate.

Two supplied-seed runs reached31 in855 and637 charts. The
[remaining seed panel](CURVE302_SEED_UNIVERSALITY_RESULTS_2026-09-09.md) brought the
combined result to13/14 tested seeds reaching31; the other stopped at29 under
its declared exposure. These are strong **conditional visibility/amplification**
results. The policy, target and/or supplied seeds were historically informed.
They neither constitute an unconditioned new rank31 discovery nor transfer a
13/14 success probability to fresh curves.

What works operationally is changing pointed charts as the known subgroup grows,
retaining coverage outside extreme metric scores, reconciling the entire returned
point cloud, and independently certifying the enlarged subgroup. There is still
no demonstrated cheap equation-only rule that selects another302-scale jump.

## 3. Rank knowledge: keep the fields and endpoints separate

| Object | Established endpoint | Not established by that endpoint |
|---|---|---|
| ICARM302 over Q | Rank at least31, torsion0, exact arithmetic fingerprint | Unconditional exact31; original construction; a32nd direction |
| Recovered X1092 parent | Arithmetic/geometric generic MW17; Picard19; saturated determinant1092 | Exceptional ranks of other fibres |
| Realized X1092 class1 | Another rational generic MW17 fibration on the same surface | A better specialization tail; a new K3 surface |
| X948 production portfolio | Six compact working MW17 presentations; two complete rootless J2 frame types | Six independent surfaces; automatic equivalence of their rational fibres |
| ICARM398 | Rank at least30; recovered A1/MW16 fibration on X948; blind displayed-quotient recovery | Unconditional exact30; two independent fibrations from its duplicate survivor labels |
| ICARM273 | Rank at least30 | A general record-producing family inferred from that lower bound |
| Explicit quadratic/base-change constructions | Independently proved MW18 or MW19 subgroup lower bounds on the specified new base | MW18 on the original Picard19 K3; new rank on an unchanged rational fibre |
| Fresh R17/curated foundry outputs | Exact independent point-subgroup lower bounds, up to27 in the cited panel | Exact ranks or worldwide novelty without a separate catalogue comparison |

Sources: [302](ICARM_CURVE302_RANK31.md),
[398](ICARM_CURVE398_RANK30_AND_CONSTRUCTION.md), [reproduction index](../REPRODUCE.md),
[rank balance and arithmetic marking](../../elkies-k3/RANK_MUTATION_AND_LIFT_THEOREMS.md).
Do not confuse `rank E(Q(t))`, `rank E(Qbar(t))`, a base-change rank, the rank of a
known finite subgroup, and `rank E_t(Q)`.

For a Jacobian fibration on either current Picard19 K3, Shioda–Tate gives
`generic geometric rank = 17 - reducible-root rank`. Base changing a surface can
leave the K3 category. A suggestion to obtain rational MW18 merely by selecting
another Picard20 K3 is not a licensed construction route; its arithmetic existence
requires a separate valid theorem, not geometric dimension counting.

## 4. Constructive toolbox: how we got objects, not just what frames exist

### A. Established fibration hopping plus project-specific inverse planning

A primitive U-embedding in NS fixes a frame. Changing the embedding changes the
fibration; removing reducible-fibre root rank releases MW rank by Shioda–Tate.
Kneser–Nishiyama, Weyl movement, linear systems and two-neighbour equation
conversion are established infrastructure. The project's useful additions are
inverse root/low-norm constraints, finite obstruction masks, marked physical
witnesses and fail-closed transport/equation certificates. See the
[theorem ledger](../../elkies-k3/RANK_MUTATION_AND_LIFT_THEOREMS.md) and
[route index](../../elkies-k3/CONSTRUCTION_ROUTES.md).

The historical H3/Kumar root-rich source and lattice corridors are useful
construction knowledge, but lattice closure and equation closure must be checked
separately. Published compact R17 is an independent usable endpoint. The practical
alternate-Q80 route is the **direct norm12/orbit11952 degree-two hop**, not the
old degree11511 transport or unfinished giant-coefficient suffixes. The noncyclic
`R17 -> 4A1/MW13 -> R17` route is another equation-explicit reusable bridge.
[Current route boundaries](../../elkies-k3/AGENTS.md).

### B. Point-cloud/lattice inversion and projective reconstruction

The curve302 pipeline in section2 is a concrete algorithmic template:
`specialized points -> candidate height form -> polarization -> lines/conics -> quartic -> pencil -> sections`.
It is valuable for discovering hidden parents, but its target-conditioned
construction must not be relabelled a prospective parameter selector.

### C. Complete frame census, then rational marked-U realization

The [X1092 census](DET1092_ROOTLESS_J2_CENSUS_2026-09-10.md) uses a rank7 auxiliary
with the opposite discriminant form,16 complete D5-anchor cases, residual ADE
chambers, exact shared norm budgets, rational LDL shell enumeration, primitive
embeddings and integral-isometry deduplication. It yields208 embeddings in the
declared Weyl cover, all in the final `2A7+2D5` anchor, and19 rootless frame types.
Curve302's type is class6; eighteen other types are nonisometric.

The first fifteen zero anchors were not evidence against the last. An old
assertion requiring every complement to match302 would abort a new type; it was
removed. Shared-budget pruning and an LDL indexing fix made the replacement
practical. Re-execution of the same enumerator and separate witness assembly
are not an independently implemented proof of completeness. J2 types are not
J1 surface-automorphism orbits and are not rational equation packets.

[Class1 realization](X1092_CLASS1_RATIONAL_MW17_REALIZATION_2026-09-10.md) then used
an existing degree-two quotient representative, mask109158 with norm12, to place
the desired frame in the rational source NS. Exact complement isometry identifies
embedding33. It was **not necessary to force the originally suggested common-core
identification to extend**. A negative-wall enumeration proves nefness; the old
rational zero remains a zero. A rank8 RR constraint matrix has kernel dimension2,
giving the pencil; a pointed quartic conversion supplies maps both ways.

The old degree-one section window spans only16 child directions. A generically
chosen norm10 bisection of intersection1 with the new fibre supplies the17th.
This is an important positive use of Euclidean bisections even though the separate
specialization-seed pilot failed. The final lattice has index1, determinant1092,
and exact rank17. Coefficient normalization reduces the recorded bit measure
from11921 to1347; it is not a proven minimal chart.

### D. A1/MW16 parent construction and inverse recognition

Minimum-norm-eight neighbours, RR pencils, pointed quartic conversion, rational
zeroes and section transport construct real A1/MW16 parents. The [398 recovery](ICARM_CURVE398_RANK30_AND_CONSTRUCTION.md)
screened63917 candidates and identified two presentations of **one** fibration;
PGL2/Weierstrass transport proves their equality. The source authors' successful
MW16 record lane is evidence against treating generic MW17 as a universal
prerequisite. The recent small **new-parent** A1 pilot's weak tail is not an
exclusion of the older productive MW16 families or of A1 constructions generally.

### E. Mestre/square-completion, Kihara and different-surface work

Polynomial/root configurations, square-completion and point-rich quartic models
produce explicit curves and seeds; projective linear systems can then produce
K3 parents. [Mestre determinant468 parents](MESTRE_DETERMINANT468_PARENTS_2026-09-07.md)
and [Kihara parent proofs](KIHARA_FIRST_PARENT_PICARD_AND_RANK_2026-09-07.md) document
actual additional surfaces, rather than new coordinates on X948. Their lower
rational Picard/generic-rank ceilings and bounded pilots remain explicit.
The different-NS foundry also requires a rational point on the **fully marked**
moduli curve before equation work; abstract lattice availability alone is not
rational realization. [Arithmetic-first scope](../../elkies-k3/AGENTS.md).

### F. Covers, bisections and genuine constructive M18 seeds

The [Euclidean formula](DET1092_EUCLIDEAN_BISECTION_FORMULA_2026-09-09.md) replaces
an RR solve/factorization with exact polynomial arithmetic once a generic trace
is obtained. Cold trace cost still counts. Constructed rational covers provide
function-field rank18 examples; two covers give a rank19 subgroup over a stated
common genus-one base, whose rational points may remain unknown.

There is already an [explicit infinite progression of certified M18 fibres](DET1092_CONIC_SEED_PROGRESSION_2026-09-08.md).
This is stronger than “a cover might split”, but its parameters/heights and later
amplification matter. The [orbit8044 pilot](ORBIT8044_SEED_FACTORY_2026-09-08.md)
produced exact seeds while its six cascades remained18 after114 charts each.
A [separate follow-up](DET1092_FUNNEL_FIRST_SEEDS_2026-09-08.md) reached21. These
positive constructions must not be erased by the four negative Euclidean cases.

### G. Arithmetic-class construction and membership certificates

[The constructed3/17 block](../rank-jump/CONSTRUCTED_CLASS_BLOCK_AND_RATIONAL_LIFTS.md)
starts with equation/generic-only principal dependencies, obtains two additional
strict classes and two independent half-ideal directions, then checks solubility
using known points **after the constructor froze**. The Artin matrix proves
class-group information, not a Cassels–Tate pairing or a rank lower bound.
The missing links remain cheap prospective generation and oracle-free lifting.

An exact finite/local homomorphism can certify nonmembership of a **supplied**
class if its joint image escapes the inherited span. It need not compute a whole
class group. Conversely, in-span local data are not a global relation. Coefficients
must be compatible across places; allowing a different inherited word at every
place discards information. A novel mod2 class can reflect saturation instead of
new rational rank. Final point/subgroup independence remains a separate gate.

### H. Pointed quartics, adaptive V3 and factor-free arithmetic

Exact point maps and bounded PARI point search turn a certified subgroup into
new point candidates. Exact-CVP checks concern the recorded rounded search
metric; finite-group certificates prove rank. Cloud reconciliation can recover
several directions from one call. One need not factor a discriminant or compute
an exact conductor to run the raw point search. See the
[pointed API](POINTED_QUARTIC_SEARCH.md), [V3](ADAPTIVE_HALF_LATTICE_V3_2026-09-07.md)
and [unattended foundry](HIGH_RANK_SEARCH_FOUNDRY_2026-09-09.md).

## 5. Experimental evidence with denominators and qualifications

| Experiment | Observed result | What it licenses |
|---|---|---|
|302 calibrated V3|17 to31 in1169 charts|A working recovery/amplification engine on this calibration, not fresh-fibre success odds|
|302 supplied-seed panel|13/14 reach31|Conditional robustness to these supplied seeds; not a seed selector|
|Fresh sixty R17 fibres|48 seeds;38 gain afterward;6090 calls; best27, three26|Strongest inspected prospective evidence for scaling this particular seed-first pipeline|
|R17 lower-height stratum|27/30 seed versus21/30 higher; all four bounds>=26 in lower stratum|Descriptive height/visibility evidence within a score-selected panel, not causal or universal|
|R17 late-gain cases|Six last gains at call75 or later; one first seed at62|Do not gate all depth on reaching20 in24 calls|
|Curated foundry snapshot|30 new internal curves22–26; one23-to25 improvement|Real prospective production beyond a calibration; not worldwide novelty|
|Retained MW16 score strata|10/1/0 added directions,860 boxes per arm|Top scores helped in this finite population; nine of the top arm's ten gains came from one curve|
|Fixed M27/M28 norm12 queues|6196 calls, no gain|Do not make deep repetition on these stalled states the default|
|Euclidean four certified-M17 UNKNOWNs|All four exact integer relations in M17|Those candidates are inherited; not all40917 bisections or all seed factories excluded|
|Class1 reported checkpoint|44 completed fibres;1056 calls; no gain;35.7 CPU minutes|A shallow detector result, not a family rank bound; larger capped panel's final endpoint not inspected here|

Sources: [sixty panel](R17_SIXTY_SEED_COMPLEMENT_PANEL_2026-09-09.md),
[foundry publication](FOUNDRY_CURVE_LEDGER_2026-09-09.md),
[score-strata comparison](RETAINED_MW16_SCORE_STRATA_2026-09-06.md),
[stalled queues](UNATTENDED_NORM12_SEARCH_2026-09-09.md),
[Euclidean relations](EUCLIDEAN_FOUR_SPLIT_ADMISSIONS_2026-09-10.md),
[class1 frozen protocol](CLASS1_PROSPECTIVE_ORDINARY_SEARCH_2026-09-10.md).
The44-fibre totals are a **user-reported local checkpoint**, not a committed
complete64+16 result. Earlier control selection at indices7 modulo8 was sign-
confounded; v2 uses a hash-selected control roster and a disjoint window. Do not
pool those versions as one randomized control group.

The60-panel rank histogram is `17:12,18:7,19:4,20:5,21:5,22:4,23:6,24:9,25:4,26:3,27:1`.
In particular103b2 at877/781 reaches27 in one seed plus100 complement calls,
with its last gain at67. These packets are useful warm-start research objects,
not “fresh discoveries” if searched again.

## 6. Latest theory reports not yet independently read from main

The user reports separately replayed local theory work:214 carrier rows in149
arithmetic cover groups,50 exceptional directions in distinct groups, a
42-pencil/173-pair smooth branch-incidence obstruction, exact inherited Kummer
image dimension17 for both generic MW17 presentations, two independent3/17
additions modulo all16 inherited classes, and inherited transfers for all four
class1 carriers. The reports distinguish generated subgroup rank1 from full
quadratic-twist rank. Their portable archives/checkers must be pushed and bound
before this ledger or the new runner treats them as available APIs.

The subsequent reported inherited-image timing is about2.5ms per parameter on
four inputs. All four ambient calculations time out after5s at discriminant
factorization. Three **synthetic** pipeline controls replay, with four points
from one soluble class. That does not demonstrate production solubility on class1.

Two corrections to the earlier discussion are essential:

* **Generic is not specialized.** A dimension-at-most3 everywhere-even generic
  quotient does not imply at most8 arithmetic classes on every rational fibre.
  Infinity in a function field is a parameter place, not automatically the real
  place of a specialized number field. Neither the bound nor the infinity-local
  injectivity shortcut transfers without a specialization theorem.
* **Membership is not construction or rank.** A fast novelty checker needs an
  actual candidate. It does not enumerate an unknown ambient quotient, solve a
  cover, or prove a new rational MW direction. Numerical/local features cannot
  replace missing certificates, and “no candidate found” is not quotient zero.

A discriminant smoothness/signature condition is retained, but it has not been
shown to predict quotient jumps. There is no general rational-polynomial matrix
in t whose minors have already been proved to detect these varying number-field
class-group phenomena. Stop treating that speculative suggestion as an API.

## 7. What has been excluded, and what remains genuinely open

Common smooth2-cover geometry, rationality of its first conic, and normalized
pencil determinant are automatic and survive rational/Sha switches. They cannot
serve as the missing global-solubility discriminator. [Exact counterexamples](../rank-jump/TWO_COVER_GEOMETRY_IS_AUTOMATIC.md).
Many low-degree relation and fitted-carrier multiplicity proposals fail specificity;
[the theory index](../rank-jump/README.md) preserves those finite scopes.

None of these results excludes singular carriers, other trace classes, deeper
pointed-chart coverage, different rational surfaces, or arithmetic methods with
an explicit new candidate source. Nor does finding19 nearly similar minimum4
frame lattices show that all their arithmetic specialization tails are similar.

The constructive frontier is still:

`generate an equation/parameter or a new arithmetic class without consulting its extra points -> find rational lifts -> certify independence -> amplify`.

Proving another obstruction is worthwhile only when it changes this pipeline or
closes a well-defined finite construction universe. It need not block an already
certified point-search family.

## 8. Executable next experiment and its interpretation

The [new runner](../cas/run_broad_rank_search.py) implements the
[broad runbook](BROAD_RANK_SEARCH_2026-09-10.md). It supersedes the shallow broad
proposal in the initial version of this ledger; the concurrently added
`run_broad_mw17_search.py` remains an intact historical alternative, not the
recommended default for this campaign.

Default exposure:

* six compact X948 MW17 presentations:256 score-selected plus64 score-independent
  slots each, using the existing complete43/49-class first-seed routine and up
  to100 complement calls **after any seed**, including M18/M19;
* recovered X1092 class6 and realized class1:64 selected plus16 controls each,
  using the existing generic-parent worker with up to198 initial calls;
*262144 fixed rational addresses per presentation, beginning at index131072:
  2097152 feature rows and2080 selected slots before exact-isomorphism attrition.

The initial native/generic detectors differ, so cross-parent comparisons measure
**parent plus detector**. Within each parent, score/control arms have the same
frozen rules. Same numerical t on two different fibrations is not the same
elliptic curve. Address sampling is not uniform over all rational numbers, and
parameter/coordinate height is chart-dependent. Report those limitations.

After every baseline stage, continue all eligible genuine gains under bounded
rank/staleness budgets. A predetermined approximately1/8 of rank17 misses gets
one deeper bank rescue independent of outcomes/scores. First-stage results remain
available separately from adaptively allocated endpoints. No iterative model-
training, scoring changes, new frame construction or LLM decisions occur inside
the run. Exact Q-isomorphism deduplication uses a j-index but never merges twists
by j alone. Exclusions apply only to the pinned inventory/provided equation lists
and this queue; no worldwide novelty is asserted.

Retain per-prime residue tables and every candidate feature row, exact selected
equations, point-call and gain timelines, complete raw clouds, both rank replays,
CPU/resource receipts, failed work and explicit UNKNOWNs. A saturated initial
subgroup is not presumed for every specialization; each must pass the existing
finite generic-seed gate. Conductor/class-group factorization is absent from the
critical path. Default four workers are configurable up to eight; the campaign
is finite, resumable, disk-guarded and source-frozen.

This implementation has Python policy/arithmetic/controller tests, including
independent small-prime point-count enumeration. **Sage/PARI are unavailable in
the implementation environment**, so no end-to-end Sage search or new rank is
claimed here. The runbook includes a small real smoke test; launch only after
its native and generic backend results replay successfully.

## 9. Decision after the run

Use baseline yield, incremental gain per CPU, first/last gain depth, score/control
comparisons within families, and UNKNOWN rates. Do not train on “miss=exact17”.
Do not infer underlying rank scarcity from a detector that fails its positive
calibrations. A useful lower bound discovered during a rescue is evidence about
premature pruning, not just another lucky parameter.

If X948 again produces gains and X1092 does not, continue the productive pool
while testing X1092 visibility separately; do not require another frame before
using compute. If X1092 becomes productive at fuller coverage, deepen it under
new frozen rules. If few gains occur anywhere, inspect calibration, coordinate
heights and score transfer before merely increasing address counts. Data from
this fixed run can train a simple later selector, but freeze it before exposing
a new address range and split validation by family/height/equation aliases.

The target remains a new certified rank32 subgroup or another explicitly chosen
arithmetic improvement. We already have useful new22–27 examples; “find a new23”
is not a new historical milestone. The next engineering milestone is a productive
unattended dataset, not another chain of prerequisite theorems.
