# Elkies--K3 process atlas

This is the chronological and reusable-process map for the Elkies--K3 work. It
does not replace a proof, certificate, canonical construction note, or
`MATH_STATUS.json`; it connects those sources and records why the successful
steps worked, why plausible steps failed, and which mechanisms can be reused.

The snapshot runs from the first recovered rank-17 lattice work on 2026-08-19
through the orbit42 output visible at 17:12 CEST on 2026-08-24. The final
orbit42 row is deliberately provisional until that concurrent lift is admitted
to the canonical note and `MATH_STATUS.json`.

## The invariant that organizes the search

For an elliptic K3 surface with a section, Shioda--Tate gives

```text
rho = 2 + rank(reducible-fibre roots) + rank(MW).
```

Thus the generic rank-19 family has a fixed budget of `17`, while its CM24
rank-20 specialization has a budget of `18`. An elliptic-neighbor move changes
the embedded hyperbolic plane `U=<F,O+F>` and therefore changes which divisor
classes look vertical and which look horizontal. On the same surface, a rank
that "appears" in MW is exactly balanced by root rank that disappears, and
conversely. No Neron--Severi class is created by the neighbor.

This distinction prevents three common category errors:

- A **same-surface neighbor** re-presents one fixed Neron--Severi lattice.
- A **specialization** can raise `rho` and genuinely add algebraic classes.
- An **equation lift** realizes a lattice-selected divisor as a rational pencil;
  it does not by itself discover a new lattice neighbor.

The state of the navigation engine is consequently not an ADE/MW label. It is
a marked elliptic fibration: the full NS lattice, an embedded `U`, fibre and
zero classes, physical components, marked horizontal classes, chamber data,
and lossless forward/inverse transports. Two nodes can have the same ADE/MW
label and still be different states.

## Chronology of discoveries and corrections

Times use, in order of preference, Git commit time, imported artifact mtime,
and dated-note evidence. The files imported into `artifacts/local/elkies-k3`
do preserve modification times, ranging back to 2026-08-20 and often resolving
the calculation order to seconds. Their filesystem birth/change times were
reset by the 2026-08-24 import, so an mtime is discovery-order evidence, not a
creation timestamp and never a mathematical certificate.

<!-- BEGIN GENERATED PROCESS CHRONOLOGY -->
| Time | Turn / behavior | How it was found | Mathematical meaning | Reusable rule | Status |
|---|---|---|---|---|---|
| 08-19 22:17 | endpoint recovery: The determinant-948 rank-17 MW lattice and its 1311 height-4 pairs entered the repository. | Recover the Gram matrix and mine additive motifs, cliques, and saturation data among short vectors. | The endpoint became an exact positive lattice constraint, but not yet an equation or a source construction. | Treat a large MW basis as a lattice fingerprint first; search for simpler fibrations before solving for all sections. | exact computation |
| 08-20 09:36 | direct reconstruction: The 167-variable three-hub and 48-variable rank-3 motif systems timed out; small residuals did not provide sections. | Export sliced polynomial section systems to msolve over GF(101). | The direct equation representation discarded too little structure and was computationally wrong-sized. | Compress with exact group-law/lattice operations before Groebner elimination; timeout is an experiment, not an obstruction theorem. | negative experiment |
| 08-20 12:22 | arithmetic source hunt: The determinant-948 transcendental lattice was matched to the (D,N)=(6,79) Eichler order and CM vectors of discriminants -3 and -24. | Even ternary-form filtering, exact quaternion/Eichler construction, Gross vectors, Atkin-Lehner orbits, and inverse-Clifford transport. | The recovered lattice was tied to Shimura arithmetic; the -3 vector cuts out the singular K3 with transcendental lattice A2. | Use discriminant forms, local ramification, and CM vectors to turn a lattice fingerprint into moduli coordinates and explicit boundary anchors. | exact computation |
| 08-20 14:51 | CM deformation: The explicit discriminant-3 K3 anchor was found; one E8 breaks to 2A2, producing E8+3A2/MW3, and an index-9 glue correction recovered regulator 316/9. | Embed the rank-19 NS as the orthogonal complement of a square -316 class and saturate the root-orthogonal section quotient. | Seventeen sections were exchanged for three sections plus root components; the missing factor 81 was a square of a saturation index, not a regulator contradiction. | Whenever regulator ratios are perfect squares, compute NS/(trivial+sections) before rejecting the model; glue is geometric data. | exact computation |
| 08-20 16:00 | reverse neighbor search: A reverse q90,4,4 path changed rootless/MW17 into E6+2A3+2A1/MW3. | Enumerate primitive isotropic classes f=(a,b,v), split a new U integrally, classify child roots, and retain lower-MW frames. | MW rank disappeared only relative to the new fibration: the same rank-19 NS budget moved into reducible fibre roots. | Search the graph of primitive U embeddings, using root rank and equation cost as independent objectives. | exact path from bounded search |
| 08-20 21:38 | equation attack: Finite-field A10/E6 section charts produced modular points and triangular reductions, but the split E6 chart and several seeds did not come from the transported fibration. | Tate coefficient gates, coordinate slicing, GCD scans with collision-factor saturation, finite-field enumeration, and Hensel probes. | An abstract ADE/MW label does not fix a marked elliptic fibration or a rational equation chart. | Execute the neighbor geometrically and transport F,O,components,sections before normalizing an equation ansatz. | mixed exact reductions and superseded experiments |
| 08-20 23:39 | Kodaira correction: The all-IV j=0 deformation was obstructed by MW-rank parity; a mixed I3/IV non-isotrivial family replaced it. | Compare the order-three automorphism/parity constraint with Shioda-Tate and derive a family with nonconstant j. | The root lattice A2 does not determine Kodaira type IV versus I3; fibre realization controls automorphisms and rank parity. | Enumerate Kodaira realizations of an ADE root system and test monodromy/automorphism constraints before section solving. | exact obstruction and viable replacement family |
| 08-21 00:23 | low-q backtrack: The bounded q25,4,4,4 backtrack reached E6+D4+2A2+A1/MW2 and was later reconstructed over QQ. | Widened beam search, exact U-splits, component-glue recovery, exhaustive GF(23) chart, CRT, and Hensel reconstruction. | The endpoint can be represented by two sections and root rank 15, but this backtrack does not identify the level-474 source polarization. | Minimize symbolic equation cost after lattice rank: number of sections, pole order, pair gates, and Weierstrass degree all matter. | exact path; explicit specialization has a Picard-rank-20 correction |
| 08-21 00:48 | Kumar comparison: Three E7+E8/MW2 Kumar frames were classified; H2 exposed q60 and q80 comparison fibrations, with Q80=E6+D5+A3/MW3. | Binary height forms, discriminant glue, exact neighbors, and bounded CM-stability scoring. | Q80 is a useful secondary polarization/compiler chart on the same lattice genus, not the source entrance. | Maintain a fibration atlas on one NS lattice; different polarizations should be scored for source meaning and compiler cost separately. | exact candidates from bounded search |
| 08-21 16:16 | H3 route selection: The H3 q6 child E8+E6/MW3 and its degree-two rank-growing continuations were found at lattice/chamber level. | Exact neighbor enumeration, root adaptation, nef chamber reduction, and old-fibre-degree scoring. | A forward source-to-R17 corridor became available independently of the reverse low-MW ancestry. | A selected route is a certified path, not a shortest-path theorem; preserve alternatives and the search boundary. | exact retained path from bounded searches |
| 08-21 17:39 | source recovery: The H21 intersection H92 component identified H3, not H2, as the level-474 source; the published genus-two model was normalized into the exact E7+E8/MW2 family. | Humbert intersection, modular solving, p-adic lifting, exact rational reconstruction, and comparison with the published rational point. | The construction direction changed from reverse discovery to forward execution from the actual QM/Kumar source. | Use moduli geometry to choose the source polarization; lattice isometry alone cannot tell which embedded U is historically or arithmetically canonical. | proved |
| 08-21 20:20 | first equation neighbor: The first q6 H3 neighbor was executed exactly, giving E8+E6/MW3. | Resolved Riemann-Roch pencil, binary quartic, Jacobian conversion, and fibre classification. | The lattice neighbor acquired an actual characteristic-zero equation and transported markings. | A neighbor certificate should have two layers: integral NS transport and an equation-level RR/Jacobian realization. | proved |
| 08-22 08:11 | q8 false frontier: Large q8 module ambients repeatedly appeared to have zero surviving kernels, generating bounded modular obstructions that were later superseded. | Principal-part matrices, local E7/E8 component modules, saturation probes, and multi-prime rank checks. | The computations were correct for the supplied wrongly marked/normalized modules but did not obstruct the true q8 divisor. | A full-rank obstruction is only as valid as its marking and module; pin the divisor identity and normalization before scaling linear algebra. | superseded bounded experiments |
| 08-22 15:13 | module saturation: Reverse saturation exposed adapted coordinates such as x-m^2 and separated local component modules from global fibre twists; several apparent obstructions vanished or changed dimension. | One-uniformizer saturation, exact/two-prime kernel comparison, and resolved chart bases. | Cancellations belong to the saturated pushforward module, not to arbitrary monomial pole floors. | Compile local line-bundle modules by saturation and quotient conditions; never count exceptional components as independent scalar rows by default. | exact local modules plus superseded global interpretations |
| 08-22 15:34 | Q80 low-q route: The generic Q80 low-q corridor reached a rootless/MW17 frame with every retained suffix divisor of old-fibre degree two. | Alternate-neighbor beam, chamber scoring, exact section/CVP walls, and stagewise integral transport. | A second complete lattice route demonstrated that the H3 selected path is not unique and supplied simpler compiler experiments. | Represent the field as a route graph, not one chain; retain multiple exact paths and optimize them by different costs. | proved lattice corridor |
| 08-22 20:31 | Q80 CM24 compiler: The Q80 specialization corridor was compiled over GF(73), showing repeated pole/intersection collapse and connected ADE quotient conditions. | Specialize actual divisors, saturate local modules, solve quotient-line conditions, compile quartics, and classify Jacobians. | At CM24, extra NS classes change roots, pole orders, and MW labels; specialized horizontals are regressions, not generic definitions. | Specialization is a typed operation with its own marking and Picard rank; recompute every divisor invariant after base change. | exact modular regression, later lifted to characteristic zero |
| 08-22 21:09 | q8 repair: Two bugs were isolated: the quartic covariant 2-cover was doubled again, and the q-normalizer omitted Dx. Collision degree fell 46 to 10; the exact 13-to-2 RR pencil produced D13/MW4. | Global height audit, direct rational-function residue derivation, two-prime regression, and exact QQ Jacobian certification. | The apparent high-degree obstruction came from a wrong group-law scale and a hidden vertical pole, not from the geometry. | Track isogeny/cover degree on every point map and clear complete rational expressions before CRT/local normalization. | proved correction and equation |
| 08-23 05:55 | Q80 characteristic-zero closeout: The terminal Q80 CM24 q6 marking was recovered as a difference of easier exact sections, closing the exact QQ(sqrt(-3)) RR pencil and 4A2+A3+A5/MW2 child. | Replay unimodular transports, reconstruct high-incidence sections, identify a group-law difference modulo 73, then certify over characteristic zero. | A singular direct lift was bypassed without changing the selected divisor; generic rootless/MW17 remained distinct from the CM24 child. | When a marked section is singular to lift, construct a small exact MW basis and recover the target by group operations plus modular identification. | proved specialization equation |
| 08-23 11:24 | q24 route recovery: The D13 q24 orbit85 child and the complete D12-to-R17 degree-two lattice suffix were recovered; route-scout also showed valid but differently marked alternatives. | Forward replay, orbit classification, chamber reduction, suffix beam search, and endpoint isometry. | Matching ADE/MW sequences is insufficient: the selected D13/D12 marking determines whether the terminal rootless frame is pinned R17. | Every route edge must carry its full U embedding and forward/inverse transport, not just q and ADE labels. | proved lattice/marking path |
| 08-23 15:56 | q24 resolved geometry: The I9* resolution distinguished blow-up centres from twelve physical components; an ordinary-jet shortcut overcounted conditions, while the resolved component quotient left a two-plane. | Explicit blow-up graph, effective-cluster transport, smooth-collision quotient, and resolved valuation RR matrices. | Cancellation and dependence occur along connected exceptional geometry; centre counts are not component counts. | Build a physical component graph and valuation atlas before imposing local conditions; validate Euler/root data after compilation. | exact modular discovery, later lifted to QQ |
| 08-23 22:17 | q24 horizontal lift: Multi-prime reconstruction and structured Hensel lifting produced the exact q24 horizontal section over QQ. | Modular section sampling, CRT/rational reconstruction, structured Newton-Hensel lift, and literal equation substitution. | The selected lattice divisor acquired an exact rational horizontal representative. | Separate section reconstruction from the full RR lift; checkpoint exact identities so expensive local geometry need not be repeated. | proved |
| 08-24 11:16 | q24 equation closeout: The characteristic-zero resolved RR lift compiled a quartic whose Jacobian is D12/MW5, closing H3-03 at equation level. | Resume checkpointed valuation elimination, certify a 56-to-48-to-8-to-6-to-2 dimension chain, compile and classify the Jacobian. | The H3 equation frontier moved from D13 to D12 while the later suffix remained lattice-only. | Checkpoint long exact eliminations by invariant subspace and hash; certify each dimension loss by its geometric condition block. | proved |
| 08-24 11:42 | orbit42 preflight: The D12-to-A11 orbit42 profile, physical I8* marking, and eighteen identity-class zero-pole sections were recovered; fast point transport and rational-halving shortcuts were rejected. | Resolved component marking, finite-field scans, Hensel/QQ section recovery, halving audits, and exact degree preflight. | The selected lattice edge survived, but easy point transports did not supply its resolved line-bundle trivialization. | Before a lift, audit actual rational-map degrees and component orientations; negative shortcut results should redirect the compiler without reopening the lattice selection. | proved prerequisites and shortcut rejections |
| 08-24 16:31 | orbit42 exact section recovery: The opposite spinor zero-pole pair completed the twenty-point shell, and group-law combinations produced four exact P.O=3 orbit42 section candidates. | Characteristic-zero elimination for the spinor pair, exact QQ(u) group law, and pinned mod-100003 shell regression. | The hard horizontal classes became exact equation points, while physical-orientation selection still belonged to the resolved RR problem. | When a marked target is hard to lift directly, complete a small exact section shell and reconstruct candidates through exact group operations plus modular marking data. | proved construction aids with a modular marking boundary |
| 08-24 17:12 | orbit42 active result: The active lift emitted an exact-QQ artifact reporting a resolved RR plane, quartic, Jacobian, and A11 fibre classification. | Resolved component-valuation Riemann--Roch elimination followed by quartic/Jacobian compilation and exact identity checks. | If admitted by the canonical note and status authority, this closes the D12/MW5-to-A11/MW6 equation edge; at this snapshot the artifact remains active concurrent work. | Record newly emitted outputs immediately but keep a provisional type until the canonical proof boundary, software lock, and status entry are updated. | provisional active artifact; not yet promoted by MATH_STATUS.json |
| 08-24 18:09 | orbit42 equation promotion: The resolved RR result was promoted to the status authority, and the exact identity-shell degree fingerprint selected equation-side orbit64/mapping7 in the physical C10 orientation. | Replay the exact QQ pencil and Jacobian certificate, restrict the exact shell at the pinned good prime 100003, and compare all A11 q6 shell mappings. | The D12/MW5-to-A11/MW6 equation edge is closed; only the shell-to-lattice orientation identification retains a stated good-reduction boundary, and the active frontier moves to A11-to-2A5. | Separate exact equation construction from finite-field marking orientation, state both proof scopes explicitly, and only then advance the canonical frontier. | proved equation edge with pinned-good-reduction marking boundary |
| 08-24 19:24 | A11 q8 construction audit: The historical q8/orbit922 divisor was recognized as O+P-2F with zero vertical-root correction and transported through the complete root/MW/glue marking. The equation-side nef targets collapsed to orbits 12 and 2162, and orbit12 was selected by a declared minimum-MW-L1 convention. | Decompose the positive-frame isometry into the A11 root chain and rank-six MW quotient, enumerate the order-16 MW automorphism group, impose integral glue, match every transported ray to the exhaustive q8 search, and recursively eliminate eight infinity coefficients in the modular section chart at the preferred prime 43. | The active q8 lift now has a construction-compatible equation-side divisor and a materially smaller discovery system; the section function and characteristic-zero pencil remain open. | Before solving a new neighbor, replay how the historical divisor was assembled. Transport formulas, vertical support, root orientation, MW quotient and integral glue are stronger fingerprints than q, ADE type or a local component label alone. | exact partial target certificate plus bounded modular benchmark |
| 08-24 20:30 | A11 q8 target-coset audit: The eighteen exact identity-shell points were found to span an index-five sublattice of the saturated rank-five MW hyperplane. The previously proposed pole-order-four missing direction lay in the wrong coset, while exhaustive exact enumeration selected the minimum-pole bridge M=(1,0,0,0,0,1) and the literal word P12=M+S6-2*S2-2*S8. | Compute the exact shell Smith index, enumerate all MW vectors through the pole-order-five height bound, filter by the transported target coset, verify the group-law coordinate identity, and compile a free-infinity finite-field chart for the bridge rather than a fixed leading branch. | The construction fingerprint extends beyond the target orbit to a precise quotient-lattice coset. The open symbolic solve drops from pole order six to pole order five without changing the marked target. | After transporting a historical divisor, compute the saturation quotient generated by known equation sections and solve in the target coset; a generic missing-rank generator may be both cheaper and irrelevant. | exact target-coset theorem plus bounded modular benchmark |
<!-- END GENERATED PROCESS CHRONOLOGY -->

## Rank-flow routes

Every row is checked against Shioda--Tate. A status of
`lattice-proved-equation-open` means that the integral neighbor and marking are
exact, but the characteristic-zero Riemann--Roch/Jacobian realization is not
yet closed. A path returned by a bounded search is a certified path, not a
shortest-path theorem.

<!-- BEGIN GENERATED RANK ROUTES -->
| Route | Stage change | q / old degree | Root rank | MW rank | Rank interpretation | Equation status |
|---|---|---:|---:|---:|---|---|
| Low-q reverse backtrack | start: rootless/MW17 | - | 0 | 17 | rho=19; budget=17 | proved |
|  | rootless/MW17 -> A3+7A1/MW7 | 25 / - | 10 | 7 | roots +10; MW -10; sum fixed | proved |
|  | A3+7A1/MW7 -> D4+A3+2A2+2A1/MW4 | 4 / - | 13 | 4 | roots +3; MW -3; sum fixed | proved |
|  | D4+A3+2A2+2A1/MW4 -> A5+D4+2A2+A1/MW3 | 4 / - | 14 | 3 | roots +1; MW -1; sum fixed | proved |
|  | A5+D4+2A2+A1/MW3 -> E6+D4+2A2+A1/MW2 | 4 / - | 15 | 2 | roots +1; MW -1; sum fixed | proved |
| H3 selected corridor | start: E7+E8/MW2 | - | 15 | 2 | rho=19; budget=17 | proved |
|  | E7+E8/MW2 -> E8+E6/MW3 | 6 / 2 | 14 | 3 | roots -1; MW +1; sum fixed | proved |
|  | E8+E6/MW3 -> D13/MW4 | 8 / 2 | 13 | 4 | roots -1; MW +1; sum fixed | proved |
|  | D13/MW4 -> D12/MW5 | 24 / 2 | 12 | 5 | roots -1; MW +1; sum fixed | proved |
|  | D12/MW5 -> A11/MW6 | 6 / 2 | 11 | 6 | roots -1; MW +1; sum fixed | proved with pinned-good-reduction marking boundary |
|  | A11/MW6 -> 2A5/MW7 | 8 / 2 | 10 | 7 | roots -1; MW +1; sum fixed | lattice-proved-equation-open |
|  | 2A5/MW7 -> 3A3/MW8 | 4 / 2 | 9 | 8 | roots -1; MW +1; sum fixed | lattice-proved-equation-open |
|  | 3A3/MW8 -> A3+2A2/MW10 | 4 / 2 | 7 | 10 | roots -2; MW +2; sum fixed | lattice-proved-equation-open |
|  | A3+2A2/MW10 -> 5A1/MW12 | 4 / 2 | 5 | 12 | roots -2; MW +2; sum fixed | lattice-proved-equation-open |
|  | 5A1/MW12 -> 4A1/MW13 | 4 / 2 | 4 | 13 | roots -1; MW +1; sum fixed | lattice-proved-equation-open |
|  | 4A1/MW13 -> 3A1/MW14 | 4 / 2 | 3 | 14 | roots -1; MW +1; sum fixed | lattice-proved-equation-open |
|  | 3A1/MW14 -> 2A1/MW15 | 4 / 2 | 2 | 15 | roots -1; MW +1; sum fixed | lattice-proved-equation-open |
|  | 2A1/MW15 -> A1/MW16 | 4 / 2 | 1 | 16 | roots -1; MW +1; sum fixed | lattice-proved-equation-open |
|  | A1/MW16 -> rootless/MW17 | 6 / 2 | 0 | 17 | roots -1; MW +1; sum fixed | lattice-proved-equation-open |
| Q80 generic corridor | start: E6+D5+A3/MW3 | - | 14 | 3 | rho=19; budget=17 | proved |
|  | E6+D5+A3/MW3 -> D9+A4/MW4 | 4 / 2 | 13 | 4 | roots -1; MW +1; sum fixed | proved |
|  | D9+A4/MW4 -> D7+D5/MW5 | 4 / 2 | 12 | 5 | roots -1; MW +1; sum fixed | proved |
|  | D7+D5/MW5 -> D7+D4/MW6 | 6 / 2 | 11 | 6 | roots -1; MW +1; sum fixed | proved |
|  | D7+D4/MW6 -> A6+A4/MW7 | 4 / 2 | 10 | 7 | roots -1; MW +1; sum fixed | proved |
|  | A6+A4/MW7 -> A6+A3/MW8 | 4 / 2 | 9 | 8 | roots -1; MW +1; sum fixed | proved |
|  | A6+A3/MW8 -> A4+A2+A1/MW10 | 6 / 2 | 7 | 10 | roots -2; MW +2; sum fixed | proved |
|  | A4+A2+A1/MW10 -> A3+A2/MW12 | 4 / 2 | 5 | 12 | roots -2; MW +2; sum fixed | proved |
|  | A3+A2/MW12 -> 4A1/MW13 | 4 / 2 | 4 | 13 | roots -1; MW +1; sum fixed | proved |
|  | 4A1/MW13 -> A1/MW16 | 4 / 2 | 1 | 16 | roots -3; MW +3; sum fixed | proved |
|  | A1/MW16 -> rootless/MW17 | 6 / 2 | 0 | 17 | roots -1; MW +1; sum fixed | proved |
| Q80 CM24 specialization shadow | start: 2A6+3A1/MW3 | - | 15 | 3 | rho=20; budget=18 | proved |
|  | 2A6+3A1/MW3 -> A5+2A4+2A1/MW3 | 6 / - | 15 | 3 | roots +0; MW +0; sum fixed | proved |
|  | A5+2A4+2A1/MW3 -> 2A4+2A3+A1/MW3 | 4 / - | 15 | 3 | roots +0; MW +0; sum fixed | proved |
|  | 2A4+2A3+A1/MW3 -> A1+2A3+2D4/MW3 | 4 / - | 15 | 3 | roots +0; MW +0; sum fixed | proved |
|  | A1+2A3+2D4/MW3 -> A1+A2+A3+A4+A5/MW3 | 4 / - | 15 | 3 | roots +0; MW +0; sum fixed | proved |
|  | A1+A2+A3+A4+A5/MW3 -> 4A2+A3+A5/MW2 | 6 / - | 16 | 2 | roots +1; MW -1; sum fixed | proved |
<!-- END GENERATED RANK ROUTES -->

The two generic forward corridors show that navigation is a graph rather than
a unique chain. The H3 route is the arithmetically correct level-474 source
route; Q80 is a valuable secondary compiler route. The CM24 rows are a typed
specialization shadow and must never be substituted for the generic Q80
endpoint.

## Where cancellation and false rank enter

"Cancellation" has several mathematically different meanings. Conflating them
caused most of the expensive false frontiers. The ledger keeps them separate.

<!-- BEGIN GENERATED CANCELLATION MECHANISMS -->
| Mechanism | Diagnostic symptom | Exact response | Meaning | General rule |
|---|---|---|---|---|
| Shioda-Tate rank exchange | MW ranks appear or disappear under a neighbor. | Hold NS and rho fixed and verify delta(root_rank)+delta(mw_rank)=0. | Sections and reducible-fibre components are two presentations of the same divisor budget. | Use root rank as a controllable storage channel for MW rank when navigating fibrations. |
| Weyl fixed-component cancellation | A primitive isotropic class has negative intersections or large raw old-fibre degree. | Reflect in supplied effective (-2)-curves until the movable nef representative is reached; retain the removed fixed component sequence. | The divisor class is unchanged up to the Weyl action, while its effective presentation becomes a fibration pencil. | Make chamber reduction a first-class, replayable operation and never infer nefness from root type alone. |
| NS glue saturation | A computed MW regulator differs from the discriminant prediction by a perfect square. | Saturate trivial lattice plus section classes inside NS and project all glue cosets. | The missing classes are finite-index glue, not new analytic sections. | Record Smith invariants and torsion/glue before reporting an MW height lattice. |
| Cover-degree normalization | Heights and collision degrees are multiplied by four. | Account for the degree-two binary-quartic covariant map before applying the elliptic group law. | The transported point was 2P rather than P. | Annotate every point map with isogeny/cover degree and test it by canonical heights. |
| Full-denominator pole cancellation | A supposedly regular CRT frame retains a high-degree vertical pole. | Clear Nx/Dx and Ny/Dy completely; the q8 residue is Ny*Dx/(h*Dy), not Ny/(h*Dy). | A missing denominator factor changed the line bundle being compiled. | Derive local residues symbolically before modular normalization and regression-test generic branch degree. |
| Resolved-module saturation | Independent component rows overconstrain an RR pencil or raw minors share a spurious factor. | Use the saturated pushforward module and connected ADE quotient lines on physical resolved components. | Exceptional conditions have syzygies and cancellations dictated by the resolution graph. | Compile valuation modules and quotient conditions, not a row-per-centre heuristic. |
| Specialization rank and pole change | Generic P.O, MW coordinates, roots, or rank do not survive at CM24. | Transport the divisor, specialize, re-chamber, and recompute Picard rank, roots, heights, poles, and vertical support. | Extra algebraic classes can move rank between roots and sections; a specialized section label need not define the generic divisor. | Treat specialization as a new typed stage linked by a specialization map, never as an in-place coefficient substitution. |
| Marking-preserving transport | Two frames have the same ADE/MW label but lead to different endpoints or enormous pullback degree. | Retain the embedded U, zero, fibre, components, marked horizontals, and determinant-one NS transport in both directions. | An NS isometry may change the elliptic fibration; ADE/MW labels are not object identifiers. | Use marked-fibration nodes and lossless transport matrices as the graph state. |
| Numerical collision false rank | A tiny residual suggests a new section but its x-coordinate approaches an existing section and target/non-target incidence errors are comparable. | Require scale-free separation, edge-gap and wrong-sign-gap tests before promoting numerical rank growth. | Near-collision can mimic an independent MW direction. | Rank discovery needs independence/collision guards before exact lifting. |
| Collision-factor saturation | Repeated factors arise from normalized fibres or lambda=mu instead of distinct reducible fibres. | Saturate known boundary/collision factors and demand squarefree residual factors before lifting. | Raw gcd degree can count a degenerate collision locus rather than the desired Kodaira configuration. | Classify and remove geometric boundary factors before interpreting finite-field hits. |
| Construction-fingerprint transport | Many low-q neighbors have the same ADE child and local component profile, while a generic section solve is prohibitively large. | Recover the historical divisor formula and vertical support, split the marking isometry into root and MW blocks, enforce integral glue, match the transported divisor literally against the equation-side chamber search, and compute its coset modulo the saturated lattice generated by known equation sections. | The way a divisor was assembled can identify its equation-side orbit and the exact quotient-lattice coset that must be solved, even when ADE and q cannot. | Always audit the construction before launching symbolic elimination; preserve O/P/F coefficients, vertical-root support, root orientation, MW coordinates, full glue transport and the shell saturation quotient as first-class fingerprints. |
<!-- END GENERATED CANCELLATION MECHANISMS -->

There are four especially durable lessons:

1. A perfect-square regulator discrepancy points first to finite-index glue.
   In the rank-3 source calculation the ratio `2844/(316/9)=81=9^2` identified
   an index-9 saturation, rather than a contradiction.
2. A binary-quartic covariant is already a 2-cover. Doubling once more turned a
   height-24 class into height 96 and inflated collision degree from 10 to 46.
3. Local regularity belongs to the complete rational expression. The q8
   normalizer must use `Ny*Dx/(h*Dy)`; dropping `Dx` changed the compiled line
   bundle by leaving a vertical pole.
4. Resolved ADE constraints live on connected physical-component modules.
   Blow-up centres, monomial pole floors, and scalar row counts are not reliable
   substitutes for saturated quotient conditions.

## A generalized navigation and compilation engine

The reusable engine is a typed graph with independently checkable layers:

```text
arithmetic/moduli source
        |
        v
marked NS state --exact neighbor--> marked NS state --...--> target frame
        |                              |
        | RR pencil                    | specialization (new typed node)
        v                              v
equation + transported markings     specialized equation/marking
        |
        v
fibre, height, discriminant, and endpoint certificates
```

For each node, record:

- `NS`, its discriminant form and saturation/glue data;
- `F`, `O`, physical reducible-fibre components and the chosen nef chamber;
- the MW basis, height matrix, torsion and component incidences;
- generic versus specialized coefficient field and Picard rank;
- an equation certificate when one exists.

For each edge, record:

- the primitive isotropic divisor and `q` factorization;
- the fixed-component/Weyl reflection replay;
- determinant-one forward and inverse NS transports;
- the old-fibre degree, orbit identity and search boundary;
- the resolved Riemann--Roch module, binary quartic and Jacobian certificate;
- exact hashes/checkpoints for long eliminations.

A practical search/compile loop is:

1. Enumerate primitive isotropic candidates in the current marked NS lattice.
2. Remove fixed components by replayable Weyl reflections and require nefness.
3. Split the new `U` integrally; saturate the frame and classify roots, glue and
   MW height data.
4. Reject rank-budget, determinant, discriminant-form or marking failures.
5. Score surviving edges separately for mathematical destination and equation
   cost: old-fibre degree, pole order, RR ambient dimension, component-module
   complexity, coefficient height and available exact sections.
6. Keep multiple Pareto-optimal routes instead of one greedy chain.
7. Compile a selected divisor on resolved charts; derive saturated local
   modules and connected-component quotient conditions.
8. Use finite fields only for discovery/regression. Lift by CRT/Hensel or exact
   group-law reconstruction, then certify by literal substitution.
9. Compile the quartic/Jacobian and verify fibres, Euler number, Shioda--Tate,
   heights, discriminant and the pinned endpoint marking.
10. On specialization, create a new node and recompute every invariant; never
    inherit generic MW coordinates, poles or root labels unchecked.

This architecture generalizes beyond the current K3: the lattice navigator is
independent of equation coordinates, while the compiler is a collection of
local resolved-module adapters. New fibre types or coefficient fields should
extend adapters and scoring, rather than fork the mathematical state model.
The precise conservation, determinant, nef-pencil, specialization, and
conditional lift statements are developed in
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md).

## Machine-readable ledger and checks

The source of this atlas is
[`data/process/elkies_k3_process_ledger.json`](data/process/elkies_k3_process_ledger.json),
validated against
[`../schemas/elkies_k3_process_ledger.schema.json`](../schemas/elkies_k3_process_ledger.schema.json).
It records stages, transitions, chronology, mechanisms, evidence and literature.
The dependency-free checker is
[`scripts/analyze_process_ledger.py`](scripts/analyze_process_ledger.py).

```bash
# Validate ADE ranks, Shioda--Tate budgets, edge rank exchange, chronology,
# evidence paths, MATH_STATUS identifiers and route connectivity.
python3 elkies-k3/scripts/analyze_process_ledger.py

# Verify that the three generated tables in this document match the ledger.
python3 elkies-k3/scripts/analyze_process_ledger.py --check-document

# Refresh the generated tables after an intentional ledger update.
python3 elkies-k3/scripts/analyze_process_ledger.py --update-document

# Write the compact, deterministic consistency audit stored in
# artifacts/generated-results/elkies-k3-process-ledger-audit.json.
python3 elkies-k3/scripts/analyze_process_ledger.py \
  --write-audit artifacts/generated-results/elkies-k3-process-ledger-audit.json

# Audit the still-untracked imported artifacts by preserved mtime. This reads
# local metadata only; it does not promote mathematical status.
python3 elkies-k3/scripts/analyze_process_ledger.py --audit-local
```

The checker deliberately does not rerun Sage, Magma, Singular or the expensive
symbolic calculations. Those remain in the cited canonical notes and
certificates, with `MATH_STATUS.json` as the sole status authority.

## Literature bridge

The engine packages standard ideas into a provenance-aware workflow rather
than claiming a new general theorem. Its mathematical backbone is:

- T. Shioda, [*On the Mordell--Weil Lattices*](https://rikkyo.repo.nii.ac.jp/records/10027):
  height pairings, trivial lattices, regulators and discriminants.
- A. Kumar, [*K3 surfaces associated with curves of genus two*](https://arxiv.org/abs/math/0701669):
  the explicit genus-two/Kumar `E8+E7` source construction.
- N. Elkies, [*Shimura curve computations via K3 surfaces of Neron--Severi rank at least 19*](https://arxiv.org/abs/0802.1301):
  rank-19 K3 families, Shimura curves and CM rank-20 points.
- A. Kumar, [*Elliptic fibrations on a generic Jacobian Kummer surface*](https://arxiv.org/abs/1105.1715):
  primitive nef isotropic classes, Weyl chambers and explicit 2-neighbor moves.
- N. Elkies and A. Kumar,
  [*K3 surfaces and equations for Hilbert modular surfaces*](https://arxiv.org/abs/1209.3527):
  moduli navigation using explicit elliptic K3 and Shioda--Inose models.
- O. Padurariu and F. Saia,
  [*Shimura curve Atkin--Lehner quotients of genus at most two*](https://arxiv.org/abs/2509.25368):
  the published quotient equations used in the level-474 source identification.

The local contribution is the exact, lossless coupling of those layers:
arithmetic source identification, marked-lattice navigation, resolved equation
compilation, specialization typing, negative-result provenance and automatic
rank/cancellation guards.
