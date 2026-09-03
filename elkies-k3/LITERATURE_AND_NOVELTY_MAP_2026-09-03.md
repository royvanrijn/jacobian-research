# Literature and novelty map for target-directed fibration hopping

Date: 2026-09-03

This is the provenance companion to
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md).
It incorporates the primary-source literature audit supplied for repository
snapshot `fb056eaa153f70da49e5268b6cefaca039b9b819`.  The theorem note is the
canonical proof source; this file classifies claims and does not duplicate
their proofs.

## Public framing

Use **target-directed fibration hopping** or **a lattice-guided construction
calculus for elliptic K3 fibrations**.  “Integral rank transfer” is the project
name, not one new foundational theorem.  The forward machinery is established:

- Shioda--Tate and the Shioda height/discriminant formula;
- Nikulin's isotropic-overlattice and primitive graph-gluing formalism;
- Kneser--Nishiyama frame classification and the `J0`/`J1`/`J2` distinction;
- Kneser `p`-neighbours, with their stated local and connectivity hypotheses;
- primitive-`U` changes, Weyl movement to nef fibre classes, and explicit
  equation-level fibration hopping.

The defensible new-looking centre is the inverse controller: prescribe a
low-norm or root-system outcome, compile finite discriminant-coset and
physical-witness conditions, direct the neighbour search, minimize one-edge
incidence over copies of a `J2` class, and preserve fail-closed boundaries
through marking, equation, and arithmetic certification.

## Source loci used below

- `SHIODA-MW`: Shioda, *On the Mordell--Weil lattices*, especially the
  Néron--Severi quotient, height pairing, torsion and discriminant formula;
  Schütt--Shioda, *Elliptic Surfaces*, Sections 6 and 11.
- `SHIODA-GAL`: Shioda, *Mordell--Weil lattices and Galois representation I,
  III*, especially III, Section 5; Schütt--Shioda, *Mordell--Weil Lattices*,
  chapter “Galois Representations and Algebraic Equations”.
- `NIKULIN`: Nikulin, Propositions 1.4.1 and 1.5.1 for even overlattices and
  primitive graph gluing; the indefinite uniqueness theorem for `H2`.
- `KN`: Nishiyama's Kneser--Nishiyama method; Braun--Kimura--Watari,
  Sections 2--3 and Proposition C' for `J0`/`J1`/`J2` and multiplicity.
- `HOP`: Brandhorst--Elkies, Section 2, especially Lemmas 2.5--2.6, Remark
  2.10, and Sections 2.3--2.4; Kumar, Appendix; Elkies--Kumar, Section 5.
- `KNEIGH`: Chenevier, equations (1.2)--(1.3) and the large-prime
  equidistribution theorem; Voight, Section 3 and Theorem 3.18 for algorithms
  and connectivity hypotheses.
- `THETA`: Bruinier--Stein for the Weil representation and Hecke operators;
  Kane--Kim for lattice-coset theta series and neighbour algorithms; Müller
  for the Weil-representation basis problem.
- `PASTEN-SALGADO`: Pasten--Salgado, Theorem 1.1.
- `TATE`: Kubert's marked-point Tate normal form; the two-point identity in
  `F0` is proved directly and has no separate antecedent located.
- `MASS`: the Minkowski--Siegel mass formula and Conway--Sloane's use of
  neighbour closure plus exact mass for definite genera.
- `ELKIES-2026`: Elkies, arXiv:2608.25406, whose abstract announces a later
  construction paper.

Full bibliographic data are in
[`references/integral-rank-transfer.bib`](references/integral-rank-transfer.bib).

## Claim ledger

The “new axes” column separates mathematical (`math`), algorithmic (`alg`),
determinant/model-specific computational (`comp`), and certification/software
(`cert`) novelty.  A dash means no novelty claim.  Exact hypotheses are those
in the linked canonical statement; a source applies only under its own stated
characteristic, definiteness, primitivity, good-prime, spinor-genus, or
descent hypotheses.

| ID | Canonical term | Class | Strongest authority / antecedent | Project adaptation | New axes |
| --- | --- | --- | --- | --- | --- |
| `A` | Shioda--Tate rank balance under change of fibration | `ESTABLISHED` | `SHIODA-MW` | subtraction in the fixed-`NS` notation | -- |
| `A1` | fixed-`NS` rank-budget corollary | `TAILORED_COROLLARY` | `SHIODA-MW` | corridor accounting and equation-warning boundary | cert |
| `A2` | Galois-equivariant Shioda--Tate quotient identity | `TAILORED_COROLLARY` | `SHIODA-GAL`, `SHIODA-MW` | representation-ring subtraction for two marked fibrations | cert |
| `A2.1` | rational divisor-span promotion gate | `TAILORED_COROLLARY` | `SHIODA-GAL` | fail-closed source-first rank promotion | cert |
| `A2.2` | section fields from stabilizers | `ESTABLISHED` | `SHIODA-GAL` | arithmetic-marking output convention | cert |
| `B` | Shioda discriminant/regulator comparison | `ESTABLISHED` | `SHIODA-MW` | two-fibration ratio in fixed `NS` | -- |
| `B1` | determinant/index-square identity | `ESTABLISHED` | `NIKULIN`; elementary Gram determinant identity | saturation diagnostic | cert |
| `B2` | Hermite obstruction specialized to rootless K3 frames | `TAILORED_COROLLARY` | `SHIODA-MW` plus classical Hermite bounds | determinant-28 cutoff and E6 control | comp |
| `C` | primitive nef isotropic class gives a genus-one pencil | `ESTABLISHED` | standard K3 linear systems; `HOP`, Theorem 2.1 | marked section criterion | -- |
| `C1` | Weyl reduction to a nef chamber | `ESTABLISHED` | `HOP`, Remark 2.10 | fail-closed declared-wall boundary | cert |
| `C2` | finite horizontal-wall enumeration at fixed degree | `TAILORED_COROLLARY` | `HOP`; no explicit antecedent located for the displayed finite norm identity | exact finite gate used by compiler | alg, cert |
| `C2.1` | score only after physical Weyl reduction | `TAILORED_COROLLARY` | `HOP` and `C2`; no separate antecedent located | q104-to-q10 exact control | comp, cert |
| `C2.2` | old-zero coefficient-swap obstruction | `TAILORED_COROLLARY` | `HOP`; no explicit antecedent located for this coefficient-swap packaging | determinant-36 complete shell control | comp |
| `C3` | certified zero-changing neighbour loop | `TAILORED_COROLLARY` | `HOP` plus primitive-`U` algebra | exact route-specific cost mechanism | comp, cert |
| `D` | unimodular integral marking transport | `ESTABLISHED` | `NIKULIN`; elementary integral lattice algebra | mandatory lossless-marking gate | cert |
| `E` | Shioda--Tate specialization balance | `TAILORED_COROLLARY` | `SHIODA-MW` | separates Picard, root and MW jumps | cert |
| `E2` | non-thin rank jumps on double elliptic K3 surfaces | `ESTABLISHED` | `PASTEN-SALGADO`, Theorem 1.1 | full hypothesis audit for published R17 | comp, cert |
| `F0` | one- and two-section Tate-chart identities | `TAILORED_COROLLARY` | `TATE`; no explicit antecedent located for the displayed two-point Bézout packaging | exact compiler ansatz and scope gate | alg, cert |
| `F` | conditional correctness of an equation-level fibration hop | `TAILORED_COROLLARY` | `HOP`, Sections 2.3--2.4 | explicit list of sufficient certificate layers | cert |
| `F0b` | quadratic parent/coercion gate | `TAILORED_COROLLARY` | no explicit antecedent located; elementary polynomial algebra | prevents a documented coefficient-tower failure | cert |
| `F0c` | unordered incidence does not imply section descent | `TAILORED_COROLLARY` | `SHIODA-GAL` | exact E6 descent counterexample | comp, cert |
| `F1` | regular chord construction from a height-ten trace | `LIKELY_NEW_ALGORITHM` | no explicit antecedent located; compare `HOP` fibrewise trace | closed-form R17 bisection compiler | alg, comp |
| `F1.1` | height-eight genus-one bisection pencil | `LIKELY_NEW_ALGORITHM` | no explicit antecedent located | exact parameter/incidence/Kummer gate | alg, comp |
| `F1.2` | height-twelve regular quartic | `LIKELY_NEW_ALGORITHM` | no explicit antecedent located | complete 43-class R17 construction | alg, comp |
| `F2` | injectivity of the 39,120-class R17 extension map | `NEW_COMPUTATION` | `ELKIES-2026` is overlap context, not this census | complete equation-level squareclass certificate | comp, cert |
| `F2.1` | translation invariance of bisection visibility | `TAILORED_COROLLARY` | `SHIODA-MW`; no explicit antecedent located for this atlas boundary | bounds what translated shells can reveal | alg |
| `F3` | distinct R17 bisections give genus-one `V4` bases | `NEW_COMPUTATION` | no explicit antecedent located for this all-pairs certificate; classical Kummer theory and Riemann--Hurwitz | all-pairs exact application | comp |
| `F4` | multiquadratic character decomposition | `ESTABLISHED` | no special antecedent needed beyond rational character idempotents and quadratic twisting | R17 height/rank specialization | comp |
| `F5` | rootless multisections as coset minima | `TAILORED_COROLLARY` | `SHIODA-MW`; no explicit antecedent located for the exact CVP packaging | finite foundry coordinate | alg |
| `F6` | translation-quotient coset metric and degree scaling | `TAILORED_COROLLARY` | no explicit antecedent located; elementary lattice-torus algebra | invariant graph/hypergraph gate | alg |
| `H-1` | cross-Gram reconstruction of two primitive `U`-embeddings | `TAILORED_COROLLARY` | `HOP`, Lemmas 2.5--2.6; `NIKULIN` for saturation | 84-presentation replay; no identity-level novelty | comp, cert |
| `H-1a` | bounded relative-marking enumeration | `TAILORED_COROLLARY` | `HOP`, `NIKULIN`, and positive-definite shell finiteness | exact box completeness with coercivity warning | alg, cert |
| `H-1b` | one-edge elliptic incidence distance | `LIKELY_NEW_ALGORITHM` | no explicit optimization antecedent located; compare `HOP` | finite classifier and exact determinant-948 value two | alg, comp |
| `H` | common-core graph-glue decomposition | `TAILORED_COROLLARY` | `NIKULIN` | complete 42-edge rank-15/cyclic-rank-2 corpus | comp, cert |
| `H0` | low-norm coset-theta convolution | `TAILORED_COROLLARY` | `THETA`, `NIKULIN` | inverse zero-support enumerator | alg, cert |
| `H0b` | discriminant-form reconstruction of the core genus | `TAILORED_COROLLARY` | `NIKULIN`, Propositions 1.4.1 and 1.5.1 | forced-genus search and 84-presentation replay | alg, comp |
| `H0c` | bounded theta-decorated completion classifier | `LIKELY_NEW_ALGORITHM` | no explicit inverse-completion antecedent located | rootless-completion acceptance procedure | alg, cert |
| `H0d` | reverse low-norm coset-theta obstruction masks | `LIKELY_NEW_ALGORITHM` | no explicit antecedent located in `THETA` sources | mask antichain and lazy CVP | alg, comp, cert |
| `H0e` | zero-orbit Weil compression preserving mask cells | `LIKELY_NEW_ALGORITHM` | `THETA` for commuting actions; no explicit zero-slice antecedent located | exact four-module compression | alg, comp |
| `H0f` | nonselectivity of the linear modular mask sieve | `NEW_COMPUTATION` | `THETA` plus invariant Riemann--Roch | exact four-control negative result | comp |
| `H0g` | forced-genus mask-aware neighbour generation | `LIKELY_NEW_ALGORITHM` | `KNEIGH` plus `H0d` | exact Golay construction, not a complexity theorem | alg, comp, cert |
| `H0h` | prospective masked-core controls | `NEW_COMPUTATION` | no separate antecedent applies to this finite certificate; method from `H0g` | exact H3/NS0024 successes and Q80 near miss | comp |
| `H0i` | witness-resolved directed neighbour search | `LIKELY_NEW_ALGORITHM` | `KNEIGH`, equations (1.2)--(1.3) | physical-witness transition state and exact Q80 path | alg, comp, cert |
| `H0i.1` | coset-resolved `p`-neighbour transition formula | `LIKELY_NEW_ALGORITHM` | bare affine layers are likely implicit in `KNEIGH`/`THETA`; no transition-oracle antecedent located | predicts all forbidden births and deaths | alg, comp, cert |
| `H0i.2` | fixed-prime directed defect traps | `NEW_COMPUTATION` | `KNEIGH` supplies unrestricted neighbours | mass-closed ternary counterexamples | comp |
| `H0j` | NS0024 completed-core rank profile | `NEW_COMPUTATION` | `SHIODA-MW`, `NIKULIN`, `KN` | exact `4,12,12,17` frame-level path | comp, cert |
| `H0k` | marked metric reconstruction of completion roots | `TAILORED_COROLLARY` | standard ADE recovery plus `NIKULIN` | physical witness coordinates recover saturation/torsion | alg, cert |
| `H0l` | target root-system constraints for a `p`-neighbour | `LIKELY_NEW_ALGORITHM` | Chenevier's visible-root principle; no antecedent for all affine/glue births located | exact inverse ADE predicate and NS0024 control | alg, comp, cert |
| `H1` | good-prime neighbour as isotropic-line swap | `ESTABLISHED` | `KNEIGH` | decorated low-norm bookkeeping | cert |
| `H2` | Kneser--Nishiyama `J2` frame-genus realization | `ESTABLISHED` | `KN`, `NIKULIN` | rank-19 length specialization | -- |
| `H2a` | rank-three Hodge rigidity and finite `J1` bound | `TAILORED_COROLLARY` | `KN`, Proposition C' | exact determinant-948 bound `2..8` | comp |
| `H3` | conditional large-prime one-click existence | `TAILORED_COROLLARY` | `KNEIGH`, Chenevier's equidistribution theorem | rootless target specialization | -- |
| `H4` | mass-closing neighbour enumeration | `TAILORED_COROLLARY` | `MASS` plus `KNEIGH` | decorated fail-closed completion criterion | cert |
| `H5` | determinant-78 frame genus is globally rootful | `NEW_COMPUTATION` | `KN`, `NIKULIN`, `MASS` | complete determinant-specific obstruction | math, comp, cert |
| `I` | involution eigensublattices joined by 2-primary graph glue | `ESTABLISHED` | `NIKULIN` plus rational idempotents | character-saturation bookkeeping | cert |
| `I1` | one character-glue type on eleven rank-28 lifts | `NEW_COMPUTATION` | `NIKULIN` and the eigenspace argument in `I` | exact repeated pattern on distinct covers | comp |
| `I2` | exact visible glue on cover `0x103b2` | `NEW_COMPUTATION` | `NIKULIN` and the eigenspace argument in `I` | exact rank-18 visible lattice and specialization | comp, cert |
| `G` | completeness in a declared lattice box | `TAILORED_COROLLARY` | no explicit antecedent needed beyond positive-definite shell finiteness | explicit bounded-search proof boundary | cert |

## What is new, and what is not

Not new as general mathematics: rank/root balance, regulator comparison,
primitive-`U` changes, Kneser neighbours, Weyl reduction, graph gluing,
frame-genus classification, equation changes between elliptic fibrations,
vector-valued theta series, or Galois actions on Mordell--Weil lattices.

New exact computations in this repository include the 42-edge bridge corpus,
the two-class determinant-948 rootless `J2` census, the distance-two/33--10
accessibility computation, the directed Q80 defect path, the NS0024 rootless
frame and `4,12,12,17` profile, the determinant-78 global obstruction, the
39,120-class extension injectivity result, and the bounded relative-`U`
obstructions.  These remain genuine results even though their infrastructure
is classical.

Likely-new algorithmic contributions, stated conservatively, are reverse
low-norm theta masks and their antichain/lazy-CVP use, mask-aware forced-genus
generation, physical-witness-resolved directed neighbour dynamics, the
complete masked birth/death oracle, target root-system constraints including
affine/glue births, and one-edge incidence optimization inside `J2` classes.
For each, the present claim is only: **no explicit antecedent was located in
the sources checked**.

## Unresolved prior-art and publication risks

Before any priority statement or submission:

1. search MathSciNet, zbMATH, theses, and software literature for low-norm
   discriminant-coset masks in glue-code enumeration;
2. ask whether the affine dual-layer formula appears explicitly in
   lattice-coset Hecke or neighbour literature, and claim novelty only for
   its complete low-norm transition-oracle use;
3. search for target ADE/root-system constraints in neighbour algorithms;
4. search for minimization of fibre-intersection degree across `J2` copies;
5. recheck Elkies's announced construction sequel (`ELKIES-2026`);
6. obtain review from an elliptic-K3 expert and a lattice-neighbour/theta
   expert.  The review should receive this claim map, not an unfiltered
   repository dump.

## Integrity boundary

This reclassification changes terminology and novelty scope, not mathematical
status or checker outputs.  Historical theorem IDs and `MATH_STATUS.json`
claim IDs remain stable.  A status entry is changed only if its mathematical
scope changes; no such change is implied by recognizing classical
infrastructure.
