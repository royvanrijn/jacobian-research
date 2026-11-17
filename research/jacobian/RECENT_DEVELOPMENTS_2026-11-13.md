# Recent Jacobian research and its relevance to this repository

Literature cutoff: **13 November 2026**. This expanded review replaces the
earlier four-paper preflight. The important additions are two November HC4
claims, Fonollosa's complete missing-divisor calculation, and Prellberg's
revised treatment of Roy van Rijn's homogeneous map. They change what is
worth auditing before further calculation. The mathematical-status authority
remains [MATH_STATUS.json](../MATH_STATUS.json).

The [48-record catalog](literature-2026-11-13/catalog.json) records authors,
versions, dates, review depth and repository relevance. The
[coverage record](literature-2026-11-13/coverage.json) pins five API snapshots
and access limitations. This is a broad literature triage, not an independent
verification of 48 proofs or a claim of bibliographic exhaustiveness.

## Current boundaries after the 14 November follow-ups

The [current problem review and literature supplement](OPEN_PROBLEMS_AND_LITERATURE_2026-11-14.md)
supersedes the continuation choices in the historical progression below.
The 48-record catalogue remains the dated literature snapshot.

- Degree-five HC4 is closed locally by reproduction, and the standard
  Meng--Yang polynomial-graph collision route is excluded in every degree.
- The degree-six gap-21 stratum with one entire ramified omitted prime
  of index two and residue degree one is excluded under its stated hypotheses.
  Residual r7 cases have additional ramification; actual B/H images, joint
  cover relations and the affine source remain the global gate.
- F2's two-face contact spectrum starts at507; affine target normalization
  leaves five forms. The nonlinear source coefficient gate remains open.
- Compatible quartic Wright presentations have pole height at least9.
  The current pole-(6,9) noncontact types are7+2 and5+4, with sixteen
  marked leading candidates. PF2RULED22 proves that all sixteen lift
  through the full first cusp jet, with affine116-dimensional spaces.
  The reduced next test is eight immutable high coefficients plus eleven
  matrix-membership conditions on46 boundary parameters; higher coefficients,
  actual finiteness and compatible completion remain unresolved.
- GAOF7FOUR2 completes the quadruple-root strata and excludes(4,3,2).
  Exactly912 normalized targets there have six-point fibres and118 have
  seven-point fibres. Only(3,2,2,2) remains for global three-point F7
  existence; GAOF7TRP1 excludes it on the completed triple-root pole charts.
- PF2D6E1 remains partial: the old normal-even ramification argument did
  not enumerate omitted unramified rows.

These statements follow the [live mathematical ledger](../MATH_STATUS.json)
and the canonical sources linked in the current review. They do not
promote the newly screened literature to proved repository claims.

## Constructive follow-up: Gao F4

The current Gao revision (4 November) explicitly leaves complete F4--F7
fibre classifications open. The [F4 follow-up](gao-f4-fibres/README.md) now
gives a written full complex classification, with exact checks: the possible
sizes are 0,1,2,3,5, and the previously uncounted second boundary piece
contributes two points precisely off an explicit quadratic locus. The
two-double-root stratum is an elliptic curve minus thirteen points; the
full nonproper locus has two irreducible hypersurface components. These
deductions are registered as GAOF4FIB1, GAOF4ELL1 and GAOF4NP1, without
external review or established novelty priority. F5, F7 and JC2 remain open;
the subsequent F6 classification is linked below.

A supplementary search found [Korea Superintelligence Labs, *Exact fibers
of the Keller counterexamples F4, F5 and F6*](https://ideosphere.ai/papers/keller-fibers).
Its inspected 11-page manuscript and KF1--KF16 list include selected rational
fibres, real-root constraints and an earlier restricted root correspondence.
They are credited in the follow-up; no complete F4 boundary or multiplicity
stratification was found there. This supplementary primary source is pinned
in the follow-up manifest and is additional to the earlier 48-record catalog.

**14 November scoped follow-up:** the [F6 classification](gao-f6-fibres/README.md)
now gives every complex fibre as well: precisely 0,1,2,3,4,6 points occur.
The two-double-root coefficient stratum is rational with thirteen punctures,
contrasting with the elliptic F4 stratum. The F6 nonproper locus again has
two irreducible hypersurface components. GAOF6FIB1/GAOF6NP1 record the
written proofs and exact checks. The latest arXiv page still identifies v2
of 4 November. A [same-agent F4 audit](gao-f4-fibres/AUDIT_2026-11-14.md)
checks its geometric gates; independent review and novelty confirmation
remain outstanding for both maps. F5, F7 and JC2 remain unresolved.

**Further 14 November F7 progress:** the [F7 follow-up](gao-f7-fibres/README.md)
counts every boundary fibre (one or two points), certifies the occurrence of
1,2,4,5,6,7,8,9,10,12-point fibres, and excludes11. Only three now remains
undecided after the surjectivity proof below; this is not a complete F7 stratum classification. The four/five/six
occurrences have exact nonsingular finite-field coefficient-factorization
certificates and written characteristic-zero lifting proofs. The full
nonproper locus has two irreducible components, with fifth-order intersection
on the target boundary. GAOF7BND1/GAOF7CARD1/GAOF7NP1 retain the scopes,
source credit, incomplete sixfold-root attempt and assurance limits.

**F7 surjectivity, 14 November:** [GAOF7SURJ1](gao-f7-surjectivity/README.md)
now proves every complex target is attained. Every hypothetical empty
eliminant reduces to Q2^3*T3^2 or Q4^3. Centering the quadratic reduces
the first coefficient system to two variables; exact rational identities
exclude both families. Singular generated the identities and FLINT replays
them independently, with SymPy reconstructing their necessary equations.
This closes the empty-fibre gate without proving a uniform gcd-degree bound.
Three-point fibres, full F7 strata, F5 and JC2 remain open. Surjectivity in
other Keller families was already reported; no first-example, external
review or novelty-priority claim is made.

## November claims that could change the HC4 programme

[Qiyue Tang, *Equivalence of the Four-Variable Hessian and Planar Jacobian
Conjectures*, Hexagon 2611.00156v1, 6 November](https://hexagonmath.org/pdf/2611.00156v1)
claims the converse JC2 implies HC4. Theorem 1.1 says that a noninvertible
four-variable constant-Hessian potential has, after a constant linear change,
the form

\[
 f=R(x,y)+zA(x,y)+wB(x,y),\qquad \det D(A,B)\in\mathbb C^\times.
\]

The proof runs through a translation direction for singular gradient images
and weighted leading forms. Its central structural argument has now received a
[complete written follow-up review](../hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md),
including Theorems 2.1/2.3, the De Poi classification hypotheses and
Lemma 3.2 through Proposition 4.1. Checking only the final cotangent
identity would miss that burden of proof.

**Quantitative consequence, now registered as HCTANG3:** a potential of
degree d gives A and B of degree at most d-1. Combining the reproduced
structural theorem with Moh's plane result through degree 100 gives HC4
through potential degree 101. This also excludes every v2 Meng--Yang graph
of degree at most 48 retaining the marked collision. Equivalence alone,
without the degree-preserving normal form, would not imply these bounds.

The [subsequent structural audit](../hessian/tang-hc4-audit/README.md) reproduces the
local Frobenius and phase-cofactor identities and proves a sharp denominator
extension of characteristic transport. It also gives explicit rational
examples showing why polynomiality cannot be dropped. HCTANG1 records those scoped additions. HCTANG2 now reproduces the full
structural theorem by the separate written audit, supported by exact
symbolic identities. Geometric classification and case exhaustion are
written arguments, not finite-checker outputs; no formal verification or
external review is claimed.

[David Niedbala Giraudin, *The four-dimensional Hessian conjecture for
potentials of degree at most five*, Hexagon 2611.00161v1,
6 November](https://hexagonmath.org/pdf/2611.00161v1) claims all complex
degree-at-most-five cases. Its reductions use Schur complements, partial
Legendre transforms, a constant direction and explicit inverses. Five
verifiers account for 247 author-reported checks, with an additional
dependency on the preceding ternary calculation. The paper explicitly says
it is unrefereed. The initial literature pass did not execute its programs;
the [subsequent degree-five audit](../hessian/hc4-degree-five-audit/README.md) did: all
247 main checks and all 29 ternary-dependency checks passed. Its written
normalization review accounts for the complete branch partition, with an
additional independent full 35-coefficient block exclusion.

This directly overlaps [OP-HC4-D5](../foundations/HC4_PROJECTIVE_POLAR_GEOMETRY.md),
including the binary and unary top forms left in our archived programme.
The [5 November ternary-top note](https://hexagonmath.org/pdf/2611.00134v1)
has narrower coverage and is a dependency of the next-day claim. The follow-up reproduction is now registered as HC4D5EXT1 and closes
OP-HC4-D5; HC4D5Z1 records the extra full-block exclusion. The external theorem
retains its original attribution. Full independent reimplementation, formal
verification and external review are not claimed. Neither paper closes
arbitrary-degree HC4.

**Hexagon identifiers are not arXiv identifiers.** The versioned primary PDFs
and source archives were accessible. The initial abstract-route failures
were later bypassed with a browser user agent; all source bytes are pinned
in the follow-up audit.

The later [all-degree polynomial-graph application](../hessian/hc4-polynomial-graphs/README.md)
goes beyond the quantitative transfer: an elementary fixed-slice and
polynomial-degree argument excludes every affine fibre two-plane in the
standard Meng--Yang graph chart. HCTANG2 therefore makes any unit-Hessian
graph gradient invertible in every degree. HC4MYGALL1 closes the complete
standard graph collision route; HC4MYGAP1 separately proves the no-plane
obstruction and an explicit nonrational constant-determinant formal model.
Arbitrary nonlinear reparameterizations and mixed/coisotropic mechanisms
remain separate.

The [dual-substitution extension](../hessian/hc4-triangular-graphs/README.md) then
closes one triangular change in either order with arbitrary base-dependent
coefficients and every base-independent plane Keller dual substitution,
including arbitrary polynomial dual automorphisms. HC4MYGT1 is the
elementary mixed-jet, two-divisor and directional-ideal theorem; HC4MYGT2
uses HCTANG2 for global gradient invertibility. A unit Jacobian of the
parameterization does not make the Hessian transform covariantly, and
marked collision transfer requires equality of the pulled-back covectors.
Longer alternating words remain a gate only when their coefficients depend
on the base variables; nonlinear base changes and other reductions remain.

A separate [relative-linear row analysis](../hessian/hc4-relative-row-orbit/README.md)
uses the retained HC5T1 normalizer to classify the full SL2 output orbit
and its stabilizer. HC4MYOR1 identifies the output rows with oriented
ordered generators of (x,y). HC4MYOR2 reduces the unrestricted Keller
question in that class to a singleton-zero plane map, with a concrete
properness obligation: its finite normalization must have no omitted
fibre over zero. A nonautomorphic candidate would have omitted scheme
length d-1 there. The known marked-root cubic supplies an exact source-A3
control showing why singleton affine fibres do not establish properness.
This does not exclude the unrestricted relative-linear class or improve
a plane degree frontier.

The [pointed-pencil follow-up](../hessian/hc4-pointed-pencil/README.md) refines that
missing length into the radial divisor budget sum(e_i*f_i)=d-1, and adds
an Euler budget relating nonzero fibre deficits to generic pencil genus
and punctures. HC4MYPP1 includes a complete rational radial slice of the
known F1 three-variable control. HC4MYPP2 exhibits every nonnegative
pencil genus inside the same oriented row orbit, even with a fixed
nonconstant Jacobian polynomial. These rows are not Keller maps. They
show why the rationality and irreducibility hypotheses in Nguyen Van
Chau's published pencil theorem cannot be obtained from ideal equality
or determinant-one row operations alone. No unrestricted Keller gate is
closed. Primary de Bondt and Gao version listings were rechecked and
remain v1 (17 October) and v2 (4 November), respectively.

The [radial-intersection continuation](../hessian/hc4-radial-intersection/README.md)
adds HC4MYRB1's pointwise equations Qe=-f and
ell_q=f^T*(-Q)^(-1)*f, the source canonical-order check, and an exact
cyclic-chain congruence. A resolved chain with all weights at most -2
cannot have only one radial horizontal of degree one. HC4MYRB2 applies
the test to the previously retained uv=w^n family: even n is already
radially resolved on the minimal surface resolution, whereas odd n
requires another blowup and a degree-one horizontal -1 curve. This
excludes the declared complete-chain stratum and prevents extending it
to an unresolved pencil. Local controls do not supply a global Keller
map or close the unrestricted properness problem.

The [global attachment audit](jc2-boundary-attachment/README.md) imports
Borisov's published plane Stein-boundary structure, distinguishing its
existential ramified-dicritical statement from an all-dicritical claim.
PF2BATT1 excludes every collided uv=w^n node projection from a plane
Keller finite normalization and proves that generically unramified omitted
curves have immersed residue-degree-one images. PF2R7BND1 then forces
two or three additional unramified boundary primes for either classified
normal-r7 cusp row, under its explicit single ramified-row hypothesis.
Their images differ from the cusp curve; a one-component nonproper-value
locus is excluded, but additional immersed images remain an open global
gate. The count comes from complete finite-normalization Euler strata,
without assuming a complete one-row missing-boundary ledger.

## Results closest to our current questions

The statements below describe each paper's asserted result. The catalog
distinguishes abstract screening from inspection of theorem statements and
proof sections.

| Source and latest inspected version | Relevant statement and consequence for our work |
| --- | --- |
| [Fonollosa, 2609.33491v1, 27 September](https://arxiv.org/html/2609.33491v1) | Sections 5.1–5.4 explicitly construct the finite normalization and every omitted divisor of marked-root maps. Besides one ramified divisor, there are omitted unramified divisors. This is directly relevant to our missing-boundary issue, in dimension at least three. |
| [Van Dobben de Bruyn, 2608.27341v2, 14 September](https://arxiv.org/html/2608.27341v2) | Recognizes affine-space complements for projective bundles over P1 when each boundary component restricts to a hyperplane on every fibre. The criterion requires two components, adjacent divisor classes and a unique coincidence fibre. Version 2 corrects the main theorem. |
| [Prellberg, 2608.12543v3, 15 September](https://arxiv.org/html/2608.12543v3) | Directly uses our 12-variable map and 19-variable homogenization. The displayed collision generates the full 19-dimensional polarization algebra, with growth 2,4,11,19,19. Therefore invariant linear restriction retaining that collision cannot improve the unchanged 38-variable symmetric lift. This is a limitation of that procedure, not global minimality. |
| [De Bondt, 2610.24889v1, 17 October](https://arxiv.org/abs/2610.24889v1) | Economic stabilization gives smaller homogeneous and vanishing-conjecture witnesses, including the 12/18/90 chain examined in [our exact audit](de-bondt-2610.24889/README.md). It cites and develops the repository's work. Adding variables does not give a plane descent. |
| [Wang and Yan, 2610.36495v1, 27 October](https://arxiv.org/html/2610.36495v1) | Claims invertibility of every three-variable Keller map of ordinary total degree at most four, plus tameness over C. This is distinct from geometric degree, which counts generic inverse sheets. Its concluding questions concern total degrees five and six. |
| [Adamus and Hajto, 2610.36736v1, 27 October](https://arxiv.org/abs/2610.36736v1) | Reduces the gradient differential criterion to membership of the Legendre potential in the gradient subring, or termination of a derivative tower. It expressly notes that doubling gives no independent Nambu conditions. Useful as a proof target, but does not supply termination. |
| [Kanel-Belov and Fedorov, 2610.33919v1, 24 October](https://arxiv.org/html/2610.33919v1) | JC2 implies a generalized first-Weyl Dixmier assertion. The converse requires existence of characteristic-zero graph-module lifts; an alternative criterion is the specified uniform Berest–Wilson boundedness. Uniqueness of a lift does not establish existence. |
| [Jelonek, 2607.20597v2, 25 October](https://arxiv.org/html/2607.20597v2) | Bounded-degree Jacobian-one automorphisms form a closed locus. An irreducible component either consists of automorphisms or generically contains counterexamples. This does not let us extend one computed braid witness to every r7 curve. |
| [Gao, 2608.00222v2, 4 November](https://arxiv.org/abs/2608.00222v2) | Develops tangent sweeps and higher-dimensional inequivalence. Complete fibre descriptions for examples F4–F7 remain open. Our explicit normalization methods might contribute there after identifying which cases existing work already covers. |
| [Peretz, 2610.02961v2, 26 October](https://arxiv.org/html/2610.02961v2) | Studies integral points on asymptotic paths for integer-coefficient plane maps with Jacobian one, and theta-function criteria. These hypotheses do not automatically hold for a complex candidate in our coefficient locus. |
| [Orevkov, *On germs of mappings C2 to C2*, 2026](https://www.math.univ-toulouse.fr/~orevkov/k-en.pdf) | Published addendum to [Kulikov's germ paper](https://arxiv.org/abs/2503.07099v1). Its initial hypotheses require a smooth source germ and bijectivity on every component over the branch curve. The local classification does not specify which unramified components the affine source omits. |

Fonollosa gives the most immediate lesson for
[PF2R7CUSP1](jc2-normal-r7/README.md): a ramification calculation can miss
deleted unramified sheets. His complete divisor ledger is an example to
learn from, not a theorem that a plane companion is either present or absent.
Van Dobben de Bruyn offers a possible global recognition tool, but only after
proving that the compactification and its boundary satisfy his hypotheses;
an arbitrary sequence of surface blowups need not be a projective bundle.

The [14 November application audit](affine-complements/README.md) makes
that condition explicit. A completion preserving the entire quartic affine
normalization has no two-omitted-prime realization in this class. The
one-prime case forces F1 with boundary squares3,5 and contact4. The boundary
pair exists; an actual finite quartic map would still need a contracted
infinity section and a base-point cluster satisfying the stated budgets.

Prellberg's obstruction concerns our older 19-variable map and a fixed
collision. It does not rule out de Bondt's different five-cubic-output
12-to-18 construction, nonlinear coordinate changes, or a different
collision. Any new compression claim should specify exactly which of those
choices it changes.

## Supporting results and distinctions

[Ni's quartic HC4 paper](https://arxiv.org/html/2608.14217v1), 14 August,
overlaps our registered quartic result HC4CQ1. Proposition 6.1 also handles
an arbitrary-degree P in the special form
P(x1,x2,x3)+x4 Q(x1,x2,x3)+a x4^2 with degree Q at most two.
[Liu's v2](https://arxiv.org/abs/2608.19112v2), 2 September, treats constant
pivots and degree-five potentials under real Lorentzian signature
hypotheses. Those hypotheses cannot be removed when comparing with complex
HC4.

[Zhang](https://arxiv.org/html/2609.10180v1), 9 September, studies differential
operator extensions and holonomicity in characteristic p.
[Borisov, Gabber and Vasiu](https://arxiv.org/html/2609.05746v1), 4 September,
give positive-characteristic plane examples and special factorization
criteria. Neither supplies the unrestricted characteristic-zero lifting
step needed above. The characteristic-two plane examples of
[Mondello](https://arxiv.org/abs/2608.02634v1) and
[Huq-Kuruvilla](https://arxiv.org/abs/2607.20968v1) likewise do not settle
complex JC2.

There are smaller witnesses for related conjectures, but their definitions
matter. [Meng and Yang](https://arxiv.org/abs/2607.22198v2) give a
five-variable constant-Hessian counterexample;
[Dvorsky](https://arxiv.org/abs/2608.07338v1) gives a five-variable generalized
vanishing counterexample using a third-order operator and an auxiliary
multiplier. Neither is a five-variable homogeneous quartic counterexample
for the Laplacian statement. [Long's Gaussian-moment witnesses](https://arxiv.org/abs/2607.18186v1)
are already explicit in three variables and were obtained separately from
the large Keller-reduction chain. His
[compact-group Mathieu paper](https://arxiv.org/abs/2609.16178v1) has a
31 July submission date despite its September identifier.

[Kaplan, Shpilka and Tal](https://arxiv.org/abs/2611.07106v1), 4 November,
study identity testing and reconstruction for sums of powers. At low outer
degree some algorithms reconstruct a polynomial in another representation;
that does not promise a smaller shared cube decomposition for our matrix.
The [lineability paper](https://arxiv.org/abs/2611.03572v1), 2 November,
concerns families of higher-dimensional counterexamples. Both are useful
background rather than a current plane exclusion.

## Bounds and proof claims that need separate review

The existing [plane provenance audit](../foundations/plane-jc/PROVENANCE.md)
pins [Helali's (72,108) exclusion](https://doi.org/10.5281/zenodo.21479814)
and its dependence on the GGHV reduction. This scan found no verified new
degree frontier beyond the repository's retained 125. Live access to the
Zenodo record failed, so no claim about a fresh deposited revision is made.
The retained PJ1 scope is a **conditional lower bound of 125 on a possible
counterexample's larger coordinate degree**; it is not a theorem that all
plane maps through degree 125 are invertible. The all-degree graph exclusion
uses no such frontier input.
The [2026 institutional record for Guccione, Guccione and Valqui's Newton-polygon paper](https://cris.pucp.edu.pe/entities/publication/4ac4fae4-49ed-421c-a246-4d9534888bd6)
also requires comparison with its older manuscript before treating the
journal publication date as a new mathematical advance.

Recent revisions are not automatically established results.
[Oda v40](https://arxiv.org/abs/1203.1691v40), revised 18 October, claims a
plane open-embedding theorem and JC2; it is not imported here.
[Charbonnel v5](https://arxiv.org/abs/2304.14675v5), revised 6 August, still
asserts all-dimensional invertibility, which is incompatible with the
explicit Keller collisions verified in this repository.
[Zheglov v5](https://arxiv.org/abs/2410.06959v5), revised 19 January, claims
ordinary first-Weyl Dixmier. Kanel-Belov and Fedorov explicitly say in their
October introduction that they could not fully understand that proof.
Ordinary DC1, generalized DC1 and the missing lifting hypothesis must remain
separate. [Moskowicz's prime extension-degree claim](https://arxiv.org/abs/2407.13795v1)
has no newer revision listed and has not been independently audited here.

## Next research gates after this lookup

1. **Use the completed structural audit before more HC4 computation.**
   HCTANG2 supplies the reduction; HCTANG3 gives the quantitative bounds.
   HC4MYGALL1 now excludes the standard marked-collision graph chart in
   **every degree**, replacing its former degree-48 continuation boundary. For a new
   higher-degree candidate, test the projective constant-covector scheme
   before developing a moving-pivot search. Emptiness certifies invertibility;
   nonemptiness does not certify a counterexample.
2. **Retain the completed degree-five closure.** HC4D5EXT1 accounts for
   every binary/unary/ternary stratum and all 276 source checks, with the
   independent HC4D5Z1 block exclusion. Further direct-quintic collision
   searches are superseded unless a specific proof defect is identified.
   The higher-degree Meng--Yang graph problem remains separate.
3. **Complete the r7 missing-boundary ledger.** Use the distinction between
   ramification and omission explicitly. Our two cusp classes and conductor
   calculations remain useful inputs, but the class-II exclusion stays
   conditional until the global ledger is proved.
4. **Choose a concrete contribution to higher-dimensional examples.** Test
   whether our normalization calculations resolve an unclassified Gao fibre,
   or apply collision-generated closure to the newer 18-variable map. These
   are proposed next checks, not results or novelty claims from this pass.

The initial literature review promoted no external theorem. The subsequent
degree-five reproduction now registers HC4D5EXT1 and closes OP-HC4-D5 in
MATH_STATUS.json. The later complete written Tang audit additionally
registers HCTANG2 and HCTANG3. HC4 and JC2 themselves remain unresolved.

## 14 November: F7 coefficient obstruction and elliptic incidence

The [F7 follow-up](gao-f7-three-point/PROOF.md) proves an exact obstruction
on the Delta=0 cubic-cube slice of a putative three-point fibre. All eight
three-point multiplicity patterns are covered by two factorization families.
At that stage both full families were unresolved. A 1397-term rational unit identity
is replayed in FLINT and the retained Fraction implementation. Separately,
the auxiliary triple-root incidence has an elliptic completion with fibres
4 I4 + 64 I1 + IV, Euler number84, chi(O)=7 and geometric genus6, proving
nonrationality. These are GAOF7THREE1 and GAOF7CUSP1; they do not close JC2
or the F7 three-point question. Gao's arXiv record was rechecked on
14 November and still lists v2 of 4 November as latest. The model uses
Fisher's binary-quartic invariants and Schuett--Shioda's elliptic-surface
formulas, with sources pinned in the follow-up. No independent end-to-end
review or novelty priority is claimed.

The subsequent [section-height calculation](gao-f7-three-point-structure/PROOF.md)
registers GAOF7MW1: geometric Mordell-Weil torsion is trivial; an explicit
Q(s)-section has Shioda height58/3 and is primitive even over C(s), since
every nonzero section has height at least28/3. The incidence surface has
Zariski-dense rational points, while the actual three-point fibre gate
remains open. Modular lift failures and finite-field subsystem reductions
are preserved without characteristic-zero conclusions.

The subsequent [complete cubic-cube proof](gao-f7-three-point-certificate/PROOF.md)
registers GAOF7CUBE1. No F7 eliminant is Q3^3*S3 for arbitrary monic cubics,
including repeated or shared roots. Thus sum floor(m_r/3)<=2, at most two
distinct roots are triple-or-higher, and each multiplicity is at most eight,
without a sharpness claim. This excludes three-point nonsimple lists
(9), (6,3), (3,3,3); the five lists (7,2), (5,4), (5,2,2), (4,3,2),
(3,2,2,2) were still open at that stage. A separate streaming checker passed761 checks in
724.77seconds:734 rational predecessor nodes after integer denominator
clearing, followed by a full402x402 integer coefficient matrix containing
every actual monomial including1. Its nonzero determinant modulo a proved
prime establishes full span over Q; this is stronger than a modular
unit-ideal report. Local F4 identities were motivated by the retained
[Groebner.jl research](https://arxiv.org/abs/2607.06372v3), while the
characteristic-zero conclusion rests on the separately replayed certificate.
The [original seven-two reduction](gao-f7-linear-cube/README.md) gave seven
necessary equations in four variables; its modular unit output and capped
rational lift gave no further exclusion. No independent end-to-end or
formal review, novelty priority or JC2 conclusion is claimed.

The subsequent [sharp maximum-five proof](gao-f7-linear-cube/HIGH_ROOT_PROOF.md)
registers GAOF7MULT1: every F7 eliminant root has multiplicity at most five,
and the known exact Y0 polynomial t^5*(t^7-2t^4+2t-1) attains it. It removes
(7,2), leaving four three-point lists:(5,4),(5,2,2),(4,3,2),(3,2,2,2).
The738473byte certificate passed2426 checks in28.94seconds. A different
representation cancels the shared quadratic coefficient in b, then uses
three univariate resultants, verified by1157 exact integer evaluations
with both FLINT determinants and plain Python Bareiss. A full actual-degree
623x623 integer Sylvester determinant proves the final coprimality over C;
the exceptional loci s0 and3s^3=2 have separately replayed34/355-term unit
identities. Nominal degrees are retained at every specialization, so no
leading-coefficient locus is lost. The old raw lifts remain inconclusive.
No independent end-to-end review, formal verification, novelty priority or
JC2 conclusion is claimed. The three relevant arXiv version records were
rechecked on14November and remain Gao v2, de Bondt v1 and Groebner.jl v3.

The next [complete fivefold-root classification](gao-f7-quintuple/PROOF.md)
registers GAOF7FIVE1: exactly118 distinct normalized F7 targets have a root
of multiplicity at least five, and each has a single fivefold root and
seven other simple roots. They are the known origin plus117 targets
parametrized by a degree117 integer polynomial and rational recovery
formulas. On c!=0 they give118 disjoint C* curves of seven-point fibres.
At that stage this excluded(5,4),(5,2,2), leaving(4,3,2),(3,2,2,2).
GAOF7FOUR2 subsequently excludes the former, as recorded below. The5490-check replay passed in10.15seconds. Exact source coverage,
integer resultant/cofactor identities, eight finite FREE good-reduction
unit certificates, all exceptional cases and forward homogeneous source
identities prove both exclusion and the exact count. Finite modular unit
reports alone would not justify the conclusion. The complete geometric
point set is classified; no scheme-structure, independent end-to-end review,
formal proof, novelty priority or JC2 conclusion is claimed. Gao's earlier
axis fibre and isolated-quintuple existence retain source credit.

The next [quadruple-root incidence theorem](gao-f7-quadruple/PROOF.md)
registers GAOF7FOUR1. The marked-root projection has finite fibres at every
s-value and generic degree8; generic targets in this locus have ordinary
eight-point F7 fibres. A full54x54 integer coefficient resultant excludes
vertical components. A checked rank8 fibre at s1 proves nonzero GENERIC
norms, whose zeros/poles form a finite set; it does not prove units at all
complex s. Thus the nonordinary normalized quadruple-root locus, including
the then-undecided(4,3,2) case, is finite. This predecessor theorem did
not compute its complete finite projector or list. All24 points on s*(3s^3-2)=0 are classified:
one seven-point origin and23 eight-point quadruple targets. The1535-check
replay passed in10.30seconds. Bounded global lifts and two partial DAG
pilots were inconclusive. This predecessor left both(4,3,2) and(3,2,2,2);
the complete follow-up below supersedes that continuation boundary.
Full F7 strata,F5 and JC2 remain open, with no independent end-to-end or
formal review, external review or novelty-priority claim.


The [complete quadruple-root classification](gao-f7-quadruple/COMPLETE_PROOF.md)
now registers GAOF7FOUR2: its only multiplicity types are4+1^8,
4+2+1^6 and5+1^7. The nonordinary normalized locus has912 six-point
and118 seven-point targets. Thus(4,3,2) is excluded globally.
The [triple-root pole-chart classification](gao-f7-triple-doubles/README.md),
GAOF7TRP1, additionally excludes(3,2,2,2) on those charts and counts288
seven-point and23 eight-point targets there; the generic fibres have nine
points. The global D0!=0 case of(3,2,2,2) is the remaining three-point
question. These are scoped exact classifications, not a full F7
stratification or a JC2 result.
