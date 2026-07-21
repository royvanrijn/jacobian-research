# Research status

Status reconciled after the C04, C12, C14, and C24 strengthening pass. The
evidence labels below follow the canonical [claim ledger](CLAIMS.md); workflow
state is separate from proof status.

## Programme overview

| Claims | Evidence state | Workflow state | Current meaning |
|---|---|---|---|
| C01--C04 | C01 and C04 independently verified; C02--C03 proved internally | Core write-up | Stable theorem chain; no theorem extension planned |
| C05--C11 | Proved internally | Awaiting review | Uniform geometry is frozen pending specialist audit |
| C12 | Proved internally, decisive local algebra independently checked | Awaiting review | Full unrelated-CAS tangent-matrix reproduction remains desirable |
| C13 | Proved internally | Awaiting review | Independent arithmetic-geometry review remains desirable |
| C14 | Proved internally, central quartic algebra independently checked | Awaiting review | Full unrelated-CAS radical/path reproduction remains desirable |
| C15--C16 | Proved internally | Awaiting review | C15 minimization is open; C16 needs an independent compactification audit |
| Archived transfer programme | Historical only | Archived | Finite examples, failed filtrations, and scoped obstruction results are not active claims or dependencies |
| C24 | Proved construction with incomplete global classification | **Active** | Sole active constructive and classification programme |

## Strengthening completed

### C04

C01 now also has a
[positive invariance regression](verified/C01_INVARIANCE_REGRESSION.md):
nonlinear source/target conjugation, a different primitive element, a
Möbius root chart, and stabilization all transport the same canonical
boundary object.  The primitive-element change deliberately creates an
extra square in the raw resolvent discriminant, confirming that the test
detects and discards presentation artifacts.

C04 is now **verified independently**. The
[clean-room audit](verified/C04_INDEPENDENT_AUDIT.md) rederives polynomiality,
constant Jacobian, marked-root reconstruction, and normalization, and gives
an alternative vertical-line branch-cycle proof of `S_n`. Its checker uses
only standard-library sparse-polynomial arithmetic and does not import the
project weighted-model implementation.

### C12 and C14

C12 and C14 are now **proved internally**, rather than computational evidence.
For C12, an elementary triangular ideal argument replaces the decisive
sixfold-block Groebner step and independently recovers the two dual-number
blocks and transverse length four. For C14, a dependency-free implementation
reconstructs the quartic map, determinant, collision, incidence,
discriminant, repeated-root normalization, boundary relation, and omitted
node. These targeted audits do not yet justify labeling either entire package
independently verified.

## C24: exact current boundary

| Subproblem | Status | Result or remaining obligation |
|---|---|---|
| Polynomial cancellation construction | Proved internally for all `m,r>=1` | Constant Jacobian, finite cancellation operator, Hensel recurrence, reconstruction, collision, inverse degree, and monodromy |
| Comparison with generic weighted seeds | Settled uniformly | For every noncubic C24 member, the full intersection of its two divisorial boundary images is nonreduced, while the generic weighted intersection is reduced; this gives one stable proof for all `r`; `(1,1)` is C01 |
| Canonical boundary-incidence invariant | Proved internally with hypothesis audit | Every quasi-finite map in the stated normal affine setup has stable diagrams of both scheme-theoretic and reduced intersections of divisorial boundary images; Keller maps satisfy the setup automatically |
| Boundary-invariant ladder | Proved abstractly; C04/C24 reduced target layer explicit | The stable hierarchy `I^formal -> I^sch -> I^red` records successively reduced incidence, full intersection schemes, and completed maps/conductors/differents/valuation filtrations; C24 intersection strata upstairs remain to be computed |
| Boundary exhaustion and thick intersections | Exhaustion independently audited; thick intersection proved internally | A clean-room normalization and two-chart audit confirms the C04/C24 boundary lists are exhaustive; for every noncubic `(m,r)`, the full canonical C24 intersection has nilradical index `mr(m+1)`; the generic weighted intersection is reduced |
| Tail deformations | Settled | Every allowed `A^(r+1)` tail is removed by a polynomial source automorphism |
| Generalized two-weight mechanism | Settled inside the stated ansatz | Spectral coprimality forces the monomial C24 branch; no nonmonomial polynomial branch survives |
| Parameter separability and discriminant | Proved uniformly | Closed discriminant formula, exact square criterion, complete even-degree square locus, and an infinite square family for each fixed `r` |
| Parameter irreducibility | Partially solved | Full `m=1` column; three uniform arithmetic criteria; every pair `mr<=30`; arbitrary remaining `(m,r)` open |
| Parameter Galois groups | Exact through `mr<=30` | Symmetric, alternating, dihedral, and exceptional degree-six actions occur; all-degree classification open |
| Distinct parameter-root equivalence | Partially settled | Deck rigidity settles the target-fixed case. Every labelled target automorphism restricts on `P=0` to weighted scaling, which fixes all parameter roots; unrestricted equivalence is reduced to the cover-lifting congruence kernel acting trivially on the boundary plane |
| Minimal collision fields | Open | Determine the least fields over which the full collisions are defined |
| Universal three-weight monomial-triangular class | Settled relative to explicit axioms | All integer weights, arbitrary `f,g`, and every one-factor polynomial derivative are exhausted: the Jacobian forces `c=a-1`, arbitrary target-dependent derivatives then collapse, and polynomiality yields only C24. An assumption/failure-mode table identifies the exact boundary of the theorem |
| Finite factorized-resolvent relaxation | Settled | Every normalized factor `1-tf_i(Q-Pt)` is forced either to equal the original factor or to be identically one; all nontrivial factors coalesce and polynomiality again yields only C24 |
| Arbitrary target-dependent polynomial derivative | Settled inside the skeleton | Algebraic independence of `(s,P,Q)` forces `H(T,P,Q)=lambda(1-Tf(Q-PT))^n`; polynomiality again yields only C24 |
| Beyond the completed skeleton | Open | Change the reconstruction skeleton itself, for example with an additional source function, source variable, or inverse variable |
| Independent verification | Partial | Boundary normalization and exhaustion now have a separate clean-room audit; polynomial cancellation, monodromy, and the thick scheme-intersection formula still need genuinely separate all-parameter audits |

The boundary comparison is proved in
[BOUNDARY_INTERSECTION_OBSTRUCTION.md](experimental/transfer-and-cancellation/BOUNDARY_INTERSECTION_OBSTRUCTION.md).
Its reusable invariant-theoretic form is proved in
[CANONICAL_BOUNDARY_INTERSECTION_INVARIANT.md](experimental/transfer-and-cancellation/CANONICAL_BOUNDARY_INTERSECTION_INVARIANT.md).
The arithmetic results are collected in
[PARAMETER_IRREDUCIBILITY.md](experimental/transfer-and-cancellation/PARAMETER_IRREDUCIBILITY.md),
[PARAMETER_DISCRIMINANT.md](experimental/transfer-and-cancellation/PARAMETER_DISCRIMINANT.md),
and
[PARAMETER_GALOIS_GROUPS.md](experimental/transfer-and-cancellation/PARAMETER_GALOIS_GROUPS.md).
The first parameter-equivalence theorem is
[TARGET_FIXED_PARAMETER_RIGIDITY.md](experimental/transfer-and-cancellation/TARGET_FIXED_PARAMETER_RIGIDITY.md).

## Other open or review-level work

1. Independently audit the remaining all-parameter C24 inputs beyond the now
   clean-room-verified boundary exhaustion: polynomial cancellation,
   monodromy, and the thick scheme-intersection formula.
2. Record a primary-source hypothesis audit for the Mason--Stothers uses in
   C07 and C10.
3. Independently review C13's effective Chebotarev and twist argument; its
   present uniform constants are deliberately coarse.
4. Reproduce C12's complete tangent matrices and C14's full singular-radical
   and escape-path suite in an unrelated CAS.
5. Determine the least exceptional exponent and minimal dimension in C15's
   stable normal-form consequences.
6. Add proof-assistant artifacts for the remaining all-degree results if
   formalization becomes a project goal.

## Archived transfer programme

The former transfer route, including its finite examples, conditional
all-`k` arguments, two refuted filtrations, and scoped cubic obstruction, is
preserved under [archive/transfer-program](archive/transfer-program/).  It is
not an active claim family or a dependency of C24.
