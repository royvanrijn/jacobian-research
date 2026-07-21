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
| C17 | Proved internally | Frozen example | Exact two-transfer local algebra |
| C19--C20 | Computational evidence | Frozen examples | Exact three- and four-transfer presentations only |
| C18, C21 | Conjectural / incomplete | Archived | Their unrecovered all-`k` input is not an active dependency |
| C22 | Conjectural / incomplete; proposed mechanisms refuted | Archived | Neither failed filtration may support an all-`k` theorem |
| C23 | Proved internally in its stated cubic scope | Scoped | Higher affine-factorization complements remain open |
| C24 | Proved construction with incomplete global classification | **Active** | Sole active constructive and classification programme |

## Strengthening completed

### C04

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
| Comparison with generic weighted seeds | Settled | Ramification separates `r>=2`; the reduced boundary intersection separates `r=1,m>1`, even stably; `(1,1)` is the known C01 identification |
| Tail deformations | Settled | Every allowed `A^(r+1)` tail is removed by a polynomial source automorphism |
| Generalized two-weight mechanism | Settled inside the stated ansatz | Spectral coprimality forces the monomial C24 branch; no nonmonomial polynomial branch survives |
| Parameter separability and discriminant | Proved uniformly | Closed discriminant formula, exact square criterion, complete even-degree square locus, and an infinite square family for each fixed `r` |
| Parameter irreducibility | Partially solved | Full `m=1` column; three uniform arithmetic criteria; every pair `mr<=30`; arbitrary remaining `(m,r)` open |
| Parameter Galois groups | Exact through `mr<=30` | Symmetric, alternating, dihedral, and exceptional degree-six actions occur; all-degree classification open |
| Distinct parameter-root equivalence | Open | Determine whether maps attached to conjugate or otherwise distinct roots `q` are polynomially left--right equivalent over an algebraic closure |
| Minimal collision fields | Open | Determine the least fields over which the full collisions are defined |
| Beyond the two-weight skeleton | Open | Find the smallest coordinate relaxation admitting a new cancellation family, or classify all polynomial cancellation branches |
| Independent verification | Pending | The all-`m,r` cancellation, normalization, and boundary package still needs a genuinely separate audit |

The boundary comparison is proved in
[BOUNDARY_INTERSECTION_OBSTRUCTION.md](experimental/transfer-and-cancellation/BOUNDARY_INTERSECTION_OBSTRUCTION.md).
The arithmetic results are collected in
[PARAMETER_IRREDUCIBILITY.md](experimental/transfer-and-cancellation/PARAMETER_IRREDUCIBILITY.md),
[PARAMETER_DISCRIMINANT.md](experimental/transfer-and-cancellation/PARAMETER_DISCRIMINANT.md),
and
[PARAMETER_GALOIS_GROUPS.md](experimental/transfer-and-cancellation/PARAMETER_GALOIS_GROUPS.md).

## Other open or review-level work

1. Independently audit the global compactification and boundary exhaustion in
   C05 and C16.
2. Record a primary-source hypothesis audit for the Mason--Stothers uses in
   C07 and C10.
3. Independently review C13's effective Chebotarev and twist argument; its
   present uniform constants are deliberately coarse.
4. Reproduce C12's complete tangent matrices and C14's full singular-radical
   and escape-path suite in an unrelated CAS.
5. Extend C23 from the cubic and first `(2,3)` obstruction to a uniform
   classification of higher affine factorization complements.
6. Determine the least exceptional exponent and minimal dimension in C15's
   stable normal-form consequences.
7. Add proof-assistant artifacts for the remaining all-degree results if
   formalization becomes a project goal.

## Archived transfer conclusion

The divided-power ribbon norm misses `X^2`, and the later
quadratic-remainder proposal misses `X^3` already at `k=2`. The latter
proposed special fiber is positive-dimensional for every `k>=2`, so it
cannot be the tangent cone of the finite low-degree transfer algebras. C18,
C21, and C22 therefore remain archived, and no C22-to-C18-to-C21 dependency
chain is active. See the exact
[counterexample](experimental/transfer-and-cancellation/QUADRATIC_REMAINDER_ALGEBRA.md)
and the preserved [archive](archive/transfer-all-k/).
