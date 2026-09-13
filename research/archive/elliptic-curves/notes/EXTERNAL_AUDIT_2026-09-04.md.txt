# External mathematical and machinery audit — 2026-09-04

Audit begun September 4; report completed September 5 (Europe/Amsterdam).

The principal record lower bounds survived the bounded replays. The audit
found a real quotient-linear-algebra bug, an unsafe Cassels–Tate certification
handoff, a missing rigorous primality step, a false CRT monotonicity statement,
and an unsupported generic index claim. Small corrections are applied below.
This is an audit of the current working tree, which already contained extensive
uncommitted work; it is not a clean-release or formal-verification certificate.

`../../MATH_STATUS.json` remains the mathematical authority. In particular,
curve 302 has certified rank **at least 31**, not proved exact rank 31;
rank at least 32 and the outstanding global Selmer computations remain open.

## Scope and evidence

The initial inventory covers **1,411 files** in the active elliptic programme,
its archive, and its generated-results directory, together with 70 directly
associated status entries. It includes the 74 active Markdown documents,
Python/Sage machinery, historical producers, saved certificates, and tests.
The audit also follows the K3 dependencies needed for specialization, height
geometry, generic MW16/MW17 subgroups, and class-group pressure.

The depth of examination is deliberately unequal: theorem and algorithm
arguments were examined in the central proof paths; historical searches and
large K3 enumerations received scope, dependency, and evidence triage. File
inventory, passing fixture tests, and source hashes do **not** mean that every
historical calculation was independently rerun or every line formally checked.

The [execution manifest](../../artifacts/generated-results/elliptic-curves/external_audit_20260904.json)
records commands, limits, final module outcomes, and unresolved failures. The
linked compressed evidence bundle preserves the initial file hashes, status
entries, complete command logs, and the bounded runner. No record search,
large descent, neighbour enumeration, or parameter campaign was launched.
Individual replay limits were 8–25 seconds; the longest successful recorded
execution took about six seconds. A timed-out test was killed with its process
group and was not counted as passed.

## Findings and corrections

### 1. High: quotient reduction could invent residual dimensions — fixed

In [residual_selmer_quotient.py](../cas/residual_selmer_quotient.py),
`EarlyQuotient.reduce` stopped when the highest occupied bit was not a known
pivot. Lower pivot bits then remained. For example, modulo the known span
of `01`, the vectors `10` and `11` represent the same class; the old routine
returned two independent ambient vectors and could report residual rank two
instead of one.

Reduction now eliminates **every** known pivot, in descending order. The
regression exhausts all two-generator known spaces and candidate pairs in
`F_2^3`, comparing projected rank with
`rank(known + candidates) - rank(known)`. All 4,096 cases pass.

No generated or archived certificate with the affected
`candidate_residual_rank` output field was found. This therefore identifies a
machinery error without evidence that a recorded theorem depended on its
incorrect output. Historical local outputs from that helper should be rerun
before reuse. Early stopping in a *membership* reducer is a different case:
it can correctly detect non-membership without producing a canonical quotient
representative.

### 2. High: a supplied alternating matrix could become an exact-rank claim — fixed handoff, arithmetic still required

[audit_residual_cassels_tate.py](../cas/audit_residual_cassels_tate.py)
previously accepted a truthy basis-certification value and an arbitrary
alternating matrix as enough for an exact-rank classification. It also emitted
an upper bound when the input basis was uncertified. Alternating linear
algebra alone does not show that the matrix is the arithmetic pairing.

The handoff now requires literal `True` attestations for the complete residual
Selmer basis, actual pairing entries, known point rank, and unconditionality,
plus a named algorithm and evidence. Otherwise the unconditional upper-bound
field is `null` and the result is explicitly uncertified. Tests include a
nondegenerate matrix with no entry certificate and a string `"false"`.

This remains an **attestation and linear-algebra checker**, not an independent
Cassels–Tate implementation. A production proof still needs a curve-bound
complete Selmer basis and replayable local pairing contributions. The numerical
formula is valid only under those premises; a JSON flag cannot establish them.
The global arithmetic bottleneck has not been solved by this repair.

### 3. High for exact factor claims: probable primes were labelled as proved — repaired by a supplement

The historical curve-302 and curve-356 producers use SymPy `isprime` when their
optional primality switch is enabled. Above `2^64`, that is a strong BPSW
probable-prime test, not a primality proof. See the
[SymPy documentation](https://docs.sympy.org/latest/modules/ntheory.html#sympy.ntheory.primetest.isprime).

The new [factor verifier](../cas/verify_record_prime_factors.py) reconstructs
all recorded `c4`, discriminant, and conductor factor products and calls
PARI `isprime` on all **35 distinct factors**. All pass. PARI distinguishes this
rigorous test from probable-prime factorization defaults in its
[FAQ](https://pari.math.u-bordeaux.fr/faq.html).

The [supplemental certificate](../../artifacts/generated-results/elliptic-curves/record_prime_factor_proofs_20260904.json)
is source-bound and replayable. The original producers and their historical
artifacts were preserved. Curve 302's pinned checker now invokes the rigorous
supplement, curve 356's reproduction commands include it, and both status
entries reference it. The independent point-rank arguments did not rely on
these large factors and were unaffected.

Some exploratory factorization scripts still use probable primes for scoring.
That is acceptable as heuristic triage, but fields such as
`factorization_complete` in `build_conductor_first_family_pilot.py` must not be
reused as exact primality evidence without a rigorous promotion step.

### 4. Medium: generic index 12,288 lacked a replayable proof — claim narrowed

The [hidden-section note](NEWFAMILY_HIDDEN_GENERIC_BASIS_AND_HALF_LATTICE.md)
asserted that the recovered generic sections enlarge the automatic subgroup
by index 12,288. The committed checker verifies eleven rational-function
curve identities; it does not replay the claimed relations to the automatic
sections. A specialized regulator ratio, specialized relations, and
interpolation samples do not prove those generic relations.

The generic index is now **UNKNOWN** pending an exact relation matrix,
coefficientwise group-law verification over `Q(T)`, and an integral inclusion
and determinant/Smith certificate. The specialized index calculation remains
historical evidence, not a replacement for this missing generic replay.
The eleven hidden sections and their rank-11 specialization, and the separate
exact rank-14 result at `T=83/6`, remain valid. The latter replay completed.

The rank-14 note also incorrectly implied that every later point must already
belong to the displayed subgroup. Exact rank proves rational-span containment
modulo torsion, not full saturation. That sentence is corrected.

### 5. Medium: the CRT monotonicity explanation was false — corrected

For a fixed CRT branch, adding congruences shrinks its feasible set. Its true
minimum admissible height therefore cannot decrease. What can change is the
**ordering of different branches**. A greedy beam can discard the eventual
best branch, and a bounded representative enumeration does not certify a
minimum in the first place.

The utility docstring, [CRT note](CRT_LATTICE_PIPELINE.md), `THEORY.md`, the
K3 sieve note, and the corresponding status scope now say this. The existing
160-choice counterexample still verifies the beam's false negative. Its
historical status identifier is retained rather than creating a duplicate
status narrative.

### 6. Medium: exact point addition could use floating-point division — fixed

Both `alternate_quartic_covers.short_add` and
`ecsearch.rank_certification.add_rational_points` accepted integer points but
could divide two Python integers with `/` for unequal abscissas. Coordinates
are now converted to `Fraction` before the slope calculation.

On `y^2=x^3-36x`, the regression checks exactly
`(-3,9) + (12,36) = (-144/25,-504/125)` through both implementations.
Existing rational-point certificate replays pass after the change.

### 7. Medium: the documented Fermigier scale missed a factor three — fixed

For the raw model `Y^2=X^3-27I*X-27J`, one has `c4=1296I`.
Combining that with the stated covariant normalization gives coordinate scale
`607392*u`, not `202464*u`, and discriminant multiplier `(607392*u)^12`.
The [reproduction note](FERMIGIER_REPRODUCTION.md) and family metadata are
corrected; exact invariant regressions check two rational parameters. The
canonical curve coefficients and recorded ranks did not change.

### 8. Medium: residual dimensions need the actual Kummer image

The reusable residual gate subtracts `rank(M)` from `dim Sel_2`. Calling this
the quotient by the known points assumes their Kummer image has that rank.
Rational independence alone is insufficient: `M=2E(Q)` has full free rank
and zero image modulo 2. The API documentation now states the precondition.
Current R17 applications supply separate full-image certificates.

The numerical difference remains a conservative rank allowance even without
that image premise. With rational 2-torsion of dimension `t`, the sharper
total rank bound is `dim Sel_2-t`; the current helper retains torsion and can
therefore miss a possible rejection. A future schema should distinguish
known free rank, certified Kummer-image dimension, and total rank allowance
explicitly rather than overloading the historical residual field.

### 9. Smaller correctness and reproducibility repairs

- Both production and legacy bounded-search validators now reject NaN and
  infinite limits. These had evaded ordinary positivity/comparison checks.
- The BNF-free explanation now distinguishes even valuations from a unit
  modulo squares: a nonprincipal square-root ideal can carry class-group
  2-torsion. It also correctly calls `-alpha^2` the negative of a square;
  its norm is negative in the odd-degree field.
- NumPy was absent from both the dependency declaration and the repository
  environment, breaking three test modules. Version 1.26.4, already used by
  the system interpreter, is now pinned and installed in the virtualenv.
- The unrun zero-gain rescue freeze depended on the budget validator's source
  hash. Its old bytes are [archived](../../archive/elliptic-curves/external-audit-2026-09-04/README.md).
  Only that implementation hash and the derived protocol hash were amended;
  population, treatment assignments, order, and budgets are identical. No
  worker was launched.
- The historical Fermigier generic-rank artifact retains its old finite-group
  helper hash. That source is archived as a replay input; every mathematical
  field was compared with a fresh replay under the corrected point-addition
  code and agrees exactly. Historical provenance was not silently rewritten.
- Stale baseline/navigation wording was corrected, and `STATUS.md` was
  regenerated through the status renderer.

## Mathematical assessment of the central machinery

The finite-reduction independence argument is sound: full column rank forces
every integral relation to be divisible by the chosen prime; excluding
rational torsion at that prime permits iteration. Failure to attain full rank
is inconclusive. The new
[Nagao replay](../scripts/verify_pinned_nagao_lower_bounds.py) checks the saved
point lists, torsion witnesses, and every saved finite-reduction row without
repeating discovery or saturation.

The two exact generic-rank-12 arguments correctly use the rational
Néron–Severi classes, good-reduction Frobenius bounds, and Shioda–Tate. Their
full small-field replays passed; a geometric upper allowance of 13 is not
mistaken for arithmetic generic rank 13.

The specialization sandwich, quotient Schur-complement metric, pointed-quartic
midpoint identity, and ideal-square-root class-group lower-bound argument are
consistent with their stated hypotheses. The fixed-cubic-field pilot takes a
kernel on the **whole** known span; it does not incorrectly discard individual
basis elements instead of testing their combinations. Its local completeness
argument uses the independently known local Kummer dimension. The full local
campaign was not rerun here.

The higher-Selmer note correctly uses images of `Sel_(2^j)` in `Sel_2`, not
literal inclusions between different Selmer groups. Known rational classes
survive every stage. A first Cassels–Tate drop has even dimension; surviving
classes need not come from rational points. The complete ambient Selmer space,
all bad places including 2 and infinity, global class/unit completeness, and
the actual pairing remain necessary proof inputs.

No audited bounded point-search miss supplies a rank upper bound. A finite
quotient escape proves escape from the current **subgroup**, which may be a
same-rank finite-index enlargement. Numerical height/LLL/CVP results remain
triage unless followed by exact identities or certified height bounds.

## Verification outcome and remaining reproduction failures

All **115 test modules** were attempted separately under wall limits. The
final outcomes are **113 passing modules / 615 tests**, one version-sensitive
failure, and one capped module. This is not a claim that the unrestricted
`make verify-elliptic-curves` target passes.

Successful arithmetic replays include:

| Result | Audit replay |
| --- | --- |
| ICARM 302, rank at least 31 | Both finite-quotient implementations, pinned artifact, rigorous factor supplement |
| ICARM 273, rank at least 30 | Independent Sage finite groups, minimality, conductor, root number, torsion |
| ICARM 398, rank at least 30 | Exact record checker including local arithmetic |
| ICARM 356, rank at least 29 | Exact record replay plus rigorous factor supplement |
| ICARM 285/286 and 394, rank at least 21 | Exact pinned point/local-conductor replays |
| Fermigier rank 22, rank-20 near miss; Kihara rank 14; E29 baseline | Exact point-independence replays and relevant existing checks |
| Nagao four rank-17 fibres and rank-20 fibre | Every saved finite-reduction row, point equation, and 2-torsion witness |
| Mestre generic rank at least 13 and the `u=197` rank-17 fibre | Exact construction/independence replays |
| Nagao section-7 and Fermigier generic rank exactly 12 | Full small-field Picard-bound replays |
| Six-root family at `T=83/6`, exact rank 14 | Exact eclib subgroup growth and PARI interval `[14,14]` |
| CRT calibration and Fermigier CRT seed | Existing exact bounded calibration replays |

The remaining reproduction issues are explicit:

1. `artifacts/generated-results/elliptic-curves/CATALOG.tsv` is missing.
   The advertised catalogue audit, and therefore the Make target's first
   command, fail before reaching the tests. Evidence categories were not
   invented to manufacture a passing catalogue.
2. The strict curve-273 JSON replay differs only in the recorded PARI version:
   local 2.15.4 versus producer 2.17.4. The Fermigier candidate-record test
   similarly differs in version metadata, embedded verifier output containing
   that version, and the derived digest. These old artifacts were not repinned.
   Use the original environment for byte reproduction, or implement an
   explicitly separate mathematical-payload comparison; do not silently
   discard arbitrary fields. The independent curve-273 replay passes.
3. `test_nagao_section7_linear_sections` exceeded its eight-second cap. It is
   **not replayed**, not failed mathematically. The separate Picard-bound
   replay passed and does not substitute for every test in that module.
4. All 450 inspected JSON/gzip files parsed. All 517 historical-manifest
   hashes and 52 snapshot hashes match. However, the archive audit also objects
   to the revived, different active `test_rank_jump_benchmark.py`; its
   relocation-only policy needs an explicit restoration record.

## Improvements worth doing next, without a larger search

First, finish the missing proof interfaces: export generic relation matrices;
make Kummer-image dimensions explicit; and require reproducible arithmetic
pairing evidence. Keep exact certificate replay separate from point discovery,
class-group completion, and numerical height estimation. The new Nagao replay
illustrates that separation.

Second, account explicitly for projective boundary points after quartic
reduction. `half_lattice_direct_reduction.py` processes the affine points from
`hyperellratpoints` but does not enumerate rational points at infinity of the
**reduced** model. A Möbius inverse can send those to finite points of the raw
quartic. There are at most two to check per nonsingular genus-one chart. This
is a possible search-completeness improvement, not evidence that a particular
record run missed a new independent point. PARI documents the returned points
as [affine rational points](https://pari.math.u-bordeaux.fr/dochtml/ref-stable/Elliptic_curves.html#hyperellratpoints).
Amend future frozen protocols before changing that search policy.

Third, retain branch diversity in CRT selection and exploit certified lower
bounds only when pruning is meant to be exhaustive. In higher descent, compute
on a basis of the residual complement and its kernels rather than enumerating
all nonzero classes. Keep local image dimensions, class-group lower bounds,
and full Selmer upper bounds distinct.

Finally, preserve outcome blindness and count independent arithmetic families,
not coordinate presentations, in detector evaluation. The quotient-height
tables correctly disclose their floating-point boundary; agreement at two
precisions is not an interval certificate. The promising engineering target
is coordinate compression and exact subgroup classification, while any
prospective claim of rank prediction still needs held-out arithmetic outcomes.
