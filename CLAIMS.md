# Claim ledger

Only four evidence statuses are used:

- **Verified independently** — reproduced by a genuinely separate proof or
  implementation.
- **Proved internally** — a complete written proof is present in this
  repository, without an independent external audit.
- **Computational evidence** — exact or bounded computation supports the
  statement, but the repository does not claim a complete general proof.
- **Conjectural / incomplete** — a material classification or proof step is
  unresolved.

| Claim | Description | Status | Main proof | Independent check |
|---|---|---|---|---|
| C01 | Explicit Keller map and rational collision | Verified independently | [Foundation](verified/FOUNDATIONAL_GEOMETRY.md) | [Lean and exact implementations](verified/LEAN_C01.md) |
| C02 | Cubic marked-root realization | Proved internally | [Marked-root model](verified/MARKED_ROOT_MODEL.md) | Exact two-chart scripts |
| C03 | Exact cubic image, fibers, and nonproperness | Proved internally | [Image theorem](verified/IMAGE_AND_NONPROPERNESS.md) | Exact boundary scripts |
| C04 | Weighted marked-root family and `S_n` monodromy | Proved internally | [Weighted theorem](verified/WEIGHTED_SEED_THEOREM.md) | Partial external overlap only |
| C05 | Generic discriminant geometry | Proved internally | [Proof](experimental/geometry/GENERIC_DISCRIMINANT_CURVE.md) | None |
| C06 | Weighted image and boundary theorems | Proved internally | [Canonical family](experimental/geometry/CANONICAL_FAMILY_IMAGE.md) | None |
| C07 | Full-contact omission and uniqueness | Proved internally | [Classification](experimental/geometry/OMITTED_VALUE_CLASSIFICATION.md) | None |
| C08 | Contact strata and dimensions | Proved internally | [Rebuild](experimental/geometry/COINCIDENT_ROOT_REBUILD.md) | None |
| C09 | Contact-atom principle | Proved internally | [Proof](experimental/geometry/CONTACT_ATOM_PRINCIPLE.md) | None |
| C10 | Exceptional components and closure order | Proved internally | [Rebuild](experimental/geometry/COINCIDENT_ROOT_REBUILD.md) | None |
| C11 | Component normalizations | Proved internally | [Proof](experimental/geometry/COMPONENT_NORMALIZATION.md) | None |
| C12 | Degree-twelve local singularity | Computational evidence | [Argument](experimental/geometry/DEGREE12_LOCAL_SINGULARITY.md) | Exact local-algebra computation |
| C13 | Effective finite-field Chebotarev law | Proved internally | [Proof](experimental/geometry/FINITE_FIELD_CHEBOTAREV.md) | None |
| C14 | Explicit quartic weighted model | Computational evidence | [Model](experimental/geometry/QUARTIC_WEIGHTED_GEOMETRY.md) | Exact script suite |
| C15 | Stable normal-form consequences | Proved internally | [Audit](experimental/geometry/C15_INDEPENDENT_AUDIT.md) | Second exact implementation |
| C16 | Dicritical compactification | Proved internally | [Proof](experimental/geometry/DICRITICAL_COMPACTIFICATION.md) | None |
| C17 | Two-transfer block | Proved internally | [Transfer blocks](experimental/transfer-and-cancellation/TRANSFER_BLOCKS.md) | Exact Groebner calculation |
| C18 | Archived global allocation equalizer | Conjectural / incomplete | [Archived conditional argument](archive/transfer-all-k/ALLOCATION_BRANCH_INTERSECTIONS.md) | C22 input not recovered |
| C19 | Three-transfer block | Computational evidence | [Transfer blocks](experimental/transfer-and-cancellation/TRANSFER_BLOCKS.md) | Exact Groebner calculation |
| C20 | Four-transfer block | Computational evidence | [Transfer blocks](experimental/transfer-and-cancellation/TRANSFER_BLOCKS.md) | Exact Groebner calculation |
| C21 | Archived master quotient theorem | Conjectural / incomplete | [Archived conditional argument](archive/transfer-all-k/MASTER_QUOTIENT_THEOREM.md) | C22 input not recovered |
| C22 | All-`k` transfer block; proposed filtration refuted | Conjectural / incomplete | [Quadratic-remainder counterexample](experimental/transfer-and-cancellation/QUADRATIC_REMAINDER_ALGEBRA.md) | Exact low-`k` examples only |
| C23 | Cubic factorization-slice classification and `(2,3)` obstruction | Proved internally | [Cubic classification and exact obstruction](experimental/transfer-and-cancellation/UNIVERSAL_FACTORIZATION_GEOMETRY.md) | Exact finite-field checks; general higher-degree classification remains open |
| C24 | Master cancellation family and left-right classification | Conjectural / incomplete | [Proved construction](experimental/transfer-and-cancellation/MASTER_CANCELLATION_CONSTRUCTION.md) | Small exact cases only |

Workflow state is tracked separately in [RESEARCH_STATUS.md](RESEARCH_STATUS.md).
