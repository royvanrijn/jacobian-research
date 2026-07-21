# Mathematical status

This is the repository's single source of truth for mathematical completion
and external review.  “Proved” means that the repository contains a written
argument for the stated result.  Repository clean-room checks are named in the
external-review column but are not presented as external review.

| Result | Mathematical status | External review | Location |
|---|---|---|---|
| C01 — explicit Keller map and rational collision | Proved | External formal certificate by Dean Cureton; local exact implementations | [Foundational geometry](verified/FOUNDATIONAL_GEOMETRY.md) |
| C02 — cubic marked-root realization | Proved | None recorded | [Marked-root model](verified/MARKED_ROOT_MODEL.md) |
| C03 — exact cubic image, fibers, and nonproperness | Proved | None recorded | [Image theorem](verified/IMAGE_AND_NONPROPERNESS.md) |
| C04 — weighted marked-root family and `S_n` monodromy | Proved | None recorded; clean-room repository proof and checker | [Weighted theorem](verified/WEIGHTED_SEED_THEOREM.md) |
| C05 — generic discriminant geometry | Proved | None recorded | [Generic discriminant curve](experimental/geometry/GENERIC_DISCRIMINANT_CURVE.md) |
| C06 — weighted image and boundary theorems | Proved | None recorded | [Canonical family image](experimental/geometry/CANONICAL_FAMILY_IMAGE.md) |
| C07 — full-contact omission and uniqueness | Proved | None recorded | [Omitted-value classification](experimental/geometry/OMITTED_VALUE_CLASSIFICATION.md) |
| C08 — contact strata and dimensions | Proved | None recorded | [Coincident-root rebuild](experimental/geometry/COINCIDENT_ROOT_REBUILD.md) |
| C09 — contact-atom principle | Proved | None recorded | [Contact-atom principle](experimental/geometry/CONTACT_ATOM_PRINCIPLE.md) |
| C10 — exceptional components and closure order | Proved | None recorded | [Coincident-root rebuild](experimental/geometry/COINCIDENT_ROOT_REBUILD.md) |
| C11 — component normalizations | Proved | None recorded | [Component normalization](experimental/geometry/COMPONENT_NORMALIZATION.md) |
| C12 — degree-twelve local singularity | Proved | None recorded; targeted clean-room repository check of the decisive local algebra | [Degree-twelve singularity](experimental/geometry/DEGREE12_LOCAL_SINGULARITY.md) |
| C13 — effective finite-field Chebotarev law | Proved | None recorded | [Finite-field Chebotarev](experimental/geometry/FINITE_FIELD_CHEBOTAREV.md) |
| C14 — explicit quartic weighted model | Proved | None recorded; targeted dependency-free repository check of the central quartic algebra | [Quartic weighted geometry](experimental/geometry/QUARTIC_WEIGHTED_GEOMETRY.md) |
| C15 — stable normal-form consequences | Proved | None recorded; second local exact implementation | [Independent audit](experimental/geometry/C15_INDEPENDENT_AUDIT.md) |
| C16 — dicritical compactification | Proved | None recorded | [Dicritical compactification](experimental/geometry/DICRITICAL_COMPACTIFICATION.md) |
| Cancellation construction | Proved for every `m,r>=1` | None recorded; bounded local exact regressions | [Construction](experimental/cancellation/CONSTRUCTION.md) |
| Cancellation-parameter arithmetic | Partial: uniform discriminant and several irreducibility criteria; exact irreducibility and Galois groups for `mr<=30`; general classification open | None recorded; local exact certificate checks in the stated range | [Arithmetic](experimental/cancellation/ARITHMETIC.md) |
| Boundary distinction | Proved for every noncubic cancellation pair against generic weighted seeds | None recorded; boundary exhaustion has an independent repository audit | [Boundary geometry](experimental/cancellation/BOUNDARY_GEOMETRY.md) |
| Rigidity within the current ansatz | Proved under the explicitly stated monomial-triangular, one-inverse-variable axioms; not a universal classification | None recorded; bounded local exact regressions | [Rigidity](experimental/cancellation/RIGIDITY.md) |

The historical label **C24** is only a family name for the cancellation maps.
It is not one result and does not assert complete arithmetic, unrestricted
parameter equivalence, or classification beyond the current reconstruction
skeleton.  Those boundaries are listed in
[the cancellation open problems](experimental/cancellation/OPEN_PROBLEMS.md).

The transfer programme formerly numbered C17--C23 is archived under
[archive/transfer-program](archive/transfer-program/) and is not part of the
active result ledger.
