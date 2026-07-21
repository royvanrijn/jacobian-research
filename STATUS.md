# Mathematical status

This is the repository's single source of truth for mathematical completion
and external review.  “Proved” means that the repository contains a written
argument for the stated result.  Repository clean-room checks are named in the
external-review column but are not presented as external review.

| Result | Mathematical status | External review | Location |
|---|---|---|---|
| Foundational Keller map and rational collision | Proved | External formal certificate by Dean Cureton; local exact implementations | [Foundational geometry](verified/FOUNDATIONAL_GEOMETRY.md) |
| Exact symplectic cotangent lift and Weyl quantization | Proved generally for every Keller map; the foundational lift is an exact symplectic noninjective etale map of `A^6`, and its Weyl lift is an injective non-surjective endomorphism of `A_3`; stable-class lower bounds transfer unchanged | Independent explicit Weyl write-up by Omniscience Research Agent and Jeff Pickhardt; local exact checker; symplectic stable-class transfer not externally reviewed | [Classical-quantum bridge](extended-geometry/SYMPLECTIC_WEYL_LIFT.md) |
| Cubic marked-root realization | Proved | None recorded | [Marked-root model](verified/MARKED_ROOT_MODEL.md) |
| Exact cubic image, fibers, and nonproperness | Proved | None recorded | [Image theorem](verified/IMAGE_AND_NONPROPERNESS.md) |
| Weighted marked-root family and `S_n` monodromy | Proved | None recorded; clean-room repository proof and checker | [Weighted theorem](verified/WEIGHTED_SEED_THEOREM.md) |
| Generic discriminant geometry | Proved | None recorded | [Generic discriminant curve](extended-geometry/GENERIC_DISCRIMINANT_CURVE.md) |
| Weighted image and boundary theorems | Proved | None recorded | [Canonical family image](extended-geometry/CANONICAL_FAMILY_IMAGE.md) |
| Full-contact omission and uniqueness | Proved | None recorded | [Omitted-value classification](extended-geometry/OMITTED_VALUE_CLASSIFICATION.md) |
| Contact strata and dimensions | Proved | None recorded | [Coincident-root rebuild](extended-geometry/COINCIDENT_ROOT_REBUILD.md) |
| Contact-atom principle | Proved | None recorded | [Contact-atom principle](extended-geometry/CONTACT_ATOM_PRINCIPLE.md) |
| Exceptional components and closure order | Proved | None recorded | [Coincident-root rebuild](extended-geometry/COINCIDENT_ROOT_REBUILD.md) |
| Component normalizations | Proved | None recorded | [Component normalization](extended-geometry/COMPONENT_NORMALIZATION.md) |
| Degree-twelve local singularity | Proved | None recorded; targeted clean-room repository check of the decisive local algebra | [Degree-twelve singularity](extended-geometry/DEGREE12_LOCAL_SINGULARITY.md) |
| Effective finite-field Chebotarev law | Proved | None recorded | [Finite-field Chebotarev](extended-geometry/FINITE_FIELD_CHEBOTAREV.md) |
| Explicit quartic weighted model | Proved | None recorded; targeted dependency-free repository check of the central quartic algebra | [Quartic weighted geometry](extended-geometry/QUARTIC_WEIGHTED_GEOMETRY.md) |
| External quartic-island classification | Proved for Zhuang's pinned `F4a`, `F4b`, and `F4c`: weighted-seed resolvents, normalized boundary profiles, thick intersections, `S_4` monodromy, and exclusion from every cancellation normal form | Juntang Zhuang's external expanded-polynomial and collision certificates; independent local compact reconstruction and boundary audit | [Quartic islands](extended-geometry/EXTERNAL_QUARTIC_ISLANDS.md) |
| Stable normal-form consequences | Proved | None recorded; second local exact implementation | [Independent audit](extended-geometry/STABLE_NORMAL_FORM_AUDIT.md) |
| Dicritical compactification | Proved | None recorded | [Dicritical compactification](extended-geometry/DICRITICAL_COMPACTIFICATION.md) |
| Marked-point dimension barrier | Proved: in the universal marked-point/hyperplane-complement family, an affine-space source can occur only in dimension three; the entire two-dimensional family has a nonconstant unit or nontrivial Picard group | None recorded; exact divisor-localization and boundary-lattice proof | [Marked-point dimension barrier](extended-geometry/MARKED_POINT_DIMENSION_BARRIER.md) |
| Cancellation construction | Proved for every `m,r>=1` | None recorded; bounded local exact regressions | [Construction](cancellation/CONSTRUCTION.md) |
| Cancellation-parameter arithmetic | Partial: uniform discriminant and several irreducibility criteria; exact irreducibility and Galois groups for `mr<=30`; general classification open | None recorded; local exact certificate checks in the stated range | [Arithmetic](cancellation/ARITHMETIC.md) |
| Boundary distinction | Proved for every noncubic cancellation pair against generic weighted seeds | None recorded; the canonical paper has an explicit local-DVR boundary theorem, two independent all-parameter thick-contact derivations, and repository audits of both inputs | [Boundary geometry](cancellation/BOUNDARY_GEOMETRY.md) |
| Rigidity within the current ansatz | Proved under the explicitly stated monomial-triangular, one-inverse-variable axioms; not a universal classification | None recorded; bounded local exact regressions | [Rigidity](cancellation/RIGIDITY.md) |
| Degreewise stable-multiplicity theorem | Proved for every `N>=4`: at least `tau(N-1)` pairwise stably inequivalent maps of generic degree `N`, comprising one cancellation type per proper divisor of `N-1` plus one weighted type | None recorded; the [five-lemma audit](DEGREEWISE_MULTIPLICITY_AUDIT.md) is a repository verification companion, with [local-DVR exhaustion](papers/marked-root-multiplicity/boundary-exhaustion.tex), [stabilization of normalization and boundary data](papers/marked-root-multiplicity/stable-functoriality.tex), and [dual thick-contact proofs](papers/marked-root-multiplicity/thick-intersection.tex) isolated explicitly | [Canonical paper](papers/marked-root-multiplicity/main.tex) |
| Degree-five stable-moduli theorem | Proved on an explicit nonempty parameter open: every stable class meets the family in at most six parameters, hence there are uncountably many stable classes of generic degree five | None recorded; exact symbolic checker for the seed, Keller map, boundary exclusions, discriminant saturation, and orbit invariant | [Degree-five stable moduli](extended-geometry/DEGREE_FIVE_STABLE_MODULI.md) |

The cancellation family is deliberately split into construction, arithmetic,
boundary, and rigidity results.  The construction theorem does not assert complete arithmetic, unrestricted
parameter equivalence, or classification beyond the current reconstruction
skeleton.  Those boundaries are listed in
[the cancellation open problems](cancellation/OPEN_PROBLEMS.md).

The abandoned transfer programme is archived under
[archive/transfer-program](archive/transfer-program/) and is not part of the
active result ledger.