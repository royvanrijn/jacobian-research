# Adversarial proof audit for the generic discriminant theorem

## Verdict

The generic discriminant result survives as a standalone theorem over an algebraically closed field of
characteristic zero, provided “discriminant curve” means the reduced
projective closure of the repeated-root locus. The proof does not need the
later weighted-image through finite-field results
or any assertion about the image or boundary of the motivating Jacobian map.

The audit led to two proof-hardening changes in the standalone paper:

1. The codimension argument now uses the dimension of the closure of a
   finite-type image directly. It does not depend on the more elaborate
   homogeneous graph presentation in the research note.
2. The infinity argument now proves a relative formal-inverse lemma. This
   supplies the missing justification for the claim that an off-diagonal
   equal-image pair cannot escape to infinity in an exact-degree family.

## Audit matrix

| Risk | Test | Result |
|---|---|---|
| Singularities at infinity | Homogenize the normalization; inspect its sole point over `Z=0`; complete the relative map there; inspect the relative equal-image fiber product. | Pass. In the `T=1` chart, `S/T=(n/(n-1))q+O(q^2)` has a formal inverse. The image is smooth, has one branch, and the off-diagonal equal-image locus has closure disjoint from the infinity section on `h_n != 0`. |
| Tacnodes and higher contacts | Classify one-branch contacts by the order of `H-ell_r`; compare tangent directions at two normalization points. | Pass. Contact 3 is an ordinary cusp; contact at least 4 lies in `(4)`. Two unramified branches at `r != u` have directions `(1,r)` and `(1,u)` with determinant `u-r`, so a tacnode is impossible. |
| Transversality of tangent–chord coincidences | Differentiate the equal-image equations in the two normalization parameters. | Pass. The Jacobian columns are the two branch tangents, so its determinant is a nonzero scalar times `u-r`. Every retained bitangent solution is transverse; nontransverse endpoints are exactly cusp/diagonal bad strata. |
| Codimension of bad strata | Factor `H-ell=M_mu Q` for `(4)`, `(3,2)`, and `(2,2,2)` and count the marked-root and residual coefficients before taking the coefficient-space closure. | Pass. Each possible source has dimension `n-2` in the `(n-1)`-dimensional space `B_n`. Closure cannot increase image dimension. Root collisions and residual collisions stay in the same closure. |
| Admissible slice meets the generic open | Use the divided tangent–chord incidence and normalize the two contact points to `0,1`. | Pass. The incidence is irreducible because its equation solves uniquely for the quadratic coefficient, and it surjects onto the exact-degree coefficient space. The good inverse image and the endpoint-admissible locus are nonempty opens, hence intersect. |
| Genus used to exclude singularities | Inspect the logical order of the proof. | Pass after reorganization. The three contact strata and the infinity lemma first exclude every non-nodal/non-cuspidal singularity. Squarefreeness then gives exactly `n-2` cusps. The genus formula is invoked only to count the remaining ordinary nodes. |

## Exhaustiveness check

For a finite normalization parameter, every singular branch is controlled by
the contact order of its tangent line:

- order `2`: a smooth immersed branch;
- order `3`: an ordinary cusp when `H''' != 0`;
- order at least `4`: the bad pattern `(4)`.

For a singular image with several normalization preimages:

- a ramified branch plus another branch dominates `(3,2)`;
- three or more preimages dominate `(2,2,2)`;
- exactly two unramified preimages give an ordinary transverse node.

These alternatives also cover cusp–cusp collisions, tritangents, higher
contacts, fourfold fibers, and two node pairs sharing an image. The relative
infinity lemma closes the only projective escape route.

## Scope and residual caveats

- The statement is geometric: the base field is algebraically closed and has
  characteristic zero. No positive-characteristic version is claimed.
- The theorem concerns the reduced discriminant curve. If the resultant is
  retained with scheme multiplicity, that scheme-theoretic formulation needs
  a separate sentence identifying its generic multiplicity; it is irrelevant
  to the singularities of the reduced curve.
- The symbolic checks through degree ten remain regression evidence only. The
  all-degree proof rests on the contact classification and dimension argument,
  not on those computations.
- This is an internal adversarial audit, not independent external peer review.

## Standalone output

The paper is [main.tex](main.tex). Its introduction contains the only
Jacobian-map motivation; all subsequent statements and proofs concern the
pencils `H(W)-sW+t` alone.
