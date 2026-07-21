# Verification and research checklist

The shortest certificate comes first. Later geometry and consequences must not
be used to support the determinant-and-collision certificate on which they
depend. A checked item has an executable or written certificate in this
repository; it does not by itself claim independent replication.  For a
quantified all-degree item, a check normally records completion of a uniform
written proof.  Finite-degree scripts are regression tests for that proof, not
formal derivations of its universal quantifier.

## 1. Minimal certificate

- [x] Verify `det DF = -2` coefficient by coefficient with SymPy.
- [x] Repeat the determinant computation with an independent standard-library
  sparse-polynomial implementation.
- [x] Substitute the three rational source points exactly and check that they
  are distinct with a common image.
- [x] Check coordinate degrees `(7,6,4)`.
- [x] Audit every stored restatement of the map for signs and coordinate order.
- [x] Give the hand-checkable structural determinant proof through the
  `(t,r,c)` inverse chart.
- [x] Locate, scope, credit, pin, and reproduce Dean Cureton's external Lean 4
  formalization of the determinant and collision (`make verify-lean-c01`).
- [ ] Reproduce the certificate in an unrelated computer algebra system.

## 2. Cubic inverse and boundary geometry

- [x] Verify the primitive cubic and rational reconstruction identities.
- [x] Verify the discriminant formula and representative fiber cardinalities.
- [x] Classify every exceptional inverse-cubic fiber.
- [x] Prove both image and nonproperness inclusions with explicit charts and
  escaping paths.
- [x] Check explicitly that denominator clearing loses no `C=0` or
  discriminant-boundary stratum.
- [x] Certify the discriminant normalization, singular locus, and local
  transposition monodromy algebraically.
- [ ] Recheck irreducibility, normalization, and `S_3` monodromy in a second
  proof environment.
- [x] Construct the normalized inverse-graph compactification and classify
  every dicritical prime responsible for nonproperness and omission.
- [ ] Compute an explicit minimal blow-up sequence and discrepancies for that
  normalized graph model.

## 3. Weighted-seed theorem

- [x] Implement the generalized weighted-seed builder and constant-Jacobian
  identity.
- [x] Derive the inverse pencil `E(W)=H(W)-sW+t` and repeated-root
  normalization for every admissible seed.
- [x] Prove irreducibility of the generic inverse pencil and geometric and
  arithmetic monodromy `S_n` in characteristic zero.
- [x] Prove the canonical-family image theorem, including direct `C=0` fibers,
  boundary saturation, nonproperness, and surjectivity from inverse degree
  five onward.
- [x] Prove the boundary-clean deformed-seed theorem and its quartic omitted
  double-double curve.
- [x] Implement the complete omitted-value classifier by root-multiplicity
  partitions.
- [x] Prove the repeated-root boundary theorem, including exact saturation
  exponents and nonsplit factors.
- [x] Interpret omitted partitions as singularities of the dual discriminant
  curve.

## 4. Quartic geometry

- [x] Verify the exact quartic map, determinant, inverse, and reconstruction.
- [x] Classify the quartic discriminant as two ordinary cusps and one ordinary
  node.
- [x] Compute the direct `C=0` fibers and every escaping boundary sequence.
- [x] Prove both inclusions for the nonproperness set and the converse
  properness statement.
- [x] Certify the radical decomposition of the singular locus.
- [x] Prove the exact quartic image and all fiber sizes.
- [x] Prove geometric and arithmetic monodromy `S_4`.

## 5. Generic discriminant geometry

- [x] Implement the tangent-line normalization
  `nu_H(r)=(H'(r),rH'(r)-H(r))`.
- [x] Certify that simple roots of `H''` give ordinary cusps.
- [x] Implement diagonal-divided bitangent equations and account explicitly
  for saturation by `(r-u)H''(r)H''(u)`.
- [x] Check the unique point at infinity and its smooth local parameter.
- [x] Verify exact rational admissible examples in degrees `3,...,10`, with
  `n-2` cusps and `(n-2)(n-3)/2` nodes.
- [x] Use the global singular-scheme length to exclude hidden tritangents,
  cusp-branch collisions, coincident nodes, and higher singularities in those
  examples.
- [x] Define the compactified `(4)`, `(3,2)`, and `(2,2,2)` incidence
  varieties as saturated projective graph closures with explicit homogeneous
  coefficient equations.
- [x] Prove from normalization branches and common-tangent contact orders that
  those three patterns exhaust every bad affine singularity.
- [x] Classify all finite-collision, residual, marked-root-at-infinity,
  residual-at-infinity, coefficient-at-infinity, and degree-drop boundaries.
- [x] Prove that each compactified graph remains irreducible of dimension
  `n-2` and that proper projection keeps its closed image at dimension at
  most `n-2`.
- [x] Prove the divided tangent-chord incidence is an irreducible affine bundle
  dominating coefficient space, and use its explicit admissible witness to
  show every normalized weighted slice meets the good open.
- [x] Deduce generic surjectivity for every inverse degree `n>=5`.
- [x] Make the closed bad loci and their saturations into reusable symbolic
  certificates rather than relying on the global Tjurina check for examples.
- [x] Generalize the incidence API to arbitrary contact partitions, including
  residual factors and equal-part permutation quotients.
- [x] Recover the degree-five exceptional polynomial from the uniform `(3,2)`
  root curve.
- [x] Classify the degree-six full-contact strata and their closure relations.
- [x] Prove that the degree-seven nonsurjective locus has codimension two.
- [x] Prove the uniform exceptional-locus union and the formula
  `dim E_lambda=ell(lambda)-1` using weighted Newton sums and a weighted
  Vandermonde determinant.
- [x] Implement simultaneous multiple-omission incidences with separate
  common-value and distinct-value ideals.
- [x] Settle the degree-six `(2,2,2)`/`(3,3)` intersection: the exact strata
  are disjoint and their closures meet on `(6)`.
- [x] Use degree eight as a theorem test: verify the `(2,2,2,2)` and `(3,3,2)`
  dimensions, their shared `(6,2)` and `(8)` boundaries, and the absence of
  exact off-diagonal intersections.
- [x] Formalize the collision partial order `lambda<=mu` by merging parts.
- [x] Construct a uniform tangent-chord root-splitting deformation proving
  `E_mu` is contained in the closure of `E_lambda` whenever `lambda<=mu`.
- [x] Build the uniform affine-difference incidence
  `a M_lambda-b M_mu=alpha W+beta`.
- [x] Use Mason--Stothers to exclude every off-collision two-omission solution
  for distinct maximal partitions with parts only two and three.
- [x] Prove uniform irreducibility of every maximal
  `Phi_(2^a 3^b)` hypersurface.
- [x] Deduce that maximal 2/3 partitions index the irreducible components of
  the exceptional-locus closure.
- [x] Derive the component dimensions, total exceptional codimension
  `ceil(n/2)-2`, component-count generating function, and exact
  common-coarsening intersection formula.
- [x] Prove the contact-atom principle: `{2,3}` are the indecomposable allowed
  contact orders, and every higher multiplicity is their collision boundary.
- [x] Strengthen Mason separation to every pair of distinct full-contact
  partitions using the excess identity.
- [x] Close the even-degree all-double equality case and prove that every seed
  has at most one omitted inverse-pencil value.
- [x] Upgrade the exceptional-locus union to the disjoint exact-partition
  stratification `N_n = disjoint_union E_lambda`.
- [x] Prove that every admissible maximal quotient hypersurface is smooth,
  including all collision diagonals.
- [x] Identify its finite generically degree-one seed map as the normalization
  of the corresponding irreducible component.
- [x] Count geometric branches over the generic point of each collision
  stratum by the local allocations `2i+3j=m`.
- [x] Verify that degree twelve, `C_(3,2)` over `E_(6,6)`, is the first
  geometric two-branch collision, and distinguish it from the first abstract
  equal-multiplicity orbit-type ambiguity in degree fourteen.
- [x] Compute the generic completed-local fiber product, conductor, and
  transverse multiplicity four for `C_(3,2)` along `E_(6,6)`.
- [x] Classify the affine-difference two-transfer block `Z_2`, including its
  rank-four flat structure and non-Gorenstein coincident-root fiber.
- [x] Prove all global equalizer theorems for compensating transfers by the
  universal Wronskian, killing the two shared affine coefficients at once.
- [x] Audit the first `Z_3/Z_4` global configurations, including two `Z_3`
  blocks, elementary compensation of `Z_3` and `Z_4`, mixed
  `Z_3 completed-tensor Z_2`, and a three-block transfer.
- [x] Prove the strong Hensel-product decomposition and reduce affine rigidity
  to a transverse length upper bound.
- [x] Prove the quadratic-cubic colength-sixteen upper bound in the `(2,-2)`
  case and identify its affine equalizer with `Z_2 completed-tensor Z_2`.
- [x] Prove the mixed `(2,-1,-1)` affine equalizer theorem with transverse
  length sixteen.
- [x] Classify `Z_3` as a rank-eight finite-flat affine-difference block and
  compute its non-Gorenstein coincident-root fiber.
- [x] Classify `Z_4` as a rank-sixteen finite-flat affine-difference block and
  compute its Hilbert function and four-dimensional socle.
- [x] Prove the all-`k` Boolean-quotient model, finite flatness of rank `2^k`,
  affine/strong equality, and collided Hilbert series `(1+t)^k`; retain exact
  Groebner regressions through `k=6`.
- [x] Independently reconstruct every `Z_k` as the divided-power symmetric
  product of the cusp conductor ribbon, using confluent divided differences
  and a split-surjection/Nakayama comparison instead of invariant theory.
- [ ] Determine ramification divisors, conductors, and scheme-theoretic branch
  multiplicities of the remaining component normalizations.
- [x] Verify the threshold-`r` atom set `{r,...,2r-1}` and automatic abc bound
  for `r>=3`.
- [ ] Determine the scheme-theoretic component intersections, embedded
  components, and collision multiplicities.
- [ ] Reconcile the proof fully with the precise hypotheses of the classical
  projective-duality literature.

## 6. Finite-field statistics

- [x] Derive the fixed-point distribution of `S_n` and its factorial moments.
- [x] State the good-reduction Chebotarev law
  `N_j(q)=p_{n,j}q^3+O_H(q^(5/2))`.
- [x] Deduce limiting image density `1-D_n/n!`.
- [x] Check representative degree-three, four, and five models over a finite
  field against the predicted square-root scaling.
- [x] Make the excluded primes effective from denominators, Jacobian units,
  `n!`, and the squarefree-factor discriminant/resultant certificate; audit
  monodromy preservation and every Chebotarev twist hypothesis.
- [x] Separate the exact discriminant and `C=0` contributions and prove the
  universal boundary total and first-moment identities.
- [ ] Compute exact lower-order corrections for selected higher-degree and
  small-characteristic models.
- [ ] Determine how special discriminant singularities alter finite-field
  error terms without changing generic monodromy.

## 7. Master cancellation construction

- [x] Prove the localized determinant `det J(P,Q,R)=-C` for arbitrary `m,r`.
- [x] Derive the finite operator `L_(m,r)` and prove that its vanishing is
  equivalent to polynomiality.
- [x] Replace finite parameter/coefficient tables by the truncated-binomial
  polynomial `M_(m,r)` and the exact Hensel recurrence for `h(A)`.
- [x] Prove generic degree `r(m+1)+1`, irreducibility, separability, exact
  reconstruction, and a full symbolic collision for every `m,r`.
- [x] Identify the `(1,1)` member with the original cubic map and separate
  higher ramification from generic weighted seeds.
- [ ] Classify all polynomial cancellation branches up to field descent,
  Galois conjugacy, tail deformation, and polynomial coordinate equivalence.

## 8. Derived constructions and implications

- [ ] Reproduce the 95D and 451D generated artifacts from a clean checkout.
- [ ] Audit each stable-equivalence and pairing step against the cited theorem.
- [x] Audit the four recorded downstream conjectural implications, including
  hypotheses, fields, dimensions, and quantifiers.
- [ ] Independently reproduce that implication audit from the primary sources.
- [ ] Determine whether affine, triangular, or tame changes reduce degree,
  support, or coefficient height.

## 9. Reproducibility and provenance

- [ ] Record a pinned commit, platform, dependency versions, commands,
  runtimes, and hashes in an archival verification run.
- [ ] Confirm that generators and verifiers do not share hidden failure modes.
- [ ] Recover or document the absence of original prompts, search code, seeds,
  intermediate candidates, and a primary timestamped announcement.
- [ ] Run the complete executable identity and regression suite independently
  from a clean checkout and archive its logs.

## Recommended next order

1. Classify the polynomial cancellation branches of `L_(m,r)`.
2. Make good-prime exclusions effective seed by seed.
3. Continue the exceptional-strata analysis in degrees eight and higher.
4. Perform a clean-checkout archival reproduction; the independent exact
   implementation is now part of `make verify-derived`.
5. Audit compactifications, stable equivalence, and external implications.
6. Compute compact coefficient ideals for selected lower-dimensional strata.

The routine executable identity and regression audit is:

```bash
make verify
```

Large generated normal forms and archival checks remain separate from this
routine target.  Passing it does not replace the written proofs of the
all-degree theorems.
