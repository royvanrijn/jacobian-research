# Verification and research checklist

The shortest certificate comes first. Later geometry and consequences must not
be used to support the determinant-and-collision certificate on which they
depend. A checked item has an executable or written certificate in this
repository; it does not by itself claim independent replication.

## 1. Minimal certificate

- [x] Verify `det DF = -2` coefficient by coefficient with SymPy.
- [x] Repeat the determinant computation with an independent standard-library
  sparse-polynomial implementation.
- [x] Substitute the three rational source points exactly and check that they
  are distinct with a common image.
- [x] Check coordinate degrees `(7,6,4)`.
- [x] Audit every stored restatement of the map for signs and coordinate order.
- [ ] Write a short structural determinant proof from the weighted chart.
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
- [ ] Resolve the rational extension between compactifications and inventory
  every boundary divisor and dicritical component.

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
- [x] Prove nonemptiness of the good parameter locus uniformly in every degree
  using contact-incidence dimensions and tangent-chord normalization.
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
- [ ] Reconcile the proof fully with the precise hypotheses of the classical
  projective-duality literature.

## 6. Finite-field statistics

- [x] Derive the fixed-point distribution of `S_n` and its factorial moments.
- [x] State the good-reduction Chebotarev law
  `N_j(q)=p_{n,j}q^3+O_H(q^(5/2))`.
- [x] Deduce limiting image density `1-D_n/n!`.
- [x] Check representative degree-three, four, and five models over a finite
  field against the predicted square-root scaling.
- [ ] Make the finite set of excluded primes effective for an individual
  seed.
- [ ] Compute exact lower-order corrections for selected higher-degree and
  small-characteristic models.
- [ ] Determine how special discriminant singularities alter finite-field
  error terms without changing generic monodromy.

## 7. Derived constructions and implications

- [ ] Reproduce the 95D and 510D generated artifacts from a clean checkout.
- [ ] Audit each stable-equivalence and pairing step against the cited theorem.
- [x] Audit the four recorded downstream conjectural implications, including
  hypotheses, fields, dimensions, and quantifiers.
- [ ] Independently reproduce that implication audit from the primary sources.
- [ ] Determine whether affine, triangular, or tame changes reduce degree,
  support, or coefficient height.

## 8. Reproducibility and provenance

- [ ] Record a pinned commit, platform, dependency versions, commands,
  runtimes, and hashes in an archival verification run.
- [ ] Confirm that generators and verifiers do not share hidden failure modes.
- [ ] Recover or document the absence of original prompts, search code, seeds,
  intermediate candidates, and a primary timestamped announcement.
- [ ] Run the complete certificate independently from a clean checkout and
  archive its logs.

## Recommended next order

1. Make good-prime exclusions effective seed by seed.
2. Continue the exceptional-strata analysis in degrees eight and higher.
3. Perform the second-CAS and clean-checkout archival reproduction.
4. Audit compactifications, stable equivalence, and external implications.
5. Compute compact coefficient ideals for selected lower-dimensional strata.

The routine executable audit is:

```bash
make verify
```

Large generated normal forms and archival checks remain separate from this
routine target.
