# Verification and research checklist

The shortest certificate comes first. Later geometry and consequences should
never be used to support the determinant-and-collision argument on which they
already depend.

## 1. Minimal certificate

- [x] Verify `det DF = -2` coefficient by coefficient with SymPy.
- [x] Repeat the determinant computation with an independent standard-library
  sparse-polynomial implementation.
- [x] Substitute the three rational source points exactly and check that they
  are distinct with a common image.
- [x] Check coordinate degrees `(7,6,4)`.
- [ ] Write a short structural determinant proof from the weighted chart.
- [ ] Reproduce the certificate in an unrelated computer algebra system.

## 2. Inverse and boundary geometry

- [x] Verify the primitive cubic and rational reconstruction identities.
- [x] Verify the discriminant formula and representative fiber cardinalities.
- [x] Check exceptional inverse-cubic cases and denominator-safe image and
  nonproperness inclusions.
- [ ] Independently certify irreducibility of the discriminant hypersurface,
  its normalization, and the `S_3` monodromy argument.
- [ ] Resolve the rational extension between compactifications and list all
  boundary divisors and dicritical components.

## 3. Derived constructions

- [ ] Reproduce the 95D and 510D generated artifacts from a clean checkout.
- [ ] Audit each stable-equivalence and pairing step against the cited theorem.
- [ ] Audit every downstream conjectural implication, including hypotheses,
  fields, dimensions, and quantifiers.
- [ ] Determine whether affine, triangular, or tame changes reduce degree,
  support, or coefficient height.

## 4. Weighted family

- [x] Verify the general seed construction and its constant Jacobian formula.
- [x] Complete the exact quartic-sheet image, fiber, singularity, and
  nonproperness checks.
- [x] Record the bounded seed scan separately as exploratory evidence.
- [ ] Formulate and prove the general theorem with explicit assumptions on
  primitive zeros, critical behavior, monodromy, and boundary branches.

## 5. Reproducibility and provenance

- [ ] Record a pinned commit, platform, dependency versions, commands,
  runtimes, and hashes in an archival verification run.
- [ ] Confirm that generators and verifiers do not share hidden failure modes.
- [ ] Recover or document the absence of original prompts, search code, seeds,
  intermediate candidates, and a primary timestamped announcement.
