# Paper audit

This audit covers only the C01--C04 paper. It does not treat later repository
claims.

## Remaining proof obligations

1. **Root-one algebraization at `C=0`.** The weighted notes use the completed
   expansion at `W=1` to show formally that `z` has no pole. The paper marks
   the missing citation or lemma identifying this formal regularity with
   regularity in the relevant algebraic local ring.

2. **Zero-cluster special points.** The chart `W=CR` proves the stated generic
   special-fiber equations. It does not classify repeated points of that chart,
   intersections with the discriminant, or every boundary divisor over
   `C=0`. The paper makes no such claim.

3. **Uniform boundary enumeration.** The global theorem identifies the affine
   source with the simultaneous regularity open and describes its complement
   as the support of the negative valuations of `x,y,z`. C01--C04 do not
   enumerate those valuations for every admissible root profile of `H`.

4. **Independent audit of the normalized-incidence argument.** The use of the
   universal property of normalization and Zariski's Main Theorem is written
   out, but it has not been independently reviewed in an algebraic-geometry
   proof environment.

5. **External monodromy inputs.** The proof depends on Zariski--Nagata purity
   and the triviality of finite etale covers of affine space over an
   algebraically closed characteristic-zero field. The paper states the exact
   reduction and cites these theorems; it does not independently reprove them.

6. **Provenance record.** The public formula and contemporaneous attribution
   are corroborated, but the original model conversation, discovery search,
   and an unambiguous archived announcement timestamp have not been located.
   This affects historical attribution, not the algebraic proofs.

## Claims deliberately removed or weakened from the previous draft

- The weighted affine source is no longer called the global simple-root locus.
- Normality is attributed to taking the normalization; the raw incidence is
  only proved integral and finite flat.
- The formal root-one expansion is not presented as a complete algebraic chart
  without a check marker.
- The zero-cluster equations are stated generically, not as a complete `C=0`
  fiber classification.
- Exact images, omitted loci, and boundary decompositions for weighted seeds
  are not asserted.
- Bounded symbolic verification is described as regression evidence, never as
  proof of an all-degree statement.
