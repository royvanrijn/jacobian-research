# Gao F7: all boundary fibres and one remaining cardinality question

14 November 2026. Scoped additions to [Gao v2](https://arxiv.org/html/2608.00222v2):

- **GAOF7BND1:** a boundary target `(0,A,B,D,Z)` has two source points
  when `A!=0`, and one when `A=0`. Both source pieces are exhausted.
- **GAOF7CARD1:** precisely ten cardinalities are certified to occur:
  **1,2,4,5,6,7,8,9,10,12**. Eleven is impossible; the
  [separate surjectivity proof](../gao-f7-surjectivity/README.md) excludes zero.
  **Only three remains undecided.** New algebraic witnesses for four, five and six points use
  nonsingular finite-field factorization certificates and Hensel lifting.
  The eliminant is never a square, excluding all-even multiplicity patterns.
- **GAOF7NP1:** the complete nonproper-value locus has two irreducible
  hypersurface components. The nonzero-chart discriminant has normalization
  A3. Its saturated pullback meets `c=0` with equation `-182528*A^5`.

This is not a complete F7 multiplicity-stratum classification. Gao's
simple-root criterion, twelve-point and seven-point witnesses are prior
results. F5, the remaining F7 gates, JC2 and unrestricted HC4 stay open.

[Written proof](PROOF.md) · [checker](verify.py) · [receipt](receipt.json) ·
[fixed lifting inputs](lifting-points.json) · [source pins](sources/manifest.json).

```sh
timeout 60s .venv/bin/python research/jacobian-followup/gao-f7-fibres/verify.py
```

The checker repeats no finite-field search or Groebner calculation. The
incomplete sixfold-root attempt is preserved and proves no exclusion.
Independent review and novelty confirmation remain outstanding.

An empty fibre would require `deg gcd(R,R')>=7`; GAOF7SURJ1 excludes it
by power-factorization identities. Three points require gcd degree at
least five. A uniform bound at most four would exclude that last size,
but remains unproved. The full multiplicity stratification is still open.
