# Gao F6: complete fibres and the nonproper locus

14 November 2026. This follows the F4 classification and addresses another
full-fibre problem left open in [Gao v2](https://arxiv.org/html/2608.00222v2).

- **GAOF6FIB1:** the possible complex fibre sizes are exactly 0,1,2,3,4,6.
  Five is impossible. On the target boundary `(0,A,B,D,Z)`, there are
  three points if `A^2+4B != 0`, and one otherwise. All eleven sextic
  multiplicity partitions are classified; seven occur.
- The exact two-double-root coefficient stratum is A1 minus thirteen
  points; the triple-root stratum with three simple roots is A1 minus
  nine points. There are six triple-double-simple coefficient points,
  three quadruple-simple-simple points, and one three-double point.
  Each nonzero-chart target stratum is its product with C* times A1.
- **GAOF6NP1:** the full nonproper-value locus has exactly two irreducible
  hypersurface components. One is the entire target boundary; the other
  is the saturated sextic discriminant. The discriminant in coefficient
  space has an explicit A2 normalization.

[Proof and prior-credit boundary](PROOF.md) · [checker](verify.py) ·
[receipt](receipt.json) · [source pins](sources/manifest.json).

The simple-root correspondence, omitted locus, boundary isomorphism and
identity factor are Gao's prior results. Earlier selected rational fibres
are also credited. This written proof and exact arithmetic support the
scoped claims; independent review and novelty confirmation remain open.
F5, F7, JC2 and unrestricted HC4 remain unresolved.

```sh
timeout 60s .venv/bin/python research/jacobian-followup/gao-f6-fibres/verify.py
```
