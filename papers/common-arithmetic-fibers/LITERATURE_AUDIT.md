# Literature audit: finite étale Keller fibers

**Search date:** 24 July 2026

## Claim audited

The manuscript proves that every finite étale algebra over a
characteristic-zero field, except one of rank two, occurs as a full fiber of a
polynomial Keller map. In rank at least three the realization is explicit in
`A^3` and has determinant `-2`.

The priority language in the manuscript is deliberately qualified:

> To our knowledge, no earlier source states this finite-étale realization
> theorem or the resulting rank classification.

This file records the search supporting that sentence. It is not a substitute
for MathSciNet/zbMATH review or direct expert consultation.

## Queries

The following exact and near-exact queries were run against general web,
arXiv, MathOverflow, and bibliographic search results, with spelling variants:

```text
"finite étale algebra" "Keller map"
"finite etale algebra" "Keller map"
"Keller fiber" polynomial map Jacobian
"full fiber" "Keller map" finite etale
"every finite étale" Keller
"prescribed finite étale" polynomial map fiber constant Jacobian
finite étale scheme complete fiber polynomial étale endomorphism affine space
arbitrary number field fiber Keller map
```

Repository and coefficient-string searches were also used for the phrases
`finite étale Keller fiber`, `prescribed fiber`, and `full fiber`.

## Degree-two source audit

The audit corrected one initially over-specific citation. Chapter I,
Theorem 2.1 of Bass--Connell--Wright concerns the implication from bijectivity
to a polynomial inverse; it is not the primary source for the Galois case used
here.

The relevant historical sources are:

- L. A. Campbell,
  *A condition for a polynomial map to be invertible*, Math. Ann. 205
  (1973), 243--248,
  [DOI](https://doi.org/10.1007/BF01349234), for the complex Galois case;
- M. Razar,
  *Polynomial maps with constant Jacobian*, Israel J. Math. 32 (1979),
  97--106,
  [DOI](https://doi.org/10.1007/BF02764906);
- D. Wright,
  *On the Jacobian conjecture*, Illinois J. Math. 25 (1981), 423--440,
  [DOI](https://doi.org/10.1215/ijm/1256047158).

The manuscript does not rely on an unverified field-general formulation of
one of these sources. Instead it derives the arbitrary characteristic-zero
field statement from Campbell's complex theorem: descend the coefficients to
a finitely generated subfield, embed it into `C`, use stability of generic
degree under scalar extension, and then descend the unique formal/polynomial
inverse. This supplies exactly the ground-field generality required for the
rank-two exclusion.

## Closest located sources to the realization theorem

### Gallagher: split complete fibers

A. Gallagher,
[Counterexample atlas: generic fiber degrees 3 through 100](https://jacobianfun.org/counterexamples),
2026.

The atlas gives explicit Keller maps in every generic degree and complete
split rational fibers. It realizes the split algebra `K^N`, not arbitrary
finite étale `K`-algebras, and does not state the rank classification.

### Miranda--Neto: ideals of existing fibers

C. B. Miranda--Neto,
*An ideal-theoretic approach to Keller maps*, Proc. Edinburgh Math. Soc. 62
(2019), 1033--1044,
[DOI](https://doi.org/10.1017/S0013091519000099).

This work studies radicality and maximal-ideal decomposition for fibers of a
Keller map already given. It does not construct a Keller map from a prescribed
finite étale algebra.

### Lipton--Markakis: rational images and Hilbert irreducibility

R. J. Lipton and E. Markakis,
[Some remarks on the Jacobian conjecture and connections with Hilbert's
irreducibility theorem](https://arxiv.org/abs/math/0507525), 2005.

This studies rational images and Hilbert subsets under Keller maps. Its
quantifier order starts with a map and studies its values; the realization
theorem starts with a finite étale algebra and constructs the map.

### Shaska: graded arithmetic structure

T. Shaska,
[Graded Keller maps and the Jacobian conjecture](https://arxiv.org/abs/2607.20210),
2026.

This studies graded Keller maps, quotient geometry, and arithmetic thinness.
No prescribed finite-étale fiber theorem was located there.

### General specialization literature

General work on specialization of finite covers, Hilbert irreducibility, and
prescribed local behavior can realize finite étale algebras inside suitable
covers. Those results do not impose that the cover is a polynomial
constant-Jacobian self-map of affine space. They are therefore adjacent
arithmetic inputs, not prior versions of the Keller-fiber theorem.

## Negative-search conclusion

No located source states any of the following:

1. every squarefree polynomial of degree at least three is the exact full
   inverse polynomial of an explicit Keller fiber after translation;
2. every finite étale algebra of rank at least three occurs as a full fiber of
   a polynomial Keller map of `A^3`;
3. the possible ranks of nonzero Keller fibers are exactly
   `1,3,4,5,...`.

The manuscript should retain the phrase **“to our knowledge”** and should not
claim that the search proves absolute priority. Before journal submission,
the statement should also be sent to specialists in Keller maps and finite
cover arithmetic, and checked in MathSciNet and zbMATH by a reader with full
access.
