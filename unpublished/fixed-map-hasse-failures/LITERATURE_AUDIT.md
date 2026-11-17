# Literature and priority audit

This note records the external interfaces used by
*Quantitative Hasse-Principle Failures in the Fibers of a Fixed Keller Map*.
It is a source audit, not a claim of peer review or an exhaustive priority
search.

## Load-bearing sources

### Berend–Bilu

Daniel Berend and Yuri Bilu, “Polynomials with roots modulo every integer,”
*Proceedings of the American Mathematical Society* **124** (1996),
1663–1671, DOI
[`10.1090/S0002-9939-96-03210-8`](https://doi.org/10.1090/S0002-9939-96-03210-8).

The paper supplies the fixed-point/Galois criterion for an integer
polynomial to have a root modulo every nonzero integer. The manuscript uses
its low-degree consequence only. Section 5 now spells out the complete
degree-at-most-four factorization analysis: one irreducible factor of degree
`2`, `3`, or `4`, or two irreducible quadratics. This makes the precise
interface visible rather than hiding the lower bound behind a citation.

### de la Bretèche–Tenenbaum

Régis de la Bretèche and Gérald Tenenbaum, “Remarks on the
Selberg–Delange method,” *Acta Arithmetica* **200** (2021), 349–369, DOI
[`10.4064/aa201024-26-5`](https://doi.org/10.4064/aa201024-26-5);
[author preprint](https://arxiv.org/abs/2010.12929).

This source treats multiplicative functions whose Dirichlet series has the
form `ζ(s)^ρ G(s)` and the corresponding asymptotic estimates. The
manuscript now computes the exponent for every character of
`{1,4,7} ⊂ (Z/9Z)×`: `ρ = 1/2` for the trivial character and `ρ = 0` for
the two nontrivial characters. Character orthogonality then yields the
factor `1/3` in the full-family count. The clean family uses the standard
prime-density exponent `ρ = 1/6`.

## Contextual sources

Jack Sonn, “Polynomials with roots in `Q_p` for all `p`,” *Proceedings of
the American Mathematical Society* **136** (2008), 1955–1960, is cited for
the broader intersective-polynomial context. The paper’s explicit local
argument does not depend on Sonn’s theorem.

The prescribed-fiber mechanism is taken from the companion manuscript
`papers/common-arithmetic-fibers`; the present paper rederives the
specialized inverse calculation needed for its fixed map.

## Priority boundary

The checked sources support the intersective-polynomial and analytic inputs.
They do not appear to state the fixed-map conclusion: one explicit
Jacobian-one polynomial map carrying a height-counted family of full
Hasse-failing fibers. That comparison supports the manuscript’s stated
novelty boundary, but it is not a substitute for MathSciNet/zbMATH searches
or specialist review.

Before journal submission, the recommended human checks are:

1. an analytic-number-theory reading of the two character-twisted
   Selberg–Delange paragraphs;
2. an arithmetic-geometry reading of the rank-versus-geometric-degree and
   Berend–Bilu interfaces;
3. a database search for fixed morphisms with quantitatively many
   zero-dimensional Hasse failures.
