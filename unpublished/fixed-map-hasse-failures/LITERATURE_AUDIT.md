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

The rational count uses **Theorem 1.2**, formulae (1.10)–(1.12), and
the proof in Section 4 (author PDF pages 3 and 8–12). The inspected
[author PDF](https://tenenb.perso.math.cnrs.fr/PPP/On-SD.pdf) notes minor
corrections to the published version. Its SHA-256 on 17 November 2026 was
`cc20403e28ec1513664bc10c528473f8b20744135e3582580a769a588fee688f`.

For the paper's indicator and its cubic-character twists, the theorem's
class condition holds with `r=1/2, σ=3/4`; the prime means are `ρ=1/2,0,0`.
Fixed-modulus prime distribution supplies any fixed logarithmic-power
error, and prime-power values are bounded. Formula (1.10) gives zero
coefficients for `ρ=0`. Section 5 proves the coprimality, region and
valuation-stratum passages separately. The sufficient integer-family
application has the same hypotheses; the clean prime-support family has
prime mean `1/6`.

## Contextual sources

Jack Sonn, “Polynomials with roots in `Q_p` for all `p`,” *Proceedings of
the American Mathematical Society* **136** (2008), 1955–1960, is cited for
the broader intersective-polynomial context. The paper’s explicit local
argument does not depend on Sonn’s theorem.

The prescribed-fiber mechanism is taken from the companion manuscript
[Common Arithmetic Fibers](../common-arithmetic-fibers/); the present paper
rederives the specialized inverse calculation needed for its fixed map.

## Priority boundary

The checked sources support the intersective-polynomial and analytic inputs.
They do not appear to state the fixed-map conclusion: one explicit
Jacobian-one polynomial map carrying a height-counted family of full
Hasse-failing fibers. That comparison supports the manuscript’s stated
novelty boundary, but it is not a substitute for MathSciNet/zbMATH searches
or specialist review.

Before journal submission, the recommended human checks are:

1. an analytic-number-theory reading of the complete rational count,
   including the coprimality convolution and uniform valuation tails;
2. an arithmetic-geometry reading of the rank-versus-geometric-degree and
   Berend–Bilu interfaces;
3. a database search for fixed morphisms with quantitatively many
   zero-dimensional Hasse failures.
