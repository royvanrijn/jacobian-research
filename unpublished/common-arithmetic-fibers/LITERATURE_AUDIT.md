# Literature comparison: finite etale Keller fibers

## Current contribution and the closest theorem

The characteristic-zero existence statement is a consequence of the earlier
coefficient family in the unbylined Ulam technical note, *A counterexample to
the Jacobian conjecture*, Theorem 5.2, pp. 6–7. The primary PDF at
[Ulam](https://www.ulam.ai/research/jacobian.pdf) was inspected on
17 November 2026. Its printed date is 20 July 2026; an archived July copy
of the exact current formulas has not been inspected. This is a comparison
of mathematical statements, not a certification of exact version priority.

For `P(a+s)=sum c_k s^k`, choose `c_1 != 0` and the constant parameters
`eta_1=3(c_3/c_1-1)`, `eta_(k-2)=k c_k/c_1` for `4<=k<=n`, with target
`(p,q,r)=(1,-2c_2/c_1,-2c_0/c_1)`. Equation (10) in that theorem then
becomes `Omega=(2/c_1) P(a+s)`. Separability makes its derivative a unit
in the quotient; `p=1` forces all source points into the reconstruction
chart. The reconstruction identities work over commutative test algebras,
so they identify the whole fiber scheme, rather than just its complex
points. The theorem's degree calculation and monogenicity give the
characteristic-zero full finite-etale realization conclusion. The formulas
are algebraic over the characteristic-zero coefficient field even though
the source states its theorem over the complex numbers.

The revised paper therefore presents a quadratic-gauge proof and Lean
formalization of this consequence, and a separate supplied-presentation
statement in characteristic different from two. No novelty or priority
claim for these refinements is inferred from the comparison. A new order
of quantifiers does not establish a new existence theorem.

The specialization is checked coefficientwise by `verify_comparison.py`.
That exact algebra check supports the comparison; it does not verify either
paper's full reconstruction or function-field proof.

## All-degree constructions

[Alexis Gallagher, *The Jacobian counterexample, explained*, Section 7](https://jacobianfun.org/jacobian-explained)
was inspected on 17 November 2026. The general weighted-lift recipe, rather
than only the finite atlas, gives every generic degree at least three.
The [atlas](https://jacobianfun.org/counterexamples) provides exact split
fibers in degrees 3 through 100. The general construction must be credited
separately from those finite certificates.

[Piotr Migus, arXiv:2607.21572v3](https://arxiv.org/html/2607.21572v3),
Introduction and Proposition 2, expressly use the preceding Ulam and
Gallagher constructions. Version 3 (27 August 2026) identifies the Ulam
note by URL and cites Theorem 5.2 and Corollary 5.3. This supports the
current construction comparison; exact July version priority remains
qualified as above.

## Classical rank-two exclusion

L. A. Campbell, *A condition for a polynomial map to be invertible*,
Math. Ann. 205 (1973), 243–248,
[DOI](https://doi.org/10.1007/BF01349234), unnumbered theorem on p. 244,
is the existing primary interface for the Galois case over the complex
numbers. Razar (1979) and Wright (1981) are later algebraic treatments.
The manuscript descends to a finitely generated characteristic-zero field,
uses finite locally free rank to preserve geometric degree under extension
to the complex numbers, applies the normal quadratic case, and descends
invertibility by faithful flatness. This is a written external interface;
it is not included in the focused formal certificate.

## Other attribution and scope

The cubic seed announced by Alpöge and its Isabelle/HOL certificate are
credited where used. Shaska's grading and quotient geometry explain the
cubic construction but are not substituted for the explicit reconstruction.
Miranda–Neto's ideal-theoretic study of existing Keller fibers and
Lipton–Markakis's Hilbert-irreducibility perspective are adjacent topics;
neither is used to claim an advance based on quantifier order.

Generic monodromy, Hilbertian specialization, exact nonproperness and stable
atomicity are not results of this focused paper. Older audit descriptions
of a generic-monodromy theorem in this manuscript have been removed.

An exhaustive MathSciNet/zbMATH search and external specialist priority
review remain unperformed. The manuscript's contribution is deliberately
stated without a claim of first realization.
