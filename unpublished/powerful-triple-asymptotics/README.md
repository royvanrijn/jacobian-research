# Asymptotic counts for powerful triples with arbitrary common divisor

[PDF](paper.pdf) · [LaTeX source](main.tex) · [Metadata](metadata.json) ·
[Verification record](verification.json)

The paper studies ordered positive triples with large comparable fixed
fullness orders. Under an explicit dyadic bound for primitive quotients,
it gives a positive main term at exponent \(1/M\), where \(M\) is the
largest prescribed order. It also gives strictly alternating lower-order
terms and quantitative convergence of the primitive quotient law. The
portion with deficient-prime product \(D>w^{20/M}\), where \(w\) is the
primitive height, has an unconditional expansion in the same powers.

The full expansion uses the small-budget count and an endpoint restriction on
the scaled valuation measures. The paper proves the transfer from those inputs
to the required primitive quotient count. The supporting research contains
written arguments for the [limiting one-third count](../../research/arithmetic/ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)
and [nonatomic endpoint mass](../../research/arithmetic/ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees).
The [common-divisor expansion](../../research/arithmetic/POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#asymptotic-expansion-for-comparable-orders)
and [primitive quotient law](../../research/arithmetic/POWERFUL_TRIPLE_PRIMITIVE_QUOTIENT_LAW.md)
give the corresponding research statements and evidence.

## Build

Run the following command from this directory with TeX Live and `latexmk`.

```sh
latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex
```

The bibliography is included in the LaTeX source.

## Review

The full expansion is stated as a conditional result in this edition. To
prepare an unconditional edition, verify the uniformity in growing valuation
degrees in the signed monomial counts and the limiting passage to the
nonatomic endpoint restriction, and include the complete checked proof or a
precise citable source for those inputs.
