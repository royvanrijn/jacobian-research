# Pointwise Binary Correlations in Abelian Splitting Semigroups

The flagship opens with concrete arithmetic consequences: `C_12 + C_12`,
Smilansky's five quadratic-form pair predictions, consecutive sums of two
squares with coefficient `1/2`, and prime-`3 mod 4` self-sums on every
sufficiently large even target. It then develops the shared comparison and the
general pointwise theorem for fixed subgroup pairs and prime-coset self-sums.
The arithmetic application now treats primitive values of arbitrary fixed
positive definite forms, mixed pairs and genus invariance. Principal norms
remain the prototype; the short basis-order-thirteen corollary remains in
the main text.

The shared arithmetic-model comparison theorem isolates exact invariance under
a common prime class, centered one-variable rough means and positive sparse
sieve bounds. Its assumptions contain no centered binary cancellation or
major-arc energy estimate. The reflected-target and fixed-shift frequency
estimates are proved separately from the same imported analytic inputs.
Subgroups and prime cosets then require explicit arithmetic application checks.
The `3 mod 10` example retains all four quartic-character branches, including
the complex pair. The ordinary quadratic-indicator application uses the
fixed-shift comparison after checking the inert valuations and exceptional-prime
conditions. It includes Connors–Keating and the five forms treated by Smilansky.

The editorial decision is to include Connors–Keating in this article.
Its five-page quadratic-character application uses the shared fixed-shift
theorem and introduces no new graph, energy or prime-distribution input.
It proves Smilansky's five ordinary value-pair predictions for `x^2+d*y^2`,
`d=1,2,3,4,7`, after checking fixed valuation masks and the coupled
exceptional-prime factor. The consecutive-value coefficients for `d=1,2`
are `1/2` and `1/sqrt(2)`. This application claims the stated little-o
remainder, with no growing-shift uniformity or arbitrary principal-form
identification.

The coprime `C_12` corollary follows by Möbius inversion from the uniform
count. It multiplies the leading term by
`prod_{p|N, p=1 mod 12}(1-1/p)` and retains relative error
`O((log N)^(-1/65536))`. The proof controls both the large-divisor tail
and the loss from dividing by this finite product. It directly recovers
the coprime scope in the 1993 comparison, using the article's existing
analytic dependencies.

The arithmetic application retains two proof routes. Its order-of-magnitude
theorem uses older weighted-correlation and positive-restoration inputs
independently of the recent zero-free and graph theorems. The leading asymptotic
uses the new general binary theorem. The all-order leading-term transfer gives
uniform `o(1)`; an explicit logarithmic remainder for the full leading term is
asserted for the fixed-data semigroup theorem. Its modulus is even and fixed, and all primes dividing
it are excluded from the core.
The independent order-bound theorem and dependency distinction stay prominent
in the main text; its full classical core proof is in Appendix A. Numerical
saving exponents and localization budgets are collected in the technical
reconstruction section, after the early statements with positive savings.

The class-group error now has a saving independent of the class number.
For `G = Pic(O)` and `B = 2G != 0`, every exponent
`eta < min_{ell prime, ell | |B|} (ell-1)/(2ell)` is allowed; this
threshold is at least `1/4`. In particular every fixed order has
`R_Lambda(N)-R_P(N) <<_O N J_O(N)/(log N)^(6/5)` for all `N >= 3`.
When `2G = 0`, the value sets are exactly equal. The constants depend on
the fixed order. The proof charges failure to an entire deficient collection
outside a maximal subgroup and uses the existing marking sieve. It adds no
analytic input. This estimate does not supply a quantitative remainder for
the full conductor-prefix leading-term transfer, which still asserts `o(1)`.

For fixed forms `Q_i,Q'_i` of one discriminant, replacing each form within
its own genus gives `R_{Q1,Q2}/R_{Q'1,Q'2} = 1+O((log N)^(-eta))` on
sufficiently large common locally admissible targets. The two slots may
belong to different genera. This quantitative comparison and the mixed
order bounds use the classical route, independently of the recent inputs.
The recent binary theorem supplies the explicit leading coefficient, using
the mixed local **value** densities. No `1/h` factor is inserted.

The form/ideal dictionary targets `-[A]`, and its local shadow means that
some primitive invertible ideal of norm `n` lies in `-[A]+2Pic(O)`.
This remains existential at conductor primes: a number can occur in more
than one genus. Separate local genus biases are retained in the two prefix
filters. A common-field corollary then allows two different fixed conductors,
marking one order's deficient prime families at a time. It needs no joint
ring-class prime-distribution input. These statements concern numerical
values at coprime coordinates, not ordinary unrestricted coordinate values
or forms over different quadratic fields. The exact companion remains frozen.

[PDF](paper.pdf) · [Source](main.tex) · [Metadata](metadata.json) ·
[Verification record](verification.json) · [Standalone source archive](source.zip)

The [exact companion](../eisenstein-restricted-prime-sums/) evolves the existing
Eisenstein manuscript and retains its proofs and certificates while integrating
discriminants -7, -8 and -11. The separately
[deposited window paper](../../results/eisenstein-pointwise-counts/) addresses
where the summands can be placed. These are different statements and proof routes.

## Build

From this directory, with TeX Live and `latexmk`, run:

```sh
SOURCE_DATE_EPOCH=1795046400 FORCE_SOURCE_DATE=1 latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex
```

All essential application arguments are included in the article. No research
directory is needed to compile the source archive. The imported primary
literature is cited with exact versions and statements. In particular the
general correlation depends on the September 2026 OpenAI zero-free and
weighted-dilation preprints, and adapts Walsh's November 2026 framework.
The audit of this article does not independently revalidate those entire
preprints. Their source editions and inspected application scopes are recorded
in the verification record.

## Review scope

The verification record separates mathematical application review, source
inspection, isolated PDF reading, rendered-page inspection and source-bundle
builds. None implies external human review, independent external replay,
formal verification or established historical priority. There is no numerical
onset or complete small-exception assertion in this paper.

Before submission, obtain specialist review of the shared comparison hypotheses,
both frequency regimes, the complex-character coset model, the older
positive-restoration route, the form-to-ideal dictionary at conductor primes,
the mixed genus filters and their conductor/class-group transfer.
The common invariant prime class is essential: arbitrary disjoint prime
prescriptions and growing sublinear shifts are outside the comparison theorem.
Broader historical-priority checking remains
separate from the exact comparisons in the introduction.
