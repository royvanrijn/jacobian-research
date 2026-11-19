# Pointwise Binary Correlations in Abelian Splitting Semigroups

The flagship article develops pointwise binary counting for two fixed abelian
subgroup prime prescriptions and unfiltered self-sums for a fixed prime coset.
It opens with the explicit `C_12 + C_12`
asymptotic, proves the general correlation and joint congruence criterion, and
applies them to primitive norm values in every fixed imaginary quadratic order.
The short basis-order corollary gives thirteen for `C_12` with zero adjoined.

The shared arithmetic-model comparison theorem isolates exact invariance under
a common prime class, centered one-variable rough means and positive sparse
sieve bounds. Its assumptions contain no centered binary cancellation or
major-arc energy estimate. The reflected-target and fixed-shift frequency
estimates are proved separately from the same imported analytic inputs.
Subgroups and prime cosets then require explicit arithmetic application checks.
The `3 mod 10` example retains all four quartic-character branches, including
the complex pair. The Connors–Keating supporting note also uses the fixed-shift
comparison; its inert-valuation and dyadic model calculations are separate
arithmetic work and are not included as an application theorem in this article.

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
uses the new general binary theorem. The all-order transfer gives uniform
`o(1)`; the explicit logarithmic saving belongs to the fixed-data semigroup
theorem. The modulus in that theorem is even and fixed, and all primes dividing
it are excluded from the core.

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
positive-restoration route and the conductor/class-group transfer.
The common invariant prime class is essential: arbitrary disjoint prime
prescriptions and growing sublinear shifts are outside the comparison theorem.
Broader historical-priority checking remains
separate from the exact comparisons in the introduction.
