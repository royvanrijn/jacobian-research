# Pointwise representation counts for Eisenstein restricted sums

**Published preprint · [Hexagon v1](https://hexagonmath.org/2611.00593v1) · 19 November 2026 · [Deposited files](published/)**

For every sufficiently large admissible \(n\), the number of distinct ordered
Eisenstein restricted-prime decompositions in a prescribed summand interval
has order \(HJ(n)/\log n\) under the stated window condition. The bounds cover
ordinary sums and both orientations of the mixed sum with one ramified factor
of \(3\). They apply to central widths above the power \(2/3\) and prescribed
small summands at scales above \(1/2\).

[PDF](paper.pdf) · [LaTeX](main.tex) · [Source archive](source.zip) ·
[Metadata](metadata.json) · [Verification](verification.json)

Build here or in the extracted archive:

~~~sh
SOURCE_DATE_EPOCH=1794960000 FORCE_SOURCE_DATE=1 latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex
~~~

The proof combines [moving-window estimates](../../research/arithmetic/EISENSTEIN_MOVING_WINDOWS.md),
local density recurrences and a cofactor sieve. The
[exact companion](../eisenstein-restricted-prime-sums/) gives complete primitive
norm sumsets at four discriminants. The [pointwise companion](../pointwise-binary-splitting-semigroups/)
develops broader correlation asymptotics and common-field genus comparison.
Both companions are published preprints on Hexagon; their [pointwise](../pointwise-binary-splitting-semigroups/published/) and [exact-sumset](../eisenstein-restricted-prime-sums/published/) deposited editions preserve the public files.
