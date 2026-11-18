# Pointwise representation counts for Eisenstein restricted sums

For every sufficiently large admissible target, this paper gives matching
bounds of order **H J(n) / log n** for the number of distinct ordered
Eisenstein restricted-prime decompositions in a prescribed summand interval.
The theorem covers ordinary sums and each orientation of the mixed sums
with one ramified factor of 3. Its window condition gives central widths
above the power \(2/3\) and prescribed small summands at scales above \(1/2\).

[PDF](paper.pdf) · [LaTeX](main.tex) · [Standalone source archive](source.zip) ·
[Verification record](verification.json)

The [moving-window result](../../research/arithmetic/EISENSTEIN_MOVING_WINDOWS.md)
is the mathematical basis of the paper. The manuscript includes the angular
Hecke estimate, the singular local recurrences, the cofactor sieve and the
conversion to distinct value pairs. The
[exact sumset paper](../eisenstein-restricted-prime-sums/) treats existence
at every admissible target and the primitive-norm consequences.

## Build

From this directory or the extracted source archive, run:

~~~sh
SOURCE_DATE_EPOCH=1794960000 FORCE_SOURCE_DATE=1 latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex
~~~

The argument is analytic. The retained exact local computation checks
selected split-prime recurrences and the prime-square density underlying
the cofactor sieve. The proof of the pointwise theorem uses the written
estimates. The comparison constants and onset are existential.

A fresh internal mathematical audit reconstructed the Hecke operator,
its stabilizers and measure, and the scalar exceptional spectrum. The
paper includes the prime-power normalization and a worked split-prime
recurrence at \(7^2\). The literature comparison checks the hypotheses
of the Frobenian correlation theorem and the half-dimensional sieve
against the fixed-target problem. External specialist review and broader
historical priority assessment remain open.
