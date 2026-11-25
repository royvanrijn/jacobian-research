# The Strong Factorial Conjecture for Univariate Binomials

**Unpublished manuscript for mathematical review**

[PDF](strong-factorial-binomials.pdf) · [LaTeX source](main.tex) ·
[Electronic supplement](supplement/README.md) ·
[Metadata](metadata.json) · [Verification record](verification.json)

The paper proves the Strong Factorial Conjecture for univariate
binomials over every field of characteristic zero. The proof compares
a complex bound for a Padé resultant with divisibility from short
factorial blocks. At the remaining orders, Newton polygons restrict
the degrees of common factors and exact arithmetic rules them out.
The paper includes a worked example of the prime sieve.

## Build

Run the following from this directory with TeX Live and latexmk.

    latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=strong-factorial-binomials main.tex

The LaTeX source bundle consists of the eleven TeX files in this
directory. It includes the bibliography. The electronic supplement is
part of the proof bundle and is read separately.

## Review

The verification record identifies the reviewed edition, exact
arithmetic replays, and findings from the PDF-only internal reader.
External specialist review and independent external replay are open.
