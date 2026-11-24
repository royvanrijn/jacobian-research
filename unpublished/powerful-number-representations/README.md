# Additive representations by powerful integers

Roy van Rijn · Draft in preparation

[PDF](paper.pdf) · [LaTeX source](main.tex) · [Metadata](metadata.json) · [Verification record](verification.json)

The manuscript proves pointwise bounds for sums of two positive \(k\)-full
integers and for large coprime differences. The signed graph estimates also bound primitive \(k\)-full triples \(a+b=c\le X\) by \(O_{k,\varepsilon}(X^{2/(k+1)+\varepsilon})\) for fixed \(k\ge3\).
For fixed offsets between unequal orders, the mixed-exponent graph gives a saving of leading size \(2/r^2\) under the stated gcd condition.
For consecutive orders the leading saving is \(2r^{-2}\) for all sufficiently large \(r\), using a grouped square variable in the even case.
It covers the square-full case of Erdős #943 and every fixed order \(k\ge3\).

The written proofs and imported theorem hypotheses were checked internally.
Independent PDF-only review, external specialist review and an exhaustive
priority search remain pending.

Build from this directory with

~~~sh
latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex
~~~
