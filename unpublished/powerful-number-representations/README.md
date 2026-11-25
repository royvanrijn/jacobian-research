# Additive representations by powerful integers

Roy van Rijn · Draft in preparation

[PDF](paper.pdf) · [LaTeX source](main.tex) · [Metadata](metadata.json) · [Verification record](verification.json)

The headline theorem saves \(1/(153k)\) over \(2/k\) for
primitive diagonal \(k\)-full triples. For fixed offsets
\(p=r+u,\ q=r+v\), it saves \(1/(306r)\) over \(1/p+1/q\)
for sufficiently large \(r\). Section 9.4 combines
minimum-energy subset selection with adaptive amplification.
Fixing \(c\) gives a saving of \(1/(300k)\).
Choosing the fixed summand by residual-factor size proves
the headline bound. Section 9.3 gives the fixed-depth
divided-difference estimate, and Section 9.2 gives the
compression and the earlier \(k^{-3/2}\) saving.

The [parameter checker](../../research/arithmetic/verify_powerful_adaptive_determinant.py)
replays the strict inequalities at \(1/320\), \(1/300\)
and \(1/153\) and the prime-power ceiling calculation.
Run it from this directory with

~~~sh
python3 ../../research/arithmetic/verify_powerful_adaptive_determinant.py
~~~

The draft retains the fixed-offset graph estimates and the two-power
weighted graph bound of Section 9.1. Its exact certificates give
leading coefficients \(6\) for even \(k\) and \(16/3\) for odd \(k\)
at scale \(k^{-2}\). Section 6.3 gives the earlier bounds from
common-divisor grouping.

The manuscript also proves pointwise bounds for powerful sums and
large coprime differences, including the \(N^{3/10+\varepsilon}\)
bound for square-full representations in Erdős #943, and a bound
\(O_{k,\varepsilon}(X^{2/(k+1)+\varepsilon})\) for primitive
diagonal triples at every fixed \(k\ge3\).

The written proofs and imported theorem hypotheses were checked internally.
Independent PDF-only review, external specialist review and an exhaustive
priority search remain pending.

Build from this directory with

~~~sh
latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex
~~~

## Rational certificates

The [certificate tool](../../research/arithmetic/verify_powerful_fixed_offset_constraints.py)
proposes certificates with SciPy and verifies them with standard-library
rational arithmetic. From this directory, generate the first grouped-graph
obstruction and verify it with

~~~sh
python3 ../../research/arithmetic/verify_powerful_fixed_offset_constraints.py generate 3 1 62 > certificate.json
python3 ../../research/arithmetic/verify_powerful_fixed_offset_constraints.py verify certificate.json
~~~

The witness is supported on explicit coordinates. The mathematical argument
in the paper checks the entire residue class.

Verify the [certificate for multiples of six](../../research/arithmetic/powerful_even_order_divisible_by_six_certificate.json.gz)
from this directory with

~~~sh
python3 ../../research/arithmetic/verify_powerful_fixed_offset_constraints.py verify \
  ../../research/arithmetic/powerful_even_order_divisible_by_six_certificate.json.gz
~~~

To generate the even-diagonal certificates with SciPy, run

~~~sh
python3 ../../research/arithmetic/verify_powerful_fixed_offset_constraints.py \
  generate 0 0 62 --even-diagonal --cutoff 4 --saving 8/3
python3 ../../research/arithmetic/verify_powerful_fixed_offset_constraints.py \
  generate 0 0 60 --even-diagonal --cutoff 5 --saving 4
~~~

The representative order selects the permitted congruence constraints.
The proof applies throughout the stated even orders or multiples of six.

For the two-power bounds, replay the
[even certificate](../../research/arithmetic/powerful_two_power_even_certificate.json.gz)
and the [odd certificate](../../research/arithmetic/powerful_two_power_odd_certificate.json.gz)
from this directory with

~~~sh
python3 ../../research/arithmetic/verify_powerful_fixed_offset_constraints.py verify \
  ../../research/arithmetic/powerful_two_power_even_certificate.json.gz
python3 ../../research/arithmetic/verify_powerful_fixed_offset_constraints.py verify \
  ../../research/arithmetic/powerful_two_power_odd_certificate.json.gz
~~~

Generation uses the same tool's `--two-power-diagonal` option
with `--cutoff 6` and `--saving 6` at representative order 62,
or `--saving 16/3` at representative order 61.
