# Additive representations by powerful integers

Roy van Rijn · Draft in preparation

[PDF](paper.pdf) · [LaTeX source](main.tex) · [Metadata](metadata.json) · [Verification record](verification.json)

The headline theorem gives a leading saving of \(2/r^2\) over
\(1/(r+u)+1/(r+v)\) for primitive triples of orders \((r+u,r+v,r)\)
at every fixed pair of offsets \(u\ge v\ge0\). It uses a uniform
short rational congruence bound over powerful moduli. The exponent
error in the manuscript is \(O_{u,v}(r^{-21/10})\).

For all sufficiently large diagonal orders \(k\), the leading
saving over \(2/k\) is \(6/k^2\) for even \(k\) and \(16/(3k^2)\)
for odd \(k\), with an exponent error \(O(k^{-21/10})\).
Section 9.1 proves the weighted graph bound. It gives a leading
coefficient of \(5\) in both parities before the exact rational
optimization. Section 6.3 gives the earlier bounds obtainable
from common-divisor grouping alone.

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
