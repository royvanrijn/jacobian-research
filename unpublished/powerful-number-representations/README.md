# Primitive sums of powerful integers and the abc endpoint

Roy van Rijn · Draft in preparation

[PDF](paper.pdf) · [LaTeX source](main.tex) · [Metadata](metadata.json) · [Verification record](verification.json)

For ordered primitive \(p,q,r\)-full triples, Section 11 derives
the exponent \((1/3-\delta)(1/p+1/q+1/r)+\varepsilon\) for
sufficiently large fixed orders, with an absolute \(\delta>0\).
For equal orders this is \((1-3\delta)/k+\varepsilon\).
The transfer is proved in the manuscript under two explicit
small-budget box hypotheses. The analytic dependencies are the
[limiting one-third box bound](../../research/arithmetic/ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)
and the [balance and nonatomic-mass theorem](../../research/arithmetic/ABC_WIDE_COFACTOR_CONGRUENCES.md#necessary-conditions-at-normalized-exponent-13).
Review of the full analytic proofs is still needed for an
unconditional presentation of this bound.

Section 12 gives the [reverse transfer](../../research/arithmetic/POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target).
The least common multiplier makes the target \(M\)-full and
preserves its primitive quotient. A uniform power saving below
\(1/M\) for the resulting balanced count with \(M=\lceil K/\lambda\rceil\)
is equivalent under the endpoint hypotheses and the theorem's
parameter inequalities to a uniform coefficient gap below \(1/3\)
for the small-budget \(abc\) count. Proving either saving remains open.

The independent determinant theorem saves \(1/(153k)\) over \(2/k\)
for sufficiently large diagonal orders and \(1/(306r)\) over
\(1/p+1/q\) for fixed offsets \(p=r+u,\ q=r+v\) above sufficiently
large \(r\). Section 9.4 combines minimum-energy subset selection
with adaptive amplification. Fixing \(c\) gives a saving of
\(1/(300k)\), and choosing the fixed term by residual-factor size
gives the stated \(1/(153k)\). Sections 9.2 and 9.3 give the
compression and fixed-depth estimates.

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

The transfer deductions and their quantifiers were checked internally
against the cited and stated inputs. The determinant parameter check,
the three retained graph certificates and the PDF layout were checked
for this edition. An internal PDF-only reader pass prompted clearer
conditional wording and direct links to the rational verification files.
External specialist review and a broader literature comparison remain pending.

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
