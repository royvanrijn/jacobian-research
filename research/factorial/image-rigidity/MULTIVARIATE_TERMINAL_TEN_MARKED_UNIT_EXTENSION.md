# Lower outer-mean thresholds for two marked terminal-ten rows

Let two distinct ordinary exponent columns have only shared unequal
rows, and let \(u\) be the order-nine outer mean after positive
endpoint normalization and removal of binomial factors. The ninth
and tenth moments have disjoint root arguments under independent
positive rescalings in either of the following cases.

- A row has entries \((3,4)\) in either orientation and \(u\ge79/100\).
- A row has entries \((4,5)\) in either orientation and \(u\ge21/25\).

All remaining shared exponents and the dimension are arbitrary.
Rational Bézout gives coprimality in characteristic zero.

The same normalized marked cores, exact circle tests, and
centered Taylor ray projections apply with these thresholds.
For a marked row with smaller exponent \(a\), write
\(v=u/u_a\), where \(u_a^2=(9a+1)/(9a+9)\). The new lower
bounds \(v\ge447/500\) for \(a=3\) and \(v\ge463/500\) for
\(a=4\) follow by exact comparison of the squares.
The unit-row residual condition gives reciprocal caps
\(1/37\) and \(1/55\), respectively. For nonunit residual
rows the bound \((1-v)/(9v-1)\) gives smaller caps, and
the product \(v^{-2}\) further bounds each order-nine
reciprocal by \((v^{-2}-1)/8\). The checker verifies the
thresholds and cap comparisons. These bounds hold along
the same real shift of residual exponents towards the
marked core. The circle and ray exclusions transfer the
core sector counts through the deformation.
[Marked-core proof and exact interval method](MULTIVARIATE_TERMINAL_TEN_MARKED_UNIT_ROWS.md).

The extended tables certify 58,930 ray cells across 736
order, ray, and subinterval checks. The replay verifies
the four marked-core anchors and all the expected counts.
[Exact replay wrapper](verify_multivariate_terminal_ten_marked_unit_extension.py).

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_marked_unit_extension.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_marked_unit_extension.py 3 10
~~~
