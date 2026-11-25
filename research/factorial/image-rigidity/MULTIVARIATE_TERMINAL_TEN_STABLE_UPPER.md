# Shared terminal-ten separation with outer mean at least .93

Let two distinct ordinary exponent columns have only shared unequal
rows. Let \(u\) be their order-nine outer geometric mean after positive
endpoint normalization and removal of the binomial factors. If
\(u\ge93/100\), the ninth and tenth factorial moment polynomials have
disjoint root arguments under independent positive rescalings.
The dimension and positive exponents are arbitrary. Rational Bézout
gives coprimality in characteristic zero.

Use the same grouped Taylor enclosure and Cayley contours through
deficit weight seven as the terminal-ten near-flat proof. Put
\(p=1-u\), and put \(q=1-w\) for the order-ten outer mean \(w\).
For \(p>0\), write \(r=q/p\). At the new lower endpoint \(u=.93\),
the matched reciprocal bound \(c=(u^{-2}-1)/8\) and the outer product
inequalities give
\[
  1+\frac{u(1+u)}{160(1+(1-u^2)/80)^{3/2}}
       \le r\le\frac{81}{8(10-c)} .
\]
The checker encloses the resulting ratios by
\([20223/20000,5073/5000]\) and uses
\(\sqrt p\le2645752/10^7\), \(q\le q_0=9/125\).
The log moment majorants and every grouped raw tail bound use
\(C=1/(1-q_0)\). As before, cancellation is checked before
division by the Cayley scale.
[Near-flat proof and exact coefficient enclosure](MULTIVARIATE_TERMINAL_TEN_NEAR_FLAT.md).

A new ray table covers every deficit from zero to \(.07\).
The checker verifies the nine upper sectors at the Hermite anchor
uniformly in the enlarged ratio interval. Leading dominance keeps
the transformed roots in \(|Z|<2\), which maps into the physical
left half-plane since \(4q_0<1\). The same centered complex Taylor
tests exclude the rays on the larger domain. Thus the reciprocal
deformation preserves the sector counts, and the positivity of
the even factorial moment separates the arguments in the sector
joined through the negative real axis.

The extended parameters, rays and expected counts are supplied to
the near-flat checker by the [replay wrapper](verify_multivariate_terminal_ten_stable_upper.py).
Its six exact partitions certify 170,066 projection cells across
224 order, ray, and subband checks. Run all partitions or select
an order and a range of deficit knots in thousandths.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_stable_upper.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_stable_upper.py 10 -1 50 70
~~~
