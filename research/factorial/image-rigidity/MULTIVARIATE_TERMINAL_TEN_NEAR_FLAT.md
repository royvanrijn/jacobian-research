# Shared terminal-ten separation with outer mean at least .97

Let two ordinary exponent columns be distinct and have only shared
unequal rows. Write \(u\) for their order-nine outer geometric mean
after positive endpoint normalization and removal of the binomial
factors. If \(u\ge97/100\), the ninth and tenth factorial moment
polynomials have disjoint root arguments under independent positive
rescalings. The dimension and positive exponents are arbitrary.
Rational Bézout gives coprimality in characteristic zero.

## Outer deficits and common angular scale

Put \(p=1-u\), and let \(q=1-w\) be the order-ten deficit.
Write \(r=q/p\) for \(p>0\). For lower-order reciprocals
\(x=d/(9\ell+h)\), their matches are \(x'=9x/(10-ax)\),
where \(a=h/d\in(0,1]\). Because \(u^{-2}\) is the product of
\(1+8x\), we have \(x\le c=(u^{-2}-1)/8\).
The matched ratios \(9x'/(8x)\) lie in
\[
               [81/80,\;81/(8(10-c))].
\]
Put \(b_0=81/(8(10-c))\). Expanding the outer products gives
\(w\le u(1+(1-u^2)/80)^{-1/2}\). For the other direction,
\((1+s)^{b_0}\ge1+b_0s\) gives \(w\ge u^{b_0}\).
Integrating the derivative of \((1+s)^{-1/2}\) gives
\[
       1+\frac{u(1+u)}{160(1+(1-u^2)/80)^{3/2}}
                  \le r\le b_0 .
\]
Both bounds hold uniformly at the worst endpoint \(u=.97\).
The checker uses the rational interval
\(10119/10000\le r\le5067/5000\), and verifies these inequalities
by rational comparison. It also checks that the upper-order
deficit stays below \(q_0=31/1000\).
These bounds hold along the deformation obtained by multiplying
all matched reciprocals at both orders by a common factor in
\((0,1]\). The individual matched ratios stay fixed, and the
lower bound on \(u\) improves as the factor decreases.

At order \(j\), use its own deficit \(h=p\) or \(q\), scale
\(\lambda_j=\sqrt h\), and transform
\[
 t=\frac{\lambda_j Z-1}{\lambda_j Z+1},\qquad
 C_j(Z)=\frac{(1+\lambda_jZ)^j}{2^j\lambda_j^j}P_j(t).
\]
The transformed polynomial has a continuous limit at zero
deficit as verified below.

## A grouped Taylor enclosure through deficit weight seven

For a factor of reciprocal \(x_i\) in an order \(j\) with
deficit \(h\), put
\[
 s_i=\tfrac12\log(1+(j-1)x_i),\qquad T_i=\tanh s_i,\qquad
 e_{j,k}=k(j-k)/(j-1).
\]
The reference multipliers are \((1-h)^{e_{j,k}}\).
The exact logarithmic correction to this reference is a sum
of signed even-power skew terms and unsigned odd-power
mean-defect terms. For \(a_k=\min(k,j-k)\), put
\(p_l=1-2l/(j-1)\) and
\(\omega_{k,l}=({\bf1}_{a_k\le l<j-a_k}-(j-2a_k)/j)/2\).
The even-power direction at power \(d\) is
\(-\sum_l \omega_{k,l}p_l^d\) for \(k\le j/2\), with opposite
sign at the reflected index. Its scalar is
\(\sum_i\epsilon_i T_i^d/d\), with the sign \(\epsilon_i\)
fixed by the row orientation. The odd-power direction at
\(d\ge3\) is \(\sum_{l<a_k}(p_l-p_l^d)\), with scalar
\(\sum_i T_i^d/d\). These follow by expanding the centered
Gauss logarithms and their reflected means.
[Grouped logarithmic coefficient formulas](MULTIVARIATE_TERMINAL_NINE_PRIVATE_LINEAR.md#grouping-the-shared-logarithmic-corrections).

Write the scalar at degree \(d\) as \(h^d\eta_d\).
Since \(\sum_i T_i\le\sum_i s_i=-\log(1-h)\le C h\) with
\(C=1/(1-q_0)\), one has \(|\eta_d|\le C^d/d\), and
\(\eta_d\ge0\) at odd \(d\). Retain the six directions
\(d=2,\ldots,7\), and write their sum in the log as \(A_k(h)\).
The omitted log term is bounded by \(h^8L_k\), with
\[
 L_k=\frac{(\sum_l|\omega_{k,l}|) C^8/8+a_k C^9q_0/9}
                 {1-(Cq_0)^2}.
\]
This sums the skew tail from \(d=8\) and the mean tail from
\(d=9\). It also applies to the endpoint indices, where both
terms vanish.

In the exponential of \(A_k\), keep every monomial in the
six directions whose weights sum to at most seven. A repeated
direction with multiplicity \(m_d\) contributes the
multinomial factor \(1/\prod_d m_d!\).
Write \(H_k(h)\) for the polynomial obtained from \(A_k(h)\)
by replacing each coefficient by its absolute upper bound.
If \(H=H_k(q_0)\), exponential powers of order at least four
are bounded at \(q_0\) by \(H^4/(24(1-H))\).
Add the absolute monomials of the first three powers whose
weights are at least eight, and divide by \(q_0^8\).
Their contribution is bounded by this constant times \(h^8\)
because every omitted monomial has nonnegative majorant and
weight at least eight. The effect of the log tail is bounded
by \(L_k h^8/(1-(H+L_kq_0^8))\). The checker verifies every
denominator positive. This gives an independent raw remainder
of order \(h^8\) only after all lower-weight monomials have
been kept as grouped directions.

For each grouped direction, the raw Cayley coefficient of
degree \(a\) is a rational linear combination of
\((1-h)^{e_{j,k}}\). Expand this combination through degree
eleven in \(h\) with exact rational binomial coefficients.
After its initial zeros have been removed, evaluate the remaining
polynomial by interval arithmetic. The tail is bounded by the sum
of the absolute twelfth coefficients times a geometric series,
since successive absolute binomial coefficients decrease from
that index onward. Include the direction's raw power of \(h\)
before dividing by \(\lambda_j^j\). The checker verifies that
all surviving powers of \(\lambda_j\) are nonnegative. Also,
the independent raw errors \(O(h^8)\) contribute only powers
\(\lambda_j^{16+a-j}\), which vanish at zero for these two orders.
The reference gives the only nonzero family at zero, with
imaginary Hermite variance \(1/(2(j-1))\) in each respective
scale. The checker evaluates the resulting rational coefficients
directly.

## Exact rays and sector counts

Leading-coefficient dominance confines the transformed roots
to \(|Z|<2\). This disk maps into the physical left half-plane
because \(4q_0<1\). The eight separating rays have angles
\(\pi-2\arctan(\sqrt p\,b)\), with rational levels \(b\)
piecewise affine in \(p\). In the order-nine scale the level is
\(b\), and in the order-ten scale it is \(b/\sqrt r\).
For deficit \(h\) and the corresponding level \(b_*\), the
ray has equation
\[
        (Y-b_*)(1+h b_*Y)=-h b_*x^2,\qquad 0<Y\le b_*.
\]
Three interval updates enclose \(Y\). Reflection of the rows
sends \(Z\) to \(-\overline Z\), so the nonnegative range
\(0\le x\le2\) suffices. Subdivide the scale, the real coordinate,
or the order-ten ratio until a strict projected sign holds.

To evaluate the projections, first try the coarse coefficient
intervals. The sharper grouped test expands each polynomial family
on a complex ray rectangle at its center \(Z_c\). It evaluates the
constant and first derivative there with interval arithmetic,
and bounds the remaining term by
\(\sup |F''|\,|Z-Z_c|^2/2\).
The supremum uses the coefficient absolute bounds and the
maximum modulus on the rectangle. The direction of a projection
may be proposed from the centers of the resulting complex
rectangles. The acceptance test only uses exact interval signs
after each of the grouped contributions is multiplied by its
scalar interval.

At zero deficit the checker evaluates the two Hermite anchors
on the imaginary axis, dividing the odd polynomial by its
central factor. The strict signs on the eight rays and the
endpoints place the nine upper roots in sectors of orders
\(10,9,10,9,10,9,10,9,10\), with one real nonic root.
The simultaneous reciprocal deformation preserves the
root disk and ray exclusions, hence preserves those counts.
Join the last upper and lower sectors through the negative
real axis. The sole nonic root there stays real by conjugacy.
The actual decic has no real root because its even factorial
moment is strictly positive for a nonzero real linear
combination of distinct monomials. This separates all arguments.

The [exact checker](verify_multivariate_terminal_ten_near_flat.py)
certifies 17,029 projection cells across 96 order, ray, and
subband checks. It compares each count with the exact replay
receipt. Run the full check or a partition from the repository
root. Deficit bounds are in thousandths.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_near_flat.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_near_flat.py 10 -1 20 30
~~~
