# Common-factor degree bounds from deep factorial blocks

Fix two ordinary factorial moment polynomials at consecutive
orders. Write \(m\) for a marked order and \(m'\) for the other.
All valuations below are normalized by \(v_p(p)=1\).
For each row strip the scalar factorial of its smaller
exponent at each order. This preserves the roots.

## Marks in arbitrary blocks

In an increasing row of smaller exponent \(A\) and gap \(d\),
a selected factor \(mA+h\), \(1\le h\le md\), first occurs at
coefficient index \(r=\lceil h/d\rceil\), and occurs at every
later index. In a decreasing row a factor in block \(s\)
occurs at indices \(k\le m-s\).

Select one or two factors, with at most one of each orientation,
at the same prime \(p\). Suppose the indices at which none of
the selected factors occurs form a nonempty interval
\(I=[a,b]\cap\mathbb Z\). Let \(E\) be the least valuation
of the selected factors. At order \(m\) delete the selected
factors where they occur and let \(U\) be an upper bound
for the range of the resulting coefficient valuations.
Let \(L\ge0\) bound the absolute values of all slopes in the
other Newton polygon. One may take \(L\) to be the range of
the coefficient valuations at order \(m'\).

**Degree bound.** If
\[
                          E>U+L\max(a,m-b),
\]
then every rational common factor has degree at most \(b-a\).
In particular, width zero or one certifies physical coprimality.
An increasing factor in the second block can therefore suffice.
With opposite marks in blocks \(r,s\) and \(r+s\ge m+2\), the
width is \(r+s-m-2\).

Indeed for any slope \(|\xi|\le L\), compare the marked ordinate
at \(k\notin I\) with that at the nearest index \(i\in I\).
After subtracting \(\xi k\) and \(\xi i\), respectively, the
difference is at least \(E-U-L|k-i|>0\).
Therefore every marked edge with a slope in \([-L,L]\) has both
endpoints in \(I\). Their total widths are at most \(b-a\).
Newton polygons count root valuations by edge widths, and a common
factor can only use these roots. This is the ordinary Newton polygon
product rule, applied over the algebraic closure.
[Guàrdia, Montes and Nart, Newton polygons of higher order](https://arxiv.org/pdf/0807.2620)

The even ordinary moment is
strictly positive at every real argument by the exponential
integral of an even power. Hence rational common factors have
even degree, proving the coprimality assertion.

## A divisor condition for an existing mirrored pair

Suppose the original profile contains
\((A,A+d),(A+d,A)\), and take \(m\ge4\).
Put \(l=\lfloor m/2\rfloor\) and \(r=\lfloor(m+3)/2\rfloor\).
Choose any integer \(h_0\) with
\[
                     d(r-1)<h_0\le dr,\qquad m\nmid h_0 .
\]
Such an offset exists. If \(d=1\), take \(h_0=r<m\).
If \(d\ge2\), the block contains consecutive integers.

For each order \(j=m,m'\), and for every row other than the
selected pair, write \(\ell_i,d_i\) for its smaller exponent
and absolute gap and put
\[
 \begin{split}
 B_j&=j!\prod_{i\ {\rm retained}}\prod_{h=1}^{j d_i}(j\ell_i+h),\\
 C_m&=B_m\prod_{\substack{1\le h\le md\\h\ne h_0}}|h-h_0|^2,\\
 C_{m'}&=B_{m'}\prod_{h=1}^{m'd}|mh-m'h_0|^2.
 \end{split}
\]
The last product has no zero, since \(\gcd(m,m')=1\) and
\(m\nmid h_0\).

**Divisor condition.** A physical common root forces
\[
                           mA+h_0\mid C_m C_{m'}^{\,l}.
\]
The right side depends on the retained profile, \(d\), and
the two orders. It is independent of \(A\). This test
permits primes dividing the marked order.

To prove the condition, take a prime \(p\mid mA+h_0\) and set
\(E=v_p(mA+h_0)\), \(U=v_p(C_m)\), \(L=v_p(C_{m'})\).
Suppose \(E>U+lL\). For \(h\ne h_0\) in the marked ramp
one has \(v_p(mA+h)=v_p(h-h_0)\). The selected factor is used
exactly once outside \(I=[m-r+1,r-1]\cap\mathbb Z\),
whereas it is absent from both rows at every index in \(I\). After
removing it wherever it occurs, the other valuations lie in \([0,U]\),
since the retained
coefficient factors divide \(B_m\).

In the other order use
\[
                m(m'A+h)=m'(mA+h_0)+(mh-m'h_0).
\]
The valuation of the first term is strictly larger than that
of the second. Thus
\(v_p(m'A+h)=v_p(mh-m'h_0)-v_p(m)\) for every sampled \(h\).
Each is nonnegative because \(m'A+h\) is an integer.
In particular the entire other valuation array lies in
\([0,L]\). The interval \(I\) has width zero at even \(m\)
and width one at odd \(m\). Its farthest outside index is
at distance \(l\). The degree bound excludes a common factor.
Hence a common root requires \(E\le U+lL\) at each prime,
which is the stated divisibility.

At \(N\ge5\), one may mark either order. If offsets \(h_m,h_{m'}\)
give respective proposed divisors \(M_m,M_{m'}\), then a common
root also forces
\[
        (mA+h_m)(m'A+h_{m'})
        \ \le\ \operatorname{lcm}(M_m,M_{m'})
                           |mh_{m'}-m'h_m|.
\]
Indeed the gcd of the left factors divides the nonzero determinant
on the right. The lcm of the left factors divides the lcm of
the proposed divisors.

For a single marked row the same calculation gives an
arbitrary-block divisor test. With \(r=\lceil h_0/d\rceil\),
\(m\nmid h_0\), take the definitions above with that one row
omitted from \(B_j\) and with both offset products to the
first power. Reverse the columns if necessary so the row
increases. The unmarked interval is \([0,r-1]\), so a rational
common factor of degree at least \(r\) forces
\[
                         mA+h_0\mid C_m C_{m'}^{m-r+1}.
\]
For the first two blocks even-degree parity makes this condition
necessary for any physical common root.

For an exact small-prime illustration, take the original pair
\((2,3),(3,2)\) at orders \(5,6\) and \(p=2\). The marked
factor is \(6\cdot2+4=16\). The orders five and six have
stripped valuation arrays, respectively,
\[
                     (3,3,5,5,3,3),\qquad
                     (6,6,6,4,6,6,6).
\]
The order-five polygon is horizontal and the order-six
polygon has no horizontal edge. Here the local comparison
uses the exact arrays and improves the sufficient uniform
divisor bound.

For a fixed retained profile and fixed \(d\), choose one allowed
offset in the even order of every pair with \(4\le N<J\), where
\(J\) is any integer at least the
[absolute ordinary terminal cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Taking \(A\) large enough that \(mA+h_0>C_m C_{m'}^{\,l}\)
at each of these finitely many pairs certifies the original
profile in every degree. If the retained bases also vary and are
at most \(B\), with absolute gap sum \(G\), then
\(C_m C_{m'}^l\) is bounded by an explicit constant times
\((1+B)^{(m+l m')G}\), by its product definition.
This gives a native all-degree region for independently varying
retained bases. Terminal three follows from
[quadratic nondivisibility](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#an-inequality-using-the-endpoint-coefficients),
terminal two from even positivity, and higher orders from the cutoff.
All physical conclusions transfer to each characteristic-zero
field by a rational Bézout identity.

The [single-row deep marker test](MULTIVARIATE_NATIVE_DEEP_PRIME_DIVISOR_GATE.md)
uses the first block of one selected row. The
[opposite large-prime test](MULTIVARIATE_OPPOSITE_MARKED_ROW_PRIME.md)
uses a narrow interval where both marks are absent and the other
array is a unit array. The valuation range comparison here
allows nonzero slopes in the other array and simultaneously
marks the two reflected rows even when their factors coincide.
