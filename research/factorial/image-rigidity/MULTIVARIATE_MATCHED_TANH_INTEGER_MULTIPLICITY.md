# Integer row counts in the joint matched tanh data

The exact joint tanh data at one consecutive pair determine all
unequal ordinary rows up to permutation. Finite divisibility and
integrality constraints characterize which joint data come from
such rows. Equal rows supply only scalar factors and are omitted.
The characterization applies to private and shared rows.

## Recovering the Gauss grid

Fix \(n\ge2\) and \(N=n+1\). In a row of gap \(d\ge1\) and smaller
exponent \(\ell\ge0\), write \(u=\ell/d\) and \(a=h/d\) at the
Gauss factor \(1\le h\le d\). The exact order-\(j\) tanh value is
\[
              T_j=\frac{j-1}{j(2u+1)+2a-1}.
\]
It follows from the
[matched tanh formula](MULTIVARIATE_MATCHED_TANH_MOMENTS.md)
by substituting the reciprocal \(d/(j\ell+h)\).
Consequently a pair \((T,T')=(T_n,T_N)\) determines the
two grid coordinates:
\[
 \begin{aligned}
  u&=\frac12\left(\frac{n}{T'}-\frac{n-1}{T}-1\right),\\
  a&=\frac12\left(1+\frac{n^2-1}{T}-\frac{n^2}{T'}\right).
 \end{aligned}
\]
These coordinates are rational and satisfy \(u\ge0,\ 0<a\le1\).
Conversely the displayed formula for \(T_j\) recovers the pair
from its coordinates.

The data consist of an unsigned counting measure \(\mathcal U\)
and a signed counting measure \(\mathcal M\) on pairs
\((T,T')\in(0,1)^2\), with the signed weights given by orientation.
The two orientation measures are
\((\mathcal U+\mathcal M)/2\) and \((\mathcal U-\mathcal M)/2\).
The following criterion applies separately to each of them.
If the data are expressed as all mixed moments of both measures,
they determine the measures uniquely among finitely supported
measures. Indeed, for a point in the union of the two candidate
supports, a product of squared distances to the other points,
divided by its value at that point, gives a polynomial that
extracts the weight at that point.

## Exact multiplicity test

Fix one orientation and one value of \(u=L/D\) in lowest terms,
with \(u=0\) written as \(0/1\). For each \(e\ge1\), consider all
values \(a=b/e\) with \(1\le b\le e\) and \(\gcd(b,e)=1\).
Missing weights are zero. An ordinary profile has the same
multiplicity \(C_e(u)\) at each of these values. In fact a row
of gap \(d\) at that ratio contributes the complete grid
\(a=1/d,\ldots,d/d\), so it contributes one to each of the
listed weights precisely when \(e\mid d\). If \(\nu_d(u)\)
is the number of rows of gap \(d\) at this ratio, then
\[
               C_e(u)=\sum_{d:\ e\mid d}\nu_d(u).
\]
Thus the data satisfy the additional conditions
\[
    \nu_d(u)=\sum_{k\ge1}\mu_{\mathrm{Mob}}(k)C_{dk}(u)
                        \ \in\ \mathbb Z_{\ge0},\qquad
    \nu_d(u)=0\quad\text{if }D\nmid d,
\]
where \(\mu_{\mathrm{Mob}}\) denotes the ordinary Möbius function.
The sums are finite because the support is finite. The formula
follows from \(\sum_{k\mid r}\mu_{\mathrm{Mob}}(k)=\mathbf1_{r=1}\).

Conversely, suppose both orientation measures have rational grid
coordinates in the ranges above, with equal weights at equal
reduced denominators for each \(u\), and the inverse counts satisfy
these conditions. Insert \(\nu_d(u)\) copies of the row of gap
\(d\) and smaller exponent \(Ld/D\) in the chosen orientation.
This smaller exponent is an integer by the divisibility condition.
Inverting the finite divisor sums shows that the constructed
rows have exactly the given weights at every point. This also
proves uniqueness of the row multisets.

In particular, the weight at \(a=1\) for a fixed \(u\) equals the
number of rows of that ratio and orientation, so it is an
integer. The full Möbius condition is stronger. For example,
weights one at each of \(1,1/2,1/3,2/3\) give
\(C_1=C_2=C_3=1\) and inverse count \(\nu_1=-1\).

## Constraints from just the inverse moments

For one row the Gauss grid gives the exact marginal identity
\[
               (j-1)\sum_{h=1}^d T_{j,h}^{-1}
                         =j(2\ell+d)+1.
\]
For one orientation put \(A=(n-1)\sum T^{-1}\),
\(B=n\sum (T')^{-1}\), and let \(F=\sum 1\) count its
Gauss factors. Then
\[
    B-A=\sum_{\rm rows}(2\ell+d),\qquad
    N A-n B=\#\{\text{rows}\},\qquad
    B-A-F=2\sum_{\rm rows}\ell .
\]
In particular the middle expression is an integer and the last
is a nonnegative even integer. These aggregate tests do not
require knowing which atoms are matched.

In the [weighted angular collision](MULTIVARIATE_WEIGHTED_TANH_ANGLE_OBSTRUCTION.md)
there is one ratio \(u=1/30\) with a complete gap-\(30\) grid of
weight \(\mu=\kappa/30\) per atom in each orientation.
Its \(a=1\) atom has noninteger multiplicity \(\mu\). The
marginal expression \(N A-n B\) also equals \(\mu\) in each
orientation. Adding any number of ordinary rows cannot turn
that expression into an integer. Even adding fractional weights
at other ratios cannot change the defective \(a=1\) count at
\(u=1/30\).

This supplies an integer constraint that the homogeneous matched
tanh inequalities miss. It does not give a degree-uniform angular
gap for the integer lists themselves.

