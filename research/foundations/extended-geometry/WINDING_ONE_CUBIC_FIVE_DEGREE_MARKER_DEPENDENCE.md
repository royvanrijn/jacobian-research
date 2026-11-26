# A marked-moment dependence in the degree-five cubic boundary

Let \(K\) have characteristic zero. Use the sphere functional \(\Lambda\)
on the regular winding-one phase
\[
 f=x^{-1}(C+Bt+Dt^2+Et^3),\qquad E=e(x-r)^4(x-s),
\]
where \(e,r,s\ne0\), \(r\ne s\), and the least sphere degree of \(f\)
is seven. We may extend \(K\) to name the two roots. As in the
[boundary theorem](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md), put
\(v=x-r\). If \(f\) is eventually pure, the lower coefficients allow
the regular polynomial presentation
\[
 f=\frac v x R(x,vt),\qquad
 R(x,T)=U_0+U_1T+U_2T^2+UT^3,
 \quad U=e(x-s),\quad \deg U_j\le 7-2j.
\]

**Proposition.** If \(f\) is eventually pure, then there exist scalars
\(a,b\) in the coefficient field containing the roots such that, for
every \(l\ge0\),
\[
 \Lambda\bigl(t(x-r)^2 f^l\bigr)
  =a\Lambda\bigl(x(x-r)f^l\bigr)
   +b\Lambda\bigl(x^2(x-r)f^l\bigr).                     \tag{1}
\]
In particular, three linearly independent moment columns for these
markers disprove eventual purity. The statement is an exact sequence
dependence and includes the cases in which the two columns on the
right are dependent.

## The second prime progression

Work first over a number field. Take a sufficiently large good prime
of residue characteristic \(p\equiv2\pmod3\), and set
\(n=(2p-1)/3\). The height of \(f^n\) is \(2p-1\). For any fixed
\(l\ge0\), reduce the identity for \(p\Lambda(f^{lp+n})\).
Frobenius in the residue field replaces \(f^{lp}\) by the \(l\)-th
power of \(f(x^p,t^p)\) with all coefficients raised to the
\(p\)-th power. The only heights of \(f^n\) that can contribute
are \(p-1\) and \(2p-1\), since the denominator divided by \(p\)
must be integral. For \(p>3l+2\) its scaled denominator is a unit.

If the height in the other factor before Frobenius is \(j\), then
height \(p-1\) contributes the interval value of \(t^j\), with
denominator \(j+1\) and the parity test for \(j\). Height \(2p-1\)
contributes the interval value of \(t^{j+1}\), with denominator
\(j+2\) and its parity test. Pure angular matching uses just the
multiples of \(p\) in the Laurent exponents of the two height rows.

At height \(p-1\), the numerator in \(f^n\) is
\[
 v^{n+p-1}[T^{p-1}]R^n=(x^p-r^p)
              v^{n-1}[T^{p-1}]R^n\pmod p.
\]
Put \(S_1=v^{n-1}[T^{p-1}]R^n\). Its degree is at most
\(8n-1-2(p-1)=(10p-5)/3\). The selections
\(c_i=[x^{n+ip}]S_1\) can occur only for \(i=0,1,2\), since
\(n<p\). Hence the angular polynomial for this row, in the
variable whose powers are the selected multiples of \(p\), is
\((X-r^p)(c_0+c_1 X+c_2X^2)\).

At height \(2p-1\), the numerator is \(E^n\). The power of \(v\)
satisfies \(4n=2p+n-1\), so
\[
 E^n=(x^p-r^p)^2v^{n-1}U^n\pmod p.
\]
The degree of \(v^{n-1}U^n\) is \(2n-1<p+n\), so the
angular polynomial for this row is \(d_p(X-r^p)^2\), where
\(d_p=[x^n]v^{n-1}U^n\).

At \(l=0\), eventual purity kills the scaled moment. The high
row integrates an odd power and only \(-r^pc_0\) survives from
the low row, so \(c_0=0\). At each fixed \(l>0\) and all
sufficiently large good primes, eventual purity gives
\[
 c_1\Lambda(xv f^l)^p+c_2\Lambda(x^2v f^l)^p
              +d_p\Lambda(tv^2 f^l)^p=0                     \tag{2}
\]
in the residue field. The powers on the moments denote residue-field
Frobenius. This accounts for primes of any absolute residue degree.

## Exactness at primes and the characteristic-zero obstruction

Suppose (1) fails. Some finite selection of rows then witnesses that
the third column is outside the span of the first two. Indeed the
kernels of the successive finite row matrices stabilize in the
three-dimensional vector space of column coefficients. For all but
finitely many primes the same finite matrix witnesses noncontainment,
by its nonzero minors and the identities of its vanishing minors.
Frobenius preserves this noncontainment. Apply (2) to these fixed
rows to obtain \(d_p=0\) for all sufficiently large good primes with
\(p\equiv2\pmod3\).

Consider the smooth projective curve with affine function field
\[
 W^3=F_0=x^2vU,qquad \eta=\frac{dx}{vW}.
\]
The four branch points \(0,r,s,\infty\) are fully ramified. The curve
has genus two by Riemann–Hurwitz. The differential \(\eta\) has a
unique pole, of order two, at the point above \(r\). For example
there \(v\) has order three, \(W\) order one, and \(dx\) order two.
At zero these respective orders for \(v,W,dx\) are zero, two, two;
at \(s\) they are zero, one, two; and at infinity they are
\(-3,-4,-4\). Therefore \(\eta\) cannot be exact in characteristic
zero. A rational primitive would have a unique simple pole and hence
would give a degree-one map of a genus-two curve to the line.

For comparison, take any sufficiently large good prime. In
characteristic \(p\) the Cartier operator on a polynomial differential
in \(x\) selects the coefficients at exponents \(ip+p-1\).
Multiplication by a \(p\)-th power factors out of Cartier as its
\(p\)-th root. These facts follow by writing a rational differential
as \(\sum_{j=0}^{p-1}a_j^p x^j dx\), since \(x\) is separating.
In particular a differential with zero Cartier has a rational
primitive by integrating the terms with \(j<p-1\).

For \(p=3k+1\), write
\(\eta=(W^{-1})^p F_0^k\,dx/v\). The only possible Cartier
coefficient is \([x^{p-1}]F_0^k/v=[x^k]v^{k-1}U^k\), by the
degree bound. This is zero by the unmarked prime-height test in the
[fixed-deficit calculation](WINDING_ONE_CUBIC_FIVE_DEGREE_BOUNDARY_MOMENTS.md#the-selected-prime-coefficient)
at deficit zero. For \(p\equiv2\pmod3\), the identity
\(2p=3n+1\) instead gives \(\eta=(W^{-2})^p F_0^n dx/v\).
Here the numerator polynomial has degrees between \(2n\) and
\(4n-1\), so only the coefficient at \(2p-1\) contributes. It is
\([x^n]v^{n-1}U^n=d_p=0\). Thus \(\eta\) is exact on the
reduction of the curve at all but finitely many prime ideals.

On the open curve where \(\eta\) is regular, take the flat additive
principal bundle given by the equation \(dz=\eta\). An exact reduction
has zero \(p\)-curvature. Bost's Theorem 2.9 for flat principal
bundles with solvable neutral structural group therefore makes this
bundle trivial on a finite étale cover. The group here is
\(\mathbb G_a\). A horizontal section gives a rational primitive
of the pullback. Its field trace divided by the covering degree gives
a rational primitive of \(\eta\) on the original curve. This contradicts
the pole computation. [Bost, Theorem 2.9](https://numdam.org/item/10.1007/s10240-001-8191-3.pdf)

For an arbitrary characteristic-zero field, failure of (1) is again
detected by a finite matrix. Specialize the finitely generated
coefficient algebra to a number field, preserving its nonzero minors
and the requisite branch parameters. The vanishing minors and all
eventual pure moment identities persist because they are rational
polynomial identities in the coefficients. The number-field case
gives the same contradiction. This proves (1) over the stated field.
