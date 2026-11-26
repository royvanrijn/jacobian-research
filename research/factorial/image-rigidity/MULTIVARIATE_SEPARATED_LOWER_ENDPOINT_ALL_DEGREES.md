# All-degree coprimality from a separated lower endpoint

Fix integers \(A\ge1\), \(B\ge0\) and a rational
\(0<\gamma\le1\). In an original ordinary two-monomial profile,
omit equal rows and choose one unequal row. Orient the columns
so that the chosen row is \((a,a+D)\), and write \(M=a+D\).
Select any group containing this row whose gaps are at most \(AD\)
and whose larger exponents are in \([M-BD,M+BD]\).
Suppose every decreasing row in the group has smaller exponent
at least \(a+\gamma D\). There is no count, multiplicity or
minimum-gap requirement in the group.

Take \(J\ge4\) at least the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
and put \(H=J^2(A+B+2)\). For rows outside the group, let
\(M_i\) be the larger exponent and \(d_i\) the absolute gap.
Sum over those with \(JM_i>HD\), putting
\[
  \mathcal B=\sum_{\rm increasing}d_i\log(JM_i)
                  +\frac1J\sum_{\rm decreasing}d_i\log(JM_i).
\]

**Theorem.** For
\[
 D\ge\max\{\lceil2J/\gamma\rceil,\ 144,\ 2^{\lceil512JH/\gamma\rceil}\},
 \qquad
                 \mathcal B\le\frac{\gamma D\log M}{32J^3},
\]
the original factorial moment polynomials are coprime in every
consecutive pair over every characteristic-zero field. The
endpoint \(M\) need not be the global maximum.
Rows outside the group with upper exponent at most \(HD/J\) have
no count or gap restriction. The numerical value of \(J\) is unevaluated.

## A first-block progression at the chosen lower endpoint

The [curvature twelve theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
allows total curvature below \(12\) throughout. In particular
\(M>D^2/12\). Fix a pair with \(4\le m=N<J\), \(m'=m-1\).
After stripping each row's smaller-end factorial, take the marks
\[
                         x=ma+h,\qquad
        1\le h\le\gamma D,\qquad h\equiv1\pmod m .
\]
The progression has count \(k\ge\gamma D/m-1\ge\gamma D/(2J)\),
span less than \(D\), and step \(m\) coprime to every term.
Every \(x\) is greater than \(m(M-D)>M\) under the stated bounds.

Consider a prime \(p>HD\) dividing a mark \(x\), and put
\(E=v_p(x)\). A full stripped factorial ramp in group row \(i\)
at order \(j\) consists of integers \(jM_i-t\),
\(0\le t<j d_i\). At order \(m\) use the offset
\[
                       mM_i-t-x=m(M_i-a)-t-h .
\]
At the other order use
\[
                    m(m'M_i-t)-m'x
                         =mm'(M_i-a)-mt-m'h .
\]
Both magnitudes are bounded by \(HD\). The latter is nonzero
modulo \(m\) since \(h\equiv1\pmod m\), so all group factors
at order \(m'\) are \(p\)-units. The former shows that only
actual copies of \(x\) in the marked ramps can contribute.
All increasing rows have an empty ramp at coefficient index zero.
Every decreasing group row has all factors strictly above
\(m(a+\gamma D)\ge x\). Thus every group contribution at
index zero is a \(p\)-unit. At every nonzero index the chosen
row samples \(x\), because \(h\le D\) is in its first block.
The binomial coefficients at both orders are \(p\)-units.

Let \(B_m^-\) be the product of full order-\(m\) ramps for
the decreasing rows counted in \(\mathcal B\), and let \(B_{m'}^+\)
be the corresponding full order-\(m'\) product for the increasing
ones. Uncounted outside factors are at most \(HD\).
In the marked polynomial the constant coefficient has valuation
\(U=v_p(B_m^-)\), and all other coefficients have valuation
at least \(E\). In the other polynomial all ordinates are
nonnegative and the leading ordinate is \(V=v_p(B_{m'}^+)\).
Every slope in that Newton polygon is at most \(V\), by the
last-edge slope formula and the monotonicity of the slopes.

If \(E>U+mV\), then for every slope \(s\) of the other polygon,
index zero is the unique minimizer of the marked ordinate
minus \(s\) times its index. Indeed each nonzero index \(k_1\)
exceeds the zero index after this subtraction by at least
\(E-U-k_1s\ge E-U-mV>0\). There is no common slope.
The ordinary Newton polygon root count then rules out a common
root. Consequently a possible common root requires
\[
                         v_p(x)\le v_p\!\left(B_m^-(B_{m'}^+)^m\right)
\]
for each \(p>HD\) dividing a mark. The bound charges only the
decreasing outside rows at the marked order and the increasing
ones at the other order.

## The product comparison

A prime above \(HD\) divides at most one mark, since the span is
less than \(D\). Thus the large-prime part of the product of
the marks divides \(B_m^-(B_{m'}^+)^m\). For each \(p\le HD\)
dividing a mark, remove one mark of maximal \(p\)-valuation.
As in the [retained-prime product proof](MULTIVARIATE_LARGE_FACTORIAL_MARK_PRIME_BUDGET.md#a-product-bound-at-a-possible-common-root),
the remaining small-prime part divides \((k-1)!\).
For clarity, if \(p\mid m\) then all marks are units. Otherwise
the indices divisible by \(p^\nu\) form a single residue class
modulo \(p^\nu\). After removing a maximal-valuation term at
most \(\lfloor(k-1)/p^\nu\rfloor\) such indices remain.

Since \(\log D\ge144JH/\gamma\), the elementary prime count
\(\pi(HD)\le9HD/\log(HD)\) gives at most \(k/8\) removed
marks. At a common root the logarithm of the remaining product
would therefore exceed \((7/8)k\log M\). The product divisibility
bounds that logarithm above by
\[
             k\log D+J^2\mathcal B
                   \le \left(\frac23+\frac1{16}\right)k\log M .
\]
In this bound \(k\le D\), \(D\ge144\) and \(M>D^2/12\) give
\(D\le M^{2/3}\). Also \(J^2\mathcal B\le\gamma D\log M/(32J)
\le k\log M/16\). The incompatible bounds exclude a common root.

Terminal two follows from even positivity, terminal three from
[quadratic nondivisibility](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#an-inequality-using-the-endpoint-coefficients),
and the absolute cutoff covers all \(N\ge J\). A rational Bézout
identity gives the stated characteristic-zero conclusion.

## A balanced example with a dense gap spectrum

Take \(D\) large, \(M=2D^3\) and \(q=\lfloor3D/4\rfloor\).
Use increasing rows \((M-d,M)\) for every \(1\le d\le D\).
Use their decreasing reflections for \(1\le d\le q\) and add
\(\sum_{d=q+1}^D d\) copies of \((M,M-1)\).
The signed gap is zero. All the rows have common upper endpoint,
and every gap \(1,\ldots,D\) occurs. The chosen row of gap
\(D\) has smaller exponent \(M-D\), whereas all decreasing rows
have smaller exponent at least \(M-3D/4\). Thus \(A=1\), \(B=0\)
and \(\gamma=1/4\) suffice. One may keep the fixed balanced hard
pair \((1,2),(2,1)\) outside the group, or any fixed balanced
ordinary core. Its factors eventually lie below the prime cutoff.
The criterion therefore applies in every degree for all
sufficiently large \(D\), despite a complete gap spectrum and
no positive signed degree difference.
