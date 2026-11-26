# All-degree coprimality from simultaneous midpoint marks

Take an integer \(J\ge4\) at least the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md),
and fix integers \(A\ge1\), \(B\ge0\) and a rational \(\xi>0\).
In an original ordinary two-monomial profile, select \(r>0\)
increasing and \(r\) decreasing unequal rows.
Write \(D\) for their minimum gap, and \(M\) for their maximum
upper exponent. Their gaps \(d_i\) belong to \([D,AD]\), and
their upper exponents \(M_i=M-e_i\) satisfy \(0\le e_i\le BD\).
Suppose an interval \([u_0,u_1)\) of length at least \(\xi D\)
is contained in
\[
                 \bigcap_{i\ {\rm selected}}
            \left[e_i+\left(\frac12-\frac1{4J}\right)d_i,\
                       e_i+\frac12d_i\right).
\]
The rows need not form reflected pairs.

Among the retained rows, call a row passive if its larger
exponent is in \([M-BD,M]\), its gap is at most \(AD\), and either
its larger exponent is at most \(M-u_1\) or its smaller
exponent is at least \(M-u_0\). Put \(H=J(A+B+1)\), and let
\(\mathcal L=\sum_{\rm retained}d_i\log(JM_i)\) where \(M_i\)
is the retained upper exponent. In that sum omit passive rows
and all rows with \(JM_i\le HD\). There are no count or gap
bounds on the omitted rows beyond these requirements.

**Theorem.** For
\[
 \begin{split}
 D&\ge\max\{\lceil2/\xi\rceil,\ 24(A+B),\ 12\cdot2^8,\
                           2^{\lceil512H/\xi\rceil}\},\\
 \mathcal L&\le\frac{\xi rD\log M}{32J^2},
 \end{split}
\]
the original factorial moments are consecutively coprime at
every degree over every characteristic-zero field.
In particular the retained prime budget increases linearly
with the number of selected rows per orientation.
The cutoff \(J\) is unevaluated.

## The same mark at both orders

Fix adjacent orders with \(4\le N<J\); write \(m=2l\) for the
even order and \(m'\) for the odd order.
The [curvature twelve theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
allows us to assume total curvature below twelve, so
\(M>D^2/12\). Take the integers \(u\in[u_0,u_1)\) and mark
the integer \(x_u=M-u\). The literal factor at order \(j\) is
\(jx_u\). Unlike a mark used at only one order, this factor
is allowed to be sampled in both orders.

For row \(i\), put \(t_i=(u-e_i)/d_i\).
The selected rows satisfy
\(\frac12-\frac1{4J}\le t_i<\frac12\).
After stripping the smaller-end factorial,
an increasing row samples \(jx_u\) exactly at indices
\(k\ge j-\lfloor jt_i\rfloor\).
A decreasing row samples it exactly at
\(k\le\lfloor jt_i\rfloor\).
Indeed the factor lies in the full stripped ramp since
\(0\le t_i<1\), and the inequalities check its inclusion in
the sampled prefix or reversed prefix.

At even order \(m=2l\), \(\lfloor mt_i\rfloor=l-1\).
Thus there are \(r\) marks at every index except \(l\), which
has none. At odd order \(m'=2s+1\), \(\lfloor m't_i\rfloor=s\).
Every index of that order has exactly \(r\) marks, using the
equal counts in the two orientations.

## Native prime comparison and the multiplicity gain

Let \(p>HD\) divide \(x_u\), and put \(E=v_p(x_u)\).
For any selected or passive row, the difference between a
stripped factor \(jM_i-h\) and \(jx_u\) is
\[
                      j(u-e_i)-h,\qquad 0\le h<j d_i .
\]
Its magnitude is at most \(HD\): both \(ju\) and \(je_i+h\)
are nonnegative and at most \(J(A+B)D\), using
\(u<(B+A/2)D\). Thus any such factor different from \(jx_u\)
is a \(p\)-unit. A passive row samples no \(jx_u\): its
full ramp lies entirely below the candidate interval of
factors or starts strictly above it. Rows omitted because
\(JM_i\le HD\), and the binomial coefficients, also give units.

Let \(B_j\) be the product of the full stripped ramps at
order \(j\) of the rows counted in \(\mathcal L\).
At the even order the middle ordinate is at most
\(U=v_p(B_m)\), and every other ordinate is at least \(rE\).
At the odd order all ordinates after subtracting \(rE\)
are in \([0,V]\) with \(V=v_p(B_{m'})\).
Every slope of the odd Newton polygon lies in \([-V,V]\).
If \(rE>U+lV\), the middle index \(l\) is the unique
minimizer of the marked even ordinate minus any such slope
times its index. Hence there is no shared slope and no
common root. A possible common root therefore forces
\[
                         r\,v_p(x_u)
                                    \le v_p(B_m B_{m'}^l).
\]

The marks form a step-one progression of count \(k\ge\xi D/2\),
span at most \(D/(4J)\) and count at most \(D\). Each mark
exceeds \(M/2\). A prime above \(HD\) divides at most one mark,
so the \(r\)-th power of the product of all large-prime parts
divides \((B_m B_{m'}^l)_{>HD}\). Remove one maximal-valuation
term for every smaller prime that occurs, as in the
[retained-prime product proof](MULTIVARIATE_LARGE_FACTORIAL_MARK_PRIME_BUDGET.md#a-product-bound-at-a-possible-common-root).
The remaining small-prime part divides \((k-1)!\).
From \(\log D\ge144H/\xi\) we get
\(\pi(HD)\le\xi D/16\le k/8\). Hence the \(r\)-th power
of the product of the remaining marks has logarithm greater
than \((7/8)(15/16)rk\log M\), using \(M>2^{16}\).
Its logarithm is bounded above by
\[
       r k\log D+J^2\mathcal L
                           \le\left(\frac23+\frac1{16}\right)rk\log M .
\]
Here \(D\le M^{2/3}\), and each retained row contributes at
most \((m+l m')d_i\le J^2d_i\) factors of size at most \(JM_i\).
These bounds contradict one another.

Terminal two follows from even positivity, terminal three from
[quadratic nondivisibility](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#an-inequality-using-the-endpoint-coefficients),
and the absolute cutoff covers \(N\ge J\).
Rational Bézout identities transfer the conclusion to every
characteristic-zero field.

## A balanced example with a passive lower band

Take large even \(D\) and \(M=4D^2\). Use one reflected pair of
gap \(2D\) and upper exponent \(M\), and one reflected pair of
gap \(D\) and upper exponent \(M-D/2\).
The midpoint window is \([D-D/(4J),D)\).
At every integer \(s\) in
\([\lceil2D(1-1/J)\rceil,2D)\) one may add \(D\) copies of
each orientation of the unit-gap row with upper exponent \(M-s\).
They are all passive at the midpoint marks, yet their total
gap is of order \(D^2\). A fixed balanced hard core such as
\((1,2),(2,1)\) is permitted, with factors eventually below
the prime cutoff. All profiles are original ordinary profiles
of signed gap zero. The marked upper endpoints have spread
\(D/2\), and the midpoint argument applies independently of
the added lower-band multiplicities.

For comparison, a bottom mark at \(u=s\) in that lower band
samples the uppermost factor of its unit-gap rows. That factor
lies in their last block. The midpoint window instead avoids
all that sampling.

## A general bound from the nominal slopes

Suppose a simultaneous mark gives order-\(j\) coefficient
ordinates \(E C_{j,k}+u_{j,k}\), where \(C_{j,k}\) are the
integer mark counts and the range of the residual ordinates
\(u_{j,k}\) is at most \(U_j\). One may also add an arbitrary
constant in each order and the same affine slope in both orders.
These changes preserve all shared slopes. After discarding the
constants and the common linear term, assume
\(0\le u_{j,k}\le U_j\).

The lower convex envelope of the coefficient ordinates is between
\(E\) times the lower convex envelope of the counts and that
same function plus \(U_j\). This follows by taking the minimum
over convex combinations of points with the prescribed horizontal
coordinate. The increments of either envelope on the integer grid
give its successive slopes, each repeated by its width.
Subtracting the two endpoint values of an increment therefore shows
that each actual slope differs from \(E\) times its nominal
count slope by at most \(U_j\).

Let \(w_j(s)\) be the width assigned to slope \(s\) by the count
polygon in order \(j\), with missing slopes assigned width zero.
If \(E>nN(U_n+U_N)\), different nominal slopes cannot
produce the same actual slope at adjacent orders \(n,N\).
In fact distinct nominal rational slopes have denominators
at most \(n,N\) and so differ by at least \(1/(nN)\).
The degree of every common factor is consequently at most
\[
                              \sum_s\min(w_n(s),w_N(s)).
\]
A bound of one suffices by even-degree parity.
A larger known separation between distinct nominal slopes may be
used in place of \(1/(nN)\).

For instance, at orders six and seven take one group in each
orientation with common gap \(D\) and upper exponents
\(M-3D/10,M\), respectively, with \(D\) divisible by ten.
For \(59D/100\le u\le61D/100\), the increasing and decreasing
threshold floors are \(1,3\) at order six and \(2,4\) at order
seven. Write the orientation counts as \(r_+,r_->0\).
The order-six count polygon has slopes \(-r_-/4,r_+/2\).
At order seven it has a horizontal edge and, if the counts
differ, an edge of slope \((r_+-r_-)/3\) for \(r_+>r_-\),
or \(-(r_--r_+)/5\) for \(r_->r_+\).
The two slope sets are disjoint and separated by at least
\(\min(r_+,r_-)/4\). This gives a native prime test on the
indicated whole mark interval outside a common midpoint window:
a common factor at that pair requires
\(E\min(r_+,r_-)\le4(U_6+U_7)\).

## Unequal counts in the midpoint window

The all-degree theorem above also holds for unequal positive
orientation counts \(r_+,r_-\), with \(r=\min(r_+,r_-)\)
in the retained budget and every other hypothesis unchanged.

At even degree \(m=2l\) the nominal count slopes are
\(-r_-/l,r_+/l\). At odd degree \(m'=2s+1\), the count
is \(r_-\) up to index \(s\) and \(r_+\) thereafter.
Its lower convex envelope has slope zero on the lower side
and, if the two counts differ, a slope of signed magnitude
\(|r_+-r_-|/(s+1)\) on the other side.
Because \(s+1\in\{l,l+1\}\), the odd and even slope sets are
disjoint with distance at least \(r/(l+1)\).
For example, for \(r_+>r_-\) the positive-side distance is
\[
               \frac{r_+}{l}-\frac{r_+-r_-}{s+1}
                              \ge\frac{r_-}{s+1}.
\]
The other orientation is symmetric and the distance to zero
is at least \(r/l\).

At a possible common root the nominal-slope bound therefore
forces \(rE\le(l+1)(U+V)\), with \(U=v_p(B_m)\) and
\(V=v_p(B_{m'})\) as above.
The same product proof uses
\((B_mB_{m'})^{l+1}\) in the large-prime bound.
Each retained row contributes at most
\((l+1)(m+m')d_i\le J^2d_i\) factors to this product,
so the stated budget and constants still apply.
