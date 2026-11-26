# Exact radial separation under ordered integer gap matching

Start with any ordinary mirrored core and add a finite shared
tail with zero total signed gap. For a tail row of absolute
gap \(d\ge1\) and midpoint \(m\), put \(r=m/d\).
Give an increasing row \(d\) positive units and a decreasing
row \(d\) negative units, each unit carrying the same pair
\((r,d)\) as its row. Suppose the positive and negative units
can be matched so that in every matched unit pair,
\[
                         r_+\ge r_-,\qquad d_+\ge d_- ,
\]
and at least one pair has a strict inequality. Define
\[
 \Sigma=\max_{\rm matched\ pairs}
       \left\{ r_-^{-2},\
               \mathbf 1_{r_+>r_-}\frac{d_-}{r_-}\right\}.
\]

**Theorem.** There is a computable absolute \(\Sigma_0>0\)
such that \(\Sigma\le\Sigma_0\) gives physical ordinary
SFC(2) in every degree. There is no bound on the number
of rows, their individual exponents, or the sizes of the
inflation ratios \(d_+/d_-\). The matching is in integer gap
units, and need not pair whole rows or preserve their
primitive exponent ratios.

For same-gap pairs this recovers the
[ordered midpoint comparison](MULTIVARIATE_ORDERED_OPPOSITE_PAIR_RADIAL_CLOSURE.md).
For equally normalized midpoints \(r_+=r_-\), the hypothesis
uses the squared relative gap \(r_-^{-2}\) and allows
arbitrary curvature in an individual inflated row. In
particular, let \(a<b\) be positive integers and \(t\ge2\).
Add one increasing row \((ta,tb)\) and \(t\) copies of the
decreasing row \((b,a)\) to any mirrored core. If
\(2(b-a)/(a+b)\) is sufficiently small absolutely, this
profile is coprime at every consecutive pair of orders,
uniformly in \(t\) and the entire core. This is an exact
analytic row inflation and identification comparison.

## Exact per-unit functions

Fix \(n,N=n+1\). For an increasing row with gap \(d\)
and midpoint \(dr\), its shared exponents are
\(d(r-1/2),d(r+1/2)\). In the formulas below extend
\(r,d\) to real positive values with \(r\ge2,d\ge1\).
Let \(g(x)=\log\Gamma(1+x)\).
The row's contribution to the adjacent log center,
per unit of gap and after removal of \(\log(N/n)\), is
\[
 E(r,d)=\int_{-1/2}^{1/2}
          [g'(Nd(r+y))-g'(nd(r+y))-\log(N/n)]\,dy .
\]
This follows by writing the exact endpoint log
ratio per degree as an integral of \(g'\) across the
row interval. Reversal negates the contribution.

At order \(j\) define the normalized chord per degree
and per unit gap by
\[
 L_{j,r,d}(v)=\frac1{jd}
       \left[g(jd(r+v-1/2))
       -(1-v)g(jd(r-1/2))-v g(jd(r+1/2))\right].
\]
Thus the row's log multiplier at index \(k\) is
\(jd L_{j,r,d}(k/j)\).
Write \(L^o_{j,r,d}(v)\) for its odd part under
\(v\mapsto1-v\). Reversal also negates this odd part.

The kernel \(K(x)=N g''(Nx)-n g''(nx)\) is positive.
The [paired midpoint comparison](MULTIVARIATE_ORDERED_OPPOSITE_PAIR_RADIAL_CLOSURE.md#same-gap-pairs-give-a-center-larger-than-their-skew)
gives from its Laplace integral
\(K(x)\ge1/(108nN x^2)\) for \(x\ge1\).
Differentiation of \(E\) gives
\[
 \partial_r E=\int_{-1/2}^{1/2}dK(d(r+y))\,dy,\qquad
 \partial_d E=\int_{-1/2}^{1/2}(r+y)K(d(r+y))\,dy .
\]
Since \(d(r+y)\ge3/2\) and \(r+y\le5r/4\),
\[
       \partial_r E\ge\frac{1}{200nN d r^2},\qquad
       \partial_d E\ge\frac{1}{200nN d^2 r}.                     \tag{1}
\]
In particular the center per unit is strictly increasing
in both parameters.

## Two exact derivatives of the odd chord

Put \(y=v-1/2\). The odd curvature is
\[
 (L^o_{j,r,d})''(v)
   =\tfrac12 jd\,[g''(jd(r+y))-g''(jd(r-y))].
\]
Differentiating in \(r\), then applying the mean value bound
\(g^{(4)}(x)\le2/x^3\) to the difference of the two \(g'''\)
terms, bounds the result by \((r-1/2)^{-3}\).
The zero endpoint values of the chord give its supremum
bound by integrating with the interval Green kernel of
total mass at most \(1/8\). Thus
\[
                     \|\partial_r(jL^o_{j,r,d})\|_\infty
                                      \le\frac{j}{3r^3}.        \tag{2}
\]

For the scale derivative use the exact trigamma kernel
with \(q(s)=s/(e^s-1)\):
\[
 (L^o_{j,r,d})''(v)
      =\tfrac12\int_0^\infty
         [e^{-(r+y)u}-e^{-(r-y)u}]\,q(u/(jd))\,du .
\]
It follows from the series for \(g''\) by Laplace
integration and the substitution at degree \(jd\).
The bound \(0\le-q'\le1/2\) follows since
\(q\) decreases and \(q(s)+s/2=(s/2)\coth(s/2)\)
increases. The difference of exponentials has modulus
at most \(u e^{-(r-1/2)u}\).
Differentiation in \(d\) therefore bounds the derivative
of the curvature by
\[
      \frac1{4 j d^2}\int_0^\infty u^2 e^{-(r-1/2)u}\,du
                    =\frac1{2 j d^2(r-1/2)^3}.
\]
The same Green kernel bound, and \(r-1/2\ge3r/4\), give
\[
                       \|\partial_d(jL^o_{j,r,d})\|_\infty
                                        \le\frac1{6d^2r^3}.     \tag{3}
\]
All differentiations are on the positive argument domain,
and the displayed bounds justify the integral derivative
there.

## Integration along the ordered matching

For a matched pair of units move from \((r_-,d_-)\) to
\((r_+,d_-)\), and then to \((r_+,d_+)\).
On the first segment the ratio of the skew derivative
bound (2) to the center derivative bound (1) is at most
\((200/3)j nN\,d_-/r_-\).
This segment is absent if \(r_+=r_-\). On the second
segment the ratio from (3) and (1) is at most
\((200/6)nN/r_+^2\le(200/6)nN/r_-^2\).
Both center increments are nonnegative. Consequently
\[
 \big\|jL^o_{j,r_+,d_+}-jL^o_{j,r_-,d_-}\big\|_\infty
      \le70J^3\Sigma
                       [E(r_+,d_+)-E(r_-,d_-)]                 \tag{4}
\]
for \(j=n,N<J\), where \(J\) is the absolute eventual
angular cutoff from the ordered midpoint comparison.

Every unit is used once. The signed \(\log(N/n)\)
terms cancel by zero total signed gap, and the paired
core has zero center. Therefore the exact adjacent
log center \(C\) is the sum of the rightmost bracket
in (4), and \(C>0\). Summing (4) gives
\(\delta_j\le70J^3\Sigma |C|\) for the exact tail skew
norm of that comparison. Also the maximum tail relative
gap, measured against the midpoint, is at most
\(\sqrt\Sigma\), because every \(r\) is at least the
smaller \(r_-\) in its matched pair.

Take \(B,\eta_1\) from the
[uniform exact skew test](MULTIVARIATE_ORDERED_OPPOSITE_PAIR_RADIAL_CLOSURE.md#an-exact-skew-test-for-any-diffuse-tail-of-a-mirrored-core)
and choose
\[
              \Sigma_0=\min\left(\eta_1^2,
                                         \frac1{560 B J^3}\right).
\]
This ensures \(r\ge2\) and all hypotheses of that test
at a possible collision. It would force
\[
                  |C|\le2B(\delta_n+\delta_N)
                            \le280B J^3\Sigma |C|\le |C|/2 .
\]
This contradicts \(C>0\). The absolute curvature
and eventual angular theorems handle the orders or
curvatures outside the skew test, and terminal two
is separated by even real positivity. Rational Bézout
gives the result over every characteristic-zero field.

For an example where both parameters move, take the
increasing row \((2A+2,2A+4)\) and two decreasing rows
\((A+1,A)\), with any positive multiplicity of this block.
The increasing units have \(d_+=2,r_+=A+3/2\);
the decreasing units have \(d_-=1,r_-=A+1/2\).
For sufficiently large \(A\) the conclusion holds
uniformly in the block multiplicity and the mirrored core.
The method needs a one-sided matching of the gap masses.
It gives no exclusion if the center contributions cancel.

## Inflation compared by spread at a fixed normalized midpoint

There is also an all-degree test allowing both increases and
decreases of scale within a fixed normalized midpoint.
Group the tail by its values of \(r\). At each \(r\) form the
signed measure on gaps
\[
           \nu_r=\sum_{i:r_i=r}\operatorname{sgn}(d_i)
                                      |d_i|\,\delta_{|d_i|},
\]
where \(\delta_x\) is a point mass. Suppose
\(\int 1\,d\nu_r=\int d\,d\nu_r(d)=0\), and
\[
                       A_r(t)=\int (d-t)_+\,d\nu_r(d)\ge0
                                  \qquad\text{for all }t\ge1.
\]
Assume at least one \(A_r\) is not identically zero.
The same inequalities with every sign reversed are also
allowed. If the tail relative gap \(\eta=\max_r 1/r\) is
sufficiently small absolutely, the enlarged mirrored core
has ordinary SFC(2) in every degree.
The support of \(A_r\) is bounded, and its piecewise linear
form makes the hypotheses a finite rational test.
They say that the positive gap-unit measure is more spread
than the negative one in the usual convex order.

For the proof, the Laplace lower bound on \(-K'\) is
\[
                         -K'(x)\ge \frac1{162 nN x^3}
                                  \qquad(x\ge1).
\]
Indeed restrict its integral to \(0<u<1/x\) as in (1),
using the extra factor \(u\) from differentiation.
It follows that
\[
                         -\partial_d^2 E(r,d)
                                      \ge\frac1{250nN d^3r}.    \tag{5}
\]
Here the integrand for \(-\partial_d^2 E\) contains
\((r+y)^2[-K'(d(r+y))]\).

The odd log chord has the upper bound
\[
                   \|\partial_d^2(jL^o_{j,r,d})\|_\infty
                                      \le\frac6{d^3 r^3}
                             \qquad(r\ge2,\ d\ge1).             \tag{6}
\]
For details, in the Laplace kernel of (3) the second
derivative of \(q(u/(jd))\) has absolute value at most
\(u/(j d^3)+18u^2/(j^2d^4)\).
The first term uses \(|q'|\le1/2\) and the second uses
\(|q''|\le18\).
In fact \(q''(t)=e^t[(t-2)e^t+t+2]/(e^t-1)^3\).
The inner bracket and its derivative vanish at zero
and its second derivative is \(t e^t\). Thus on \(0<t\le1\)
the formula gives \(0\le q''(t)\le e^2/6\).
On \(t\ge1\), divide numerator and denominator by \(e^{3t}\)
and bound by
\((t+2)(e^{-t}+e^{-2t})/(1-e^{-t})^3<18\).
Use the same exponential difference as in (3).
With \(a=r-1/2\), integrating and taking the chord
supremum gives
\[
      \frac1{8d^3a^3}+\frac{27}{4j d^4a^4}
                    \le\frac6{d^3 r^3},
\]
since \(j\ge2,\ a\ge3/2\), and \(a\ge3r/4\).
This proves (6).

For a twice differentiable function \(F\) on the gap range,
two integrations give
\[
                     \int F(d)\,d\nu_r(d)
                            =\int_1^\infty A_r(t)F''(t)\,dt,
\]
because both affine terms vanish under \(\nu_r\).
By (5) and \(A_r\ge0\), each group contributes a nonpositive
exact center and at least one contributes a strictly
negative center. Integrating (6) against the same
\(A_r\) bounds the skew by \(1500 nN/r^2\) times the
magnitude of that center. Summing the groups gives
\(\delta_j\le1500J^2\eta^2 |C|\).
With the same \(B,\eta_1\) as above, it suffices to require
\(\eta\le\min(\eta_1,(12000 B J^2)^{-1/2})\).
The radial test (1) of the ordered midpoint comparison
then contradicts \(C\ne0\) at all remaining orders.
The other orders are covered as before.

For example, at primitive ratio \((a,a+1)\), take three
increasing copies at scale one and one at scale three,
and three decreasing copies at scale two. The positive
gap-unit measure is \(3\delta_1+3\delta_3\), and the
negative one is \(6\delta_2\). Its stop-loss difference
\(A_r\) is nonnegative and nonzero by convexity.
For sufficiently large \(a\), this block, with arbitrary
multiplicity, gives an all-degree ordinary extension of
any mirrored core. The scale lists cross, so a
coordinatewise ordering of the gap units cannot give
this comparison.

## Any fixed higher order of integer scale matching

Fix an integer \(p\ge3\). Keep the same gap-mass grouping
by \(r\), and suppose that at every \(r\)
\[
  \int d^k\,d\nu_r(d)=0\quad(0\le k<p),\qquad
  A_{r,p}(t)=\int\frac{(d-t)_+^{p-1}}{(p-1)!}\,d\nu_r(d)\ge0
                                              \quad(t\ge1).
\]
Assume at least one of these functions is nonzero.
Reversing all the inequalities is also allowed. This gives
ordinary SFC(2) in every degree when the tail relative
gap satisfies the effective absolute bound
\[
             \eta\le\min\left(
                     \eta_1,\ [4000\cdot2^p(p+1)!\,B J^2]^{-1/2}
                                      \right).                  \tag{7}
\]
The bound depends on the chosen order \(p\), not on the
number or sizes of the ordinary rows.

For \(x\ge1\), differentiating the Laplace formula for
\(K(x)\) gives a strict alternating sign in every order.
Using the same interval \(0<u<1/x\) as in (1) gives
\[
        (-1)^{p-1} K^{(p-1)}(x)
                   \ge \frac1{54(p+1)nN x^{p+1}} .
\]
Thus
\[
     (-1)^{p-1}\partial_d^p E(r,d)
                        \ge\frac1{100(p+1)nN d^{p+1}r}.          \tag{8}
\]
Indeed the \(p\)-th derivative integrates
\((r+y)^p K^{(p-1)}(d(r+y))\), and \(r+y\le5r/4\).

The odd chord derivative satisfies
\[
             \|\partial_d^p(jL^o_{j,r,d})\|_\infty
                           \le \frac{5\cdot2^p p!}{d^{p+1}r^3}.
                                                                    \tag{9}
\]
To prove it, extend the scale variable to \(z\) on
\(|z-d|\le d/2\). For \(u>0\), the point
\(w=u/(jz)\) has \(\Re w\ge |w|/3>0\), so
\(q(w)=w/(e^w-1)\) is holomorphic there. In this sector
\(|q(w)-1|\le8|w|\). For \(|w|\ge1/2\) use
\(|e^w-1|\ge e^{\Re w}-1\ge\Re w\), which gives
\(|q(w)|\le3\). For \(|w|\le1/2\) put
\(h=(e^w-1)/w-1\). The exponential remainder gives
\(|h|\le |w|e^{1/2}/2<1/2\). Then
\(|q(w)-1|=|h/(1+h)|\le2|h|\le2|w|\).
Since \(|w|\le2u/(jd)\) on the scale disk, Cauchy's
estimate on \(q(u/(jz))-1\) gives
\[
             |\partial_d^p q(u/(jd))|
                          \le \frac{16\cdot2^p p!\,u}{j d^{p+1}}.
\]
Apply this inside the odd kernel integral, with the
same exponential difference \(u e^{-(r-1/2)u}\)
and Green kernel as in (3). It gives
\(2\cdot2^p p!/[d^{p+1}(r-1/2)^3]\), bounded by (9).

The function \(A_{r,p}\) vanishes above the maximum gap.
It also vanishes below the smallest gap because the
truncated polynomial there reduces to a polynomial of
degree \(p-1\) in \(d\). Taylor's formula from \(d=1\)
and the moment equalities give
\[
                     \int F(d)\,d\nu_r(d)
                            =\int_1^\infty A_{r,p}(t)F^{(p)}(t)\,dt .
\]
It follows from (8) that all nonzero group centers
have sign \((-1)^{p-1}\). There is a nonzero center,
since a nonzero continuous nonnegative \(A_{r,p}\)
is positive on an interval. Comparing (9) with (8)
and summing the groups gives
\(\delta_j\le500\cdot2^p(p+1)! J^2\eta^2 |C|\).
The exact skew test gives a contradiction under (7),
as in the lower-order comparisons.

For explicit ordinary blocks, fix positive integers
\(D,p\). Take the gap masses
\((-1)^{p-k}\binom pk\) at gaps \(D+k\),
\(0\le k\le p\), multiplied by a common multiple
of \(D,\ldots,D+p\). Dividing each absolute gap
mass by its gap gives the integer number of rows
of the indicated orientation. Use a common
normalized midpoint \(r=a+1/2\) with \(a\) a
positive integer, so the row at gap \(d\) is
\((da,d(a+1))\) or its reverse. The moment
equalities are the usual forward differences of
polynomials of degree less than \(p\).
For each \(t\), the \(p\)-th forward difference
of \((d-t)_+^{p-1}/(p-1)!\) is nonnegative.
Indeed the \((p-1)\)-st forward difference is
the integral of \(\mathbf 1_{D+u_1+\cdots+u_{p-1}>t}\)
over the unit cube. Shifting \(D\) by one and
subtracting gives the nonnegative \(p\)-th difference,
which is positive for \(D<t<D+p\).
The same identity follows by smoothing the truncated
power before applying the fundamental theorem of calculus.
Taking \(a\) sufficiently large yields exact all-degree
coprimality for every multiplicity of this block and
any mirrored core, despite cancellation of all the
specified gap-mass moments.
