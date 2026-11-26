# Radial and angular separation by matching linked Gauss factors

Fix an ordinary two-monomial profile with distinct exponent
columns. Delete equal rows. Each remaining row has smaller
exponent \(\ell\), positive gap \(d\), and one orientation.
Represent it by the \(d\) Gauss slots
\[
                        (u,x)=(\ell/d,h/d),\qquad 1\le h\le d
\]
of that orientation. A slot is shared when \(u>0\) and
private when \(u=0\). Its base at order \(m\) is \(mu+x\).
The slope \(u\) and offset \(x\) are determined jointly by the
[matched tanh data](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid).

At any fixed matching of opposite slots, the estimates below
use the same matching across the entire consecutive pair.
They give all-degree certificates for original mixed profiles,
including near matches between different ratios and profiles
with positive and negative surplus in different ratio groups.
Let \(J\ge4\) be at least the
[computable absolute eventual angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Its numerical value is unevaluated.

## A radial criterion with arbitrary slopes and private rows

A matching may pair any opposite shared slots and any opposite
private slots. All slots not used in pairs are unmatched.
Orient the columns so the signed gap sum \(\Delta_{\rm sh}\)
of the shared rows is positive. Define
\[
\begin{split}
 C_0={}&\sum_{\rm paired\ shared}|\log(u/v)|
        +\sum_{\rm unmatched\ shared} u^{-1}\\
      &+\sum_{\rm paired\ private}|\log(x/y)|
        +\sum_{\rm unmatched\ private}\log(J/x),\\
 C_1={}&\sum_{\rm paired\ shared}|x/u-y/v|,
\end{split}
\]
where the paired slots are denoted \((u,x),(v,y)\).
If a matching gives
\[
                   \Delta_{\rm sh}\ge (J-1)C_0+\frac32 C_1,    \tag{1}
\]
the profile satisfies SFC(2) in every degree. More precisely,
\(Q_n(z)\) and \(Q_{n+1}(\lambda z)\) have no common root for
\(\lambda\ge1\) in the indicated orientation. No curvature
restriction is required for this radial criterion.

Use the
[paired Gamma annuli](MULTIVARIATE_PAIRED_GAMMA_ROOT_ANNULUS.md#matched-factor-bounds-for-the-ordinary-profiles)
for the Gauss factors. Their variable differs from the physical
one by a positive common scale independent of the order, which
cancels in comparisons of consecutive root moduli.
A matched pair of bases \(A_m,B_m\) contributes the product
annulus with endpoints \(\min(1,B_m/A_m),\max(1,B_m/A_m)\).
We use this interval without the optional tighter intersection
in that theorem. For a shared pair \(A_m=mu+x\), \(B_m=mv+y\),
put
\[
 L_m=\log(A_m/B_m)=\log(u/v)+\delta_m,\qquad
 \delta_m=\log(1+x/(mu))-\log(1+y/(mv)).
\]
The sign of \(\delta_m\) is fixed and its modulus decreases
with \(m\), by differentiation in \(1/m\).
The loss in comparing the lower radial endpoint at \(n\)
with the upper radial endpoint at \(N=n+1\) is
\[
 (L_n)_+ +(-L_N)_+
       \le |L_N|+|\delta_n-\delta_N|
       \le |\log(u/v)|+\frac{|x/u-y/v|}{n}.
\]
The last bound follows from the monotonicity of
\(|\delta_m|\) and the derivative bound for \(\log(1+s)\)
at \(s\ge0\). A private pair has bases \(x,y\) at both
orders, so its loss is at most \(|\log(x/y)|\).

An unmatched increasing shared slot contributes a radial
interval \([1/(mu+x+m-1),1/(mu+x)]\).
In the comparison it yields at least
\(\log(N/n)-1/u\), since
\(\log(1+(x+n-1)/(nu))\le1/u\).
The corresponding decreasing slot contributes at least
\(-\log(N/n)-1/u\) by the reciprocal interval.
An unmatched private slot of either orientation loses
at most \(\log(N/x)\le\log(J/x)\) for \(N<J\).
Adding the terms gives the lower bound
\[
          \Delta_{\rm sh}\log(N/n)-C_0-C_1/n
\]
for the log of the lower endpoint at order \(n\) divided
by the upper endpoint at \(N\). For \(2\le n<N<J\), it is
strictly positive by (1), since
\(\log(N/n)>1/N\) and \(N/n\le3/2\).
The roots at the upper order thus have smaller moduli in the
same physical variable, which also excludes any common root
under the stated directional rescaling. The first and second
moments already separate by strict real positivity of the
quadratic moment. At \(N\ge J\) use the absolute angular theorem.

## Offset matching within one slope

For a shared group of common slope \(u\), let \(r_u\) be its
number of rows and \(D_u\) its signed sum of gaps. Pair slots
within each shared group as follows. If the two orientations
have different slot counts, add to the deficient side one
dummy grid of size \(|D_u|\), with offsets
\(h/|D_u|\) for \(1\le h\le |D_u|\). Sort and pair the
equal-length lists and discard all dummy pairs. The actual
cost of the retained pairs is bounded by
\[
                     \sum |x-y|
                      \le \frac{r_u+\mathbf1_{D_u\ne0}}2 .      \tag{2}
\]
In fact, for two equally sized sorted lists, the sum of
the absolute pair distances is the integral of the absolute
difference of their counting functions on \([0,1]\).
This follows by counting the pairs whose intervals cross each
point. Sorted matching gives intervals through that point
only in one direction. For a grid of size \(d\), the
counting function is \(\lfloor d v\rfloor\). In the
balanced pair of lists all terms \(dv\) cancel, and
the integral of each fractional part \(\{dv\}\) is \(1/2\).
The triangle inequality gives (2).

For a fully shared profile this supplies the concrete
sufficient condition
\[
 |\Delta_{\rm sh}|\ge (J-1)\sum_u |D_u|/u+
                \frac34\sum_u(r_u+\mathbf1_{D_u\ne0})/u ,
              \qquad \Delta_{\rm sh}\ne0.
\]
In the single-slope case there is a sharper criterion.
After orientation, all \(\Delta_{\rm sh}\) unmatched slots
increase. At \(u\ge n-1\) the margin per slot is at least
\[
          \log\frac{Nu+x}{nu+x+n-1}
                       \ge\frac{u-n+1}{Nu+1}.
\]
The sum of paired losses is at most \((r_u+1)/(2nu)\).
It follows that all degrees are covered whenever
\[
                     u\ge J-3+(r_u+1)/|\Delta_{\rm sh}|.
\]
Indeed in a remaining pair \(n\le J-2\) the unmatched
margin is at least \((r_u+1)/(Nu+1)\), strictly larger
than the paired loss since \((n-1)u>1\).
This criterion does not bound the total magnitudes of the
opposed gaps.

## Angular stability near a mirrored profile

Suppose the private rows, if present, are already mirrored.
Choose one orientation of the shared rows and form an ordinary
mirrored reference by copying those rows with their reversals
and keeping the private pairs. Denote the reference curvature
by \(S_{\rm ref}\). Match any shared slot in the copied
orientation to an opposite shared slot of the original
profile, leaving any excess slots in either list unmatched.
For positive \(u,v\), let
\[
   \varepsilon=\frac18
      \sum_{\rm paired\ shared}
          \frac{(J-1)|u-v|+|x-y|}{\min(u,v)^2}
               +\frac{J-1}{8}\sum_{\rm unmatched\ shared}u^{-1}.
\]
If \(S_{\rm ref}\ge12\) and
\[
                4(1+2^J)\varepsilon^{1/J}
                              <\frac1{2000J^2},                 \tag{3}
\]
the original consecutive polynomials have disjoint root
arguments in every degree. The groups need not be balanced,
and the signed total degree difference may be zero.

For the coefficient estimate, let \(A>0\) be the base of a
rising slot and put
\[
        L_{m,k}(A)=\log(A)_k-\frac{k}{m}\log(A)_m.
\]
For \(0\le k\le m\), its negative and its derivative in \(A\)
are sums over \(0\le i<k\le j<m\), divided by \(m\), of
\(\log(A+j)-\log(A+i)\) and
\((A+i)^{-1}-(A+j)^{-1}\), respectively.
They lie between zero and \(m^2/(8A)\) and \(m^2/(8A^2)\),
respectively. Indeed
\(\sum_{i<k\le j}(j-i)=m k(m-k)/2\).
A decreasing slot uses \(m-k\) instead of \(k\).
Replacing a base \(mu+x\) by \(mv+y\) therefore costs
at most
\[
                    \frac{m|u-v|+|x-y|}{8\min(u,v)^2}
\]
in the centered log coefficient. Adding or dropping a slot
costs at most \(m/(8u)\). These estimates apply when
replacing the originally unmatched orientation by the mirror
copy. The other orientation and the private factors are exact.
Thus the log coefficient error is at most \(\varepsilon\)
for \(m<J\), excluding the common binomial factor.

The mirrored reference has a uniform angular reserve from
the
[curvature twelve](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_TWELVE.md),
[curvature thirteen](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_THIRTEEN.md),
[curvature fourteen](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_FOURTEEN.md),
and [curvature sixteen](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md)
phase bounds. Here are explicit bounds for use in (3).
With half phase \(\phi_m\), put
\(\psi_m(\theta)=m\theta/2-\phi_m(\theta)\).
In those proofs \(0<\psi_m'\le m\), and the estimates give
\[
             |\phi_N-\phi_n|
                     \le (1/2-1/2000)\theta,\qquad
                                  0\le\theta\le\pi .           \tag{4}
\]
For the curvature-twelve short halves one may replace the
first coefficient cap \(1/10\) by \(249/2500\) in its
scalar endpoint tests. With the tail ratio \(1/5\), the
individual derivative is then at most \(2121/8504<999/4000\).
The short caps at thirteen use \(39/500,3/14\), and at
fourteen use \(1/14,1/6\). They give still smaller
individual derivative bounds, as does the private factor
cap \(10^{-k}\) valid for mirrored private rows at \(S\ge12\).
At sixteen the short and medium bounds are at most
\(7/30\) and \(201/817\) for each derivative. For the other
medium and long bands the paired derivative bounds from
the cited proofs are below \(999/2000\).
In particular the curvature-twelve long bound is at most
\(5[(67/5)/400+2(7/250)]<999/2000\).
Integration from zero gives (4).

At a lower positive root angle \(\theta\), the upper phase
lies between the same and next cosine-zero levels, at
phase distance at least \(\theta/2000\) from both by (4).
Also \(\theta\ge\pi/(2n)\), since the first lower level
is \(\pi/2\) and \(\psi_n'\le n\).
The upper phase derivative is at most \(N\), so every upper
angle in \([0,\pi]\) is separated from that lower angle
by at least \(\pi/(4000nN)\). Reflection gives the same
bound across the cut. In fact opposite-half angles near
zero are each bounded away from zero by the first phase
level. Near \(\pi\), one of the degrees is even. For the
even degree, the last root is at least \(\pi/(2N)\) from
\(\pi\), because its phase at \(\pi\) is an integral multiple
of \(\pi\). The minimum chord gap of reference arguments
is therefore at least \(1/(2000nN)\).

All endpoint-normalized log coefficients apart from the
common binomial contribution are nonpositive. Exponentiating
the error bound therefore gives coefficient errors at most
\(\binom mk\varepsilon\) in degree \(m\), since the exponential
is \(1\)-Lipschitz on the negative half-line.
A root of the actual monic polynomial has modulus in
\([2^{-m},2^m]\), by the coefficient Cauchy bound and reversal.
At that root \(z\), the reference polynomial has modulus at
most \(\varepsilon(1+|z|)^m\), so some reference root \(z_0\)
is within \(\varepsilon^{1/m}(1+|z|)\). The chord distance
of their unit projections is at most \(2|z-z_0|/|z|\).
If the two actual degrees \(n,N<J\) share an argument, the
corresponding reference arguments are thus within chord
distance \(4(1+2^J)\varepsilon^{1/J}\).
This contradicts the reference lower gap under (3).
For \(N\ge J\) the absolute angular theorem again applies
to the original profile.

As an example with no exact mirror rows or per-ratio balance,
take \(u\) a positive integer and \(d=6(u+1)\). Use an
increasing row of gap \(d\) and a decreasing row of gap
\(d+1\) at ratio \(u\), and an increasing row of gap one
at ratio \(u+1\). The total signed gap is zero. The mirrored
copy of the increasing orientation has curvature
\(12+2/(u+2)\). Match the \(d,d+1\) grids by (2), leaving
one unmatched slot at each ratio. The resulting
\(\varepsilon\) tends to zero as \(u\) grows, so the criterion
applies for all sufficiently large integer \(u\).

## Additional exact scalar checks

~~~python
from fractions import Fraction as F
from math import factorial

def lo(x):
    x = F(x)
    return sum(x**k/factorial(k) for k in range(65))

def hi(x):
    x = F(x)
    return lo(x)+x**65/factorial(65)/(1-x/66)

cap = F(249,2500)
assert 2/lo(3) < cap
for V in (F(12), F(221,2)):
    assert 28*(13/(1+V))**2*hi((V/4+F(13,2))/28)/lo(6) < cap
assert max(hi(2)/1024,F(7,100)*hi(F(6,25))) < cap

for c,r in ((cap,F(1,5)),(F(39,500),F(3,14)),
            (F(1,14),F(1,6)),(F(1,7),F(1,7)),(F(1,10),F(1,10))):
    assert (c+r/(1-r)**2-r)/(1-c-r*r/(1-r)) < F(999,4000)
assert F(201,817) < F(999,4000)
assert 5*(F(67,5)/400+2*F(7,250)) < F(999,2000)
~~~
