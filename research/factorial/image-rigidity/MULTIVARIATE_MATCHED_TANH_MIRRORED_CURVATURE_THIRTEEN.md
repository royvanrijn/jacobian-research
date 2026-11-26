# All-degree angular interlacing for mirrored curvature at least thirteen

Let the unequal rows of an ordinary two-monomial profile occur in
mirrored pairs. With an optional common physical Gaussian
addition \(\tau\ge0\), let \(S\) be the absolute curvature
including \(\tau\). If \(S\ge13\), the consecutive normalized
moments strictly interlace on the unit circle in the ordering
cut at \(1\), in every degree. At \(\tau=0\), this proves
ordinary SFC(2) over every characteristic-zero field.

This strengthens the
[matched tanh phase comparison at curvature fourteen](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_FOURTEEN.md).
Two pairwise estimates give the extension. In the medium
range, the linked coefficient derivative retains the
endpoint mass in the factor multiplying the curvature.
In the long range, the largest relative gap costs at least
its square in the paired Gamma curvature, giving a smaller
connected-tree majorant.

## The same linked endpoint mass

For a mirrored pair let \(v_i=d/(2\ell+d)\) in each of its
\(2d\) matched Gauss slots, as obtained from the linked
\((T_i,T_i')\) by
\(v_i=(n/T_i'-(n-1)/T_i)^{-1}\).
The curvature is exactly \(S=\tau+2\sum_i v_i/(1+v_i)\).

If there is a private pair, the
[private-start geometric estimate](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md#coefficient-halves-and-endpoint-interpolation)
gives \(c_{m,k}\le e^{-(S-5)k/4}\le e^{-2k}\)
up to the midpoint in either half. Here the Gaussian factor
contributes at most \(e^{-\tau k/4}\) in this range.
Since \(e^{-2}<1/7\), each half on the unit disk has
logarithmic derivative of modulus below \(1/4\), by
summing the geometric series. The phase comparison below
then gives the conclusion. We may thus assume all pairs
are shared, and \(13\le S<14\) by the existing theorem.

Use \(V,W=1+V,L_0,L_*=V/2,Q,q,D=Nq\) as in
the [curvature fourteen endpoint estimates](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_FOURTEEN.md#exact-matched-tanh-endpoint-mass).
In particular
\[
  -\log q=\tau/2+\sum_i\operatorname{artanh}v_i,\qquad
  q\le e^{-S/2},\quad Wq\le(1+S)e^{-S/2},\quad L_0q\le e^{-S/2}.
                                                               \tag{1}
\]
The \(Wq\) proof at \(C=1+S\) still applies: at \(C=14\)
its comparison polynomial is \(r(r-2)^2+8>0\).
The \(L_0q\) proof there depends only on integer shared
mirroring, so it is unchanged. We have \(V<126\) by
\(V\le2S+S^2/2\).
The coefficient bounds, including \(Q(1/2)\le e^{-S/4}\)
and \(Q(w)\le qe^{L_*w}\), are the same as there.

Write \(F_m\) for the full normalized polynomial and \(H_m\)
for its symmetric half. Whenever \(H_m\) has an analytic
logarithm on the unit disk, normalized to zero at zero,
put \(\phi_m(\theta)=\Im\log H_m(e^{i\theta})\).
The centered polynomial is
\(2|H_m(e^{i\theta})|\cos(m\theta/2-\phi_m(\theta))\).
If the two centered phases increase and
\(|\phi_N-\phi_n|<\theta/2\) on \(0<\theta\le\pi\),
each upper phase at a lower root lies between the same
and the next root levels. The two cut endpoints and
reflection give strict interlacing, as in the earlier proof.

## Short and medium halves

The case \(m=1\) has \(H_m=1\). For \(N\le44\) use
the first-coefficient endpoint checks at \(m=2,44\),
and the \(k\ge2\) checks at \(m=4,44\) and \(w=1/2\).
At \(m=44\) for \(j=1,2\), use
\[
   q e^{jV/88}\le e^{-13/2}\frac{14}{1+V}e^{jV/88}.
\]
The right side is maximal at \(V=13\) or \(126\) on
this interval. This gives \(c_{m,1}<39/500\) and
\(c_{m,k}<(3/14)^k\) for \(k\ge2\).
The latter uses \(2eQ(1/2)\le 2e^{-9/4}<3/14\).
The two geometric sums of coefficients and derivatives
give \(|zH_m'/H_m|<1/4\) on the unit circle.

For \(N>44\), \(n\ge44\) and
\[
 L_*/n\le63/44,\qquad L_*D/n<27/2500.
\]
At \(D\le1/20\), each bound on \(D e^{jL_*/n}\), \(j=1,2\),
has its maximum at the switch \(D=(27/2500)/(63/44)\)
or at \(1/20\). This gives \(c_{m,1}<63/1000\)
and the same cap \((3/14)^k\) for \(k\ge2\) by the
\(k=2\) and midpoint endpoint tests. The half derivative
is again below \(1/4\) on the unit circle.

Next take \(1/20<D\le2/5\). Put \(a=27/2500\), \(D_*=2/5\).
For \(1\le k\le8\), the maximum of \(D e^{ak/D}\) is at
\(D=D_*\) on this interval. For \(2\le k\le7\), the maximum
of \(D^{k-1}e^{a k^2/D}\) is also there: both maxima
are at endpoints, and the scalar checks compare them.
Thus the first seven coefficients on radius \(6/5\)
are bounded by
\[
                  b_k=\frac{((6/5)D_*)^k e^{ak^2/D_*}}{k!}.
\]
The \(k=8\) and midpoint checks bound the remaining
coefficients by \((3/14)^k\). The sums give
\(\|H_m-1\|_{6/5}<13/20\),
\(\sum_{k=1}^7 k b_k<21/25\) and
\(\sum_{k=1}^7 k(k-1)b_k<9/20\).

For the adjacent difference use the same exact first-factor
and Gamma derivatives as in the curvature-fourteen medium
range. Indeed \(2L_0 7/n<1/2\), by (1), \(D>1/20\), and
\(N/n\le45/44\). Also \(V/n<1\), since \(WD/n<27/1250\)
and \(D>1/20\). Thus \(|(\log B_m)'|\le1/n\) and
\[
 |\partial_m c_{m,k}|\le \frac{c_{m,k}}n
                    \left(k+\frac{8W}{n}k(k-1)\right),\qquad k\le7.
\]
In the second term, keep the dependence between \(W/n\)
and the coefficient mass: use \(W/n\le(27/1250)/D\)
and maximize \(D^{k-1}e^{ak^2/D}\) at \(D_*\).
Consequently the low-coefficient difference is at most
\[
     \frac{21/25+8(27/1250)(9/20)/D_*}{44}.
\]
After adding the two geometric tails, the difference
\(\|H_N-H_n\|_{6/5}\) is below \(3/125\).
The straight line between the two halves stays within
distance \(13/20\) of one, so the difference of the
logarithms is below \((20/7)(3/125)<7/100\).
Cauchy's estimate therefore bounds the paired angular
derivative by \(7/20<1/2\). The individual angular
derivatives are below \(11/2\), by the log norm
\(\log(20/7)<11/10\). This gives the phase inequalities.

## Long halves

Suppose \(D>2/5\). Then \(n\ge266\) and
\[
  D/n<151/100000,\qquad WD/n<2113/100000,\qquad
  L_0D/n<151/100000 .
\]
Put \(K=\lceil14D\rceil+7\le34D\).
The checks at \(k=K\) and the midpoint give
\(c_{m,k}\le3^{-k}\) from \(K\) to the midpoint.
For the first endpoint use
\((e/14)e^{(34/2)(2113/100000)}<1/3\).
The exact mirrored first-factor derivative gives
\(|(\log B_m)'|\le1/n\) and \(B_m<(103/100)D\).

Use the [sharper Gamma estimate](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md#the-long-range-with-sharper-gamma-derivative-estimates)
with \(\kappa=13/250\) and \(\alpha=L_0/[(1-\kappa)n]\).
Both \(L_0K/n,K/n\) are less than \(\kappa\), and
\(K+1<n/2\).
It gives \(H_2<(28/25)W/n<3/50\) and \(\alpha K<11/200\).
With \(\beta=e^{H_2+1/2}H_2\) the bounds are
\(\beta K<71/50\), \(\beta B_m<43/1000\).
For shared mirrored rows \(L_0\le\sqrt{V/2}\):
if \(L_0>1\), the pair attaining it contributes
\(2\ell L_0^2\ge2 L_0^2\) to \(V\), and if \(L_0=1\)
use \(V\ge13\).
Thus the tree parameter satisfies
\[
 \frac{\alpha}{\beta}\le e^{-1/2}\frac{L_0}{W}
      \le e^{-1/2}\frac{\sqrt{13/2}}{14}<1/9 .
\]
Here \(\sqrt{V/2}/(1+V)\) decreases for \(V>1\).
This replaces the tree base in the
[differentiated connected comparison](MULTIVARIATE_ABSOLUTE_CURVATURE_20.md#a-smaller-differentiated-connected-logarithm-bound)
by \((10/9)^{10}<29/10\).
With the same additional-edge bound it gives \(p_m=\Xi B_m<7/25\)
and \(|\ell_k'(m)|\le(5/2)B_m p_m^{k-1}/n\),
using the first-factor derivative bound.

For the truncated logarithm \(A_m\), the coefficient norms
on \(9/5,6/5\) are below \((13/5)B_m,(3/2)B_m\),
respectively, by summing \(B_m p_m^{k-1}/k\).
The paired difference on \(6/5\) is below \(5D/n\).
The half and \(e^{A_m}\) agree through degree \(K\).
On \(6/5\) the relative exponential and actual tails
are bounded by
\[
 (2/3)^{K+1} e^{(41/10)B_m}<e^{-15/4},
 \qquad (2/5)^K e^{(3/2)B_m}<e^{-4K/5}.
\]
The first uses \(K+1\ge14D+8,\ D>2/5\).
The second uses \(K\ge14D\) and
\(e^{4/5+(3/2)(103/100)/14}<5/2\).
Since \(K\ge13\), the two tails sum to less than \(3/125\).
This gives \(\|\log H_m-A_m\|_{6/5}<1/40\).
The paired log difference is at most \(5D/n+2/40<3/50\);
the individual norm is below \((3/2)B_m+1/40\).
Cauchy's estimates give paired angular derivative below
\(3/10<1/2\) and individual derivatives below \(m/2\).
The phase comparison therefore proves interlacing in all orders.

## Exact scalar checks

~~~python
from fractions import Fraction as F
from math import factorial

def lo(x):
    x=F(x)
    return sum(x**j/factorial(j) for j in range(81))

def hi(x):
    x=F(x)
    return lo(x)+x**81/factorial(81)/(1-x/82)

e=hi(1)
U=1/lo(F(13,2)); r=F(3,14)
assert 2/lo(F(13,4))<F(39,500)
assert 2/lo(F(9,4))<r
for V in [13,126]:
    q=U*F(14,1+V)
    assert 44*q*hi(F(V,88))<F(39,500)
    assert 22*q*hi(1+F(V,44))<r

def half_ratio(c):
    return (c+r/(1-r)**2-r)/(1-c-r*r/(1-r))

assert half_ratio(F(39,500))<F(1,4)
A=F(27,2500); d=F(1,20); D=F(2,5)
p=F(63,44)
d0=A/p
assert F(45,44)*14*U/2<A
assert max(d0*hi(p),d*hi(A/d))<F(63,1000)
assert e/2*max(d0*hi(2*p),d*hi(2*A/d))<r
assert half_ratio(F(63,1000))<F(1,4)
assert (d/D)*hi(8*A*(1/d-1/D))<1
assert all((d/D)**(k-1)*hi(k*k*A*(1/d-1/D))<1
           for k in range(2,8))
assert 2*7*F(45,44)*U/d<F(1,2)
assert D*hi(1+8*A/D)/8<r
a=[(F(6,5)*D)**k*hi(A*k*k/D)/factorial(k) for k in range(1,8)]
s=(F(9,35)**8)/(1-F(9,35))
assert sum(a)+s<F(13,20)
assert sum(k*b for k,b in enumerate(a,1))<F(21,25)
assert sum(k*(k-1)*b for k,b in enumerate(a,1))<F(9,20)
assert F(45,44)*14*U<F(27,1250)
assert (F(21,25)+8*F(27,1250)/D*F(9,20))/44+2*s<F(3,125)
assert F(3,125)/F(7,20)<F(7,100)
assert lo(F(11,10))>F(20,7)

assert D*lo(F(13,2))>266
assert F(267,266)*U<F(151,100000)
assert F(267,266)*14*U<F(2113,100000)
q=F(151,100000); Wq=F(2113,100000)
assert 34*q<F(13,250)
assert (1+34*q/(1-F(13,250)))/(1-F(13,250))<F(28,25)
assert F(13,250)/(1-F(13,250))<F(11,200)
assert F(28,25)*Wq/D<F(3,50)
assert hi(Wq/(2*D))<F(103,100)
beta=F(28,25)*hi(F(14,25))
assert beta*34*Wq<F(71,50)
assert beta*F(103,100)*Wq<F(43,1000)
assert F(13,2)*81/F(14)**2<lo(1)
assert F(10,9)**10<F(29,10)
alpha=F(11,200)
p=F(29,10)*F(43,1000)/(1-alpha)*hi(F(71,50)*(F(1,2)+alpha/(6*(1-alpha))))
assert p<F(7,25)
assert hi(1+34*Wq/2)/14<F(1,3)
assert lo(F(13,5)*F(7,25))>F(125,62)
assert lo(F(3,2)*F(7,25))>F(125,83)
assert hi(F(2,5))<F(3,2)
assert F(2,5)*(14*D+8)-F(41,10)*F(103,100)*D>F(15,4)
assert hi(F(4,5)+F(3,2)*F(103,100)/14)<F(5,2)
assert 1/lo(F(15,4))+1/lo(F(52,5))<F(3,125)
assert lo(F(1,40))>F(125,122)
assert 1+4/lo(1)<F(5,2)
assert F(5,2)*F(103,100)*F(6,5)/(1-F(42,125))<5
assert 5*q+F(2,40)<F(3,50)
assert 5*(F(3,2)*F(103,100)*q+F(1,40*266))<F(1,2)
assert lo(2)>7 and F(1,7)/(1-F(1,7))/(1-F(2,7))<F(1,4)
~~~
