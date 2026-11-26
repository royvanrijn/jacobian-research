# All-degree angular interlacing for mirrored curvature at least fourteen

Let the unequal rows of an ordinary two-monomial profile occur in
mirrored pairs \((\ell,\ell+d),(\ell+d,\ell)\), with
\(\ell\ge0,d\ge1\). An optional common Gaussian strength \(\tau\ge0\)
multiplies the normalized coefficient at \(k\) of order \(m\)
by \(e^{-\tau k(m-k)/(2m)}\). Write \(S\) for the ordinary absolute
curvature plus \(\tau\).

**Theorem.** If \(S\ge14\), the consecutive endpoint-normalized
polynomials strictly interlace in the angular ordering from \(1\)
on the unit circle. At \(\tau=0\) they satisfy SFC(2) in every
degree and over every characteristic-zero field.

The exact matched tanh rates give the endpoint mass of the
ordinary pair without separating the lower and upper arrays.
The mirrored row pairing gives a sharper bound on the product
of that mass with the largest relative gap. In the intermediate
range the proof compares the actual half coefficients at the
two orders. It uses no tests of individual moment polynomials.

## Exact matched tanh endpoint mass

If there is a private pair, the
[private mirrored geometric estimate](MULTIVARIATE_MIRRORED_CORE_GAUSSIAN_ENDPOINT_TRANSFER.md#endpoint-comparison-with-the-exact-gaussian-clock)
already gives the assertion at \(S\ge14\). Hence assume
\(\ell\ge1\) on all pairs. Write
\[
 V=2\sum_{\rm pairs}d^2/\ell+\tau,\qquad W=1+V,\qquad
 L_0=\max(1,\max_{\rm pairs}d/\ell).
\]
At any adjacent orders \(n,N=n+1\), assign to each matched
Gauss factor the rate
\[
             v_i=\left(\frac{n}{T_i'}-\frac{n-1}{T_i}\right)^{-1}.
\]
By the [matched Gauss grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid),
a mirrored pair has \(2d\) rates equal to \(d/(2\ell+d)\).
The exact common endpoint envelope for the two orientations is
\[
 Q(w)=e^{-\tau(1-w)/2}
      \prod_{\rm pairs}\left(\frac{\ell+dw}{\ell+d}\right)^d,
 \qquad q=Q(0),\qquad L_*=\frac V2 .
\]
In the matched rates this gives
\[
 \log Q(w)=-\frac{\tau(1-w)}2+
        \frac12\sum_i
           \log\frac{1-(1-2w)v_i}{1+v_i},
 \qquad -\log q=\frac{\tau}{2}+\sum_i\operatorname{artanh}v_i .
                                                               \tag{1}
\]
Thus the endpoint mass uses the complete odd positive tanh
moment series at the same linked rates at both orders. The
curvature in the theorem is likewise exactly
\(S=\tau+2\sum_i v_i/(1+v_i)\). This identity includes private
pairs, since a private row has \(v_i=1\) in each of its Gauss slots.

As in the
[retained mirrored endpoint estimate](MULTIVARIATE_RETAINED_MIRRORED_ENDPOINT_CRITERION.md#the-two-coefficient-estimates),
for \(1\le k\le m/2,\ w=k/m\), the normalized coefficients obey
\[
 c_{m,k}^{1/k}\le eQ(w)/w,\qquad
 c_{m,k}\le (mq)^k e^{L_* k^2/m}/k!.                           \tag{2}
\]
The second inequality holds also for real \(m\) and integer \(k\)
in this range, using the Gamma interpolation and the same
finite row products. Also \(Q(w)\le q e^{L_*w}\).
The maximum of \(eQ(w)/w\) on any closed \(w\)-interval
is at an endpoint. For each fixed positive \(j\), the same
holds for \(mQ(j/m)\) on every closed bounded subinterval
of \(m\ge2j\).
Both properties follow by logarithmic differentiation.
The midpoint satisfies \(Q(1/2)\le e^{-S/4}\) by
\(\log(1-r/(2(1+r)))\le-r/(2(1+r))\), \(r=d/\ell\).

The endpoint mass obeys
\[
             q\le e^{-S/2},\qquad Wq\le(1+S)e^{-S/2},
                     \qquad L_0q\le e^{-S/2}.                \tag{3}
\]
For the first two inequalities, the \(q\) in (1)
is also the zero-endpoint value of the fully retained
oriented entropy envelopes, because the two signed
chord slopes in a pair add to \(-d\log(1+d/\ell)\).
The
[entropy absorption proof](MULTIVARIATE_ABSOLUTE_CURVATURE_20.md#the-endpoint-cost-absorbs-the-excess-curvature)
therefore applies. Its derivative test holds already at
\(C=1+S\ge15\): at \(C=15\) the difference polynomial is
\(r^3-4r^2+5r+9=r[(r-2)^2+1]+9>0\).

For the third bound in (3), the cost of a pair beyond its
contribution to \(S/2\) is
\(\ell[r\log(1+r)-r^2/(1+r)]\).
For \(r\ge1\) the expression in brackets is at least
\(\log r\). To verify this, use
\(\log(1+r)\ge 2(u+u^3/3)\), \(u=r/(r+2)\), and
\(\log r\le(r-1)-(r-1)^2/(2r)\), the latter by integration.
The difference of the resulting bounds has denominator
\(6r(1+r)(2+r)^3\) and numerator, with \(t=r-1\ge0\),
\[
 7t^6+49t^5+113t^4+63t^3-74t^2-28t+62.
\]
Here \(113t^4-74t^2+13>0\) by its quadratic discriminant,
and \(63t^3-28t+49>0\) separately on \(t\le1\) and \(t\ge1\).
The row attaining \(L_0>1\) thus supplies cost at least
\(\log L_0\) beyond \(S/2\), proving (3).
The case \(L_0=1\) is immediate.

We may now assume \(14\le S<16\) by
[the absolute curvature sixteen theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md).
For its use here, absolute angular separation implies strict
interlacing on the mirrored class at \(S\ge16\): increase the
Gaussian strength to infinity, preserving symmetry and
angular separation. At infinity the endpoint-normalized
polynomials tend to \(1+z^m\), whose consecutive roots
strictly interlace. The intervening polynomials have simple
unit-circle roots by
[circular composition](MULTIVARIATE_MIRRORED_CORE_SIGNED_DIFFUSE_RADIAL_SEPARATION.md#a-comparison-at-exactly-the-same-signed-center),
and no cross-degree coincidence or root through the positive
cut occurs. Also \(V\le2S+S^2/2<160\), by \(d/\ell\le1+S_i/2\)
for the curvature \(S_i\) of each mirrored pair.

## Half phases for small endpoint mass

Write \(N=n+1\), \(D=Nq\), and \(F_m(z)=\sum_k c_{m,k}z^k\).
Use its symmetric half \(H_m\) with half weight on the middle
coefficient when present. At \(m=1\) it is \(H_1=1\).
Whenever \(H_m\) is nonvanishing on the unit disk,
\[
 e^{-im\theta/2}F_m(e^{i\theta})
  =2|H_m(e^{i\theta})|\cos\psi_m(\theta),\qquad
 \psi_m(\theta)=m\theta/2-\phi_m(\theta),\quad
 \phi_m(\theta)=\Im\log H_m(e^{i\theta}).
\]
Here the logarithm is zero at zero and is real on the real
diameter. In particular \(\phi_m(0)=\phi_m(\pi)=0\).
It suffices to prove \(\psi_m'>0\) and
\(|\phi_N-\phi_n|<\theta/2\) on \(0<\theta\le\pi\).
Then at a lower root \(\psi_N-\psi_n\) is strictly between
zero and \(\theta\le\pi\). Monotonicity, the endpoint phases,
and reflection give strict interlacing.

First suppose \(N\le64\). The \(k=1\) bound is maximized
by the checks at \(m=2,64\). For \(k\ge2\) use (2) at
\(w=2/m\) and \(w=1/2\), and test the first at \(m=4,64\).
At the \(64\) endpoints we use
\[
 q e^{jV/(2\cdot64)}
   \le e^{-7}\frac{15}{1+V}e^{jV/128},\qquad j=1,2 .
\]
The right side attains its maximum on \(14\le V\le160\)
at \(14\) or \(160\), by its increasing logarithmic
derivative. The scalar tests below give \(c_{m,1}<1/14\),
\(c_{m,k}<6^{-k}\) for \(k\ge2\).
Summing these bounds and their derivatives on the unit circle
gives \(|zH_m'/H_m|<1/4\). This proves the two phase inequalities.

Suppose \(N>64\), so \(n\ge64\). Equations (2)–(3) give
\(L_*/n\le5/4\) and \(L_*D/n<7/1000\).
If \(D\le1/20\), then for \(j=1,2\) the bound on
\(D e^{j L_*/n}\) attains its maximum at
\(D=(7/1000)/(5/4)\) or \(1/20\).
Below the switch it increases, and above the switch the
logarithmic derivative \((D-7j/1000)/D^2\) is first negative
and then positive. These bounds give \(c_{m,1}<3/50\),
\(c_{m,k}<6^{-k}\) for \(k\ge2\), using the \(k=2\)
and midpoint tests in (2). Again \(|zH_m'/H_m|<1/4\)
on the unit circle.

Now take \(1/20<D\le8/25\) and \(n\ge64\).
For \(1\le k\le8\), the bound \(D e^{(7/1000)k/D}\)
is maximized at \(D=8/25\) on this interval, by comparing
the two endpoints. Hence at radius \(6/5\) the first seven
terms are bounded by
\[
  b_k=\frac{((6/5)(8/25))^k}{k!}
              \exp\left(\frac{(7/1000)k^2}{8/25}\right).
\]
The bound holds for the actual coefficients at real
\(m\in[n,N]\) as well. The \(k=8\) and midpoint
checks in (2) give \(6^{-k}\) for the remainder of
either integer half. In particular
\(\|H_m-1\|_{6/5}<1/2\), with
\(\sum_{k=1}^7 k b_k<3/5\) and
\(\sum_{k=1}^7 k(k-1)b_k<1/4\).

To match the two halves, interpolate their first coefficient
\(B_m\) and \(f_m(x)=\log(c_{m,x}x!/B_m^x)\) with positive
Gamma values. The exact mirrored factorial ratios give
\[
            0\le1/m-(\log B_m)'\le V/m^2 .
\]
For each selected pair, differentiate
\((m\ell+h)/(m(\ell+d)-d+h)\); the absolute logarithmic
derivative is \(d(\ell+h)\) divided by the product of
the two arguments, so summing is at most the pair's
\(2d^2/(m^2\ell)\). The Gaussian adds
\(-\tau/(2m^2)\). Thus \(|(\log B_m)'|\le1/n\) here,
since \(V/n\le WD/(nD)<(14/1000)/D\le7/25\).

The [Gamma interpolation estimate](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md#coefficient-halves-and-endpoint-interpolation)
gives \(|\partial_m\partial_x^2 f_m(x)|\le16W/n^2\)
on \(0\le x\le7\).
Its argument lower bounds apply because
\(2L_0\,7/n<27/100\), by (3) and \(D>1/20\),
and \(n\ge64\). The Gaussian is absorbed in this
bound by the added term \(\tau x(x-1)/(2m)\).
Since \(\partial_m f_m\) vanishes at zero and one,
\[
 |\partial_m c_{m,k}|
       \le\frac{c_{m,k}}{n}
                \left(k+\frac{8W}{n}k(k-1)\right),\qquad k\le7.
\]
Summation using the two bounds on \(b_k\) gives the
low-coefficient difference on \(6/5\) at most
\((3/5+8(7/25)/4)/64\).
The two tails have total norm at most
\(2(1/5)^8/(1-1/5)\).
Thus \(\|H_N-H_n\|_{6/5}<1/50\).
The straight line between the two halves stays in
\(\|H-1\|<1/2\), which bounds the difference of their logs
by \(2\|H_N-H_n\|<1/25\).
Cauchy's estimate gives paired angular derivative below \(1/5\)
and individual derivatives below \(5\), using \(\|\log H_m\|<1\).
Since \(m\ge64\), this proves the phase inequalities in the
intermediate range.

## Long halves

Suppose \(D>8/25\). Here \(n\ge350\) by \(q\le e^{-7}\).
Put \(K=\lceil16D\rceil+7\le41D\).
The budgets (3) give
\[
 D/n<915/10^6,\quad WD/n<1373/10^5,\quad
 L_0D/n<915/10^6,\quad
             L_0K/n<19/500,\quad K/n<19/500 .
\]
In particular both halves contain the coefficients through \(K\).
At \(k=K\), (2) bounds \(c_{m,k}^{1/k}\) by
\((e/16)e^{41(1373/10^5)/2}<1/4\).
At the midpoint the bound is even smaller, so the
geometric tail from \(K\) is \(4^{-k}\).
The same first-factor derivative proof gives
\(|(\log B_m)'|\le1/n\) and makes \(B_m/m\) decreasing.
Consequently \(B_m\le D e^{L_*/n}<(103/100)D\).

Apply the
[sharper Gamma and connected sum comparison](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md#the-long-range-with-sharper-gamma-derivative-estimates)
with \(\kappa=19/500\) and \(\alpha=L_0/[(1-\kappa)n]\).
The polygamma factor
\((1+K/[(1-\kappa)n])/(1-\kappa)<109/100\)
gives \(H_2<(109/100)W/n<47/1000\).
Thus \(\alpha K<1/25\), and with \(\beta=e^{H_2+1/2}H_2\)
we have \(\beta K<107/100\) and \(\beta B_m<27/1000\).
The same tree base bound \(18/5\) gives
\(p_m=\Xi B_m<9/50\). Since
\(|(\log B_m)'|\le1/n\), the derivative estimate is now
\(|\ell_k'(m)|\le(5/2)B_m p_m^{k-1}/n\).

The absolute coefficient norms of the truncated logarithm
\(A_m\) are therefore below \((11/5)B_m\) on radius \(9/5\)
and \((7/5)B_m\) on radius \(6/5\).
Its difference across orders on \(6/5\) is below \(4D/n\).
The polynomial \(H_m\) and \(e^{A_m}\) agree through degree \(K\).
The exponential tail and actual polynomial tail there,
relative to \(e^{A_m}\), are bounded by
\[
           (2/3)^{K+1}e^{(18/5)B_m},\qquad
                              (3/10)^Ke^{(7/5)B_m}.
\]
The first is below \(e^{-4}\), using \(K+1\ge16D+8\)
and \(D>8/25\); the second is below \(e^{-K}\), \(K\ge13\).
Their sum is below \(19/1000\).
Thus \(\|\log H_m-A_m\|_{6/5}<1/50\).
In particular the paired logarithm difference is below
\(4D/n+2/50<1/20\), and the individual norm is below
\((7/5)B_m+1/50\).
Cauchy's estimates give paired angular derivative below
\(1/4\) and individual derivatives below \(m/2\).
The phase proof closes for this last range too.

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

u=1/lo(7)
assert 2/lo(F(7,2))<F(1,14)
assert 2/lo(F(5,2))<F(1,6)
for x in [14,160]:
    e=u*F(15,1+x)
    assert 64*e*hi(F(x,128))<F(1,14)
    assert 32*e*hi(1+F(x,64))<F(1,6)
c=F(1,14)
assert (c+F(1,6)/(1-F(1,6))**2-F(1,6))/(1-c-F(1,36)/(1-F(1,6)))<F(1,4)

A=F(7,1000); D=F(8,25); d=F(1,20)
assert F(65,64)*15*u<F(14,1000)
assert F(65,64)*u<F(93,100000)
d0=A/F(5,4)
assert d0*hi(F(5,4))<F(3,50) and d*hi(A/d)<F(3,50)
assert hi(1)*max(d0*hi(F(5,2)),d*hi(2*A/d))/2<F(1,6)
assert (d/D)*hi(8*A*(1/d-1/D))<1
assert 2*7*F(93,100000)/d<F(27,100)
assert D*hi(1+8*A/D)/8<F(1,6)
a=[(F(6,5)*D)**k*hi(A*k*k/D)/factorial(k) for k in range(1,8)]
small=F(1,5)**8/(1-F(1,5))
assert sum(a)+small<F(1,2)
assert sum(k*x for k,x in enumerate(a,1))<F(3,5)
assert sum(k*(k-1)*x for k,x in enumerate(a,1))<F(1,4)
assert (F(3,5)+8*F(14,1000)/d*F(1,4))/64+2*small<F(1,50)

assert D*lo(7)>350
assert F(351,350)*u<F(915,10**6)
assert F(351,350)*15*u<F(1373,10**5)
qq=F(915,10**6); wq=F(1373,10**5)
assert 41*qq<F(19,500)
assert F(1)/(1-F(19,500))*(1+41*qq/(1-F(19,500)))<F(109,100)
assert F(19,500)/(1-F(19,500))<F(1,25)
assert F(109,100)*wq/D<F(47,1000)
assert hi(wq/(2*D))<F(103,100)
beta=F(109,100)*hi(F(547,1000))
assert beta*41*wq<F(107,100)
assert beta*F(103,100)*wq<F(27,1000)
alpha=F(1,25)
assert F(18,5)*F(27,1000)/(1-alpha)*hi(F(107,100)*(F(1,2)+alpha/(6*(1-alpha))))<F(9,50)
assert hi(1+41*wq/2)/16<F(1,4)
assert lo(F(11,5)*F(9,50))>F(250,169)
assert lo(F(7,5)*F(9,50))>F(125,98)
assert hi(F(2,5))<F(3,2)
assert F(2,5)*(16*D+8)-F(18,5)*F(103,100)*D>4
assert 1/lo(4)+1/lo(13)<F(19,1000)
assert lo(F(1,50))>F(1000,981)
assert hi(F(11,10))<F(10,3) and F(7,5)*F(103,100)/16<F(1,10)
assert 1+4/lo(1)<F(5,2)
assert F(5,2)*F(103,100)*F(6,5)/(1-F(27,125))<4
assert 4*qq+F(2,50)<F(1,20)
assert 5*(F(7,5)*F(103,100)*qq+F(1,50*350))<F(1,2)
~~~
