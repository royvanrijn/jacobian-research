# All-degree angular interlacing for mirrored curvature at least twelve

Let the unequal rows of an ordinary two-monomial profile occur in
mirrored pairs \((\ell,\ell+d),(\ell+d,\ell)\), with
\(\ell\ge0\) and \(d\ge1\) integers. Allow a common Gaussian
strength \(\tau\ge0\), which contributes the normalized
coefficient factor \(e^{-\tau k(m-k)/(2m)}\) in degree \(m\).
Write \(S\) for the absolute curvature including \(\tau\).

**Theorem.** If \(S\ge12\), consecutive normalized polynomials
have strictly interlacing root arguments on the unit circle in
the ordering cut at \(1\), in every degree. At \(\tau=0\) this
proves ordinary SFC(2) over every characteristic-zero field.

Use the [curvature thirteen theorem](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_THIRTEEN.md)
for \(S\ge13\). The remaining band is \(12\le S<13\).
The proof retains the exact endpoint mass from the matched
tanh rates and bounds its product with the square of the
curvature budget. It uses the sharper binomial entropy at
the midpoint and a less costly derivative estimate in the
connected expansion at large endpoint mass.

## Exact endpoint mass and short coefficient envelopes

For the normalized polynomial of degree \(m\) write \(H_m\)
for its half with half weight on the middle coefficient when
present. At \(m=1\) use \(H_1=1\).
If \(\log H_m\) is analytic on the closed unit disk with
\(\log H_m(0)=0\), put \(\phi_m(\theta)=\Im\log H_m(e^{i\theta})\).
The centered polynomial is
\(2|H_m(e^{i\theta})|\cos(m\theta/2-\phi_m(\theta))\).
For \(N=n+1\), it suffices that the two centered phases
increase and \(|\phi_N-\phi_n|<\theta/2\) for \(0<\theta\le\pi\).
At any lower root the upper phase is then strictly between
the same phase level and the next. The real coefficients
give \(\phi_m(0)=\phi_m(\pi)=0\). The endpoints and reflection
give strict interlacing.

If there is a private pair, reserve a factorial from one
private row on each side. The binomial cancels and the
[private reservation estimate](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_CURVATURE_INTERLACING.md#a-common-radial-bound)
gives \(c_{m,k}\le e^{-(S-2)k(m-k)/(2m)}\) for normalized
coefficients, with the Gaussian included.
Thus \(c_{m,k}\le e^{-5k/2}<10^{-k}\) for \(k\le m/2\).
The sum of coefficients of \(H_m-1\) on the unit disk is
below \(1/9\) and the derivative sum is below \(10/81\).
Consequently \(|zH_m'/H_m|<1/4\). The two required
phase estimates follow directly.

Assume all pairs are shared. Use \(r=d/\ell\) in each pair and
the quantities \(V,W=1+V,L_0,L_*=V/2,Q_M,q\) of the
[matched tanh endpoint comparison](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_FOURTEEN.md#exact-matched-tanh-endpoint-mass),
with \(Q_M\) denoting its rational envelope including the
Gaussian. Thus \(S_i=2d^2/(\ell+d)\), \(V_i=2d^2/\ell\),
\(12\le V<221/2\), and
\[
  S=\tau+2\sum_i\frac{v_i}{1+v_i},\qquad
  -\log q=\tau/2+\sum_i\operatorname{artanh}v_i,\qquad
  v_i=\left(\frac{n}{T_i'}-\frac{n-1}{T_i}\right)^{-1}.
\]
The sums are over matched Gauss slots, with \(2d\) rates
\(v_i=d/(2\ell+d)\) from each mirrored pair.
We have \(Q_M(w)\le q e^{L_*w}\),
\(Q_M(1/2)\le e^{-S/4}\), and \(L_0q\le e^{-S/2}\).
Also \(q\le e^{-S/2}\) and
\[
                      W^2q\le (1+S)^2 e^{-S/2}.                \tag{1}
\]
For (1) put \(C=1+S\ge13\),
\(K(r)=r\log(1+r)-r^2/(1+r)\) and \(g(r)=2r^3/(1+r)\).
The excess cost from a pair above \(S_i/2\) is \(\ell K(r)\)
and its excess curvature \(V_i-S_i\) is \(\ell g(r)\).
The inequality \(\log(1+r)\ge2r/(r+2)\) gives
\[
 K'(r)\ge \frac{r^2(3+2r)}{(r+2)(1+r)^2}
        \ge \frac{2g'(r)}{C+g(r)}.
\]
For the second step it suffices that \(C+g(r)\ge4(r+2)\).
After multiplication by \(1+r\) and setting \(C=13\), the
difference is \(2r^3-4r^2+r+5>0\). On \([0,2]\) use
\(2r^2(2-r)\le64/27\), and on \([2,\infty)\) it follows
from \(2r^3-4r^2\ge0\). Integrating gives
\(K(r)\ge2\log(1+g(r)/C)\).
Since \(\ell\) is a positive integer,
\((1+g(r)/C)^\ell\ge1+\ell g(r)/C\).
Multiplication over pairs gives
\(\exp(\sum\ell K(r)/2)\ge1+(V-S)/C=W/C\), proving (1).
It also gives \(Wq\le (1+S)e^{-S/2}\) since \(W\ge1+S\).
Both \((1+S)^2e^{-S/2}\) and \((1+S)e^{-S/2}\)
are decreasing for \(S\ge12\).

For small \(D=Nq\) use the
[oriented entropy envelope](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md)
\(Q_E(w)\) with every shared row retained and residual
curvature \(\tau\). It has \(Q_E(0)=q\), and
\(Q_E(1/2)\le e^{-S/4}\) by the retain-all estimate there.
In the notation \(A(w),G(w)\) for the row chords, the
two \(A(w)/w\) derivatives in a pair with endpoints \(a<b\)
sum to at most \(d^2(1/a+1/b)/2=(V_i+S_i)/4\).
Indeed the reciprocals at \(a+dwu,b-dwu\) have sum
at most \(1/a+1/b\) for \(0\le w\le1/2\), \(0\le u\le1\).
The sum of the \(G/2\) derivatives is at most its
tangent slope \(d^2/(2ab)\le S_i/4\).
With the Gaussian slope this gives
\[
              Q_E(w)\le q e^{L_E w},\qquad L_E=V/4+S/2.        \tag{2}
\]
For the binomial entropy factor
\(B(w)=\exp([-w\log w-(1-w)\log(1-w)]/w)\), the bound is
\(c_{m,k}^{1/k}\le B(w)Q_E(w)\), where \(w=k/m\).
We need that \(B(w)Q_E(w)\) attains its maximum on each
closed subinterval of \(0<w\le1/2\) at an endpoint. To see
this, set \(\Phi=w(\log Q_E)'\). The contribution to
\(\Phi'\) from the two \(A/w\) terms is at least
\[
 \frac{d^2}{2}\left(\frac{a}{(a+d/2)^2}+\frac1b\right),
\]
using \(d^2\int_0^1 u[a/(a+dwu)^2+b/(b-dwu)^2]\,du\).
The two \(G/2\) contributions to \(\Phi'\) sum to
\(d[a/(a+dw)^2-b/(b-dw)^2]/2\), which is decreasing
and bounded below by \(-d^2/[2(a+d/2)^2]\).
As \(a\ge1\), their combined lower bound is
\(d^2/(2b)=S_i/4\). Including the Gaussian gives
\(\Phi'\ge S/4\).
The
[binomial entropy derivative](MULTIVARIATE_ENTROPY_ENDPOINT_PHASE_CRITERION.md#combining-the-binomial-entropy-and-the-row-chord)
is \([w(\log B)']'\ge-4(1-\log2)\).
Thus \(w(\log(BQ_E))'\) is increasing.
The same sign argument gives endpoint maxima for
\(mQ_E(1/m)\) on closed intervals in \(m\ge2\) and
for \(B(2/m)Q_E(2/m)\) on those in \(m\ge4\).

## Short and medium halves

First take \(N\le28\).
At \(m=2\), \(c_{m,1}\le2e^{-3}<1/10\).
At \(m=28\) the tests
\(28Q_E(1/28)<1/10\) and \(14e Q_E(1/14)<1/5\)
follow from (1)–(2) using
\[
 q e^{jL_E/28}\le
       e^{-6}\left(\frac{13}{1+V}\right)^2
                       e^{j(V/4+13/2)/28},\qquad j=1,2.
\]
Each right side attains its maximum on
\(12\le V\le221/2\) at an endpoint since its logarithmic
derivative is increasing. The midpoint test is
\(B(1/2)Q_E(1/2)\le4e^{-3}<1/5\).
The coefficient and order endpoint maxima therefore give
\(c_{m,1}<1/10\) and \(c_{m,k}\le5^{-k}\) for \(k\ge2\)
in both halves.

Next take \(N>28\) and \(D\le7/100\).
From (1) and \(q=D/N\),
\[
                 L_*/n\le \min(2,1/(16\sqrt D)).
\]
Indeed \(V/(2n)<2\), and \((L_*/n)^2D\) is at most
\(169 e^{-6}(29/28)/(4\cdot28)<1/256\).
The maximum of \(D\exp(\min(2,1/(16\sqrt D)))\)
is at \(D=1/1024\) or \(7/100\).
On the part beyond the switch its logarithmic derivative
has at most one sign change, from negative to positive.
Thus \(D e^{L_*/n}<1/10\), and (2) gives
\[
       \frac{e}{2}D e^{2L_E/n}
        \le \frac12 e^{1+13/28}D e^{L_*/n}<1/5 .
\]
These are the tests for \(k=1,2\), using \(Q_M\) for the
first, \(B(w)\le e/w\) with \(Q_E\) for the second, and
\(2L_E=L_*+S\).
The midpoint and entropy maxima yield the same short caps.
For both short cases the geometric series with first term
bounded by \(1/10\) and remaining coefficients by \(5^{-k}\)
gives \(|zH_m'/H_m|<1/4\) on the unit disk, proving the
phase inequalities.

For \(7/100<D\le12/25\), still with \(n\ge28\), write
\(d_*=7/100\), \(D_*=12/25\), \(A_0=21/1250\).
The endpoint budgets give \(D/n,L_0D/n<257/100000\),
\(WD/n<167/5000\), and \(L_*D/n<A_0\).
The rational envelope gives \(c_{m,k}\le D^k e^{A_0k^2/D}/k!\)
for the actual inclusive coefficients and their real-degree
interpolation in the first seven. For \(1\le k\le8\) the
maximum of \(D e^{A_0k/D}\) on \([d_*,D_*]\) is at \(D_*\).
Its logarithmic derivative changes sign at most once from
negative to positive, and the two endpoints compare at \(k=8\).
The test at \(k=8\) and the midpoint test with \(eQ_M(w)/w\)
give \(c_{m,k}\le(3/10)^k\) for \(k\ge8\) in the half.

For the adjacent difference let \(B_m\) be the exact first
coefficient at real \(m\in[n,N]\), and interpolate
\(f_m(x)=\log(c_{m,x}x!/B_m^x)\) with Gamma functions.
The [exact mirrored first-factor derivative](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_FOURTEEN.md#half-phases-for-small-endpoint-mass)
gives \(0\le1/m-(\log B_m)'\le V/m^2\). In particular
\(|(\log B_m)'|\le1/n\), since \(V/n\le(167/5000)/D<1\).
All Gamma arguments in \(\partial_m\partial_x^2 f_m\)
are at least \(ma_i/2\) for a row starting at \(a_i\)
and \(m/2\) for the binomial, since \(7L_0/n<1/2\).
The bound \(|\psi_2(y)|\le2/y^2\) for \(y\ge2\)
therefore gives \(|\partial_m\partial_x^2 f_m|\le16W/n^2\)
including the binomial and the Gaussian.
As \(\partial_m f_m\) vanishes at zero and one,
\[
 |\partial_m c_{m,k}|\le\frac{c_{m,k}}n
            \left(k+\frac{8W}{n}k(k-1)\right),\qquad 1\le k\le7.
\]
In the first term use \(c_{m,k}/D\) and \(D/n\).
In the second use \(c_{m,k}/D^2\) and the product
\((D/n)(WD/n)\). For \(k\ge j\), the maximum of
\(D^{k-j}e^{A_0k^2/D}\) is at one of \(d_*,D_*\), by the
same logarithmic derivative test, including \(k=j\).
The scalar sums below give
\[
 \|H_m-1\|_{6/5}<17/20,\qquad
 \|H_N-H_n\|_{6/5}<3/250 ,
\]
where \(\|\cdot\|_R\) is the absolute coefficient norm on
radius \(R\). For both bounds the tail on \(6/5\) has ratio
\(9/25\). The straight line between the halves gives
\(\sup_{|z|\le6/5}|\log H_N-\log H_n|<2/25\),
by division by \(1-17/20\).
The individual log sup norm is below \(\log(20/3)<2\).
Cauchy's derivative estimate gives paired angular derivative
below \(2/5\) and the individual ones below \(10<m/2\).

## Long halves

Take \(D>12/25\). Then \(n\ge193\),
\(D/n,L_0D/n<1/400\), and \(WD/n<81/2500\).
Use \(K=\lceil10D\rceil+7\le(80/3)D\).
At \(k=K\) the rational bound for \(c_{m,k}^{1/k}\) is
at most \((e/10)e^{WK/(2n)}<1/2\). The midpoint bound
is \(2e^{-2}<1/2\), hence \(c_{m,k}\le 2^{-k}\) from \(K\)
to the midpoint. The first-factor derivative above gives
\(|(\log B_m)'|\le1/n\), and the rational envelope gives
\(B_m\le De^{L_*/n}<(207/200)D\).

Both \(K/n,L_0K/n\) are below \(\kappa=67/1000\), and
\(K+1<n/2\). Use the
[sharp Gamma derivative bound](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md#the-long-range-with-sharper-gamma-derivative-estimates)
with this \(\kappa\) and \(\alpha=L_0/[(1-\kappa)n]\).
It gives \(H_2<(29/25)W/n<2/25\) and
\(u=\alpha K<9/125\), so all finite differences
\(\nu_r\), \(2\le r\le K\), are bounded by
\(H_r=(r-2)!\alpha^{r-2}H_2\le H_2\),
with derivative bound \(2rH_r/n\).

In the
[differentiated connected expansion](MULTIVARIATE_ABSOLUTE_CURVATURE_20.md#a-smaller-differentiated-connected-logarithm-bound),
each edge has weight \(e^{\nu_r}-1\) and its positive bound
\(e^{H_2}H_r\). Differentiation costs \(2X/n\), where \(X\)
is the sum of edge sizes on the \(k\) labels.
Use \(X\le(8/e)k e^{X/(8k)}\).
Since each edge has size at most \(k\), inflate the edge
bound by \(e^{1/8}\), leaving a derivative factor \(16k/(en)\).
For \(\beta=e^{H_2+1/8}H_2\),
\(\beta K<31/25\), \(\beta B_m<6/125\).
The shared mirror budget \(L_0\le\sqrt{V/2}\) gives
\(\alpha/\beta\le e^{-1/8}\sqrt6/13<1/6\).
The discovery-tree base is therefore at most
\((7/6)^7<59/20\). With the additional-edge factor in the
cited estimate, the resulting connected parameter satisfies
\[
 p_m=\Xi B_m
 \le \frac{(59/20)\beta B_m}{1-u}
        \exp\!\left(\beta K\left(\frac12+\frac{u}{6(1-u)}\right)\right)
 <29/100 .
\]
The formal log coefficients through \(K\) consequently satisfy
\[
 |\ell_k(m)|\le B_m p_m^{k-1}/k,\qquad
 |\ell_k'(m)|\le7B_m p_m^{k-1}/n ,
\]
since \(1+16/e<7\). The unit-degree and second-degree
bounds are included in the same discovery-tree estimate.

For the truncated log \(A_m=\sum_{k\le K}\ell_k(m)z^k\),
summing these geometric estimates gives
\(\|A_m\|_2<3B_m\) and \(\|A_m\|_{6/5}<(3/2)B_m\).
The adjacent difference on \(6/5\) is below \((67/5)D/n\).
The halves and \(e^{A_m}\) agree through degree \(K\).
The exponential tail and actual polynomial tail on \(6/5\),
relative to \(e^{A_m}\) and measured in absolute coefficient
norm, are bounded respectively by
\[
      (3/5)^{K+1}e^{(9/2)B_m}<e^{-4},\qquad
      \tfrac32(3/5)^K e^{(3/2)B_m}<\tfrac32e^{-7K/20-1}.
\]
For the first use \(K+1\ge10D+8\) and \(\log(5/3)>1/2\).
For the second keep \(D\le(K-7)/10\).
Their sum is below \(27/1000\) since \(K\ge12\).
The logarithmic correction on \(6/5\) is below \(7/250\).
Thus the paired log difference is below \(1/10\) and the
individual norm is below \((3/2)B_m+7/250\).
Cauchy's estimates give paired angular derivative below \(1/2\)
and the individual derivatives below \(m/2\).
This completes the all-degree phase comparison.

## Exact scalar checks

The checks use rational bounds on the exponential series. They
involve only the displayed scalar estimates and no moment
polynomials or profile enumeration.

~~~python
from fractions import Fraction as F
from math import factorial

def lo(x):
    x=F(x)
    return sum(x**j/factorial(j) for j in range(65))

def hi(x):
    x=F(x)
    return lo(x)+x**65/factorial(65)/(1-x/66)

def ck(label, a,b):
    assert a < b, (label, float(a),float(b))

e=hi(1); U=1/lo(6)
ck('midQ',4/lo(3),F(1,5))
ck('m2',2/lo(3),F(1,10))
ck('entropy log2',hi(F(69,100)),2)
for V in [F(12),F(221,2)]:
    q=U*(13/(1+V))**2
    L=V/4+F(13,2)
    ck('m28k1',28*q*hi(L/28),F(1,10))
    ck('m28k2',14*e*q*hi(2*L/28),F(1,5))
ck('Vshort',F(221,112),2)
ck('squared short', F(1,4)*169*U*F(29,28)/28,F(1,16)**2)
d=F(7,100); D=F(12,25); A=F(21,1250)
small=max(hi(2)/1024,d*hi(F(6,25)))
ck('short exp',F(1,16)**2,d*F(6,25)**2)
ck('short B',small,F(1,10))
ck('short k2',small*hi(1+F(13,28))/2,F(1,5))
assert ((F(1,10)+F(1,5)/(1-F(1,5))**2-F(1,5)) /
        (1-F(1,10)-F(1,5)**2/(1-F(1,5))))==F(1,4)
q0=F(257,100000); W0=F(167,5000)
ck('medium q', F(29,28)*U,q0)
ck('medium Wq',F(29,28)*13*U,W0)
ck('medium W',W0/d,1)
ck('medium L0',7*q0/d,F(1,2))
ck('medium k8',d/D*hi(A*8*(1/d-1/D)),1)
ck('medium k8 exp',D*hi(1+A*8/D)/8,F(3,10))
ck('midQ rational',2/lo(2),F(3,10))
def b(k,j,m):
    return F(6,5)**k*m**(k-j)*hi(A*k*k/m)/factorial(k)
b0=[b(k,0,D) for k in range(1,8)]
s=F(9,25)**8/(1-F(9,25))
ck('medium norm',sum(b0)+s,F(17,20))
m1=sum(k*max(b(k,1,d),b(k,1,D)) for k in range(1,8))
m2=sum(k*(k-1)*max(b(k,2,d),b(k,2,D)) for k in range(2,8))
ck('medium m1',m1,F(27,10)); ck('medium m2',m2,F(11,2))
ck('medium diff',q0*(F(27,10)+8*W0*F(11,2))+2*s,F(3,250))
ck('medium log diff',F(3,250)/F(3,20),F(81,1000))
ck('medium log norm',F(20,3),lo(2))

q1=F(1,400); W1=F(81,2500); C=F(80,3); kap=F(67,1000)
ck('long n',193,D*lo(6))
ck('long q',F(194,193)*U,q1); ck('long Wq',F(194,193)*13*U,W1)
assert 10+8/D == C
ck('long kappa',C*q1,kap)
ck('long degree',kap+F(1,193),F(1,2))
ck('long first derivative',W1/D,1)
ck('long tail',hi(1+C*W1/2)/10,F(1,2))
ck('long midpoint',2/lo(2),F(1,2))
ck('long B',hi(W1/(2*D)),F(207,200))
ck('long H2 prefactor',(1+C*q1/(1-kap))/(1-kap),F(29,25))
ck('long H2',F(29,25)*W1/D,F(2,25))
u=F(9,125); c=F(31,25); beta=F(29,25)*hi(F(41,200))
ck('long u',kap/(1-kap),u)
ck('long betaK',beta*C*W1,c)
ck('long betaB',beta*F(207,200)*W1,F(6,125))
ck('long t',F(216,169),lo(F(1,4)))
ck('long base',F(7,6)**7,F(59,20))
p=F(29,100)
ck('long p',F(59,20)*F(6,125)/(1-u)*hi(c*(F(1,2)+u/(6*(1-u)))),p)
ck('long derivative',1+16/lo(1),7)
ck('long outer norm',1/lo(3*p),1-2*p)
ck('long inner norm',1/lo(F(3,2)*p),1-F(6,5)*p)
ck('long log diff',7*F(207,200)*F(6,5)/(1-F(6,5)*p),F(67,5))
ck('long exp tail slope',F(9,2)*F(207,200),5)
ck('long true tail slope',hi(F(7,20)+F(3,2)*F(207,200)/10),F(5,3))
ck('long true tail intercept',1,F(3,2)*F(207,200)*F(7,10))
ck('long combined tail',1/lo(4)+F(3,2)/lo(F(26,5)),F(27,1000))
ck('long log corrections',1/lo(F(7,250)),1-F(27,1000))
ck('long log diff final',F(67,5)*q1+2*F(7,250),F(1,10))
ck('long single phase',5*(F(3,2)*F(207,200)*q1+F(7,250)/193),F(1,2))
print('All mirrored S=12 scalar checks passed.')
~~~
