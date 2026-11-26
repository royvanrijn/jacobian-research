# All-degree angular separation at absolute curvature fifteen

Let \(a,b\) be distinct ordinary nonnegative exponent columns.
Omit equal rows and define the absolute curvature
\[
               S=\tau+\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)}.
\]
Here \(\tau\ge0\) is an optional common Gaussian strength,
with normalized factor \(e^{-\tau k(m-k)/(2m)}\) in degree \(m\).

**Theorem.** If \(S\ge15\), consecutive moment polynomials have
disjoint root arguments in every degree. There are no restrictions
on orientations or private rows. At \(\tau=0\) this proves
ordinary SFC(2) over every characteristic-zero field.

Use the [absolute curvature sixteen theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md)
for \(S\ge16\). On \(15\le S<16\), the proof combines a sharper
largest-gap budget with an exact first-coefficient derivative
that includes private ending rows. The latter comes directly from
the Gamma chord and applies to every orientation with positive
starting exponents. The connected logarithm then works with a
common weighted estimate for radial and angular derivatives.

## Two orientation-independent endpoint budgets

Use the
[oriented entropy envelopes](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md)
\(Q_\pm\) with all shared rows retained and private curvature
plus \(\tau\) reserved. Write \(q=\max(Q_+(0),Q_-(0))\),
\(W=1+V\), and \(L=3V/4\). Here \(V\) is the sum of \(d^2/\ell\)
over shared rows with gap \(d>0\) and smaller exponent \(\ell\),
plus private curvature and \(\tau\).
Let \(L_0\) be the maximum of \(1\) and the shared ratios \(d/\ell\).
The envelopes give, for \(w=k/m\le1/2\) and either half,
\[
 \begin{split}
 c_{m,k}^{1/k}&\le (e/w)Q_\pm(w),&
 c_{m,k}&\le(mq)^k e^{L k^2/m}/k!,\\
 Q_\pm(w)&\le q e^{Lw},&
 Q_\pm(1/2)&\le e^{-S/4}.
 \end{split}                                                    \tag{1}
\]
The maximum of \(eQ_\pm(w)/w\) is at an endpoint on each
closed subinterval of \(0<w\le1/2\).
The [endpoint cost absorption](MULTIVARIATE_ABSOLUTE_CURVATURE_20.md#the-endpoint-cost-absorbs-the-excess-curvature)
proves \(q\le e^{-S/2}\) and \(Wq\le(1+S)e^{-S/2}\)
also at \(S\ge15\). At \(C=1+S=16\) the required
comparison polynomial is \(r^3-4r^2+6r+10>0\), since its
derivative is positive on the positive half-line.

We improve the other endpoint budget to \(L_0q\le2e^{-S/2}\).
For \(L_0>1\), take a shared row realizing \(r=L_0\).
As in the curvature sixteen proof its excess cost in either
orientation is at least
\(\ell E(r)\) with \(E(r)=r/2+r/[2(1+r)]-\log(1+r)\ge0\).
Thus its excess absorbs the factor
\[
 r e^{-\ell E(r)}\le r(1+r)e^{-r/2-r/[2(1+r)]}<2 .
\]
The logarithmic derivative of the right side has the sign
of \(2+4r+2r^2-r^3\), with exactly one positive zero
between \(10/3\) and \(17/5\). At the maximum, bound
\(r(1+r)\) by its value at \(17/5\) and the exponential
by its value at \(10/3\). The exact scalar check gives the
strict inequality. For \(L_0=1\) use \(q\le e^{-S/2}\).

## Exact derivative at a positive start

At real \(m\ge2\), in any orientation with positive starting
exponents, let \(B_m\) be the inclusive first coefficient
using Gamma interpolation. The binomial factor contributes
\(1/m\) to \((\log B_m)'\). For a row with starting exponent
\(a>0\), slope \(d=b-a\), and \(w=1/m\), let
\(L_m(w)=m^{-1}[\log\Gamma(1+m(a+dw))]_{\rm chord}\)
use the normalized endpoint chord. Its exact contribution
to the derivative is
\[
                  L_m-wL_m'+m\partial_m L_m .
\]
The [Gamma integral](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md#a-gamma-chord-bound-in-both-orientations)
gives
\(L_m''=d^2\int_0^\infty e^{-(a+dw)y}h(y/m)\,dy\)
for \(h(t)=t/(e^t-1)\). We have \(0\le-h'(t)\le1/2\)
as in the [split-row derivative](MULTIVARIATE_RATIO_BALANCED_ANGULAR_SEPARATION.md#splitting-and-the-linked-gamma-derivative).
Thus both terms \(L_m-wL_m'\) and \(m\partial_m L_m\)
are nonpositive. For shared rows the magnitude of their
sum is at most \(3d^2/(4m^2\ell)\) by the
same result, whose derivative argument does not use
ratio balance.

For a private ending row \((a,0)\), the first term
\(-\int_0^w uL_m''(u)\,du\) has magnitude at most
\(a/[2m^2(1-1/m)]\). For the second term, bound
\[
 \partial_m L_m''(u)\le
                    \min\!\left(2a^2,\frac1{2m^2(1-u)^2}\right).
\]
The first bound follows by removing the exponential and
using \(\int_0^\infty t(-h'(t))\,dt=\int_0^\infty h(t)\,dt
=\sum_{j\ge1}j^{-2}<2\). Integration by parts and the positive
exponential series for \(h\) justify the identities.
The second uses \(-h'\le1/2\).
The chord Green kernel is bounded above by \(w(1-u)\).
Integrating with a split at \(1-u=1/(2ma)\) gives
\[
     m|\partial_m L_m(w)|
                   \le\frac{1/4+\tfrac12\log(2ma)}{m^2}.
\]
Therefore, with \(S_p\) the sum of the private exponents
in this orientation and including the Gaussian term,
\[
  0\le 1/m-(\log B_m)'
       \le \frac{(3/4)V+S_p\bigl(1/[2(1-1/m)]+1/4+
                                      \tfrac12\log(32m)\bigr)}{m^2}.
                                                               \tag{2}
\]
Here \(a\le S<16\) on private rows and their number is
at most \(S_p\). The estimate applies independently in
each orientation with positive start.

## The short ranges

Write \(N=n+1\), \(m=n,N\), and \(D=Nq\).
For \(N\le256\) use the discrete curvature estimate
\(c_{m,k}\le\binom mk e^{-(15/2)k(1-k/m)}\).
For \(k=1,2,3\) use the \(k\)-th power of
\(m e^{-(15/2)(1-k/m)}\) divided by \(k!\).
Its maximum on \([2k,256]\) is at an endpoint by
logarithmic differentiation. The three bounds are
\(147/1000,\,3/250,\,7/10000\). If \(m\le7\), the
first bound improves to \(6/125\), and the indices
\(k\ge4\) do not occur. For \(k\ge4\), use
\((e/w)e^{-(15/2)(1-w)}\) on \([4/m,1/2]\),
and at \(m=8,256\) for the first endpoint.
This gives \(c_{m,k}\le(13/100)^k\) in the half.

If \(N>256\) and \(D\le3/20\), \(V\le S(S+2)<288\) gives
\(L/n<17/20\) and the endpoint budget gives \(LD/n<67/10000\).
For each \(1\le k\le4\), maximize
\(D\exp[k\min(17/20,(67/10000)/D)]\)
at the switch or at \(3/20\). Above the switch its
logarithmic derivative has at most one sign change,
from negative to positive. The first three tests in (1)
give \(157/1000,\,7/500,\,9/10000\).
Use the fourth test and the midpoint test in (1)
for the geometric cap \((17/100)^k\) for \(k\ge4\).

In either short range use \(R_m=e^{2/(5m)}\).
The displayed caps bound each half's
\(|H_m-1|\) by \(A\) and \(|zH_m'|\) by \(B\) at
this radius, with
\[
          (1+A)/(1-A)<e^{2/5},\qquad B/(1-A)<1/4 .
\]
The scalar checks sum the geometric tails and bound
\(R_m\) at the smallest \(m\) in each range.
At \(m=1\) both halves equal one. Factoring
\(P_m(z)=L_m(z)+z^m U_m(1/z)\) with half weight on the
middle coefficient, Rouché in both orientations confines all
roots to \(|\log|z||<2/(5m)\).
The angular derivative of each of the four half logs
at the actual radii is below \(1/4\) in modulus.
Integrating from the positive radius gives a difference
of the two full phase corrections strictly below \(\theta\)
in modulus on \(0<\theta<\pi\).
Their common-angle equations would force the phase
difference \(\theta\) minus this correction to be an
integer multiple of \(2\pi\), which is impossible.
At positive arguments positivity excludes roots; at the
negative argument the half logs are real and only an odd
order can have a root. Reflection handles the lower half-plane.

## The long ranges

Now \(D>3/20\). From \(q\le e^{-15/2}\), \(n\ge271\).
The budgets give
\[
 D/n<556/10^6,\qquad WD/n<89/10000,\qquad
 L_0D/n<1112/10^6,\qquad K=\lceil20D\rceil+10\le94D .
\]
Thus \(K/n,L_0K/n<\kappa=21/200\) and \(K+1<n/2\).
The \(k=K\) test \((e/20)e^{(3/4)WK/n}<3/10\)
and the midpoint test give \(c_{m,k}\le(3/10)^k\)
from \(K\) to the midpoint in either half.

For a positive start (2) gives \(|(\log B_m)'|\le1/n\)
on \(n\le m\le N\). In fact \(V/m<(89/10000)/(3/20)\),
\(S_p<16\), and
\(\log(32m)/m\le\log(32\cdot271)/271<10/271\).
The rational bounds imply \(1/m-(\log B_m)'<1/m\).
Also \(B_m/m\) is decreasing. Apply (1) at the integer
order \(n\) to get \(B_m\le D e^{L/n}<(523/500)D\).

The [sharp Gamma bounds](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md#the-long-range-with-sharper-gamma-derivative-estimates)
apply with \(\alpha=L_0/[(1-\kappa)n]\), yielding
\(u=\alpha K<59/500\),
\(H_2<(119/100)W/n<71/1000\), and finite-difference
bounds \(H_r=(r-2)!\alpha^{r-2}H_2\) with order
derivative at most \(2rH_r/n\).
Use the [charged connected expansion](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_TWELVE.md#long-halves)
with \(\beta=e^{H_2+1/8}H_2\).
In its discovery-tree estimate we have
\(L_0\le\sqrt V\) since a shared row with ratio \(r\)
contributes \(\ell r^2\ge r^2\) to \(V\).
As \(V\ge15\),
\(\alpha/\beta\le e^{-1/8}\sqrt{15}/16<1/4\).
The tree base is at most \((5/4)^5\).
With \(\beta K<61/50\) and \(\beta B_m<27/2000\),
the connected parameter \(p_m=\Xi B_m\) is below \(1/10\).
The derivative constant is \(1+16/e<7\).
Thus the formal logarithm through \(K\), denoted by \(A_m\),
has coefficient bounds \(B_m p_m^{k-1}/k\)
and derivative bounds \(7B_m p_m^{k-1}/n\).

The absolute coefficient norms of \(A_m\) on radii
\(2,6/5\) are below \((9/4)B_m,(13/10)B_m\), respectively,
and \(\|A_N-A_n\|_{6/5}<10D/n\).
The exponential and the true coefficient tails beginning at
\(K+1\), relative to \(e^{A_m}\) in absolute coefficient
norm on \(6/5\), are bounded by
\[
     (3/5)^{K+1}e^{(71/20)B_m}<e^{-32/5},\qquad
     (9/25)^{K+1}e^{(13/10)B_m}/(1-9/25)<e^{-9K/10}.
\]
Here \(K+1\ge20D+11\), \(D>3/20\), and \(K\ge14\).
Their sum is below \(1/500\).
It follows that \(H_m\) has an analytic logarithm
on \(6/5\) with \(\|\log H_m-A_m\|_{6/5}<3/1000\).

For an orientation with a private starting row, the
[private-start bound](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md#uniform-comparison-across-degrees)
gives \(c_{m,k}^{1/k}\le e^{3/2}Q_\pm(k/m)\).
On radius \(6/5\) the ratio through \(K\) is then at
most \((6/5)e^{3/2+(3/4)WK/n}q<3/500\).
Beyond \(K\) use \(9/25\) as above.
The sum gives \(\|\log H_m\|_{6/5}<7/1000\).

Uniformly, the two orientations therefore have
\[
     \|\log H_m\|_{6/5}\le M=(7/5)D,\qquad
     \|\log H_N-\log H_n\|_{6/5}<\delta=3/200 .
\]
Use the [weighted radial and phase comparison](MULTIVARIATE_RATIO_BALANCED_ANGULAR_SEPARATION.md#paired-phases-at-their-actual-radii).
It requires only the given analytic half logs and norm estimates.
Rouché in both orientations gives
\(|u_m|<u_0=\log(101/100)\) since \(2M/n<u_0\).
On this radial range the individual first and mixed
log derivatives are bounded by \(P=(9/4)M\), \(C=20M\).
At a common argument the real equations give
\(n|u_N-u_n|\le(u_0+2\delta)/(1-2P/n)<1/20\).
Integrating the paired angular derivative and the
mixed derivatives gives a phase difference from
\(\theta\) bounded by \(((9/2)\delta+2C/(20n))\theta<\theta\).
This excludes every common argument and completes the proof.

## Exact scalar checks

~~~python
from fractions import Fraction as F
from math import factorial

def lo(x):
    x=F(x)
    return sum(x**j/factorial(j) for j in range(65))

def hi(x):
    x=F(x)
    return lo(x)+x**65/factorial(65)/(1-x/66)

def ck(label,a,b):
    assert a<b, (label,float(a),float(b))

S=F(15); sigma=S/2; U=1/lo(sigma)
en=lambda m,k:m/lo(sigma*(1-F(k,m)))
for name,k,ms,cap in [('first',1,[2,7],F(6,125)),('first',1,[2,256],F(147,1000)),
                      ('second',2,[4,256],F(3,250)),('third',3,[6,256],F(7,10000))]:
    ck('finite short '+name,max(en(m,k) for m in ms)**k/factorial(k),cap)
ck('finite short tail',hi(1)*max(en(m,4) for m in [8,256])/4,F(13,100))

d=F(3,20); Ac=F(67,10000); cap=F(17,20); d0=Ac/cap
ck('short L/n', F(3,4)*288/256,cap)
ck('short LD/n', F(3,4)*F(257,256)*16*U,Ac)
caps=[max(d0*hi(k*cap),d*hi(k*Ac/d)) for k in range(1,5)]
for k,b in enumerate([F(157,1000),F(7,500),F(9,10000)],1):
    ck('small D first three',caps[k-1]**k/factorial(k),b)
ck('small D tail',hi(1)*caps[3]/4,F(17,100))
ck('small D midpoint',2*hi(1)/lo(S/4),F(17,100))

for m,c1,c2,c3,r in [(2,F(6,125),F(3,250),F(7,10000),0),
                     (8,F(147,1000),F(3,250),F(7,10000),F(13,100)),
                     (256,F(157,1000),F(7,500),F(9,10000),F(17,100))]:
    R=hi(F(2,5*m))
    x=R*r
    A=c1*R+c2*R**2+c3*R**3+x**4/(1-x)
    B=c1*R+2*c2*R**2+3*c3*R**3+x**4*(4/(1-x)+x/(1-x)**2)
    ck('short amplitude',1+A,lo(F(2,5))*(1-A))
    ck('short phase',B/(1-A),F(1,4))

# The largest ratio estimate at the unique critical point.
r0=F(10,3); r1=F(17,5)
pol=lambda r:2+4*r+2*r*r-r**3
assert pol(r0)>0>pol(r1)
ck('L0q',r1*(1+r1)/lo(r0/2+r0/(2*(1+r0))),2)

q=F(556,10**6); Wq=F(89,10000); Kcap=94
kap=F(21,200); u=F(59,500)
ck('long n', 271,d*lo(sigma))
ck('long q',F(272,271)*U,q)
ck('long Wq',16*F(272,271)*U,Wq)
ck('long K cap',20+11/d,Kcap)
ck('long kappa',2*q*Kcap,kap)
ck('long range',q*Kcap+F(1,271),F(1,2))
ck('long u',kap/(1-kap),u)
ck('long derivative private log',32*271,lo(10))
ck('long derivative',F(3,4)*Wq/d+F(96,271),1)
ck('long B',hi(F(3,4)*Wq/d),F(523,500))
ck('long H2 prefactor',(1+q*Kcap/(1-kap))/(1-kap),F(119,100))
ck('long H2',F(119,100)*Wq/d,F(71,1000))
ck('long tail',hi(1+F(3,4)*Kcap*Wq)/20,F(3,10))
ck('long midpoint',2*hi(1)/lo(S/4),F(3,10))
beta=F(119,100)*hi(F(71,1000)+F(1,8))
ck('long betaK',beta*Kcap*Wq,F(61,50))
ck('long betaB',beta*Wq*F(523,500),F(27,2000))
ck('long tree t',15*16/lo(F(1,4)),16**2)
p=F(1,10)
ck('long p',F(5,4)**5*F(27,2000)/(1-u)*hi(F(61,50)*(F(1,2)+u/(6*(1-u)))),p)
ck('long derivative constant',1+16/lo(1),7)
ck('long outer norm',1/lo(F(9,4)*p),1-2*p)
ck('long inner norm',1/lo(F(13,10)*p),1-F(6,5)*p)
ck('long adjacent norm',7*F(523,500)*F(6,5)/(1-F(6,5)*p),10)
ck('long exp tail',F(32,5),F(11,2)+(10-F(71,20)*F(523,500))*d)
ck('long log radius',hi(F(1,2)),F(5,3))
real=F(9,25)
ck('long true tail slope',hi(F(9,10)+F(13,10)*F(523,500)/20),1/real)
ck('long true tail intercept',real/(1-real),1)
ck('long sum of tails',1/lo(F(32,5))+1/lo(F(9,10)*14),F(1,500))
ck('long log correction',1/lo(F(3,1000)),1-F(1,500))
ck('long positive δ',10*q+2*F(3,1000),F(3,200))
private=F(6,5)*hi(F(3,2)+F(3,4)*Kcap*Wq)*U
ck('long private ratio',private,F(3,500))
ck('long private norm',F(3,500)/(1-F(3,500))+real**15/(1-real),1-1/lo(F(7,1000)))
ck('long private δ',2*F(7,1000),F(3,200))
ck('uniform M for positive',F(13,10)*F(523,500)+F(3,1000)/d,F(7,5))
ck('uniform M for private',F(7,1000)/d,F(7,5))
M=F(7,5)*q
umin=F(199,20000); umax=F(1,100)
ck('long radii',2*M,umin)
ck('slope loss log',hi(F(17,100)),F(120,101))
ck('slope one',1/(lo(1)*F(17,100)),F(9,4))
ck('slope two',4/(lo(1)*F(17,100))**2,20)
ck('long shift',(umax+2*F(3,200))/(1-2*F(9,4)*M),F(1,20))
ck('long phase',2*F(9,4)*F(3,200)+2*20*M/20,1)
print('All unrestricted S=15 scalar checks passed.')
~~~
