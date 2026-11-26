# All-degree angular separation at absolute curvature fourteen

Let \(a,b\) be distinct ordinary nonnegative exponent columns,
with equal rows omitted. Define
\[
             S=\tau+\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)}.
\]
The optional common Gaussian strength \(\tau\ge0\) contributes
\(e^{-\tau k(m-k)/(2m)}\) to the normalized coefficient of
index \(k\) in degree \(m\).

**Theorem.** If \(S\ge14\), consecutive moment polynomials have
disjoint root arguments in every degree, for arbitrary private
rows and orientations. At \(\tau=0\) this proves ordinary SFC(2)
over every characteristic-zero field.

The [curvature fifteen theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_15.md)
covers \(S\ge15\). On \(14\le S<15\), the proof uses an
orientation-specific square budget for short coefficient halves.
The budget distinguishes rows whose exponent increases in the
chosen half from rows whose exponent decreases. The latter
contribute only their absolute curvature to the directed
quantity. Longer halves use the exact linked logarithm comparison.

## The square budget in each orientation

Use the
[oriented entropy envelopes](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md)
\(Q_\pm\) with all shared rows retained. Write \(q_\pm=Q_\pm(0)\)
and \(q=\max(q_+,q_-)\). Keep \(V,W=1+V,L_0,L=3V/4\)
from the [curvature fifteen endpoint estimates](MULTIVARIATE_ABSOLUTE_CURVATURE_15.md#two-orientation-independent-endpoint-budgets).
Their coefficient estimates and endpoint maximum property apply
here with the same proof. We have
\[
 q\le e^{-S/2},\qquad L_0q\le2e^{-S/2},\qquad
 Wq\le(1+S)e^{-S/2}.
                                                               \tag{1}
\]
For the last inequality the comparison polynomial at
\(C=1+S=15\) is \(r^3-4r^2+5r+9=r[(r-2)^2+1]+9>0\),
as required in the
[endpoint excess-curvature comparison](MULTIVARIATE_ABSOLUTE_CURVATURE_20.md#the-endpoint-cost-absorbs-the-excess-curvature).
The proof of \(L_0q\le2e^{-S/2}\) at curvature fifteen
uses only the one-row cost and applies at this threshold.

In one orientation, let \(s_0\) be its total private
curvature plus \(\tau\) and set
\[
          V_{\pm}=s_0+\sum_{\text{shared}}\frac{d_i^2}{a_i},
          \qquad W_{\pm}=1+V_{\pm},
\]
where \(a_i\) is the positive starting exponent in that
orientation and \(d_i\) is the signed difference. Then
\[
             W_\pm^2 q_\pm\le (1+S)^2e^{-S/2}.                  \tag{2}
\]
To prove this, let \(\ell\) be the smaller shared exponent
and \(r=|d|/\ell\). If the row decreases, its contribution
to \(V_\pm\) is exactly its absolute curvature \(S_i\),
and its chord cost is at least \(S_i/2\).
If it increases, the excess in \(V_\pm\) over \(S_i\) is
\(\ell g(r)\), where \(g(r)=r^3/(1+r)\).
Its endpoint log cost above \(S_i/2\) is \(\ell E_+(r)\),
where
\[
          E_+(r)=(1+r)\log(1+r)-r-\frac{r^2}{2(1+r)} .
\]
Using \(\log(1+r)\ge2r/(r+2)\) gives
\[
    E_+'(r)\ge\frac{r^2(3r+4)}{2(r+2)(1+r)^2}
                        \ge \frac{2g'(r)}{C+g(r)},\qquad C=1+S .
\]
For the second inequality it suffices to take \(C=15\).
After clearing positive denominators and cancelling \(r^2\),
it follows from
\(3r^4-4r^3+9r^2+53r+36>0\), since
\(3r^2-4r+9>0\).
Integration from zero gives
\(E_+(r)\ge2\log(1+g(r)/C)\).
For integer \(\ell\ge1\),
\((1+g(r)/C)^\ell\ge1+\ell g(r)/C\).
Multiplication over increasing rows therefore gives
\(\exp(\sum\ell E_+(r)/2)\ge1+(V_\pm-S)/C=W_\pm/C\).
Private rows have exactly their baseline contribution in
\(Q_\pm\), proving (2).

For \(0\le w\le1/40\) there is a matching directed increment
\[
              Q_\pm(w)\le q_\pm e^{\lambda V_\pm w},
                         \qquad \lambda=77/100 .                \tag{3}
\]
Indeed the \(A(w)/w\) derivative in the
[Gamma chord envelope](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md#a-gamma-chord-bound-in-both-orientations)
is \(d^2\int_0^1 u/(a+dwu)\,du\).
The denominator is at least \((39/40)a\), giving the bound
\((20/39)d^2/a\).
The concave log-chord tangent \(G(w)/2\) contributes at
most \(d^2w/(4a\ell)\le d^2w/(4a)\) to the increment.
Since \(20/39+1/4<\lambda\) and the residual increment is
\(s_0w/2\), this proves (3).

## Short halves

Write \(N=n+1\), \(m=n,N\), and \(D=Nq\).
For \(N\le160\), the discrete curvature bound is
\(c_{m,k}\le\binom mk e^{-7k(1-k/m)}\).
For \(k=1,2,3\) use
\([m e^{-7(1-k/m)}]^k/k!\), maximized in \(m\) at either
\(2k\) or \(160\). This follows from the sign change
of its logarithmic derivative from negative to positive.
The resulting caps are \(77/500,13/1000,1/1000\).
If \(m\le7\), the first cap improves to \(61/1000\),
and \(k\ge4\) does not occur in a half.
Use \((e/w)e^{-7(1-w)}\) at \(w=4/m,1/2\) and
\(m=8,160\) for the first endpoint. The same derivative
test gives the cap \((17/100)^k\) from \(k=4\)
in both halves.

For \(N>160\) and \(D\le3/20\), let \(D_\pm=N q_\pm\le D\).
Here \(n\ge160\) and \(V\le S(S+2)<255\).
The budgets (2)–(3) give
\[
       \lambda V_\pm/n
              \le \min\!\left(123/100,\frac{7}{250\sqrt{D_\pm}}\right).
\]
For the second bound use
\((\lambda V_\pm/n)^2D_\pm
   \le\lambda^2 (15)^2 e^{-7}(161/160)/160<(7/250)^2\),
and for the first use \(\lambda 255/160<123/100\).
For \(1\le k\le4\), the maximum of
\(D_\pm e^{k\lambda V_\pm/n}\) is bounded by the tests
at the switch of this upper envelope and at \(3/20\).
On the upper piece the logarithmic derivative changes
sign at most once from negative to positive.
The first three coefficient caps are \(81/500,2/125,11/10000\),
using (3) and the binomial bound by \(m^k/k!\).
The fourth uses \(eQ_\pm(w)/w\) at \(w=4/m\), with
(3) valid at this point. The midpoint bound
\(2e\, Q_\pm(1/2)\le2e e^{-7/2}<17/100\) then gives the
geometric cap \((17/100)^k\) for every \(k\ge4\).

For all of these short halves take \(R_m=e^{2/(5m)}\).
The coefficient caps give \(|H_m-1|\le A\) and \(|zH_m'|\le B\)
on this radius, with \((1+A)/(1-A)<e^{2/5}\)
and \(B/(1-A)<1/4\), as checked below.
At \(m=1\) the halves both equal one.
The [short radius and phase argument](MULTIVARIATE_ABSOLUTE_CURVATURE_15.md#the-short-ranges)
then excludes a common argument. In detail the amplitude bound
confines the roots to \(|\log|z||<2/(5m)\) by Rouché in
both orientations. The four angular slopes at the actual
radii are each below \(1/4\) in modulus. Integration from the
positive radius bounds the difference of the two full phase
corrections by strictly less than the candidate angle
\(\theta\) on \(0<\theta<\pi\). Positive coefficients exclude a
positive argument, and real half logs at the negative
argument permit only odd order. The lower half-plane follows
by reflection.

## Long halves

For \(D>3/20\) we have \(n\ge164\).
The bounds (1) give
\[
  D/n<918/10^6,\qquad WD/n<1377/10^5,\qquad
  L_0D/n<1836/10^6,\qquad K=\lceil16D\rceil+10\le90D .
\]
Both \(L_0 K/n,K/n\) are below \(\kappa=83/500\) and
\(K+1<n/2\).
Use the first endpoint test \((e/16)e^{(3/4)WK/n}<9/20\)
and the midpoint test to give \(c_{m,k}\le(9/20)^k\)
for \(k\ge K\).

For an orientation with positive start, let \(B_m\) be its
first coefficient at real \(m\in[n,N]\).
The [exact first-factor derivative including private endings](MULTIVARIATE_ABSOLUTE_CURVATURE_15.md#exact-derivative-at-a-positive-start)
gives \(0\le1/m-(\log B_m)'<1/m\).
In its bound use \(V/m<(1377/10^5)/(3/20)\), total private
curvature below \(15\), and
\(\log(2ma)/m\le\log(30m)/m\le\log(30\cdot164)/164<9/164\)
for private exponent \(a\). The estimate for the sum of
the private terms after multiplication by \(m\) is at most
\((165/2)/164\). Thus \(|(\log B_m)'|\le1/n\) and
\(B_m/m\) decreases. Use the coefficient envelope at
the integer order \(n\) to obtain \(B_m\le De^{L/n}<(134/125)D\).

The [sharp Gamma estimates](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md#the-long-range-with-sharper-gamma-derivative-estimates)
with \(\alpha=L_0/[(1-\kappa)n]\) give
\(u=\alpha K<1/5\), \(H_2<(33/25)W/n<61/500\), and
finite-difference bounds \(H_r=(r-2)!\alpha^{r-2}H_2\)
with order derivatives at most \(2rH_r/n\).
In the [charged connected expansion](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_TWELVE.md#long-halves)
put \(\beta=e^{H_2+1/8}H_2\).
As \(L_0\le\sqrt V\) and \(V\ge14\),
\(\alpha/\beta\le e^{-1/8}\sqrt{14}/15<1/4\).
The tree base is at most \((5/4)^5\).
With \(\beta K<21/10\), \(\beta B_m<251/10000\),
the connected parameter \(p_m=\Xi B_m\) is below \(31/100\).
The formal logarithm \(A_m\) through degree \(K\)
has coefficient bounds \(B_m p_m^{k-1}/k\) and derivative
bounds \(7B_m p_m^{k-1}/n\), using \(1+16/e<7\).
The norms at radii \(2,6/5\) are below
\((16/5)B_m,(31/20)B_m\), respectively, and
\(\|A_N-A_n\|_{6/5}<15D/n\).

Relative to \(e^{A_m}\) on \(6/5\), the exponential and
actual tails have absolute coefficient norm bounds
\[
 (3/5)^{K+1} e^{(19/4)B_m}<e^{-59/10},\qquad
 \frac{(27/50)^{K+1}}{1-27/50}e^{(31/20)B_m}<e^{-K/2}.
\]
For the first use \(K+1\ge16D+11\) and
\(\log(5/3)>1/2\). For the second retain \(D\le(K-10)/16\).
Since \(K\ge13\), the sum is below \(9/2000\).
Thus the half logarithm differs from \(A_m\) by less than
\(1/200\) on \(6/5\), in absolute coefficient norm.

If an orientation has a private starting row, reserve one
factorial as in the [private-start bound](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md#uniform-comparison-across-degrees).
Through \(K\) the ratio on \(6/5\) is then at most
\((6/5)e^{3/2+(3/4)WK/n}q<13/1000\).
Beyond \(K\) use \(27/50\).
Summing both tails gives \(\|\log H_m\|_{6/5}<7/500\).

In both orientations
\(\|\log H_m\|_{6/5}\le M=(17/10)D\) and
\(\|\log H_N-\log H_n\|_{6/5}<\delta=3/100\).
Use the [weighted radial and phase comparison](MULTIVARIATE_RATIO_BALANCED_ANGULAR_SEPARATION.md#paired-phases-at-their-actual-radii).
This comparison only requires these two estimates on the
analytic half logarithms. With \(u_0=\log(101/100)\) we have
\(2M/n<u_0\), so Rouché in both orientations confines the
root log radii to \(|u_m|<u_0\). The individual first and
mixed log derivatives there are bounded by
\(P=(9/4)M\) and \(C_1=20M\).
If a common argument existed, the real equations would
give \(n|u_N-u_n|\le(u_0+2\delta)/(1-2P/n)<2/25\).
The difference of the two full angular corrections would
be at most \(((9/2)\delta+2C_1(2/25)/n)\theta<\theta\).
This gives the same phase contradiction as in the short
ranges and completes the proof.

## Exact scalar checks

The checks below use only rational exponential bounds and
the displayed constants.

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

S=F(14); sigma=S/2; U=1/lo(sigma)
en=lambda m,k:m/lo(sigma*(1-F(k,m)))
for name,k,ms,cap in [('first',1,[2,7],F(61,1000)),('first',1,[2,160],F(77,500)),
                      ('second',2,[4,160],F(13,1000)),('third',3,[6,160],F(1,1000))]:
    ck('finite short '+name,max(en(m,k) for m in ms)**k/factorial(k),cap)
ck('finite short tail',hi(1)*max(en(m,4) for m in [8,160])/4,F(17,100))

d=F(3,20); A=F(7,250); cap=F(123,100); d0=(A/cap)**2
lam=F(77,100); beta_d=F(73,1000)
ck('short L/n',lam*255/160,cap)
ck('short mass comparison',lam**2*225*U*F(161,160)/160,A**2)
ck('short slope',F(20,39)+F(1,4),lam)
ck('small D radical bound',A*A,d*beta_d**2)
caps=[max(d0*hi(k*cap),d*hi(k*beta_d)) for k in range(1,5)]
for k,b in enumerate([F(81,500),F(2,125),F(11,10000)],1):
    ck('small D first three',caps[k-1]**k/factorial(k),b)
ck('small D tail',hi(1)*caps[3]/4,F(17,100))
ck('small D midpoint',2*hi(1)/lo(S/4),F(17,100))
for m,c1,c2,c3,r in [(2,F(61,1000),F(13,1000),F(1,1000),0),
                     (8,F(77,500),F(13,1000),F(1,1000),F(17,100)),
                     (160,F(81,500),F(2,125),F(11,10000),F(17,100))]:
    R=hi(F(2,5*m))
    x=R*r
    A=c1*R+c2*R**2+c3*R**3+x**4/(1-x)
    B=c1*R+2*c2*R**2+3*c3*R**3+x**4*(4/(1-x)+x/(1-x)**2)
    ck('short amplitude',1+A,lo(F(2,5))*(1-A))
    ck('short phase',B/(1-A),F(1,4))

q=F(918,10**6); Wq=F(1377,100000); Kcap=90
kap=F(83,500); u=F(1,5)
ck('long n',164,d*lo(sigma))
ck('long q',F(165,164)*U,q)
ck('long Wq',15*F(165,164)*U,Wq)
ck('long K cap',16+11/d,Kcap)
ck('long kappa',2*q*Kcap,kap)
ck('long range',q*Kcap+F(1,164),F(1,2))
ck('long u',kap/(1-kap),u)
ck('long derivative private log',30*164,lo(9))
ck('long derivative',F(3,4)*Wq/d+F(165,2*164),1)
b=F(134,125)
ck('long B',hi(F(3,4)*Wq/d),b)
ck('long H2 prefactor',(1+q*Kcap/(1-kap))/(1-kap),F(33,25))
ck('long H2',F(33,25)*Wq/d,F(61,500))
ck('long tail',hi(1+F(3,4)*Kcap*Wq)/16,F(9,20))
ck('long midpoint',2*hi(1)/lo(S/4),F(9,20))
beta=F(33,25)*hi(F(61,500)+F(1,8))
ck('long betaK',beta*Kcap*Wq,F(21,10))
ck('long betaB',beta*Wq*b,F(251,10000))
ck('long tree t',14*16/lo(F(1,4)),15**2)
p=F(31,100)
ck('long p',F(5,4)**5*F(251,10000)/(1-u)*hi(F(21,10)*(F(1,2)+u/(6*(1-u)))),p)
R=F(6,5); Aout=F(16,5); Ain=F(31,20)
ck('long outer norm',1/lo(Aout*p),1-2*p)
ck('long inner norm',1/lo(Ain*p),1-R*p)
ck('long adjacent norm',7*b*R/(1-R*p),15)
ck('long exp tail',F(59,10),F(11,2)+(8-(Aout+Ain)*b)*d)
ck('long log radius',hi(F(1,2)),F(5,3))
real=F(27,50)
ck('long true tail slope',hi(F(1,2)+Ain*b/16),1/real)
ck('long true tail intercept',real/(1-real),lo(Ain*b*10/16))
ck('long sum of tails',1/lo(F(59,10))+1/lo(F(13,2)),F(9,2000))
ck('long log correction',1/lo(F(1,200)),1-F(9,2000))
ck('long positive delta',15*q+F(1,100),F(3,100))
private=R*hi(F(3,2)+F(3,4)*Kcap*Wq)*U
ck('long private ratio',private,F(13,1000))
ck('long private norm',F(13,1000)/(1-F(13,1000))+real**14/(1-real),1-1/lo(F(7,500)))
ck('long private delta',2*F(7,500),F(3,100))
ck('uniform M for positive',Ain*b+F(1,200)/d,F(17,10))
ck('uniform M for private',F(7,500)/d,F(17,10))
M=F(17,10)*q
umin=F(199,20000); umax=F(1,100)
ck('long radii',2*M,umin)
ck('slope loss log',hi(F(17,100)),F(120,101))
ck('slope one',1/(lo(1)*F(17,100)),F(9,4))
ck('slope two',4/(lo(1)*F(17,100))**2,20)
ck('long shift',(umax+2*F(3,100))/(1-2*F(9,4)*M),F(2,25))
ck('long phase',2*F(9,4)*F(3,100)+2*20*M*F(2,25),1)
print('All unrestricted S=14 scalar checks passed.')
~~~
