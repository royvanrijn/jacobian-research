# All-degree angular separation at absolute curvature thirteen

Let \(a,b\) be distinct ordinary nonnegative exponent columns,
with equal rows omitted. Define
\[
             S=\tau+\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)} .
\]
Here a common strength \(\tau\ge0\) contributes the Gaussian
factor \(e^{-\tau k(m-k)/(2m)}\) to the normalized coefficient
of index \(k\) in degree \(m\).

**Theorem.** If \(S\ge13\), consecutive moment polynomials have
disjoint root arguments in every degree, for arbitrary private
rows and orientations. At \(\tau=0\) they satisfy ordinary SFC(2)
over every characteristic-zero field.

The [curvature fourteen theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_14.md)
covers \(S\ge14\). For \(13\le S<14\), each orientation uses its
own endpoint mass to decide how to estimate its coefficient
half. The small-mass estimate uses a cubic directed curvature
budget. The large-mass estimate uses the stronger endpoint
charge paid by a large ratio in an increasing row.
The two orientations may fall in different mass ranges.

## Directed endpoint estimates

Retain every shared row in each
[oriented entropy envelope](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md).
Fix one orientation and write its envelope as \(Q\), with
\(q=Q(0)\), \(V_d=s_0+\sum_{\rm shared}d_i^2/a_i\), \(W_d=1+V_d\).
Here \(a_i>0\) is the oriented starting exponent on a shared
row, \(d_i\) its signed gap, and \(s_0\) is the total private
curvature plus \(\tau\). Also put
\(L_d=\max(1,\max_{\rm shared}|d_i|/a_i)\).
Use \(V\) for the orientation-independent curvature budget
with the smaller exponent in each shared denominator.
We have \(S\le V_d\le V\le S(S+2)<224\) in the range of interest
and \(L_d\le\sqrt{V_d}\).
The [directed increment proof](MULTIVARIATE_ABSOLUTE_CURVATURE_14.md#the-square-budget-in-each-orientation)
gives, for \(0\le w\le w_*\le1/2\),
\[
    Q(w)\le q\exp(\lambda V_d w),\qquad
                   \lambda=\frac1{2(1-w_*)}+\frac14 .           \tag{1}
\]
The first coefficient estimate, the maximum of \(eQ(w)/w\)
on a closed subinterval, and \(Q(1/2)\le e^{-S/4}\) are as
at curvature fourteen.

The new endpoint bounds, valid already at \(S\ge12\), are
\[
 \begin{split}
 q&\le e^{-S/2},&
 W_d^3 q &\le (1+S)^3e^{-S/2},&
 L_d q&\le (9/8)e^{-S/2}. 
 \end{split}                                                    \tag{2}
\]
In particular \(W_d q\le(1+S)e^{-S/2}\).
To prove the cubic estimate let \(C=1+S\ge13\) and
\(g(r)=r^3/(1+r)\). An increasing row with smaller exponent
\(\ell\) and ratio \(r=|d|/\ell\) has excess curvature
\(V_{d,i}-S_i=\ell g(r)\) and excess cost
\(\ell E_+(r)\), where
\[
           E_+(r)=(1+r)\log(1+r)-r-\frac{r^2}{2(1+r)}.
\]
For \(u=r/(r+2)\), use \(\log(1+r)\ge2(u+u^3/3)\) to get
\[
   E_+'(r)\ge
   \frac{r^2(3r+4)}{2(r+2)(r+1)^2}
                       +\frac{2r^3}{3(r+2)^3}
                \ge \frac{3g'(r)}{C+g(r)} .
\]
The last inequality, at \(C=13\) after positive factors
and \(r^2\) are removed, is
\[
             192+400r+180r^2-81r^3-49r^4+20r^5+13r^6>0.
\]
On \(0\le r\le1\) combine the negative terms with \(180r^2\).
At \(r=1+t\) the coefficient list is
\(675,499,38,183,246,98,13\) in increasing degree, all positive.
The bound increases with \(C\).
Integrating, then using
\((1+g(r)/C)^\ell\ge1+\ell g(r)/C\), gives the budget by
multiplication over increasing rows. A decreasing row has
\(V_{d,i}=S_i\) and pays at least \(S_i/2\).
Private rows and the Gaussian enter \(Q\) at their baseline
cost, so no excess is required of them.

If \(L_d>1\), a shared increasing row attains \(L_d=r\).
Its excess cost absorbs \(r\) up to \(9/8\) since
\(r e^{-\ell E_+(r)}\le r e^{-E_+(r)}<9/8\).
Indeed \(E_+\ge0\) and
\(E_+''(r)=1/(1+r)-1/(1+r)^3>0\).
Thus \(\log r-E_+(r)\) is strictly concave and has its
maximum between \(17/10\) and \(7/4\), by evaluating its
derivative at those endpoints.
For the value use \(r\le7/4\) and
\(E_+(r)\ge E_+(17/10)\), checked by the scalar estimate
below. If \(L_d=1\), use \(q\le e^{-S/2}\).

## Very short halves

Write \(N=n+1\) and \(m=n,N\).
If \(N\le96\), the discrete coefficient bound
is \(\binom mk e^{-(13/2)k(1-k/m)}\).
For \(k=1,2,3\) the caps \(31/200,7/500,11/5000\)
follow from the bounds
\([m e^{-(13/2)(1-k/m)}]^k/k!\) at \(m=2k,96\).
On \(2\le m\le11\) the first cap improves to \(39/500\).
The first \(m\)-derivative sign of the log envelope can
only change from negative to positive, justifying
the endpoint maxima.

For \(k\ge4\) use the
[binomial entropy bound](MULTIVARIATE_ENTROPY_ENDPOINT_PHASE_CRITERION.md#combining-the-binomial-entropy-and-the-row-chord)
\(B(w)=\exp([-w\log w-(1-w)\log(1-w)]/w)\).
The maximum of \(B(w)e^{-(13/2)(1-w)}\) on each
subinterval occurs at an endpoint since its logarithmic
derivative multiplied by \(w\) is increasing.
The same property applies to the bound at \(w=4/m\)
as \(m\) runs from \(8\) to \(96\). At \(w=1/2\) use
\(4e^{-13/4}<4/25\), and at \(w=4/96\) use \(B(w)\le e/w\).
Thus \(c_{m,k}\le(4/25)^k\) for \(k\ge4\).

On \(R_m=e^{2/(5m)}\), the scalar checks give
\(|H_m-1|\le A\), \(|zH_m'|\le B\) in each half, with
\((1+A)/(1-A)<e^{2/5}\) and \(B/(1-A)<1/4\).
Use the cap intervals \(2\le m\le7\), \(8\le m\le11\),
\(12\le m\le96\) to bound \(R_m\).
At \(m=1\) the halves equal one.
The [short radius and phase argument](MULTIVARIATE_ABSOLUTE_CURVATURE_15.md#the-short-ranges)
therefore excludes a common argument throughout this range.

## Direct estimates for one orientation

Assume \(N>96\), and select an orientation with \(D=Nq\le3/10\).
Use the same \(D\) for both degrees.
On \(w\le4/96\) the bound (1) holds with \(\lambda=39/50\).
The cubic budget and \(n\ge96\) give
\[
          \lambda V_d/n\le
           \min\!\left(183/100,\frac{61}{1000D^{1/3}}\right).
                                                               \tag{3}
\]
Indeed its cube times \(D\) is bounded by
\(\lambda^3\,14^3 e^{-13/2}(97/96)/96^2\), and
\(\lambda V_d/n\le\lambda\cdot224/96<183/100\).
The functions \(D\exp(k\lambda V_d/n)\) are therefore
bounded by evaluating the right-hand envelope at its
switch \(D_0=(1/30)^3\) and at \(D=3/10\).
Above the switch the log derivative changes sign at most
once from negative to positive.
For \(k=1,2,3\) this gives \(33/100,33/500,11/1000\)
as coefficient caps. At \(k=4\) use the same test with
\(eQ(w)/w\); its midpoint test is weaker than \(3/10\).
Thus the higher coefficients are capped by \((3/10)^k\).
At radius \(6/5\) the absolute coefficient norm of
\(H_m-1\) is consequently at most \(h<1\) and
\(\|\log H_m\|_{6/5}<4/5\).

For positive starting exponents, interpolate the first
three inclusive coefficients at real \(m\in[n,N]\).
The [first-factor derivative](MULTIVARIATE_ABSOLUTE_CURVATURE_15.md#exact-derivative-at-a-positive-start)
gives \(|(\log B_m)'|\le3/n\). For its private terms
use \(\log(2ma)/m\le\log(28m)/m\le\log(28\cdot96)/96<8/96\)
and total private curvature below \(14\).
For \(f_m(x)=\log(c_{m,x}x!/B_m^x)\),
the Gamma arguments for \(0\le x\le3\) are at least
\(ma_i/2\) and \(m/2\) in the binomial, since
\(3L_d/n<1/2\).
The bound \(|\psi_2(y)|\le2/y^2\) for \(y\ge2\)
gives \(|\partial_m\partial_x^2 f_m|\le16W_d/n^2\).
The Gaussian term is included in \(W_d\).
As \(\partial_m f_m=0\) at \(0,1\), this yields
\[
   |\partial_m c_{m,k}|\le\frac{c_{m,k}}n
                 \left(3k+8(W_d/n)k(k-1)\right),\qquad k\le3.
                                                               \tag{4}
\]

For real \(m\), the real-degree Gamma chord bound on shared
rows still gives (3) for their contribution to the
coefficient envelope. For a private ending row \((a,0)\)
the Gamma curvature as a function of the index \(x\) on
\([0,m]\) is at least
\[
     a^2\psi_1(1+am)\ge a/m-1/(2m^2),
\]
by the trigamma integral with \(h(t)\ge1-t/2\).
Thus the integer-order coefficient envelope, based on
\(a/m\), needs at most \(e^{k/(4m)}\) per private row for
this interpolation. There are fewer than \(14\) private
rows, so multiply the real-degree bounds by \(e^{7k/192}\).
The coefficient estimate for the binomial continues to
hold because \(\binom mk\le m^k/k!\) for real \(m\ge k\).

In (4), retain the product \(W_d D/n<213/10000\) in the
second term and estimate \(c_{m,k}/D\) using (3).
The endpoint tests for \(D^{k-1}\exp(k^2\lambda V_d/n)\)
are again at \(D_0,3/10\) by the log derivative sign.
For the adjacent coefficient difference beyond the
first three, each pair lies in the same interval
\([0,(3/10)^k]\), so a single geometric tail suffices.
Dividing by \(1-h\) along the segment of the two halves
therefore gives \(\|\log H_N-\log H_n\|_{6/5}<3/25\)
by the scalar sums below. The inverse series on the
segment gives this bound in absolute coefficient norm.

For a private starting row the
[private reservation bound](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md#uniform-comparison-across-degrees)
gives \(c_{m,k}^{1/k}\le e^{3/2}Q(w)\).
Use (3) for \(k\le3\) to obtain the ratio \(23/1000\)
on radius \(6/5\); for the remainder \(Q(w)\le Q(1/2)\)
gives the ratio \(21/100\).
The scalar sums give \(\|\log H_m\|_{6/5}<4/125\),
which also gives the required adjacent bound.

## Connected estimates for one orientation

For \(D=Nq>3/10\), we have \(n\ge199\).
The directed budgets give
\[
   D/n<1512/10^6,\qquad W_dD/n<2118/10^5,\qquad
   L_dD/n<171/10^5,\qquad K=\lceil16D\rceil+8\le46D .
\]
Then \(K/n<7/100\), \(L_dK/n<\kappa=79/1000\),
and \(K+1<n/2\).
Use (1) through \(w=K/m\) with \(\lambda=79/100\).
The \(K\) test \((e/16)e^{\lambda W_d K/n}<19/50\)
and the midpoint test give \(c_{m,k}\le(19/50)^k\)
from \(K\) to the midpoint.

For a positive start the first-factor derivative gives
\(|(\log B_m)'|\le1/n\).
Indeed \(V/m<224/199\) and the private part after
multiplication by \(m\) is bounded by \(77/199\),
using \(\log(28m)/m<9/199\).
The total loss from \(1/m\) is strictly below \(2/m\),
which suffices for the absolute derivative bound.
Also \(B_m/m\) decreases. Using the coefficient estimate
at the integer order \(n\) gives \(B_m<(529/500)D\).

The [sharper Gamma bounds](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md#the-long-range-with-sharper-gamma-derivative-estimates)
apply with \(L_d,W_d\) in place of \(L_0,W\).
Indeed each Gamma argument on \(0\le x\le K\) is
at least \((1-\kappa)ma_i\), and the binomial argument at least
\((1-\kappa)m\). Hence the polygamma sum uses \(d_i^2/a_i\)
and each additional relative gap is bounded by \(L_d\).
A private ending contributes \(a_i\) to \(V_d\) as required.
With \(\alpha=L_d/[(1-\kappa)n]\), the bounds are
\(u=\alpha K<43/500\), \(H_2<(117/100)W_d/n<83/1000\),
\(H_r=(r-2)!\alpha^{r-2}H_2\), and derivative bound \(2rH_r/n\).
In the [connected logarithm sum](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_TWELVE.md#long-halves)
put \(\beta=e^{H_2+1/8}H_2\).
As \(L_d\le\sqrt{V_d}\),
\(\alpha/\beta\le e^{-1/8}\sqrt{13}/14<1/4\).
With tree base \((5/4)^5\), \(\beta K<141/100\),
and \(\beta B_m<13/400\), the connected parameter
is \(p_m<23/100\).
The truncated log \(A_m\) through \(K\) has coefficient
bounds \(B_m p_m^{k-1}/k\) and derivative bounds
\(7B_m p_m^{k-1}/n\).

The norms of \(A_m\) on \(2,6/5\) are below
\((27/10)B_m,(29/20)B_m\), respectively, and
\(\|A_N-A_n\|_{6/5}<13D/n\).
Relative to \(e^{A_m}\) on \(6/5\), the exponential
and actual tails have absolute coefficient norm bounds
\[
 (3/5)^{K+1} e^{(83/20)B_m}<e^{-111/20},\qquad
 \frac{(57/125)^{K+1}}{1-57/125}e^{(29/20)B_m}<e^{-17K/25}.
\]
For the first use \(K+1\ge16D+9\) and \(\log(5/3)>1/2\).
For the second keep \(D\le(K-8)/16\).
Their sum is below \(9/2000\) since \(K\ge13\).
It gives \(\|\log H_m-A_m\|_{6/5}<1/200\).

For a private start, the reservation bound and (1) give
the ratio \((6/5)e^{3/2+\lambda W_d K/n}q<9/500\)
through \(K\). Beyond use \(57/125\).
This gives \(\|\log H_m\|_{6/5}<19/1000\).
For either start the resulting bounds are
\(\|\log H_m\|_{6/5}\le(8/5)D\) and
\(\|\log H_N-\log H_n\|_{6/5}<3/25\).

## Joining the two orientations

For \(N>96\), either mass estimate gives
\(\|\log H_m\|_{6/5}/n\le M_0=1/120\)
and adjacent difference norm below \(\delta=3/25\).
The direct case uses \(n\ge96\), and the connected case
uses \((8/5)D/n<1/120\).
These estimates hold independently in the two orientations.

With \(u_0=\log(51/50)\), we have \(2M_0<u_0\), so
Rouché confines all root log radii to \(|u_m|<u_0\).
The weighted angular multipliers on radius \(51/50\)
relative to \(6/5\) are at most \(12/5\) and \(24\)
for one and two Euler derivatives. These follow from
\(\log(20/17)>4/25\) and continuous maximization of
\(x e^{-cx}\), \(x^2 e^{-cx}\).
The [paired real equations](MULTIVARIATE_RATIO_BALANCED_ANGULAR_SEPARATION.md#paired-phases-at-their-actual-radii)
therefore give at a common argument
\[
   n|u_N-u_n|\le
       \frac{u_0+2\delta}{1-2(12/5)M_0}<7/25 .
\]
Integrating the fixed-radius angular differences and the
mixed derivatives, the full phase correction difference
is bounded by
\[
                 \left(2(12/5)\delta+
                              2\cdot24 M_0(7/25)\right)\theta
                                <\theta .
\]
At \(0<\theta<\pi\) this leaves the consecutive phase
difference strictly between \(0\) and \(2\pi\).
Real half logs at the negative argument permit only odd
order, positive coefficients exclude the positive argument,
and real coefficients give the lower-half-plane result.
This completes the all-degree comparison.

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

S=F(13); sigma=S/2; U=1/lo(sigma); R=F(6,5)
en=lambda m,k:m/lo(sigma*(1-F(k,m)))
for name,k,ms,cap in [('first',1,[2,11],F(39,500)),('first',1,[2,96],F(31,200)),
                      ('second',2,[4,96],F(7,500)),('third',3,[6,96],F(11,5000))]:
    ck('finite short '+name,max(en(m,k) for m in ms)**k/factorial(k),cap)
ck('finite short midpoint',4/lo(S/4),F(4,25))
ck('finite short tail',hi(1)*en(96,4)/4,F(4,25))
for m,c1,c2,c3,r in [(2,F(39,500),F(7,500),F(11,5000),0),
                     (8,F(39,500),F(7,500),F(11,5000),F(4,25)), (12,F(31,200),F(7,500),F(11,5000),F(4,25))]:
    rad=hi(F(2,5*m))
    x=rad*r
    A=c1*rad+c2*rad**2+c3*rad**3+x**4/(1-x)
    B=c1*rad+2*c2*rad**2+3*c3*rad**3+x**4*(4/(1-x)+x/(1-x)**2)
    ck('short amplitude',1+A,lo(F(2,5))*(1-A))
    ck('short phase',B/(1-A),F(1,4))

d=F(3,10); a=F(61,1000); cap=F(183,100); d0=(a/cap)**3
lam=F(39,50); b=F(23,250)
ck('direct L/n',lam*224/96,cap)
ck('direct mass comparison',lam**3*14**3*U*F(97,96)/96**2,a**3)
ck('direct slope',F(12,23)+F(1,4),lam)
ck('direct radical',a**3,d*b**3)
caps=[max(d0*hi(k*cap),d*hi(k*b)) for k in range(1,5)]
cc=[F(33,100),F(33,500),F(11,1000)]
for k,c in enumerate(cc,1):
    ck('direct heads',caps[k-1]**k/factorial(k),c)
ck('direct tail',hi(1)*caps[3]/4,F(3,10))
ck('direct midpoint',2*hi(1)/lo(S/4),F(3,10))
rad=F(9,25); tail=rad**4/(1-rad)
h=sum(R**k*c for k,c in enumerate(cc,1))+tail
ck('direct positivity',h,1)
ck('direct norm',1/lo(F(4,5)),1-h)
ck('direct Wq',14*F(97,96)*U,F(213,10000))
ck('direct first derivative log',28*96,lo(8))
ck('direct first derivative',F(3,4)*224/96+F(70,96),3)
ck('direct arguments',3*15/96,F(1,2))
loss=F(7,192)
h1=sum(k*R**k*c*hi(loss*k) for k,c in enumerate(cc,1))
h2=sum(k*(k-1)*R**k*hi(loss*k)*max(d0**(k-1)*hi(k*k*cap),d**(k-1)*hi(k*k*b))/factorial(k)
       for k in range(2,4))
ck('direct log difference',(3*h1/96+8*F(213,10000)/96*h2+tail)/(1-h),F(3,25))
headprivate=R*hi(F(3,2))*caps[2]/96
ck('direct private heads',headprivate,F(23,1000))
ck('direct private tail ratio',R/lo(F(7,4)),F(21,100))
private=F(23,1000)
tp=F(21,100)
ck('direct private log',private/(1-private)+tp**4/(1-tp),1-1/lo(F(4,125)))

r0=F(17,10); r1=F(7,4)
slope=lambda r:1/r+r*(2+r)/(2*(1+r)**2)
ck('ascending gap critical lower',1+r0,lo(slope(r0)))
ck('ascending gap critical upper',hi(slope(r1)),1+r1)
log0=F(993,1000)
ck('ascending gap lower log',hi(log0),1+r0)
e0=(1+r0)*log0-r0-r0*r0/(2*(1+r0))
ck('ascending gap absorption',r1/lo(e0),F(9,8))

q=F(1512,10**6); Wq=F(2118,100000); Kcap=46
kap=F(79,1000); u=F(43,500); lam=F(79,100)
ck('long n',199,d*lo(sigma))
ck('long q',F(200,199)*U,q)
ck('long Wq',14*F(200,199)*U,Wq)
ck('long L0q',F(9,8)*F(200,199)*U,F(171,100000))
assert 16+9/d==Kcap
ck('long kappa',F(171,100000)*Kcap,kap)
ck('long range',q*Kcap+F(1,199),F(1,2))
ck('long local slope range',q*Kcap,F(7,100))
ck('long local slope',F(50,93)+F(1,4),lam)
ck('long u',kap/(1-kap),u)
ck('long derivative private log',28*199,lo(9))
ck('long derivative',F(3,4)*224/199+F(77,199),2)
b=F(529,500)
ck('long B',hi(lam*Wq/d),b)
ck('long H2 prefactor',(1+q*Kcap/(1-kap))/(1-kap),F(117,100))
ck('long H2',F(117,100)*Wq/d,F(83,1000))
ck('long tail',hi(1+lam*Kcap*Wq)/16,F(19,50))
ck('long midpoint',2*hi(1)/lo(S/4),F(19,50))
beta=F(117,100)*hi(F(83,1000)+F(1,8))
ck('long betaK',beta*Kcap*Wq,F(141,100))
ck('long betaB',beta*Wq*b,F(13,400))
ck('long tree t',13*16/lo(F(1,4)),14**2)
p=F(23,100)
ck('long p',F(5,4)**5*F(13,400)/(1-u)*hi(F(141,100)*(F(1,2)+u/(6*(1-u)))),p)
Aout=F(27,10); Ain=F(29,20)
ck('long outer norm',1/lo(Aout*p),1-2*p)
ck('long inner norm',1/lo(Ain*p),1-R*p)
ck('long adjacent norm',7*b*R/(1-R*p),13)
ck('long exp tail',F(111,20),F(9,2)+(8-(Aout+Ain)*b)*d)
ck('long log radius',hi(F(1,2)),F(5,3))
real=F(57,125)
ck('long true tail slope',hi(F(17,25)+Ain*b/16),1/real)
ck('long true tail intercept',real/(1-real),lo(Ain*b*8/16))
ck('long sum of tails',1/lo(F(111,20))+1/lo(F(17,25)*13),F(9,2000))
ck('long log correction',1/lo(F(1,200)),1-F(9,2000))
ck('long positive delta',13*q+F(1,100),F(3,25))
private=R*hi(F(3,2)+lam*Kcap*Wq)*U
ck('long private ratio',private,F(9,500))
ck('long private norm',F(9,500)/(1-F(9,500))+real**14/(1-real),1-1/lo(F(19,1000)))
ck('long private delta',2*F(19,1000),F(3,25))
ck('uniform M for positive',Ain*b+F(1,200)/d,F(8,5))
ck('uniform M for private',F(19,1000)/d,F(8,5))
M=F(1,120)
ck('M loss for positive long',F(8,5)*q,M)
umin=F(99,5000); umax=F(1,50)
ck('long radii',2*M,umin)
ck('slope loss log',hi(F(4,25)),F(20,17))
ck('slope one',1/(lo(1)*F(4,25)),F(12,5))
ck('slope two',4/(lo(1)*F(4,25))**2,24)
ck('long shift',(umax+2*F(3,25))/(1-2*F(12,5)*M),F(7,25))
ck('long phase',2*F(12,5)*F(3,25)+2*24*M*F(7,25),1)
print('All unrestricted S=13 scalar checks passed.')
~~~
