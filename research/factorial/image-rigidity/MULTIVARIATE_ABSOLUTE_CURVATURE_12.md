# All-degree angular separation at absolute curvature twelve

Let \(a,b\) be distinct ordinary nonnegative exponent columns,
with equal rows omitted. Put
\[
              S=\tau+\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)} .
\]
An optional common Gaussian strength \(\tau\ge0\) multiplies
the normalized coefficient of index \(k\) in degree \(m\)
by \(e^{-\tau k(m-k)/(2m)}\).

**Theorem.** If \(S\ge12\), consecutive moment polynomials have
disjoint root arguments in every degree. The exponent rows
and their orientations are arbitrary. At \(\tau=0\) this
proves ordinary SFC(2) over every characteristic-zero field.

The [curvature thirteen theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_13.md)
covers \(S\ge13\). This proof uses the same cubic directed
endpoint budget on \(12\le S<13\). It retains the sign of
the first coefficient in the early phase comparison and
the factorial prefactor in geometric coefficient tails.
For the linked Gamma estimates, a large increasing ratio
affects higher derivatives but does not shrink any Gamma
argument on the coefficient half.

## Setup and a geometric coefficient prefactor

In each orientation keep the
[directed quantities](MULTIVARIATE_ABSOLUTE_CURVATURE_13.md#directed-endpoint-estimates)
\(Q,q=Q(0),V_d,W_d=1+V_d,L_d\) from the fully retained
oriented entropy envelope. Keep \(V\) for the budget
with the smaller exponent on every shared denominator.
Thus \(S\le V_d\le V\le S(S+2)<195\) and \(L_d\le\sqrt{V_d}\).
For \(w\le w_*\le1/2\) we have
\[
 \begin{split}
 &Q(w)\le q e^{\lambda V_d w},\qquad
          \lambda=\frac{1}{2(1-w_*)}+\frac14,\qquad
          Q(1/2)\le e^{-S/4},\\
 &q\le e^{-S/2},\qquad
          W_d^3q\le (1+S)^3e^{-S/2},\qquad
          L_d q\le(9/8)e^{-S/2}.
 \end{split}
\]
In particular \(W_d q\le(1+S)e^{-S/2}\). All three
endpoint budgets hold at \(S\ge12\) by the cited proof.

For an integer degree \(m\), the normalized coefficient
at \(w=k/m\in(0,1/2]\) is bounded by
\([mQ(w)]^k/k!\). The maximum of \(eQ(w)/w\) on each
closed subinterval is attained at an endpoint.
When that maximum is \(r\), retain the factorial as
\[
              c_{m,k}\le a_k r^k,\qquad
                        a_k=\frac{k^k e^{-k}}{k!}.                \tag{1}
\]
The sequence \(a_k\) decreases
since \(a_{k+1}/a_k=(1+1/k)^k/e<1\).
In particular \(a_k<1/5\) at \(k\ge4\) and \(a_k<9/50\)
at \(k\ge5\).

## Early degrees and the sign of the first term

Write \(N=n+1\), \(m=n,N\), and first suppose \(N\le82\).
Use the discrete curvature estimate
\(c_{m,k}\le\binom mk e^{-6k(1-k/m)}\).
For each \(1\le k\le5\) in the half, the upper envelope
\([m e^{-6(1-k/m)}]^k/k!\) takes its maximum on
a closed subinterval of \(m\ge2k\) at an endpoint, by
the sign change of its logarithmic derivative.
For \(k\ge6\) use the
[binomial entropy factor](MULTIVARIATE_ENTROPY_ENDPOINT_PHASE_CRITERION.md#combining-the-binomial-entropy-and-the-row-chord)
\(B(w)=\exp([-w\log w-(1-w)\log(1-w)]/w)\).
The function \(w(\log (B(w)e^{-6(1-w)}))'\)
is increasing because its derivative is at least
\(6-4(1-\log 2)>0\).
Test \(B(w)e^{-6(1-w)}\) at \(w=1/2\) and \(w=6/m\).
For the latter test the interval of \(m\) has endpoints
\(12,82\), by the same monotonicity.
At the midpoint the bound is \(4e^{-3}<1/5\),
and at \(m=82\) it is below \(1/5\) using \(B(w)\le e/w\).
Thus \(c_{m,k}\le5^{-k}\) beyond the first five.

Set \(R_m=e^{13/(20m)}\).
For \(m\) in each of \([2,7],[8,31],[32,82]\), the
scalar checks give uniform bounds \(A\) on \(|H_m-1|\),
\(F\) on \(c_{m,1}R_m\), and \(T\) on
\(\sum_{k\ge2} k c_{m,k}R_m^k\) in either orientation.
They use the smallest order in the group to bound \(R_m\)
and the coefficient envelopes at the two endpoints of
the interval for each index in a half. They give
\((1+A)/(1-A)<e^{13/20}\).
Rouché in both orientations confines all root
log radii to \(|\log|z||<13/(20m)\).
The degree-one halves both equal one.

For \(0<\theta<\pi\), the first term gives a nonnegative
contribution to \(\Im H(re^{i\theta})\).
The remaining imaginary part in modulus is at most
\(T\theta\), using \(|\sin k\theta|\le k\theta\).
As \(\Re H\ge1-A>0\), the log phase at any half
argument of radius \(r\le R_m\) satisfies
\[
          -\frac{T\theta}{1-A}\le\Im\log H(re^{i\theta})
                    \le\frac{(F+T)\theta}{1-A}.
\]
The full phase correction at a root is the sum of the
two half log phases on the ray of angle \(\theta\), at
reciprocal radii. For two consecutive orders in the
same group the difference has modulus at most
\(2(F+2T)\theta/(1-A)<\theta\).
At a group boundary use the upper estimate from one
group and the lower estimate from the other.
The scalar checks include both consecutive group boundaries.
Consequently the difference is again strictly below \(\theta\).
The [common-angle phase contradiction](MULTIVARIATE_RATIO_BALANCED_ANGULAR_SEPARATION.md#paired-phases-at-their-actual-radii)
then excludes a common argument. At \(\theta=\pi\) the
logs are real and only odd orders can have a root there,
while positivity excludes the positive real axis.

## Direct estimates in one orientation

Assume \(N>82\) and set \(D=Nq\) separately in each
orientation. The direct range is \(D\le2/5\).
For \(w\le4/82\), the increment bound holds with
\(\lambda=39/50\). The cubic budget gives
\[
              \lambda V_d/n\le
                  \min\!\left(93/50,\frac{3}{40D^{1/3}}\right).
                                                               \tag{2}
\]
Indeed the cube times \(D\) is at most
\(\lambda^3 13^3 e^{-6}(83/82)/82^2\), and
\(\lambda V_d/n<\lambda\cdot195/82\).
The bound for \(D\exp(k\lambda V_d/n)\) at \(k\le4\)
uses the switch \(D_0=((3/40)/(93/50))^3\) and \(D=2/5\).
On the second piece the log derivative changes sign
at most once from negative to positive.
For \(k=1,2,3,4\) use these tests with the factorial
denominator to obtain the caps
\(111/250,121/1000,27/1000,7/1250\).
The \(k=4\) test with \(eQ(w)/w\) and the midpoint test
are both below \(41/100\).
By (1), \(c_{m,k}\le(9/50)(41/100)^k\) at \(k\ge5\).

All direct bounds below use the joining radius \(R=103/100\).
Write \(\|H_m-1\|_R\le h\) in absolute coefficient norm,
and use \(B_H,C_H\) for the bounds on the corresponding
first and second Euler derivatives. The four caps and
geometric tail give \(\|\log H_m\|_R<1\), and
\[
 \|z(\log H_m)'\|_R\le B_H/(1-h)<5/2,\qquad
 \|(z\partial_z)^2\log H_m\|_R
               \le C_H/(1-h)+(B_H/(1-h))^2<11.                  \tag{3}
\]

With positive starting exponents, interpolate the first
four coefficients on \(n\le m\le N\).
The [first-factor derivative](MULTIVARIATE_ABSOLUTE_CURVATURE_15.md#exact-derivative-at-a-positive-start)
gives \(|(\log B_m)'|\le2/n\).
Indeed \(0\le1/m-(\log B_m)'\), and after multiplication
by \(m\) the total loss is strictly below \(3\),
by \(V/m<195/82\), private curvature below \(13\),
and \(\log(26m)/m<8/82\).
For \(f_m(x)=\log(c_{m,x}x!/B_m^x)\), \(0\le x\le4\),
each Gamma argument is at least \((1-4/n)ma_i\),
and the binomial argument at least \((1-4/n)m\).
An increasing row has no reduction in its argument,
and a decreasing row satisfies \(|d_i|\le a_i\).
Use the exact bound \(|\psi_2(y)|\le(1+2/y)y^{-2}\)
to get \(\left|\partial_m\partial_x^2 f_m\right|
\le(6/5)W_d/n^2\).
Since \(\partial_m f_m\) vanishes at \(0,1\), this implies
\[
          |\partial_m c_{m,k}|
           \le \frac{c_{m,k}}{n}
                      \left(2k+\frac{3W_d}{5n}k(k-1)\right),
                    \qquad 1\le k\le4.                         \tag{4}
\]
For the real-degree coefficients multiply the integer-order
coefficient envelope by \(e^{13k/328}\), using the
[private ending Gamma correction](MULTIVARIATE_ABSOLUTE_CURVATURE_13.md#direct-estimates-for-one-orientation).
This factor bounds the total loss \(k/(4m)\) per private
ending row. The shared envelopes hold by the Gamma chord
bound at real degree, and the binomial bound \(m^k/k!\)
also holds there.

For the second term of (4) retain \(W_dD/n<327/10000\)
and use the envelope for \(c_{m,k}/D\).
The maximum of \(D^{k-1}\exp(k^2\lambda V_d/n)\)
uses the same two endpoint tests as (2).
The adjacent coefficient differences beginning at \(k=5\)
are each bounded by the single interval length
\((9/50)(41/100)^k\).
Let \(E_0,E_1\) bound the absolute coefficient norms
of the adjacent half difference and its Euler derivative
at \(R\), respectively.
Integration along the segment between the halves gives
\(\|\log H_N-\log H_n\|_R\le E_0/(1-h)<1/10\).
For the angular derivative the quotient rule gives
\[
  \|z(\log H_N-\log H_n)'\|_R
              \le \frac{E_1}{1-h}
                        +\frac{B_H E_0}{(1-h)^2}<2/5.          \tag{5}
\]
All numerical bounds in (3)–(5) are checked below.

For a private starting row the
[private reservation bound](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md#uniform-comparison-across-degrees)
gives \(c_{m,k}^{1/k}\le e^{3/2}Q(w)\).
Use (2) for \(k\le4\), yielding ratio \(7/200\) at \(R\).
For \(k\ge5\), use \(Q(w)\le Q(1/2)\), yielding ratio
\(231/1000\) at \(R\). These bounds give log norm below
\(1/20\), first derivative below \(1/5\), and second
derivative below one at \(R\).
Thus (3), the adjacent log bound \(1/10\), and (5) hold
also for this starting configuration.

## Connected estimates in one orientation

For \(D>2/5\) we have \(n\ge161\).
The directed budgets give
\[
 D/n<2496/10^6,\qquad W_dD/n<3245/10^5,\qquad
 L_dD/n<281/10^5,\qquad K=\lceil12D\rceil+6\le(59/2)D .
\]
Then \(K/n<\kappa=3/40\) and \(K+1<n/2\). The increment
bound through \(w=K/m\) holds with \(\lambda=4/5\).
The \(K\) and midpoint endpoint tests give
\(eQ(w)/w<1/2\) throughout \(K/m\le w\le1/2\).
Retain \(a_k<1/5\) in (1).

For a positive start the
[first-factor derivative](MULTIVARIATE_ABSOLUTE_CURVATURE_15.md#exact-derivative-at-a-positive-start)
gives \(|(\log B_m)'|\le1/n\), using \(V/m<195/161\)
and \(\log(26m)/m<9/161\) in the private bound,
so the loss from \(1/m\) is below \(2/m\).
Also \(B_m/m\) decreases, and at the integer order
\(n\) the coefficient bound gives \(B_m<(267/250)D\).

On \(0\le x\le K\), the Gamma arguments are at least
\((1-\kappa)ma_i\) and, for the binomial,
\((1-\kappa)m\). Again this uses the sign of each slope,
so the lower argument bound requires \(K/n<\kappa\).
An additional relative gap in the higher derivative
is bounded by \(L_d\). The
[polygamma finite-difference estimates](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md#the-long-range-with-sharper-gamma-derivative-estimates)
thus give, with \(\alpha=L_d/[(1-\kappa)n]\),
\[
 H_2=\frac{W_d}{(1-\kappa)n}
       \left(1+\frac{K}{(1-\kappa)n}\right)
   <(117/100)W_d/n<19/200,\qquad u=\alpha K<9/100 .
\]
The bounds are \(H_r=(r-2)!\alpha^{r-2}H_2\) and
\(2rH_r/n\) for the difference and its order derivative.
A private ending contributes \(a_i\) to \(V_d\);
the same second-difference formulas therefore apply.
The derivative bounds absorb the Gaussian contribution
using its summand in \(W_d\).

Use the [charged connected expansion](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_TWELVE.md#long-halves)
with \(\beta=e^{H_2+1/8}H_2\).
The tree base is at most \((5/4)^5\) since
\(\alpha/\beta\le e^{-1/8}\sqrt{12}/13<1/4\).
With \(\beta K<7/5\) and \(\beta B_m<51/1000\),
the connected parameter is \(p_m<9/25\).
The formal logarithm \(A_m\) through degree \(K\)
therefore has coefficient bounds
\(B_m p_m^{k-1}/k\) and derivative bounds
\(7B_m p_m^{k-1}/n\).
The norms at radii \(2,6/5\) are below
\((18/5)B_m,(8/5)B_m\), respectively, and
\(\|A_N-A_n\|_{6/5}<16D/n\).

The relative exponential and actual tails on \(6/5\)
have bounds
\[
 \begin{split}
 &(3/5)^{K+1}e^{(26/5)B_m}<e^{-379/100},\\
 &\tfrac3{10}(3/5)^K e^{(8/5)B_m}
   <\tfrac3{10}\exp[-(9/25)K-(8/5)(267/250)/2] .
 \end{split}
\]
For the first use \(K+1\ge12D+7\) and \(\log(5/3)>51/100\).
The second uses the factorial prefactor \(1/5\)
and \(D\le(K-6)/12\).
Since \(K\ge11\), the sum is below \(13/500\),
giving \(\|\log H_m-A_m\|_{6/5}<27/1000\).

For a private starting row the reservation bound gives
ratio \((6/5)e^{3/2+\lambda W_d K/n}q<29/1000\)
on \(6/5\) through \(K\).
Beyond \(K\), use the ratio \(3/5\) with coefficient
prefactor \(1/5\). The log norm is below \(4/125\).
For both types of start this proves
\[
          \|\log H_m\|_{6/5}\le (9/5)D,\qquad
          \|\log H_N-\log H_n\|_{6/5}<1/10 .                    \tag{6}
\]

## Joining direct and connected estimates

The two orientations may have different \(D=Nq\).
In both, at radius \(R=103/100\) the half log norm
divided by \(n\) is at most \(M_0=1/82\), by the
direct bound \(1\) or (6). Since \(2M_0<\log R\),
Rouché in both orientations confines root log radii
to \(|u_m|<u_0=\log R\).

The connected half bounds on \(6/5\) give first and
second Euler derivative bounds on \(R\) by multiplying
the norms by \(5/2\) and \(25\), respectively.
These constants follow from \(\log(120/103)>3/20\)
and continuous maximization of \(x e^{-cx}\) and \(x^2e^{-cx}\).
Thus the connected norm bound \((9/5)D/n\)
also gives first and mixed derivative bounds divided by \(n\)
of \((5/2)M_0\) and \(11M_0\), respectively.
The direct bounds (3) give the same estimates.
In both cases the adjacent log norm on \(R\) is below \(1/10\),
and the adjacent angular derivative at fixed radius is
below \(2/5\), by (5) or the multiplier applied to (6).

At a common angle the [real radius equations](MULTIVARIATE_RATIO_BALANCED_ANGULAR_SEPARATION.md#paired-phases-at-their-actual-radii)
therefore give
\[
      n|u_N-u_n|\le
        \frac{u_0+2/10}{1-2(5/2)M_0}<1/4 .
\]
Integrate in the angle first at a common radius, then
use the mixed derivatives for the two actual radii.
The full phase correction difference is bounded by
\((2(2/5)+2\cdot11M_0/4)\theta<\theta\).
This excludes a common angle on \(0<\theta<\pi\)
as in the early range. The two real arguments and
reflection complete the proof.

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

S=F(12); sigma=S/2; U=1/lo(sigma)
en=lambda m,k:m/lo(sigma*(1-F(k,m)))
ck('finite short tail midpoint',4/lo(3),F(1,5))
ck('finite short tail endpoint',hi(1)*en(82,6)/6,F(1,5))
short=[]
for m0,m1,kmax,r in [(2,7,3,0),(8,31,5,F(1,5)),(32,82,5,F(1,5))]:
    R=hi(F(13,20*m0)); x=R*r
    A=T=first=F(0)
    for k in range(1,kmax+1):
        start=max(2*k,m0)
        c=R**k*max(en(m,k) for m in [start,m1])**k/factorial(k)
        A+=c
        if k==1:
            first=c
        else:
            T+=k*c
    A+=x**6/(1-x)
    T+=x**6*(6/(1-x)+x/(1-x)**2)
    ck('short amplitude',1+A,lo(F(13,20))*(1-A))
    ck('short signed phase',2*(first+2*T)/(1-A),1)
    short.append(((first+T)/(1-A),T/(1-A)))
for (up,down),(next_up,next_down) in zip(short,short[1:]):
    ck('short adjacent straddle',2*max(up+next_down,next_up+down),1)

d=F(2,5); a=F(3,40); cap=F(93,50); d0=(a/cap)**3
lam=F(39,50); b=F(51,500)
ck('direct L/n',lam*195/82,cap)
ck('direct mass comparison',lam**3*13**3*U*F(83,82)/82**2,a**3)
ck('direct slope',F(41,78)+F(1,4),lam)
ck('direct radical',a**3,d*b**3)
caps=[max(d0*hi(k*cap),d*hi(k*b)) for k in range(1,5)]
cc=[F(111,250),F(121,1000),F(27,1000),F(7,1250)]
for k,c in enumerate(cc,1):
    ck('direct heads',caps[k-1]**k/factorial(k),c)
ck('direct tail',hi(1)*caps[3]/4,F(41,100))
ck('direct midpoint',2*hi(1)/lo(S/4),F(41,100))
ck('direct tail prefactor',F(5)**5/(lo(5)*factorial(5)),F(9,50))
ck('long tail prefactor',F(4)**4/(lo(4)*factorial(4)),F(1,5))
r=F(103,100); x=F(41,100)*r
z=F(9,50)
T0=z*x**5/(1-x)
T1=T0*(5+x/(1-x))
T2=T0*(25+10*x/(1-x)+x*(1+x)/(1-x)**2)
s=[sum(k**j*r**k*c for k,c in enumerate(cc,1)) for j in range(3)]
h,bh,ch=s[0]+T0,s[1]+T1,s[2]+T2
ck('direct positivity',h,1)
ck('direct norm',1/lo(1),1-h)
ck('direct P',bh/(1-h),F(5,2))
ck('direct C',ch/(1-h)+(bh/(1-h))**2,11)
ck('direct Wq',13*F(83,82)*U,F(327,10000))
ck('direct first derivative log',26*82,lo(8))
ck('direct first derivative',F(3,4)*195/82+F(65,82),3)
ck('direct Gamma derivative',(1+F(2,78))/F(78,82)**2,F(6,5))
loss=F(13,328)
Es=[]
for j,T in enumerate([T0,T1]):
    h1=sum(k**(j+1)*r**k*c*hi(loss*k) for k,c in enumerate(cc,1))
    h2=sum(k**j*k*(k-1)*r**k*hi(loss*k)*max(d0**(k-1)*hi(k*k*cap),d**(k-1)*hi(k*k*b))/factorial(k)
           for k in range(2,5))
    Es.append((2*h1+F(3,5)*F(327,10000)*h2)/82+T)
E0,E1=Es
ck('direct log difference',E0/(1-h),F(1,10))
ck('direct Euler difference',E1/(1-h)+bh*E0/(1-h)**2,F(2,5))
private=r*hi(F(3,2))*caps[3]/82
ck('direct private head ratio',private,F(7,200))
ck('direct private tail ratio',r/lo(F(3,2)),F(231,1000))
xb=F(7,200); xt=F(231,1000)
Ap=xb/(1-xb)+xt**5/(1-xt)
Bp=xb/(1-xb)**2+xt**5*(5/(1-xt)+xt/(1-xt)**2)
Cp=xb*(1+xb)/(1-xb)**3+xt**5*(25/(1-xt)+10*xt/(1-xt)**2+xt*(1+xt)/(1-xt)**3)
ck('direct private log',Ap,1-1/lo(F(1,20)))
ck('direct private P',Bp/(1-Ap),F(1,5))
ck('direct private C',Cp/(1-Ap)+(Bp/(1-Ap))**2,1)
assert 2*F(1,20)==F(1,10)

R=F(6,5); q=F(2496,10**6); Wq=F(3245,100000); Kcap=F(59,2)
kap=F(3,40); u=F(9,100); lam=F(4,5)
ck('long n',161,d*lo(sigma))
ck('long q',F(162,161)*U,q)
ck('long Wq',13*F(162,161)*U,Wq)
ck('long L0q',F(9,8)*F(162,161)*U,F(281,100000))
assert 12+7/d==Kcap
ck('long kappa',q*Kcap,kap)
ck('long range',q*Kcap+F(1,161),F(1,2))
ck('long slope',F(20,37)+F(1,4),lam)
ck('long u',F(281,100000)*Kcap/(1-kap),u)
ck('long derivative private log',26*161,lo(9))
ck('long derivative',F(3,4)*195/161+F(143,2*161),2)
b=F(267,250)
ck('long B',hi(lam*Wq/d),b)
ck('long H2 prefactor',(1+q*Kcap/(1-kap))/(1-kap),F(117,100))
ck('long H2',F(117,100)*Wq/d,F(19,200))
ck('long tail',hi(1+lam*Kcap*Wq)/12,F(1,2))
ck('long midpoint',2*hi(1)/lo(S/4),F(1,2))
beta=F(117,100)*hi(F(19,200)+F(1,8))
ck('long betaK',beta*Kcap*Wq,F(7,5))
ck('long betaB',beta*Wq*b,F(51,1000))
ck('long tree t',12*16/lo(F(1,4)),13**2)
p=F(9,25)
ck('long p',F(5,4)**5*F(51,1000)/(1-u)*hi(F(7,5)*(F(1,2)+u/(6*(1-u)))),p)
Aout=F(18,5); Ain=F(8,5)
ck('long outer norm',1/lo(Aout*p),1-2*p)
ck('long inner norm',1/lo(Ain*p),1-R*p)
ck('long adjacent norm',7*b*R/(1-R*p),16)
ck('long exp tail',F(379,100),F(357,100)+(F(153,25)-(Aout+Ain)*b)*d)
ck('long log radius',hi(F(51,100)),F(5,3))
real=F(3,5)
ck('long true tail slope',hi(F(9,25)+Ain*b/12),1/real)
ck('long sum of tails',1/lo(F(379,100))+F(3,10)/(lo(F(9,25)*11)*lo(Ain*b/2)),F(13,500))
ck('long log correction',1/lo(F(27,1000)),1-F(13,500))
ck('long positive delta',16*q+F(27,500),F(1,10))
private=R*hi(F(3,2)+lam*Kcap*Wq)*U
ck('long private ratio',private,F(29,1000))
ck('long private norm',F(29,1000)/(1-F(29,1000))+real**12/(5*(1-real)),1-1/lo(F(4,125)))
ck('long private delta',2*F(4,125),F(1,10))
ck('uniform M for positive',Ain*b+F(27,1000)/d,F(9,5))
ck('uniform M for private',F(4,125)/d,F(9,5))
M=F(1,82)
ck('M loss for positive long',F(9,5)*q,M)
umin=F(591,20000); umax=F(3,100)
ck('long radii',2*M,umin)
ck('slope loss log',hi(F(3,20)),F(120,103))
ck('slope one',1/(lo(1)*F(3,20)),F(5,2))
ck('slope two',4/(lo(1)*F(3,20))**2,25)
ck('long Euler',F(5,2)/10,F(2,5))
ck('long mixed',25*F(9,5)*q,11*M)
ck('long shift',(umax+2*F(1,10))/(1-2*F(5,2)*M),F(1,4))
ck('long phase',2*F(2,5)+2*11*M/4,1)
print('All unrestricted S=12 scalar checks passed.')
~~~
