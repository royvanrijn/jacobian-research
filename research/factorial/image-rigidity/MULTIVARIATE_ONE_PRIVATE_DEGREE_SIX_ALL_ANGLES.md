# All-degree angular separation at private degree six

Suppose an ordinary two-monomial profile has private rows on exactly
one side, of total private degree \(P=6\). The consecutive factorial
moment polynomials have disjoint root arguments in every degree, for
arbitrary shared exponents and orientations.

The [one-private entropy theorem](MULTIVARIATE_ONE_PRIVATE_ENTROPY_THIRTEEN.md)
covers \(S_{\rm sh}\ge1\), where
\(S_{\rm sh}=\sum_{\rm sh}(b_i-a_i)^2/\max(a_i,b_i)\).
Assume below that \(S_{\rm sh}<1\). In every unequal shared row,
\(d^2/(\ell+|d|)<1\) with \(\ell=\min(a,b)\), so
\(|d|\le\ell\) by integrality. Then \(V_{\rm sh}=
\sum_{\rm sh}d^2/\ell\le2S_{\rm sh}\) and the actual
private-plus-shared derivative budget is \(V=6+V_{\rm sh}<8\).
Every shared relative slope has magnitude at most one.

Orient the halves separately. The half starting at the side having
private variables is denoted \(H^{\rm e}_m\), and the half starting at
the other side \(H^{\rm s}_m\). Thus \(H^{\rm e}_m\) uses coefficients
at distance \(k\) from the endpoint with positive private exponents.
The superscript \({\rm s}\) indicates that the private exponents at
the starting endpoint are zero. Both halves use the normalization,
midpoint split and Gamma interpolation of the
[entropy proof](MULTIVARIATE_ONE_PRIVATE_ENTROPY_THIRTEEN.md#the-private-factorial-estimates).

## Short intervals with separate half bounds

For \(13\le m\le82\) the coefficients of \(H^{\rm e}_m\) obey
\[
         c_{m,k}\le\binom m k e^{-s k(1-k/m)} ,
\]
with \(s=107/20\) on \([13,31]\) and \(s=351/62\) on \([32,82]\).
In fact (3) of the entropy proof gives \(s=6(1-\delta_m)\),
with \(\delta_m\le13/120\) and \(\delta_m<7/124\) on those groups.
Use the same endpoint maxima as in its short argument for the
first five coefficients and the binomial-entropy tail from the
sixth coefficient. The corresponding tail ratios are \(7/25\)
and \(6/25\).

For \(H^{\rm s}_m\), reserve one factor \(f_{m,k}\) to cancel the
binomial numerator as in the entropy proof. The resulting
coefficient bound is
\[
 c_{m,k}\le e^k f_{m,k}^5
          \le e^k (k!)^5/L^{5k},
\]
where \(L=113/20\) and \(12\), respectively, using
\((m_0!)^{1/m_0}>L\) at each lower group endpoint \(m_0\).
The sequence \((m!)^{1/m}\) increases. For \(k\ge8\) use instead
\(\log f_{m,k}/k\le\log w+C_k\), \(w=k/m\), with \(C_k\le13/50\)
by (2) of the entropy proof. This gives the tail ratio
\(e^{23/10}/32<8/25\).

The exact replay below bounds the two half norms at
\(R_m=e^{13/(20m)}\), using the corresponding coefficient
bounds at \(R=e^{13/(20m_0)}\). Write \(A_j,F_j,T_j\) for the
bound on the coefficient norm of \(H^j_m-1\), the scaled first
coefficient, and the Euler norm beyond the first coefficient.
In each group
\[
 \frac{1+A_{\rm e}}{1-A_{\rm s}},\
 \frac{1+A_{\rm s}}{1-A_{\rm e}} < e^{13/20}.
\]
Rouché gives \(R_m^{-1}<|z|<R_m\) for every root.
At angle \(\theta\in(0,\pi)\), the upper bound on the sum of
the two half arguments, divided by \(\theta\), is
\(u=\sum_j(F_j+T_j)/(1-A_j)\), and the lower bound is
\(-v=-\sum_j T_j/(1-A_j)\). This uses the nonnegative
first sine and \(|\sin(k\theta)|\le k\theta\). The replay
checks \(u+v<1\) in each group and
\(\max(u_1+v_2,u_2+v_1)<1\) at their join.
At a common angle, the phase difference of the two half-log sums
would have to equal \(\theta\bmod 2\pi\), which is impossible
under these bounds. This covers \(14\le N=n+1\le82\).

## Long intervals and coefficient envelopes

For \(N\ge83\) use the [oriented envelope](MULTIVARIATE_ABSOLUTE_CURVATURE_13.md#directed-endpoint-estimates)
\(Q(w)\) for \(H^{\rm e}_m\), with formal private chord parameter
\(s_0=35/3\) and all shared rows retained. The formal parameter
is used in the envelope, while all derivatives involve the
original Gamma factors. Define
\[
 q=Q(0),\quad D=Nq,\quad W_d=1+s_0+\sum_{\rm sh}d_i^2/a_i ,
 \qquad S'=s_0+S_{\rm sh}.
\]
Here \(a_i\) denotes the starting exponent in this orientation and
\(d_i\) the signed gap. Since \(6\delta_m\le1/6\) for real \(m\ge82\),
the private endpoint bound (3) of the entropy proof supplies
\(c_{m,k}\le[mQ(w)]^k/k!\) also at real \(m\) and integer \(k\).
The oriented estimates for \(Q\) give
\[
 q\le e^{-S'/2}<3/1000,\quad
 Q(1/2)\le e^{-35/12},\quad W_d\le44/3,\quad
 W_d q\le (38/3)e^{-35/6}.
\tag{7}
\]
For the last bound use \(W_d\le38/3+2S_{\rm sh}\) and the
decrease of \((38/3+2s)e^{-s/2}\) for \(s\ge0\).
Also \(Q(w)\le q e^{\lambda(W_d-1)w}\) for \(w\le w_*\le1/2\)
with \(\lambda=1/[2(1-w_*)]+1/4\), and \(eQ(w)/w\) attains
its maximum on any closed subinterval at an endpoint.
These estimates use the same chord and monotonicity proof for
the nonnegative parameter \(s_0\). No cubic endpoint budget
is needed.

For \(H^{\rm s}_m\), use (6) of the entropy proof to obtain
\(c_{m,k}^{1/k}\le e^{3/2-6(1-w)}\) on the half.
For \(k\le4\), the lower bound in (2) of that proof gives
\(c_{m,k}^{1/k}\le e^6(k!)^{5/k}/82^5\).
The resulting bounds at \(6/5\) give
\[
                        \|\log H^{\rm s}_m\|_{6/5}<3/1000.
\tag{8}
\]
Throughout, \(\|\cdot\|_R\) denotes the absolute coefficient norm
at radius \(R\). The branch in (8) is real at positive arguments.

## Direct estimate for the half at positive private exponents

If \(D\le2/5\), use the increment bound with \(\lambda=39/50\)
for \(k\le4\). Since \(\lambda(W_d-1)<11\),
\[
 mQ(k/m)\le \min(2/5,3(n+1)/1000)e^{11k/n}
                     \le(2/5)e^{k/12}.
\tag{9}
\]
The function on the middle right increases to the switch
\(n=397/3\) and then decreases, for \(n\ge82,\ k\le4\).
This follows from \(n^2>11k(n+1)\) before the switch.
The endpoint test for \(eQ(w)/w\) at \(4/m\) and \(1/2\)
gives the tail from \(k=5\) with geometric ratio \(41/100\)
and factorial prefactor \(9/50\), as in (1) of the
[curvature twelve proof](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md#setup-and-a-geometric-coefficient-prefactor).

At positive starting exponents the [first-factor derivative](MULTIVARIATE_ABSOLUTE_CURVATURE_15.md#exact-derivative-at-a-positive-start)
and the [direct head derivative](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md#direct-estimates-in-one-orientation)
apply with actual \(V<8\), \(P=6\), and with \(W_d\) in the
upper quadratic budget. They give
\[
 |\partial_m c_{m,k}|
 \le \frac{c_{m,k}}n\left(2k+\frac{3W_d}{5n}k(k-1)\right),\qquad
   k\le4,\ n\le m\le N.
\]
The bounds (9) already hold for real \(m\) here. At \(R=103/100\),
the head and tail replay gives the curvature twelve direct bounds
\[
 \|\log H^{\rm e}_m\|_R<1,\qquad
 \|z\partial_z\log H^{\rm e}_m\|_R<5/2,\qquad
 \|(z\partial_z)^2\log H^{\rm e}_m\|_R<11,
\]
and the adjacent log and adjacent Euler norms are \(<1/10\)
and \(<2/5\). Indeed a tail coefficient difference is bounded
by the common upper bound on either nonnegative coefficient.
For the head-to-log step the replay uses the integral along
the line segment of the two half-polynomials and the quotient
rule for their Euler derivatives.


## Connected estimate for the half at positive private exponents

If \(D>2/5\), then \(n\ge133\), and \(D/n<q_*=(3/1000)(134/133)\)
and \(W_dD/n<19/500\). Put \(K=\lceil12D\rceil+4\).
Thus \(K\ge9\), \(K\le(49/2)D\), \(K/n<\kappa=3/40\),
and \(K+1<n/2\). The increment bound with \(\lambda=4/5\)
and the endpoint test at \(K/m,1/2\) give the actual tail
\(c_{m,k}\le(1/5)2^{-k}\) for \(k>K\) using the factorial
prefactor. The positive-start derivative calculation at actual
\(V<8,P=6\) gives \(|(\log B_m)'|\le1/n\), where
\(B_m=c_{m,1}\). Also \(B_m/m\) decreases, and (7) at \(k=1\)
gives \(B_m<(27/25)D\) on \([n,N]\).

Use the [finite-difference Gamma estimates](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md#the-long-range-with-sharper-gamma-derivative-estimates)
as in the [connected proof](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md#connected-estimates-in-one-orientation).
The actual relative slope is at most one, so set
\(\alpha=1/[(1-\kappa)n]\) and
\[
 H_2=\frac{W_d}{(1-\kappa)n}\left(1+\frac{K}{(1-\kappa)n}\right)
      <\frac{117W_d}{100n}<\frac9{80},\qquad u=\alpha K<1/12 .
\]
The bounds for differences of order \(r\ge2\) are
\(H_r=(r-2)!\alpha^{r-2}H_2\) and \(2rH_r/n\) for their
degree derivatives. The [charged connected estimate](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_TWELVE.md#long-halves)
with \(\beta=e^{H_2+1/8}H_2\) gives
\[
 \beta K<7/5,\qquad \beta B_m<31/500,\qquad
 \alpha/\beta\le e^{-1/8}/W_d<1/14.
\]
The tree factor is at most \((15/14)^{15}<17/6\); hence the
connected parameter is less than
\[
 \frac{(17/6)(31/500)}{1-u}
        \exp\!\left(\tfrac75(\tfrac12+u/[6(1-u)])\right)<2/5.
\]
With \(p=2/5\), the formal log coefficients through degree \(K\)
are bounded by \(B_m p^{k-1}/k\) and their derivatives by
\(7B_m p^{k-1}/n\). For the resulting truncated logarithm
\(A_m\), the norms at radii \(2,6/5\) are bounded by
\((403/100) B_m,(41/25) B_m\), respectively, and
\(\|A_N-A_n\|_{6/5}<18D/n\).

The exponential and actual tails relative to \(e^{A_m}\) at
\(6/5\) are bounded by
\[
 (3/5)^{K+1} e^{(403/100+41/25)B_m},\qquad
 (3/10)(3/5)^K e^{(41/25)B_m}.
\]
For the first use \(K+1\ge12D+5\) and
\(\log(5/3)>1277/2500\). For the second use \(D\le(K-4)/12\)
and \(K\ge9\). Their sum is below \(17/200\). Thus the
logarithmic correction on \(6/5\) is below \(9/100\), and
\[
       \|\log H^{\rm e}_m\|_{6/5}<2D,\qquad
       \|\log H^{\rm e}_N-\log H^{\rm e}_n\|_{6/5}<47/200.
\tag{10}
\]
The correction and norm statements follow by absolute coefficient
norm multiplication and the log series at one. Both logarithms
are real at positive arguments.

## Joining the estimates

On \(R=103/100\), the heavy and light half log norms divided
by \(n\) are at most \(1/82\) and \((3/1000)/82\), respectively.
Their sum is less than \(\log R\), so Rouché at \(R\) and in
the reciprocal variable confines the normalized root log radii
\(u_m\) to \(|u_m|<\log R\).
The direct bounds and the \(6/5\)-to-\(R\) Cauchy multipliers
\(5/2\) and \(25\) of the [curvature twelve proof](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md#joining-direct-and-connected-estimates)
give first derivative bounds divided by \(n\) of at most
\((5/2)/82\) and \((5/2)(3/1000)/82\) in the two halves.
The corresponding mixed derivative bounds divided by \(n\)
are \(16/100\) and \(25(3/1000)/82\). The adjacent log bounds
are \(47/200\) and \(6/1000\), and the adjacent Euler bounds
at fixed radius are at most \(5/2\) times those respective
bounds. In the direct case those Euler bounds follow from
the stronger direct estimate.

At a common angle, subtract the real radius equations at
degrees \(n,N\). The bounds on the adjacent logs and the
first derivatives yield
\[
 n|u_N-u_n|\le
 \frac{\log R+47/200+6/1000}{1-(5/2)(1+3/1000)/82}<3/10.
\]
Integrating the adjacent Euler difference from the positive
ray bounds the phase difference at fixed radius. To pass to
the actual radii, integrate the mixed derivative against
\(|u_N-u_n|\) and the angular displacement. The resulting
full phase difference is bounded in magnitude by
\[
 \left[\tfrac52(47/200+6/1000)
     +\tfrac3{10}\bigl(16/100+25(3/1000)/82\bigr)\right]\theta
 <\theta .
\]
The two root equations again contradict a common angle.

The ordinary real rays and small terminals \(N\le13\) are
covered as in the [one-private entropy proof](MULTIVARIATE_ONE_PRIVATE_ENTROPY_THIRTEEN.md#matching-the-long-phase-hypotheses).
Independent positive rescaling preserves angles, and rational
Bézout descent gives physical SFC(2) over every characteristic-zero
field.

## Exact scalar replay

The checks below are scalar bounds for the degree groups and
analytic estimates in the proof. They do not enumerate moment
polynomials or exponent profiles.

~~~python

from fractions import Fraction as F
from math import factorial
def lo(x):
 x=F(x);return sum(x**j/factorial(j) for j in range(65))
def hi(x):
 x=F(x);return lo(x)+x**65/factorial(65)/(1-x/66)
def c(a,b):
 assert a<b,(float(a),float(b))
def geom(x,k):
 u=x**k/(1-x)
 return (u,u*(k+x/(1-x)),u*(k*k+2*k*x/(1-x)+x*(1+x)/(1-x)**2))
early=[]
for m0,m1,s,r,den in [(13,31,F(107,20),F(7,25),F(113,20)),(32,82,F(351,62),F(6,25),F(12))]:
 R=hi(F(13,20*m0))
 en=lambda m,k:m/lo(s*(1-F(k,m)))
 c(4/lo(s/2),r);c(hi(1)*en(m1,6)/6,r); c(den**m0,F(factorial(m0)))
 hs=[(R*max(en(m0,k),en(m1,k)))**k/factorial(k) for k in range(1,6)]
 T0,T1,_=geom(R*r,6)
 Ah=sum(hs)+T0
 th=sum(j*hs[j-1] for j in range(2,6))+T1
 lcap=F(8,25)
 c(hi(F(23,10))/32,lcap)
 ls=[(R*hi(1)/den**5)**k*factorial(k)**5 for k in range(1,8)]
 T0,T1,_=geom(R*lcap,8)
 Al=sum(ls)+T0
 tl=sum(j*ls[j-1] for j in range(2,8))+T1
 c((1+Ah)/(1-Al),lo(F(13,20)));c((1+Al)/(1-Ah),lo(F(13,20)))
 upper=(hs[0]+th)/(1-Ah)+(ls[0]+tl)/(1-Al)
 lower=th/(1-Ah)+tl/(1-Al)
 c(upper+lower,1)
 early.append((upper,lower))
c(max(early[0][0]+early[1][1],early[1][0]+early[0][1]),1)
c(F(8),lo(F(52,25)))
assert lo(F(9,2))>82
c(1/lo(F(35,6)),F(3,1000))
c(F(38,3)*F(134,133)/lo(F(35,6)),F(19,500))
d=F(2,5); r=F(103,100)
cs=[(d*hi(F(k,12)))**k/factorial(k) for k in range(1,5)]
t=F(41,100)
c(hi(1)*d*hi(F(4,12))/4,t); c(2*hi(1)/lo(F(35,12)),t)
T=[F(9,50)*x for x in geom(r*t,5)]
H=[sum(c0*r**k*k**j for k,c0 in enumerate(cs,1))+T[j] for j in range(3)]
ah,bh,ch=H
c(ah,1-1/lo(1)); c(bh/(1-ah),F(5,2));c(ch/(1-ah)+(bh/(1-ah))**2,11)
Ed=[sum(c0*r**k*k**j*(2*k+F(3,5)*F(44,3*82)*k*(k-1)) for k,c0 in enumerate(cs,1))/82+T[j] for j in range(2)]
c(Ed[0]/(1-ah),F(1,10));c(Ed[1]/(1-ah)+bh*Ed[0]/(1-ah)**2,F(2,5))


q=F(3,1000)*F(134,133); Wq=F(19,500); kc=F(49,2); kap=F(3,40)
c(q*kc,kap);c(kap/(1-kap),F(1,12))
H2=F(9,80); pref=F(117,100);u=F(1,12)
c(1/(1-kap)**2,pref);c(pref*Wq/d,H2)
beta=pref*hi(H2+F(1,8));b=F(27,25)
c(hi(F(4,5)*Wq/d),b)
c(beta*kc*Wq,F(7,5));c(beta*b*Wq,F(31,500))
c(F(15,14)**15,F(17,6));c(F(14)*1/lo(F(1,8)), F(38,3))
p=F(2,5)
c(F(17,6)*F(31,500)/(1-u)*hi(F(7,5)*(F(1,2)+u/(6*(1-u)))),p)
Ao=F(403,100);Ai=F(41,25); R=F(6,5); lc=F(1277,2500)
c(1/lo(Ao*p),1-2*p);c(1/lo(Ai*p),1-R*p)
c(7*b*R/(1-R*p),18)
c(hi(lc),F(5,3));c((Ao+Ai)*b,12*lc)
c(hi(1+F(4,5)*kc*Wq)/12,F(1,2));c(2*hi(1)/lo(F(35,12)),F(1,2))
c(1/lo(5*lc)+F(3,10)*hi(Ai*b*F(5,12))*F(3,5)**9,F(17,200))
c(1/lo(F(9,100)),1-F(17,200))
c(Ai*b+F(9,100)/d,2)
c(18*q+F(18,100),F(47,200))
lcheck=hi(6)*54/82**5
assert F(24)**5 < 54**4
c(lcheck,F(6,1000000))
c(R/lo(F(3,2)),F(27,100))
ah=R*F(6,1000000)
c(ah/(1-ah)+F(27,100)**5/(1-F(27,100)),1-1/lo(F(3,1000)))
c((F(3,100)+F(47,200)+F(6,1000))/(1-F(5,2)/82-F(5,2)*F(3,1000)/82),F(3,10))
c(F(5,2)*(F(47,200)+F(6,1000))+F(3,10)*(F(16,100)+25*F(3,1000)/82),1)

assert 6*F(9,2*2*81) == F(1,6)
assert lo(8) > 32*82 and lo(9) > 32*133
for n, l, limit in [(82,8,3),(133,9,2)]:
    c((F(3,4)*8 + 6*(F(n,2*(n-1))+F(1,4)+F(l,2)))/n,limit)
c((1+F(2,78))/F(78,82)**2,F(6,5))
c(F(5)**5/(lo(5)*factorial(5)),F(9,50))
c(F(4)**4/(lo(4)*factorial(4)),F(1,5))
c(q*kc+F(1,133),F(1,2))
c(F(20,37)+F(1,4),F(4,5))
c(F(41,78)+F(1,4),F(39,50))
c(Ai*b/12,lc)
c(2*q,F(1,82))
c((1+F(3,1000))/82,F(591,20000))
c(hi(F(591,20000)),F(103,100))
c(25*2*q,F(16,100))
c(F(11,82),F(16,100))
print("Private degree six scalar checks passed.")

~~~
