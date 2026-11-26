# All-degree angular separation under exponent-ratio balance

Consider an ordinary two-monomial profile with all unequal rows
shared. Write each unordered pair of row exponents uniquely as
\((t a,t b)\) with \(1\le a<b\) coprime and \(t\ge1\) integral.
For each fixed \((a,b)\), assume the sum of \(t\) over rows
oriented from \(t a\) to \(t b\) equals the sum over rows in
the reverse orientation. Allow a common physical Gaussian
parameter \(\tau\ge0\), and let \(S\) be the usual absolute
curvature including \(\tau\).

**Theorem.** If \(S\ge12\), consecutive polynomials have disjoint
root arguments in every degree. At \(\tau=0\) they satisfy SFC(2)
over every characteristic-zero field.

The condition permits unequal rows without any exact mirror pair.
For instance \((1,5),(1,5),(10,2)\) satisfies the theorem with
\(S=64/5\) and \(\tau=0\).
We can assume \(12\le S<16\) by the
[absolute curvature sixteen theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md).
This proof uses the
[matched tanh phase estimates for mirrored curvature](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_TWELVE.md)
and an exact bound for the change of the normalized Gamma chord
under splitting. It accounts for radial motion in the original
asymmetric array.

## Splitting and the linked Gamma derivative

For a primitive row with slope \(d=b-a\) in its chosen
orientation, let \(L_m(w)=m^{-1}[\log\Gamma(1+m(a+dw))]_{\rm chord}\),
where the subscript subtracts linear interpolation of the
endpoint values. The exact
[trigamma chord formula](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md#a-gamma-chord-bound-in-both-orientations)
is
\[
 L_m''(w)=d^2\int_0^\infty e^{-(a+dw)y}h(y/m)\,dy,\qquad
 h(t)=t/(e^t-1).
\]
Here \(0\le-h'(t)\le1/2\): \(h\) is decreasing, and
\(h(t)+t/2=(t/2)\coth(t/2)\) is increasing since
\(\sinh u\cosh u>u\) for \(u>0\).
Consequently the curvature increases with \(m\); the chord
Green kernel then shows \(L_m\ge L_{tm}\).
Replacing a row of scale \(t\) by \(t\) copies of its primitive
row increases its normalized coefficient at every index,
since \(t m L_{tm}(w)\) is replaced by \(t m L_m(w)\).
As \(h(y/(tm))-h(y/m)\le (1-1/t)y/(2m)\), the increase
in the log coefficient is at most
\[
     \frac{(t-1)d^2}{4\ell^2}\,w(1-w),\qquad \ell=\min(a,b).
                                                               \tag{1}
\]
Indeed the Green kernel integrates to \(w(1-w)/2\).
The same facts apply in either orientation.

The primitive split reference is exactly mirrored by hypothesis.
Both \(S\) and \(V\) are unchanged by splitting, where
\(V=\tau+\sum d_i^2/\min(a_i,b_i)\) over the original rows.
The largest relative gap \(L_0\) is also unchanged. Use \(q,Q_M,Q_E\)
from the reference, as in the curvature twelve proof,
and \(D=Nq\), \(W=1+V\), \(N=n+1\). Both oriented entropy
envelopes of the ordinary array have mass \(q\):
the chord slopes at zero are linear under rescaling a row,
and the \(G/2\) correction vanishes there. The reference gives
\(W^2q\le(1+S)^2e^{-S/2}\), \(Wq\le(1+S)e^{-S/2}\),
and \(L_0q,q\le e^{-S/2}\).
All reference coefficient upper bounds apply to each actual
half. In particular, if the reference inclusive coefficient is
\(c^0_{m,k}\), the two oriented coefficients satisfy
\[
              |c^+_{m,k}-c^-_{m,k}|
                         \le c^0_{m,k}\,\frac{Vk}{4m}.           \tag{2}
\]
Each actual log coefficient differs from the reference
by an amount in \([-Vk/(4m),0]\), by (1) and
\(\sum(t-1)d^2/\ell^2\le V\). The exponential is
\(1\)-Lipschitz on the nonpositive half-line, which gives (2).

For the exact first coefficient \(B_m\) in either orientation,
the linked derivative is bounded by
\[
                 0\le1/m-(\log B_m)'\le\frac{3V}{4m^2}.
                                                               \tag{3}
\]
In a shared row its log contribution is \(mL_m(1/m)\),
now using the actual endpoints and slope. With \(w=1/m\),
the derivative is \(L_m-wL_m'+m\partial_m L_m\).
The first part is \(-\int_0^w uL_m''(u)\,du\), with
magnitude at most \(d^2w^2/(2\ell)\).
Also \(0\le\partial_m L_m''\le d^2/(2m^2\ell^2)\)
by the same kernel derivative bound, so
\(-m\partial_m L_m\le d^2w/(4m\ell^2)\).
Summation, the binomial first factor \(m\), and the
Gaussian contribution \(-\tau/(2m^2)\) prove (3).
These estimates hold for real \(m\) in the intervals used
to compare consecutive coefficients.

## Paired phases at their actual radii

Write the normalized polynomial as \(L_m(z)+z^m U_m(1/z)\),
with half middle weight where needed. On a disk where the
two halves have analytic logs normalized to zero at zero,
write
\[
 J_m(u,\theta)=\log L_m(e^{u+i\theta})-
                 \log U_m(e^{-u-i\theta}).
\]
A root with log radius \(u_m\) satisfies
\(m u_m=\Re J_m(u_m,\theta)\) and
\(m\theta=\Im J_m(u_m,\theta)+(2j+1)\pi\).
Thus a common \(\theta\in(0,\pi)\) is impossible if
\[
     |\Im(J_N(u_N,\theta)-J_n(u_n,\theta))|<\theta,              \tag{4}
\]
because the phase difference would be strictly between
\(0\) and \(2\pi\). Positive arguments contain no roots by
positivity. At the negative argument the logs are real,
so the root-order parity must be odd. Reflection handles
the lower half-plane.

For bounds uniform on \(|u|\le u_0\), suppose the two
logarithmic radial derivatives have modulus at most \(P\)
each, their mixed radial and angular derivatives have modulus
at most \(C\), and each difference of adjacent half logs
has norm at most \(\delta\) on the disk of radius \(6/5\).
At a common argument, the real equations give
\[
           n |u_N-u_n|\le\frac{u_0+2\delta}{1-2P/n}.            \tag{5}
\]
Use the absolute coefficient norm on \(6/5\).
When \(e^{u_0}\le101/100\), multiplication of log
coefficients by \(k(101/120)^k\) bounds the angular
derivative on the smaller disk by \(9/4\) times this norm,
and \(k^2(101/120)^k\) bounds the second derivative
by \(20\) times the norm. These constants follow by
maximizing \(x e^{-cx}\) and \(x^2 e^{-cx}\) on \(x\ge0\),
with \(c=\log(120/101)>17/100\).
At fixed \(u_N\), integrate the paired angular derivative
bound from \(\theta=0\), where the imaginary part vanishes.
For the change from \(u_N\) to \(u_n\), use the same
integration on each mixed derivative. This bounds the
left side of (4) by
\[
                       \left(\frac92\delta+
                                  2 C |u_N-u_n|\right)\theta.   \tag{6}
\]

## Short halves

Use the reference discrete chord bound
\(c_{m,1}\le m e^{-6(1-1/m)}\).
For \(N\le28\) the coefficient of degree one in the half
is below \(1/20\) at \(m=2\) and below \(11/200\)
on \(3\le m\le11\). It is below \(89/1000\) on
\(12\le m\le28\). On each interval the bound on the
inclusive coefficient takes its maximum at an endpoint,
and at \(m=2\) the middle coefficient is halved.
In the same range \(N\le28\), the higher coefficients
are at most \(5^{-k}\). Use the reference entropy endpoint
maxima and the test \(14e Q_E(1/14)<1/5\).
To check that test use \(Q_E(w)\le q e^{(V/4+S/2)w}\)
and \(W^2q\le(1+S)^2 e^{-S/2}\).
The resulting bound at \(w=1/14\) is decreasing in \(S\ge12\)
at fixed \(V\), by its logarithmic derivative.
At \(S=12\) it takes its maximum at an endpoint of
\(12\le V\le160\), by the increasing logarithmic derivative
in \(V\). The midpoint test is \(4Q_E(1/2)\le4e^{-3}<1/5\).

For \(N>28\) with \(D\le7/100\), define \(a(S)=(1+S)e^{-(S-12)/4}/208\).
The same quadratic endpoint budget gives
\[
              L_*/n\le\min(20/7,a(S)/\sqrt D),\qquad L_*=V/2.
\]
In particular \(a(S)\le1/16\). Put \(D_0=(7/320)^2\).
The two pieces about \(D_0\) give
\(D e^{L_*/n}<89/1000\), testing \(D_0\) and \(7/100\)
on the second piece. For \(k\ge2\) use the reference entropy
estimate and again the midpoint bound. At \(k=2\) it
suffices that \(D e^{1+L_*/n+S/n}/2<1/5\) because
\(2L_E=L_*+S\). For \(D\le D_0\) use \(20/7\) and \(S/n\le16/28\).
Beyond \(D_0\), use \(S/28+a(S)/\sqrt D\) for the
last two exponent terms. This is convex in \(S\) on
\([12,16]\) and takes its maximum at \(12\) or \(16\).
For either \(S\), the resulting function
\(D e^{a(S)/\sqrt D}\) takes its maximum at a
\(D\)-endpoint, since its logarithmic derivative has at
most one sign change from negative to positive.
The scalar checks give the bound.

Set \(R_m=e^{8/(25m)}\). The short caps imply both
\(|H-1|\le A\) and \(|zH'|\le B\) on this radius in
each half with \((1+A)/(1-A)<e^{8/25}\) and
\(B/(1-A)<1/4\), by two geometric sums.
The function \(z^m U_m(1/z)\) on \(|z|=R_m\)
dominates \(L_m\) and has all \(m\) zeros inside,
since \(U_m\) has no zeros on or within \(|z|=1/R_m\).
Rouché and the reverse orientation confine
\(|u_m|<\log R_m\).
At the actual radii, each of the four angular derivatives
entering (4) is below \(1/4\) in modulus.
The short comparison therefore follows by integration from
the positive radius.

## Medium halves

Assume \(N>28\), \(7/100<D\le12/25\).
All three mass budgets in the curvature twelve medium
argument still hold for \(S<16\), because
\((1+S)e^{-S/2}\) and \(e^{-S/2}\) decrease.
The real Gamma coefficient bounds from the split reference
remain valid. The derivative proof there uses only positive
Gamma arguments at least \(ma_i/2\) in each row and
\(m/2\) for the binomial, which hold for the actual
array since \(7L_0/n<1/2\). The first-factor bound
\(|(\log B_m)'|\le1/n\) follows from (3).
Thus the linked medium calculation gives
\(\|H_m-1\|_{6/5}<17/20\) and
\(\|\log H_N-\log H_n\|_{6/5}<2/25\) in each
orientation. The logarithmic bound holds in absolute
coefficient norm by the inverse series for the segment
between the two halves.

At a fixed degree use (2) and retain \(WD/n<167/5000\)
as in the medium paired derivative calculation.
The first seven coefficient differences have norm at most
\((167/5000)(27/10)/4\) on \(6/5\).
For the tail use the single coefficient envelope
\(\sum_{k\ge8}(9/25)^k\) because both coefficients are
between zero and that bound. Dividing by \(1-17/20\)
on the segment between the halves gives
\(\|\log L_m-\log U_m\|_{6/5}<31/200\).

On radius \(101/100\), the same coefficient bounds give
\(P=8/3\), \(C=12\). In detail if \(A,B_1,B_2\) bound
the coefficient norm of \(H-1\) and its first and second
Euler derivatives, respectively, the log derivative
norms are bounded by \(B_1/(1-A)\) and
\(B_2/(1-A)+(B_1/(1-A))^2\).
With \(u_0=\log(101/100)\), the modulus imbalance between
the two terms at the boundary excluding the factor \(z^m\)
is bounded in logarithm by \(31/200+2Pu_0<28u_0\).
Here use real coefficients to compare both half arguments
on the same ray and integrate the radial derivative.
Rouché in both orientations confines all roots to
\(|u_m|<u_0\).
Then (5) with \(\delta=2/25\) gives \(n|u_N-u_n|<11/50\),
and (6) gives (4).

## Long halves

Assume \(D>12/25\).
Use the curvature twelve long-half proof in both
orientations, with identical budgets and the same
continued coefficients for real \(m\). The virtual
mirror gives the tail inequality, \(B_m<(207/200)D\)
and the real interpolation envelope. The actual
derivative bound \(|(\log B_m)'|\le1/n\) follows from (3),
and the sharper Gamma bounds only require the same
\(V,L_0,n,K\) conditions. The differentiated connected
bound and the weighted tails therefore give
\(\| \log H_m\|_{6/5}\le M_m=(3/2)B_m+7/250\) and
\(\|\log H_N-\log H_n\|_{6/5}<1/10\).
Put \(M=\max M_m\) over both orientations and orders.
Then \(M/n<41/10000\).
At the radius \(101/100\) the two half log moduli
differ by at most \(2M<m\log(101/100)\).
Rouché and reversal again give \(|u_m|<\log(101/100)\).
The weighted derivative bounds now give
\(P=(9/4)M\), \(C=20M\).
With \(\delta=1/10\), (5) again gives
\(n|u_N-u_n|<11/50\), and (6) proves (4).
Together these cases establish the theorem.

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

def half(radius,c,r):
    x=radius*r
    return radius*c+x*x/(1-x), radius*c+x/(1-x)**2-x

U=1/lo(6)
for V in [12,160]:
    q=U*F(13,1+V)**2
    L=F(V,4)+6
    ck('second at 28',14*q*hi(1+L/14),F(1,5))
ck('first at 3',3/lo(4),F(11,200))
ck('first at 11',11/lo(F(60,11)),F(11,200))
ck('first at 12',12/lo(F(11,2)),F(89,1000))
ck('first at 28',28/lo(F(81,14)),F(89,1000))
limit=F(20,7)
d0=(F(1,16)/limit)**2
d=F(7,100)
ck('small cap cutoff',d0*hi(limit),F(89,1000))
ck('small cap at d',d*hi(F(6,25)),F(89,1000))
ck('small cap at d radical',F(1,16)**2,d*F(6,25)**2)
ck('small cap k2 near zero',d0*hi(1+F(16,28)+limit)/2,F(1,5))
for S,a in [(12,F(1,16)),(16,F(17,13)/(16*lo(1)))]:
    c=F(6,25) if S==12 else F(23,200)
    ck('radical',a*a,d*c*c)
    ck('small cap k2 at D0',d0*hi(1+F(S,28)+a*16*limit)/2,F(1,5))
    ck('small cap k2 at d',d*hi(1+F(S,28)+c)/2,F(1,5))

rho=F(8,25)
for m,c,r in [(2,F(1,20),0),(3,F(11,200),F(1,5)),(12,F(89,1000),F(1,5))]:
    A,B=half(hi(rho/m),c,r)
    ck('short amplitude',F(1)+A,lo(rho)*(F(1)-A))
    ck('short angular slope',B/(1-A),F(1,4))

dmin=F(7,100); dmax=F(12,25); A0=F(21,1250)
Wq=F(167,5000)
z=F(101,100)
def b(k,j,d,r=z):
    return r**k*d**(k-j)*hi(A0*k*k/d)/factorial(k)

t=F(3,10)*z
s=t**8/(1-t)
A=sum(b(k,0,dmax) for k in range(1,8))+s
B=sum(k*b(k,0,dmax) for k in range(1,8))+s*(8+t/(1-t))
C=sum(k*k*b(k,0,dmax) for k in range(1,8))+s*(64+16*t/(1-t)+t*(1+t)/(1-t)**2)
ck('medium slope',B/(1-A),F(8,3))
ck('medium mixed slope',C/(1-A)+(B/(1-A))**2,12)
t12=F(9,25)
tail=t12**8/(1-t12)
m1=sum(k*max(b(k,1,dmin,F(6,5)),b(k,1,dmax,F(6,5))) for k in range(1,8))
ck('medium first sum',m1,F(27,10))
ck('medium same order skew',(Wq/4*F(27,10)+tail)/F(3,20),F(31,200))

# Lower and upper log z by the first two terms of the power series.
umin=F(199,20000); umax=F(1,100)
ck('medium radius',F(31,200)+2*F(8,3)*umax,28*umin)
ck('medium radial difference',(umax+2*F(2,25))/(1-2*F(8,3)/28),F(11,50))
ck('medium angular comparison',2*F(9,4)*F(2,25)+2*12*F(11,50)/28,1)
ck('slope loss log',hi(F(17,100)),F(120,101))
ck('slope one',1/(lo(1)*F(17,100)),F(9,4))
ck('slope two',4/(lo(1)*F(17,100))**2,20)
ck('long norm',F(3,2)*F(207,200)/400+F(7,250)/193,F(41,10000))
ck('long radius',2*F(41,10000),umin)
ck('long radial difference',(umax+F(1,5))/(1-2*F(9,4)*F(41,10000)),F(11,50))
ck('long angular comparison',2*F(9,4)/10+2*20*F(41,10000)*F(11,50),1)
print('All ratio-balanced asymmetric scalar checks passed.')
~~~
