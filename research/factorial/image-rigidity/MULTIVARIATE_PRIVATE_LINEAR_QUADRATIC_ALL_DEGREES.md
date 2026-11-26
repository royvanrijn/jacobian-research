# All-degree coprimality with a private linear row and a private quadratic row

Suppose two ordinary monomials have exactly one private variable
each, with private exponents one and two. Their other rows are
arbitrary shared ordinary rows, with either orientation.

**Theorem.** The consecutive root arguments are disjoint in every
degree, and the profile satisfies SFC(2) in every degree over
every characteristic-zero field. There is no restriction on
the shared curvature, dimension, or signed total degree difference.

The [shared angle reserve](MULTIVARIATE_PRIVATE_ONE_TWO_SHARED_ANGLE_RESERVE.md)
handles shared curvature at least \(13/20\). For the lower range
we prove physical coprimality at every lower order \(n\ge13\).
The known smaller terminals, through
[terminal thirteen](MULTIVARIATE_TERMINAL_THIRTEEN_ALL_PROFILES.md),
then give the result. The lower-curvature proof first derives physical coprimality
from the integer total degree difference and linked half
coefficients. A phase refinement then excludes even rays
with different physical radii. The scalar checks do not
enumerate profiles.

## Exact residual and its order derivatives

Orient the degree-two private row toward index \(k\). Write
\[
 S=\sum_{\rm shared}\frac{d_i^2}{\max(a_i,b_i)}\le13/20,\qquad
 V=\sum_{\rm shared}\frac{d_i^2}{\min(a_i,b_i)},\qquad d_i=b_i-a_i.
\]
Every shared row has \(|d_i|\le\min(a_i,b_i)\). In fact if a shared
row has ratio \(|d_i|/\min(a_i,b_i)>1\), its gap is at least two
and its absolute curvature exceeds \(|d_i|/2\ge1\).
Consequently \(V\le2S\le13/10\).

Cancel the binomial against \(k!(m-k)!\), using
\((2k)!=4^k k!(1/2)_k\). After positive endpoint normalization
the remaining private coefficient factor is
\[
                       q_{m,k}=\frac{(1/2)_k}{(1/2)_m^{k/m}}.
\]
Each shared factor has normalized value at most one by log
convexity. The two oriented coefficient bounds are therefore
\(q_{m,k}\) and \(q_{m,m-k}\). In both orientations the actual
coefficients at fixed index \(k\) decrease with \(m\).
We use the following quantitative derivative bounds for real
\(m\ge k\) with \(m>0\). In the two orientations, respectively,
\[
 \begin{split}
  0&\le-\partial_m\log c^{L}_{m,k}
       \le \frac{k}{m}+ \frac{V(k^2/2+k/4)}{m^2},\\
  0&\le-\partial_m\log c^{U}_{m,k}
       \le \frac{k^2}{2m(m-k+1/2)}
                    +\frac{3k}{5m^2}+\frac{V(k^2/2+k/4)}{m^2}.
 \end{split}                                                     \tag{1}
\]

For the private part put
\(\mathcal L(m)=\int_0^m [1-s\psi'(s+1/2)]\,ds\).
The usual integral for the trigamma gives
\(\psi'(s+1/2)=\int_0^\infty e^{-st} h_*(t)\,dt\), where
\(h_*(t)=t/(2\sinh(t/2))\).
For \(t\ge0\), \(1-t^2/24\le h_*(t)\le1\).
Indeed when \(x=t/2<\sqrt6\), the terms of
\(\sinh x/x\) are bounded by those of \(1/(1-x^2/6)\),
since \((2j+1)!\ge6^j\). The bound for larger \(x\)
follows by positivity. Thus
\(0\le1-s\psi'(s+1/2)\le\min(1,1/(12s^2))\).
Its integral on the positive half-line is at most
\(1/\sqrt3<3/5\).
Integrating \(s\psi'(s+1/2)\) by parts shows
\[
 m\psi(m+1/2)-\log( (1/2)_m )=m-\mathcal L(m),
 \qquad (1/2)_m=\Gamma(m+1/2)/\Gamma(1/2).
\]
The derivative of the lower-end factor has magnitude at
most \(k/m\). The negative log derivative of the opposite
factor is
\[
  \sum_{j=0}^{k-1}\frac1{m-j-1/2}
                   -\frac{k}{m}+\frac{k\mathcal L(m)}{m^2},
\]
which is nonnegative. Subtract \(1/m\) from each summand
and use \(\sum_{j=0}^{k-1}(j+1/2)=k^2/2\) to obtain (1).

For a shared row let \(L_m(w)\) be its endpoint-normalized
log symbol per degree. Write \(x=a+(b-a)w\), \(d=b-a\),
and \(\ell=\min(a,b)\). The exact
[trigamma integral](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md#a-gamma-chord-bound-in-both-orientations)
gives \(0\le L_m''(w)\le d^2/\ell\) and
\(0\le\partial_m L_m''(w)\le d^2/(2m^2\ell^2)\), using
\(0\le-h'(t)\le1/2\) for \(h(t)=t/(e^t-1)\).
The derivative of \(mL_m(k/m)\) is
\(L_m(w)-wL_m'(w)+m\partial_m L_m(w)\).
The first two terms combine to \(-\int_0^w s L_m''(s)\,ds\).
For the last term use the zero-endpoint chord Green kernel
whose integral at \(w\) is \(w(1-w)/2\).
Both contributions are nonpositive and their magnitudes
are bounded by \(d^2 k^2/(2m^2\ell)\) and
\(d^2 k/(4m^2\ell)\). The other orientation is identical
with the columns exchanged. This proves (1).

## Exact radii and the signed degree difference

Put \(C_L=13/10\), \(C_U=17/20\). For the normalized
\(P_m\) write \(P_m(z)=L_m^*(z)+z^m U_m^*(1/z)\), with
half of the central coefficient in each half when \(m\)
is even. The scalar checks give the Rouché inequalities
\[
 \begin{split}
   1+A_L(e^{C_L/m}) &< e^{C_L}(1-A_U(e^{-C_L/m})),\\
   1+A_U(e^{C_U/m}) &< e^{C_U}(1-A_L(e^{-C_U/m})).
 \end{split}
\]
Here \(A_H(r)\) is the sum of the coefficient bounds for
\(H-1\) on radius \(r\), including the middle weights.
On each comparison circle the reciprocal term dominates
and has all \(m\) of its zeros inside, since the required
half is zero-free on the reciprocal disk by \(A_H<1\).
Apply the second inequality to the reversal. Every
normalized root has log radius in \((-C_U/m,C_L/m)\).

Let \(\Delta=\sum_i(b_i-a_i)\) include both private rows.
For a common physical root \(z\), write \(u_m=\log|G_m z|\)
with the positive endpoint normalization
\(G_m=([z^m]Q_m/Q_m(0))^{1/m}\).
We claim
\[
 \left|u_N-u_n-\Delta\log(N/n)\right|
                         \le\frac{5}{4nN},\qquad N=n+1.         \tag{2}
\]
The private factor of \(G_m\) is
\(4((1/2)_m)^{1/m}\), so its log order derivative is
\(1/m-\mathcal L(m)/m^2\). A shared row contributes
the order derivative \(\int_a^b t\psi'(1+mt)\,dt\).
Integral comparison and the trapezoidal upper sum for
\(x\mapsto(u+x)^{-2}\) give
\[
                    0\le1-u\psi'(1+u)\le 1/(2u).
\]
The error from replacing the shared integrand by \(1/m\)
is therefore at most \(|d|/(2m^2\ell)\) after integration
on the row. Sum these errors and use \(V/2\le S\).
Integration over \([n,N]\) proves (2) since
\(S+3/5\le5/4\).

## Paired estimates on a shared ray

A candidate common root is nonreal by positivity of the even
moment on the real axis. By the
[both-private angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
we may take its argument \(0<\theta<\pi/2\).
At each \(n\ge13\), in one orientation take the coefficientwise
common bounds \(b_k\) for both halves on radius
\(e^{C_H/n}\), with \(C_H=C_L\) or \(C_U\). On this radius
let \(A=\sum b_k\), \(B=\sum k b_k\), \(C=\sum k^2 b_k\),
and \(h=1-\sum_{k\ge2}b_k\).
On the sector \(0\le\theta\le\pi/2\) the first term has
nonnegative real part, so \(\Re H_m\ge h\) and
\(|\Im H_m|\le B\theta\) for both orders.
The checks give \(h>0\), and \(A<1/2\) for the left half.

In the common half indices the upper normalized coefficient
decreases by (1). The half weight increases at exactly one
possible index \(j=\lfloor N/2\rfloor\), by \(1/2\) including
a newly introduced middle term. Thus the negative part of
the adjacent difference at an index has bound \(b_k^-\)
from the lower half coefficient envelope times the
right-hand side of (1) at \(m=n\), capped at the envelope
itself. Its positive part is bounded by \(b_k^+=0\) off
\(j\) and by half of the unhalved upper coefficient bound
at \(j\). All of these bounds include the radial factors.
Let
\[
 E_0=\sum (b_k^-+b_k^+),\qquad
 E_1=\sum k b_k^+ +\frac14\sum_{k\ge3} k b_k^- .
\]
They bound the modulus of the difference and its imaginary
part from above by \(E_0\) and \(E_1\theta\), respectively.
For the second, the sine is nonnegative at \(k=1,2\).
For \(k\ge3\) use \(-\sin x\le x/4\) for \(x\ge0\):
it is immediate up to \(\pi\) and past four, and between
\(\pi\) and four use \(x-\pi\le x/4\).

The log difference is computed by integrating the
coefficientwise segment from the lower half to the upper
half. Along this segment both real-part and imaginary-part
bounds apply. Writing a log derivative as a quotient gives
\[
 \begin{split}
  |\log |H_N|-\log |H_n|| &\le E_0/h,\\
  \arg H_N-\arg H_n
       &\le (E_1/h+BE_0/h^2)\theta,\\
  |\partial_u \log |H_m||&\le B/h,\qquad
  |\partial_u\arg H_m|\le (C/h+B^2/h^2)\theta.
 \end{split}                                                     \tag{3}
\]
The derivatives are at \(z=e^{u+i\theta}\), with real
arguments normalized at the origin. They also apply to
the opposite half with argument radius \(e^{-u}\).
For the second bound, use
\(\Im(\delta H/H_t)=(\Im\delta H)\Re H_t/|H_t|^2
-(\Re\delta H)\Im H_t/|H_t|^2\) along the segment.
For the last bound replace \(\delta H\) by the Euler
derivative and use \(|\Im(zH_m')|\le C\theta\).

Write \(P,E,D,C_*\) for the sums over the two orientations
of \(B/h,E_0/h,E_1/h+BE_0/h^2,C/h+B^2/h^2\).
At either root the modulus equation is
\(m u_m=\log|L_m^*(e^{u_m+i\theta})|-\log|U_m^*(e^{-u_m+i\theta})|\).
Subtract the equations, compare the two orders at the
same radius, and use the derivative in (3). The scalar
checks give \(P<n\) and
\[
       n|u_N-u_n|
           \le\frac{C_L/n+E}{1-P/n}<1/2.
\]
If \(\Delta\ne0\), equation (2) and
\(\log(N/n)>1/N\) would instead give
\(n|u_N-u_n|>(n-5/4)/N>1/2\). Thus \(\Delta=0\).
The two radii now differ by at most \(5/(4nN)\).

The phase correction at each root is the sum
\(\phi_m=\arg L_m^*(e^{u_m+i\theta})
+\arg U_m^*(e^{-u_m+i\theta})\). Compare first at a
common radius and then move the radius using (3).
The same checks give
\[
            \phi_N-\phi_n
             \le (D+C_*(5/4)/(nN))\theta<3\theta/4 .
\]
The argument of the left half at either radius lies in
\((-\pi/4,\pi/4)\) by \(A<1/2\). The right half has
positive real part, so its argument is in
\((-\pi/2,\pi/2)\). Thus \(\phi_N-\phi_n>-3\pi/2\).
The split root identities require
\(\phi_N-\phi_n\equiv\theta\pmod{2\pi}\), impossible
between these bounds for \(0<\theta<\pi/2\).

## Angular refinement without a physical radius constraint

Suppose two roots have a common argument, allowing arbitrary
positive rescalings of the moments. The angle boxes and the
modulus-equation subtraction still apply. In the weak shared
range at \(n\ge14\), use (3) prior to any appeal to (2).
If \(R=(C_L/n+E)/(1-P/n)\), the normalized radius difference
is at most \(R/n\). At the actual radii the upper phase
difference is therefore at most
\((D+C_*R/n)\theta\). The scalar bounds give
\(D+C_*R/n<4/5\) at every \(n\ge14\).

At \(n=13\), split each of the three possible overlap intervals
\(((2j-1)\pi/13,2j\pi/15)\) into sixteen equal cells. These are
intervals of angles and require no enumeration of factorial
profiles. For a cell, bound the positive and negative parts of
the real and imaginary components of the segment \(H\) and
the numerator \(Z\) using the sine and cosine enclosures.
Write \(h>0\) for the real-part lower bound on \(H\),
\(I_+,I_-\) for its imaginary part bounds, and
\(R_+,R_-,J_+,J_-\) for the real and imaginary bounds on \(Z\).
Then
\[
 \Im(Z/H)\le J_+/h+\max(R_+ I_-,R_- I_+)/h^2 .
\]
Use this for the same-radius adjacent difference. For the
Euler derivative use \(\max(J_+,J_-)\) in the first term
and \(\max(R_+,R_-)\max(I_+,I_-)\) in the second so as to
bound the absolute radial phase derivative. The replay
checks that the two summed bounds, with the Euler bound
multiplied by \(R/13\), total less than \(4\theta/5\) on
every cell. Here the enclosure of \(\pi\) follows from
Machin's formula and the trigonometric enclosures use the
Taylor bound at each reduced midpoint with Lipschitz
constant one for the cell radius. As in the physical
proof, the real-part bounds give a full phase difference
greater than \(-3\pi/2\), excluding the other congruence
possibility on \(0<\theta<\pi/2\).

For orders below thirteen, positivity gives terminal two
and the private angle boxes give terminal three.
At terminal four the
[scale-invariant coefficient proof](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#terminal-4-for-arbitrary-supports)
also excludes a shared argument. If the cubic and quartic
had a common nonreal ray, independent positive rescalings
would place the two roots on the unit circle and give the
same conjugate quadratic. The bounds (1) and (2) in that
proof are scale invariant, as is the quartic positivity
on the real line. Its contradiction therefore applies to
these rescaled polynomials. The existing terminal-five
through terminal-thirteen angular results give the other
smaller orders.

## Uniform tails for the exact scalar replay

For \(m\ge20\) treat the first nine coefficients and
their derivative bounds at \(m=20\). Both unhalved
coefficient envelopes decrease with \(m\) at fixed
index by (1) with \(V=0\). For even \(m=2h\), the right
private factor is log-convex in \(k\) and is one at zero,
so its geometric rate on the half is bounded by its
rate at \(k=h\). The midpoint rate decreases as \(h\)
increases. Indeed the square of the midpoint value is
\((1/2)_h^2/(1/2)_{2h}\), whose consecutive ratios
\((h+1/2)^2/[(2h+1/2)(2h+3/2)]\) decrease.
Odd \(m\) is bounded by its preceding even \(m\) at each
common index. Also \(q_{m,k}\le q_{m,m-k}\) on the
half: the log difference of right minus left is
concave by the decrease of the trigamma and vanishes
at zero and at the midpoint.
At all radii at most \(e^{C_L/20}\) the resulting
geometric rate is below \(11/20\), checked at \(h=10\).
On the tail starting at index ten, bound the absolute
adjacent coefficient difference by the common cap,
since both coefficients are nonnegative. This gives
one contribution of the cap to each of \(E_0,E_1\),
regardless of where the middle weight changes.

The replay takes caps by integer power comparisons and
rational series bounds for exponentials. It needs the
small scalar orders thirteen through nineteen, plus one
uniform estimate for every higher order.

~~~python
from fractions import Fraction as F
from functools import cache
from math import factorial

def exp_lo(x):
    return sum(x**j/factorial(j) for j in range(17))

def exp_hi(x):
    return exp_lo(x) + x**17/factorial(17)/(1-x/18)

@cache
def rising(m):
    return F(factorial(2*m), 4**m * factorial(m))

@cache
def cap(m,k,rev):
    j=m-k if rev else k
    scaled = rising(j)**m * (10**8)**m / rising(m)**j
    lo,hi = 0,10**8
    while hi-lo>1:
        mid=(hi+lo)//2
        if mid**m >= scaled:
            hi=mid
        else:
            lo=mid
    assert hi**m >= scaled
    return F(hi,10**8)

def half(m,k,rev,r):
    if 2*k>m:
        return F(0)
    weight=F(1,2) if 2*k==m else F(1)
    return weight * cap(m,k,rev)*r**k

Cleft,Cright=F(13,10),F(17,20)
S,V,L=F(13,20),F(13,10),F(3,5)
center=S+L
Q,K,M=F(11,20),10,20
assert cap(M,M//2,True)*exp_hi(Cleft/2) < Q**(M//2)
tail0=Q**K/(1-Q)
tail1=tail0*(K+Q/(1-Q))
tail2=tail0*(K**2+2*K*Q/(1-Q)+Q*(1+Q)/(1-Q)**2)

def masses(m,left_radius,right_radius,large):
    sums=[]
    for rev,r in enumerate((left_radius,right_radius)):
        if large:
            sums.append(sum(cap(M,k,rev)*r**k for k in range(1,K))+tail0)
        else:
            sums.append(sum(half(m,k,rev,r) for k in range(1,m//2+1)))
    return sums

def geometry(m,large=False):
    upper=exp_hi(Cleft/m); lower=exp_hi(Cright/m)
    ul,ll=exp_lo(Cleft/m),exp_lo(Cright/m)
    l1,u1=masses(m,upper,F(1) if large else 1/ul,large)
    l2,u2=masses(m,F(1) if large else 1/ll,lower,large)
    assert u1<1 and l2<1
    assert 1+l1 < exp_lo(Cleft)*(1-u1)
    assert 1+u2 < exp_lo(Cright)*(1-l2)

def estimate(n,rev,large):
    R=exp_hi((Cright if rev else Cleft)/n)
    end=K if large else 1+(n+1)//2
    b=[]; neg=[]; pos=[]
    for k in range(1,end):
        if large:
            nxt=cur=cap(M,k,rev)*R**k
        else:
            cur=half(n,k,rev,R)
            nxt=half(n+1,k,rev,R)
        b.append(max(cur,nxt))
        if cur:
            shared=V*(F(k*k,2)+F(k,4))/n**2
            private=(F(k*k,2*n)/(n-k+F(1,2))+L*k/n**2) if rev else F(k,n)
            neg.append(cur*min(F(1),shared+private))
        else:
            neg.append(F(0))
        if not large and k==(n+1)//2:
            pos.append(F(1,2)*cap(n+1,k,rev)*R**k)
        else:
            pos.append(F(0))
    t0,t1,t2=(tail0,tail1,tail2) if large else (F(0),F(0),F(0))
    A=sum(b)+t0
    H=1-sum(b[1:])-t0
    B=sum(k*c for k,c in enumerate(b,1))+t1
    C=sum(k*k*c for k,c in enumerate(b,1))+t2
    E0=sum(neg)+sum(pos)+t0
    E1=sum(k*c for k,c in enumerate(pos,1))+sum(F(k,4)*c for k,c in enumerate(neg,1) if k>=3)+t1
    assert H>0
    if not rev:
        assert A<F(1,2)
    modulus_slope=B/H
    modulus_difference=E0/H
    phase_difference=E1/H+B*E0/H**2
    phase_radius=C/H+B**2/H**2
    return modulus_slope,modulus_difference,phase_difference,phase_radius

def check(n,large=False):
    geometry(n,large)
    P,E,D,C=map(sum,zip(estimate(n,False,large),estimate(n,True,large)))
    assert P<n
    radial=(Cleft/n+E)/(1-P/n)
    phase=D+C*center/(n*(n+1))
    angular=D+C*radial/n
    assert radial<F(1,2), (n,large,'radial',float(radial))
    assert phase<F(3,4), (n,large,'phase',float(phase))
    if n>=14:
        assert angular<F(4,5), (n,large,'angular',float(angular))
    return radial,phase,angular

for n in range(13,M):
    bounds = check(n)
    if n >= 14:
        print(n, *map(float, bounds))
    else:
        print(n, *map(float, bounds[:2]), "pointwise angular replay follows")
print('large',*map(float,check(M,True)))
assert (13-center)/14 > F(1,2)
def atan_partial(x,count):
    return sum((-1)**k*x**(2*k+1)/(2*k+1) for k in range(count))
a=F(1,5)
b=F(1,239)
pi_low=16*atan_partial(a,20)-4*atan_partial(b,21)
pi_high=16*atan_partial(a,21)-4*atan_partial(b,20)
PI=F(3141592653589793,10**15)
assert pi_low>PI-F(1,10**15) and pi_high<PI+F(1,10**15)
assert F(3)<pi_low<pi_high<F(22,7)

@cache
def trig(x,cosine):
    loc=x-2*((x+1)//2)
    y=loc*PI
    parity=0 if cosine else 1
    pol=sum((-1)**j*y**(2*j+parity)/factorial(2*j+parity) for j in range(15))
    # Bound including the argument error from PI.
    err=F(22,7)**30/factorial(30)+F(1,10**15)
    assert err<F(1,10**12)
    return pol

def trig_interval(k,left,right,cosine):
    mid=trig(k*(left+right)/2,cosine)
    error=F(22,7)*k*(right-left)/2+F(1,10**12)
    return max(-1,mid-error),min(1,mid+error)


def small_check(radial):
    n=13
    values=[]
    for rev,C in enumerate((Cleft,Cright)):
        R=exp_hi(C/n)
        coeff=[]
        for k in range(1,8):
            cur=half(n,k,rev,R)
            nxt=half(n+1,k,rev,R)
            shared=V*(F(k*k,2)+F(k,4))/n**2
            private=(F(k*k,2*n)/(n-k+F(1,2))+L*k/n**2) if rev else F(k,n)
            coeff.append((max(cur,nxt),cur*min(1,private+shared),nxt if k==7 else F(0)))
        values.append(coeff)

    def component(coeff,a,b,cosine,kind):
        pos=neg=F(0)
        for k,(bound,minus,plus) in enumerate(coeff,1):
            lo,hi=trig_interval(k,a,b,cosine)
            cp,cm=max(0,hi),max(0,-lo)
            if kind=='value':
                minus,plus=F(0),bound
            elif kind=='euler':
                minus,plus=F(0),k*bound
            pos+=max(plus*cp,minus*cm)
            neg+=max(plus*cm,minus*cp)
        return pos,neg

    for j in range(1,4):
        l,u=F(2*j-1,n),F(2*j,n+2)
        for cell in range(16):
            a=l+(u-l)*cell/16
            b=l+(u-l)*(cell+1)/16
            move=phase=F(0)
            for coeff in values:
                h=1-component(coeff,a,b,True,'value')[1]
                ip,im=component(coeff,a,b,False,'value')
                rp,rm=component(coeff,a,b,True,'difference')
                jp,jm=component(coeff,a,b,False,'difference')
                ep,em=component(coeff,a,b,True,'euler')
                fp,fm=component(coeff,a,b,False,'euler')
                assert h>0
                phase+=jp/h+max(rp*im,rm*ip)/h**2
                move+=max(fp,fm)/h+max(ep,em)*max(ip,im)/h**2
            assert phase+move*radial/n < F(4,5)*a*pi_low

small_check(check(13)[0])
~~~
