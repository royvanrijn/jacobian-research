# All-degree closure at private exponent totals one and two

**Theorem.** Every ordinary two-monomial profile with private
exponent totals one and two has disjoint consecutive root
arguments in every degree, and satisfies SFC(2) over every
characteristic-zero field. The shared rows can have arbitrary
exponents and orientations.

The angular conclusion follows from the two realizations of
the private total two. The [single quadratic angular refinement](MULTIVARIATE_PRIVATE_LINEAR_QUADRATIC_ALL_DEGREES.md#angular-refinement-without-a-physical-radius-constraint)
and the [three private units angular proof](MULTIVARIATE_THREE_PRIVATE_UNITS_ANGULAR_SEPARATION.md)
give the result. The physical coprimality argument below gives
a separate proof using the signed integer degree difference
and the cancellation of the physical coefficient drift.

The [single quadratic row theorem](MULTIVARIATE_PRIVATE_LINEAR_QUADRATIC_ALL_DEGREES.md)
covers the case where total two belongs to one private row.
Here total two comes from two private unit rows. The
[one-sided shared phase estimate](MULTIVARIATE_PRIVATE_ONE_TWO_ONE_SIDED_PHASE.md)
handles shared curvature at least \(79/100\). The
[known terminal results through thirteen](MULTIVARIATE_TERMINAL_THIRTEEN_ALL_PROFILES.md)
cover the smaller orders. We prove physical coprimality in the
remaining weak range at every lower order \(n\ge13\).

## Bounds in the physical and normalized variables

Orient the side with two private units toward \(k\).
Let \(S\le79/100\) and \(V\) be the sums of
\(d_i^2/\max(a_i,b_i)\) and \(d_i^2/\min(a_i,b_i)\) over
shared rows. Each shared ratio \(d_i/\min(a_i,b_i)\) in
absolute value is at most one, as a larger ratio would
already contribute more than one to \(S\). Therefore
\(V\le2S\le79/50\).

Cancel the binomial against one private unit on each side.
The remaining private coefficient is \(k!\), or the
endpoint-normalized factor \(f_{m,k}=k!/(m!)^{k/m}\).
Log convexity of the shared factors bounds the full
normalized coefficients by \(f_{m,k}\) and \(f_{m,m-k}\)
in the two orientations. We denote these orientations
by \(L\) and \(U\), respectively, when taking the
coefficient halves with a halved middle term.

For fixed index \(k\), both sets of normalized coefficients
decrease with \(m\). On \(n\le m\le N=n+1\), the magnitudes
of their log derivatives are bounded by
\[
 \begin{array}{c|c}
 L & \dfrac{k}{n}+\dfrac{V(k^2/2+k/4)}{n^2}\\[5pt]
 U & \dfrac{k(k-1)}{2n(n-k+1)}
           +\dfrac{L_1(N)k+V(k^2/2+k/4)}{n^2}
 \end{array}                                                    \tag{1}
\]
where \(L_1(m)=\int_0^m(1-s\psi'(1+s))\,ds\).
The shared term follows from the
[linked chord derivative](MULTIVARIATE_PRIVATE_LINEAR_QUADRATIC_ALL_DEGREES.md#exact-residual-and-its-order-derivatives).
For the private terms the log derivatives are
\(-k/m+k L_1(m)/m^2\) on the left and
\(-\sum_{h=0}^{k-1}(m-h)^{-1}+k/m-k L_1(m)/m^2\)
on the right. Integral comparison for the trigamma and its
trapezoidal bound give
\[
    0\le1-s\psi'(1+s)\le\min\{(s+1)^{-1},(2s)^{-1}\}.
\]
Integration split at one gives
\(L_1(m)\le\log 2+(\log m)/2\) for \(m\ge1\).
In particular \(L_1(N)<9/4\) if \(N\le21\).
For \(n\ge20\) we use
\(L_1(n+1)/n^2<(9/4)/20^2\) and
\((S+L_1(n+1))/(n+1)<(79/100+9/4)/21\);
the corresponding bounds with \(\log 2+(\log(n+1))/2\)
are decreasing functions of \(n\).

Let \(G_m\) be the positive endpoint normalization, so that
\(P_m(z)\) is the physical moment evaluated at \(z/G_m\)
and divided by its constant coefficient. The private factor
of \(G_m\) is \((m!)^{1/m}\). If \(\Delta\) is the signed
total degree difference, the same trigamma comparison gives
\[
 \left|\log(G_N/G_n)-\Delta\log(N/n)\right|
                  \le \frac{S+L_1(N)}{nN}.                     \tag{2}
\]
In fact the physical scale's log private derivative is
\(1/m-L_1(m)/m^2\), and each shared correction to its signed
gap divided by \(m\) has magnitude at most \(|d_i|/(2m^2\ell_i)\).
Sum and use \(V/2\le S\) and monotonicity of \(L_1\).

When \(\Delta=0\), the
[physical coefficient transport](MULTIVARIATE_BOTH_PRIVATE_PHYSICAL_COEFFICIENT_TRANSPORT.md#the-smallest-mixed-private-totals)
shows that both orientations' raw coefficients decrease
with order at a fixed index. Their negative log derivative
bounds are
\[
 \begin{array}{c|c}
 L & \dfrac{k}{n}+\dfrac{V(k^2+k)}{2n^2}\\[5pt]
 U & \dfrac{k(k-1)}{2n(n-k+1)}+\dfrac{V(k^2+k)}{2n^2}.
 \end{array}                                                     \tag{3}
\]
Thus the private center shift cancels from the right
coefficient derivative. Comparing physical coefficients
evaluated on the same physical ray will use (3) directly.

## Radii and paired sector estimates

Use \(C_L=9/10\), \(C_U=27/25\). For \(m\ge13\) all
normalized root log radii are in \((-C_U/m,C_L/m)\).
Indeed the [weighted right-half monotonicity](MULTIVARIATE_PRIVATE_ONE_TWO_ONE_SIDED_PHASE.md#monotone-weighted-coefficients-in-the-right-half)
and Young's cosine bound give
\(\Re U_m(r e^{i\theta})\ge1-u_1 r\) for \(r\le1\).
The scalar Rouché tests used below compare \(1+A_L(e^{C_L/m})\)
with \(e^{C_L}(1-b_{U,1}(e^{-C_L/m}))\).
Here \(A_H(r)\) is the mass of the corresponding half
envelope on radius \(r\) and \(b_{U,1}(r)=f_{m,m-1}r\).
For the reverse comparison use \(1+A_U(e^{C_U/m})\) and
\(e^{C_U}(1-A_L(e^{-C_U/m}))\).
Each right-hand side strictly exceeds its left-hand side;
the relevant half is zero-free on the reciprocal disk by
the same positive real-part estimates. Rouché in both
orientations gives the radii.

A common physical nonreal root has normalized arguments
equal, and by the [exact angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
we may take an upper-half-plane angle
\[
    \frac{(2j-1)\pi}{n}<\theta<\frac{2j\pi}{n+2},
                   \qquad 4j<n+2,\quad \theta<\pi/2 .          \tag{4}
\]
For each half let \(b_k\) be the maximum of its two
coefficient caps at the common radius \(e^{C_H/n}\),
including the middle weights. It bounds the coefficients
after applying any of the candidate radii in both orders
and also bounds a straight-line homotopy between the two
radially evaluated arrays. When comparing at a common
normalized radius use (1) to bound the negative part of
the difference \(H_N-H_n\) by the lower coefficient cap
times its log derivative bound. At the center one half
weight can increase by \(1/2\), so the positive part has
only that index and is bounded by \(1/2\) of its unhalved
upper cap. After proving \(\Delta=0\), the same instructions
using (3) give bounds on the difference of the arrays
evaluated at their actual corresponding radii. This uses
the same physical radius, and the compensated coefficients
are compared before endpoint scaling.

For a fixed angle or interval let
\(h=1-\sum b_k\max(0,-\inf\cos k\theta)\). It is a lower
bound on the real part of either half or the segment
between the two evaluated halves. Let \(I_+,I_-\) bound
the positive and negative parts of its imaginary part
by summing the corresponding sine bounds. For a difference
or Euler derivative \(Z\), let \(R_+,R_-,J_+,J_-\)
bound the positive and negative parts of its real and
imaginary components, from its coefficient intervals.
Whenever \(h>0\), the following quotient bounds hold along
the segment:
\[
 \begin{split}
 |\Re(Z/H)| &\le
 \frac{\max(R_+,R_-)}{h}
     +\frac{\max(J_+,J_-)\max(I_+,I_-)}{h^2},\\
 \Im(Z/H)&\le
 \frac{J_+}h+\frac{\max(R_+ I_-,R_- I_+)}{h^2}.
 \end{split}                                                     \tag{5}
\]
Both bounds follow by multiplying by \(\overline H\).
Integrating the quotient for a difference along the
segment bounds the difference of the two logarithms.
For \(13\le n\le16\), the replay applies (5) throughout
the allowed angular intervals (4). All sine and cosine
enclosures use rational Taylor bounds. The subdivision
places no restriction on the actual radius or coefficients.

For \(n\ge17\) we use a uniform version of the quotient
argument with \(h=1-\sum_{k\ge2}b_k\), since the first
real contribution is nonnegative on (4). Put
\(B=\sum k b_k\). Then \(|\Im H|\le B\theta\).
The real modulus difference has bound \(E_0/h\) where
\(E_0\) is the coefficient mass of the absolute difference;
the radial modulus derivative has bound \(B/h\).
For the upper imaginary log difference, use
\((E_1/h+B E_0/h^2)\theta\). Here \(E_1\) weights each
possible positive coefficient difference by \(k\), and
each negative one at \(k\ge3\) by \(k/4\).
Negative differences at \(k=1,2\) cannot make positive
imaginary contributions on this range. The factor \(1/4\)
uses \(-\sin x\le x/4\) for \(x\ge0\), as in the
[single quadratic proof](MULTIVARIATE_PRIVATE_LINEAR_QUADRATIC_ALL_DEGREES.md#paired-estimates-on-a-shared-ray).

## Subtracting the two root equations

Let \(u_n,u_N\) be the normalized root log radii for a
hypothetical common root. The two modulus equations
read \(m u_m=\log|L_m(e^{u_m+i\theta})|
-\log|U_m(e^{-u_m+i\theta})|\).
Let \(P\) be the sum over both orientations of the radial
log modulus bounds, and \(E\) the sum of the adjacent
log modulus difference bounds at the same normalized
radius using (1).
The replay gives \(P<n\) and
\[
          n|u_N-u_n|
          \le\frac{C_U/n+E}{1-P/n}<1/2 .
\]
To see the inequality, subtract the two equations, compare
the orders first at a common radius, and use the radial
derivative for the remaining displacement. If
\(\Delta\ne0\), then (2) and \(\log(N/n)>1/N\) instead
give \(n|u_N-u_n|>(n-S-L_1(N))/N>1/2\). The last inequality
uses the uniform loss bounds following (1).
Consequently \(\Delta=0\).

Use (3) and (5), or its uniform version, to compare
the phases at the actual radii. The replay bounds the
upper difference of their summed phase corrections by
\(3\theta/4\). In both ranges the left half has mass
strictly below \(1/2\) and arguments in
\((-\pi/4,\pi/4)\). The right half has positive real
part at the candidate rays and therefore arguments
in \((-\pi/2,\pi/2)\). The difference of full phase
corrections is strictly greater than \(-3\pi/2\).
The exact split phase equations require that difference
to be \(\theta\) modulo \(2\pi\), which is impossible
between the two bounds. This proves the claimed
coprimality.

## Exact scalar replay

For degrees \(m\ge20\), bound the first nine coefficients
by their unhalved values at \(m=20\), using the private
fixed-index monotonicity in (1) with no shared rows.
The right bound dominates the left on the half, since
the logarithm of their ratio is concave by decreasing
trigamma and vanishes at the two endpoints \(0,m/2\).
The logarithm of the right-hand cap is convex in \(k\),
so its log divided by \(k\) is bounded by its midpoint
value divided by \(m/2\). At even \(m=2h\) this rate
decreases with \(h\): \(f_{2h,h}^2=\binom{2h}{h}^{-1}\)
and the ratios of central binomial coefficients
increase. The odd degrees use the preceding even
degree's fixed-index cap. Thus for \(k\ge10\) at every
radius at most \(e^{C_U/20}\) the geometric cap
\((7/12)^k\) follows from the checked midpoint at \(m=20\).
An absolute difference of two nonnegative coefficients
with the same cap costs at most that cap also at a middle
weight change. This yields the uniform sums of the tail
masses and their first Euler derivatives in the replay.

The value of \(\pi\) is enclosed via Machin's formula and
finite alternating arctangent sums. At each angular cell,
reduce the midpoint by a period before applying the
sine or cosine Taylor bound and use their Lipschitz
constant one for the cell radius. The printed bounds
are summaries of rational tests.

~~~python
from fractions import Fraction as F
from functools import cache
from math import factorial

S,V,Loss=F(79,100),F(79,50),F(9,4)
CL,CU=F(9,10),F(27,25)

def atan_partial(x,count):
    return sum((-1)**k*x**(2*k+1)/(2*k+1) for k in range(count))
a=F(1,5)
b=F(1,239)
pi_low=16*atan_partial(a,20)-4*atan_partial(b,21)
pi_high=16*atan_partial(a,21)-4*atan_partial(b,20)
PI=F(3141592653589793,10**15)
assert pi_low>PI-F(1,10**15) and pi_high<PI+F(1,10**15)
assert F(3)<pi_low<pi_high<F(22,7)

def exp_lo(x):
    return sum(x**k/factorial(k) for k in range(17))
def exp_hi(x):
    return exp_lo(x)+x**17/(factorial(17)*(1-x/18))
assert exp_lo(F(7,10))>2
assert exp_lo(F(29,10))>17
assert F(7,10)+F(29,20)<Loss

@cache
def qcap(m,k,right):
    j=m-k if right else k
    den=10**8
    target=factorial(j)**m * den**m
    weight=factorial(m)**j
    lo,hi=0,den
    while hi-lo>1:
        mid=(lo+hi)//2
        if mid**m * weight >= target:
            hi=mid
        else:
            lo=mid
    assert hi**m*weight>=target
    return F(hi,den)

def half(m,k,right,r):
    if 2*k>m:
        return F(0)
    weight= F(1,2) if 2*k==m else F(1)
    return weight*qcap(m,k,right)*r**k

def half_envelopes(n,right):
    R=exp_hi((CU if right else CL)/n)
    coeff=[]
    for k in range(1,1+(n+1)//2):
        cur=half(n,k,right,R)
        nxt=half(n+1,k,right,R)
        center=F(1,2)*qcap(n+1,k,right)*R**k if k==(n+1)//2 else F(0)
        if cur:
            shared=V*(F(k*k,2)+F(k,4))/n**2
            h=F(k*(k-1),2*n*(n-k+1))
            norm_loss=(h+Loss*k/n**2 if right else F(k,n))+shared
            raw_loss=(h if right else F(k,n))+V*F(k*k+k,2*n**2)
        else:
            norm_loss=raw_loss=F(0)
        coeff.append((max(cur,nxt),cur*min(1,norm_loss),
                      cur*min(1,raw_loss),center))
    if not right:
        assert sum(x[0] for x in coeff)<F(1,2)
    return coeff

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

def parts(vals,mode):
    # positive upper part then negative upper part
    pos=neg=F(0)
    for k,((b,dn,dp,p),(lo,hi)) in enumerate(vals,1):
        if mode=='value':
            d,z=F(0),b
        elif mode=='euler':
            d,z=F(0),k*b
        else:
            d,z=(dn if mode=='normal' else dp),p
        cp,cm=max(0,hi),max(0,-lo)
        pos+=max(z*cp,d*cm)
        neg+=max(z*cm,d*cp)
    return pos,neg

def bounds(c,box_left,box_right):
    c_cos=list(zip(c,[trig_interval(k,box_left,box_right,True) for k in range(1,len(c)+1)]))
    c_sin=list(zip(c,[trig_interval(k,box_left,box_right,False) for k in range(1,len(c)+1)]))
    h=1-parts(c_cos,'value')[1]
    assert h>0,(box_left,box_right,float(h))
    ip,im=parts(c_sin,'value')
    ia=max(ip,im)
    def real_log(mode):
        return max(parts(c_cos,mode))/h+max(parts(c_sin,mode))*ia/h**2
    P=real_log('euler')
    E=real_log('normal')
    rp,rm=parts(c_cos,'raw')
    D=parts(c_sin,'raw')[0]/h+max(rp*im,rm*ip)/h**2
    return P,E,D

def check(n):
    values=tuple(half_envelopes(n,s) for s in [False,True])
    mx=[F(0),F(0)]
    for k in range(1,(n+1)//4+1):
        left,right=F(2*k-1,n),F(2*k,n+2)
        assert left<right
        cells=16
        for j in range(cells):
            a=left+(right-left)*j/cells
            b=left+(right-left)*(j+1)/cells
            P,E,D=map(sum,zip(*(bounds(c,a,b) for c in values)))
            assert P<n,(n,k,j,'slope',float(P))
            radial=(max(CL,CU)/n+E)/(1-P/n)
            angular=D/(a*pi_low)
            mx=[max(mx[0],radial),max(mx[1],angular)]
            assert radial<F(1,2),(n,k,j,'radial',float(radial))
            assert angular<F(1,2),(n,k,j,'angle',float(angular))
    return mx

assert (13-S-Loss)/14 > F(1,2)
for n in range(13,17):
    print(n, *map(float,check(n)))

M,K,Q=20,10,F(7,12)
assert exp_lo(F(31,10))>21
assert qcap(M,M//2,True)*exp_hi(CU/2) < Q**(M//2)
tail0=Q**K/(1-Q)
tail1=tail0*(K+Q/(1-Q))

def mass(m,right,r,large):
    if large:
        return sum(qcap(M,k,right)*r**k for k in range(1,K))+tail0
    return sum(half(m,k,right,r) for k in range(1,m//2+1))

def geometry(m,large):
    outer_left=exp_hi(CL/m)
    outer_right=exp_hi(CU/m)
    inner_left=F(1) if large else 1/exp_lo(CU/m)
    inner_right=F(1) if large else 1/exp_lo(CL/m)
    left_out=mass(m,False,outer_left,large)
    left_in=mass(m,False,inner_left,large)
    right_out=mass(m,True,outer_right,large)
    right_first=qcap(m,1,True)*inner_right
    assert left_in<1 and right_first<1
    assert 1+left_out < exp_lo(CL)*(1-right_first)
    assert 1+right_out < exp_lo(CU)*(1-left_in)

def global_bounds(n,right,large):
    if large:
        R=exp_hi((CU if right else CL)/n)
        t0,t1=tail0,tail1
        vals=[]
        for k in range(1,K):
            b=qcap(n,k,right)*R**k
            base=F(k*(k-1),2*n*(n-k+1)) if right else F(k,n)
            normalized=base+(Loss*k/n**2 if right else 0)+V*F(2*k*k+k,4*n**2)
            physical=base+V*F(k*k+k,2*n**2)
            vals.append((b,b*normalized,b*physical,F(0)))
    else:
        t0,t1=F(0),F(0)
        vals=half_envelopes(n,right)
    h=1-sum(b[0] for b in vals[1:])-t0
    A=sum(b[0] for b in vals)+t0
    B=sum(k*b[0] for k,b in enumerate(vals,1))+t1
    N0=sum(b[1]+b[3] for b in vals)+t0
    D0=sum(b[2]+b[3] for b in vals)+t0
    D1=sum(k*b[3] for k,b in enumerate(vals,1))+sum(F(k,4)*b[2] for k,b in enumerate(vals,1) if k>=3)+t1
    assert h>0
    if not right:
        assert A<F(1,2)
    return B/h,N0/h,D1/h+B*D0/h**2

def global_check(n,large):
    P,E,A=map(sum,zip(*(global_bounds(n,side,large) for side in (False,True))))
    assert P<n
    r=(max(CL,CU)/n+E)/(1-P/n)
    assert r<F(1,2),(n,large,float(r))
    assert A<F(3,4),(n,large,float(A))
    return r,A

for m in range(13,M):
    geometry(m,False)
geometry(M,True)
for n in range(17,M):
    print('global small',n,*map(float,global_check(n,False)))
print('global large',*map(float,global_check(M,True)))
~~~
