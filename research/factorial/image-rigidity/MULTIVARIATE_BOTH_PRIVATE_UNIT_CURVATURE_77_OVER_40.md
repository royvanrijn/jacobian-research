# Angular separation with a private unit pair and shared curvature at least 77/40

Suppose an ordinary two-monomial profile has exactly one private
unit variable on each side. The other unequal rows are shared
and may have either orientation. Write
\[
                 S=\sum_{\rm shared}
                        \frac{(b_i-a_i)^2}{\max(a_i,b_i)} .
\]
If \(S\ge77/40\), consecutive ordinary moments have disjoint
root arguments under independent positive rescalings in every
degree. Consequently these profiles satisfy SFC(2) over every
characteristic-zero field.

The matched tanh values give \(S=\sum 2u/(1+u)\) over shared
Gauss slots, with \(u=(n/T'-(n-1)/T)^{-1}\). This follows by
summing the \(d\) copies of \(u=d/(2\ell+d)\) for a row with
gap \(d\) and smaller exponent \(\ell\).
[Exact matched Gauss grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)

The proof below covers \(77/40\le S<7/2\) analytically for
terminals \(N\ge14\). The
[seven-halves residual theorem](MULTIVARIATE_BOTH_PRIVATE_RESIDUAL_CURVATURE_SEVEN_HALVES.md)
covers \(S\ge7/2\). Terminals five through thirteen follow
from the [complete terminal-thirteen result](MULTIVARIATE_TERMINAL_THIRTEEN_ALL_PROFILES.md)
and the preceding angular results. At terminal four an
angular contact can be positively rescaled to a common
unit-modulus root. The
[terminal-four proof](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#terminal-4-for-arbitrary-supports)
uses only scale-invariant inequalities and the resulting
common unit quadratic factor. The both-private angle boxes
separate terminal three, and even real positivity handles
terminal two.

## Normalized caps and the integer row bounds

The private factorials cancel the binomial factor up to \(m!\).
Normalize endpoints to one. Put
\(V=\sum_{\rm shared}d_i^2/\ell_i\) for smaller exponents
\(\ell_i\) and absolute gaps \(d_i\).
For a cutoff \(s\le S\) and \(k\le m/2\), in either orientation
the normalized coefficient bound is
\[
                      c_{m,k}\le e^{-s k(1-k/m)/2}.             \tag{1}
\]
In one increasing row the second coefficient log difference
is \(\sum_{j=1}^d\log[(x+j)/(x-d+j)]\ge d^2/[m(\ell+d)]\).
Reversal gives the same bound for a decreasing row.
Subtract the endpoint-zero quadratic with that second
difference from the normalized log. The resulting convex
sequence is nonpositive, proving (1) after summing.

We use the four bands below. In each, \(v\) is an upper bound
on \(V\), and \(R_m=e^{C/m}\) is the comparison radius.
\[
\begin{array}{c|c|c}
 \text{shared curvature band}&v&C\\ \hline
 [77/40,2)       &16/3&21/20\\
 [2,9/4)         &35/6&1\\
 [9/4,5/2)       &28/3&4/5\\
 [5/2,7/2)       &35/2&23/20
\end{array}
\]
For the first band a row with \(d>\ell\) must be \((\ell,d)=(1,2)\)
or \((2,3)\). These two types cannot coexist or repeat, by their
curvature costs. In the first case the remainder has curvature
below \(2/3\), and every remaining ratio \(d/\ell\le1\).
Thus \(V\le4+2(2/3)=16/3\). In the second case the remainder
has curvature below \(1/5\) and ratios at most \(1/4\),
since a row with ratio \(r\) costs \(dr/(1+r)\).
Hence \(V\le9/2+(5/4)(1/5)\).
When neither type occurs, \(V\le2S\).
For the second band, if the row \((1,2)\) occurs then
the remainder has curvature below \(11/12\) and all its ratios
are at most one, giving \(V\le4+2(11/12)=35/6\).
Without that row every ratio is at most \(3/2\), since a
ratio above \(3/2\) with \(\ell\ge2\) requires \(d\ge4\)
and costs more than \(3d/5\).
For the third band a ratio above two with \(\ell\ge2\)
requires \(d\ge5\) and costs more than \(2d/3\). At
\(\ell=1\) the only allowed ratio above two is the row
\((1,3)\). When present it consumes \(9/4\), leaving
curvature at most \(1/4\) with ratios at most \(1/3\),
so \(V\le9+(4/3)(1/4)=28/3\).
In the last band every ratio is at most four, so \(V\le5S\).

At fixed index \(k\) the normalized coefficient decreases with
real degree \(m\ge k\). Its negative log derivative is bounded by
\[
                         V(k^2/2+k/4)/m^2 .                    \tag{2}
\]
Indeed for a shared normalized log per degree \(L_m(w)\),
the trigamma integral gives
\(0\le L_m''\le d^2/\ell\) and
\(0\le\partial_m L_m''\le d^2/(2m^2\ell^2)\).
These follow from \(h(t)=t/(e^t-1)\) and
\(0\le-h'(t)\le1/2\). Differentiate \(mL_m(k/m)\).
The term \(L_m-wL_m'\) equals \(-\int_0^w tL_m''(t)\,dt\);
the remaining \(m\partial_m L_m(w)\) is bounded by the
zero-endpoint Green kernel, of integral \(w(1-w)/2\).
Both terms are nonpositive and give (2).
We also use the stronger
[Beta quotient statement](MULTIVARIATE_BOTH_PRIVATE_LINKED_COEFFICIENT_ORDER.md#a-positive-mellin-representation-for-the-quotient)
that \(c_{N,k}/c_{n,k}\) is nonincreasing in \(k\), in
each orientation, for \(N=n+1\).

## The actual root radii

Split \(P_m(z)=L_m(z)+z^m U_m(1/z)\), with half weight at
an even middle index. On \(|z|=R_m\), let \(o_k\) be the
caps from (1) with the outer radial factor \(R_m^k\) and
let \(i_k\) be those with the inner factor \(R_m^{-k}\).
At a ray angle \(\theta\) put \(\cos^\pm x=\max(0,\pm\cos x)\)
and likewise for sine. We have
\[
 \begin{split}
 \Re U_m(R_m^{-1}e^{i\theta})&\ge
                       1-\sum_k i_k\cos^-(k\theta),\\
 |L_m(R_m e^{i\theta})|^2&\le
   \left(1+\sum_k o_k\cos^+(k\theta)\right)^2+
              \max_\pm\left(\sum_k o_k\sin^\pm(k\theta)\right)^2 .
 \end{split}
\]
For the second estimate the signed real part is bounded by
\(1+\sum o_k\cos^+\), and its negative part by
\(\max(0,\sum o_k-1)<1\), using the scalar bound
\(\sum o_k<2\).
A positive lower bound in the first line holds throughout
its disk by using the same caps at all smaller radii.
The scalar replay makes the second line strictly less than
\(e^{2C}\) times the square of that positive lower bound.
This verifies Rouché domination in both orientations and gives
\(|\log|z||<C/m\) for all normalized roots.

## Comparing the two root phases on a common ray

The [both-private angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
give candidate common angles in \(0<\theta<\pi/2\), up to
conjugation. They must lie in a same-label overlap
\[
                   (2k-1)\pi/n<\theta<2k\pi/(n+2),
                                 \qquad 4k<n+2 .               \tag{3}
\]
At a common angle use common radiused caps at radius
\(e^{C/n}\) for both halves. The negative coefficient changes
are bounded by the lower caps times (2) at \(m=n\), capped
by the lower coefficient caps. Positive changes occur only
at the middle index \(\lfloor (n+1)/2\rfloor\) due to the
additional half weight, and are bounded by a half cap there.

The [paired quotient bounds](MULTIVARIATE_THREE_PRIVATE_UNITS_ANGULAR_SEPARATION.md#the-paired-equations-on-a-candidate-angle)
apply to these signed changes on a scalar angle interval.
For each half the real-part lower bound \(h\) is one minus
the cap sum weighted by the negative cosine parts.
Both signed imaginary parts are bounded by summing their
corresponding sine parts. This also bounds either part for
the Euler derivative after multiplying the caps by the index.
Let \(P_H,C_H\) bound the respective log modulus and argument
radius derivatives, \(E_H\) the adjacent log modulus change,
and \(D_H\) the upper adjacent argument change at fixed radius.
Use \(P,E,D,C^*\) for the two-orientation sums.
The real root equations give
\[
 n|u_N-u_n|\le r_*=\frac{C/n+E}{1-P/n},
                  \qquad u_m=\log|z_m|.
\]
Thus the full phase correction difference is at most
\(D+C^*r_*/n\), with the first-box refinement below.

For the first label of (3), temporarily remove the additional
half weight of the upper midpoint. The two half arrays then
have identical index support and weights. By the Beta quotient
statement their coefficient ratios, upper to lower, decrease
in the index. On pairing the indices \(i<j\), the imaginary
part of the product of the adjusted upper half with the
conjugated lower half has a nonpositive summand, because
\(\sin((j-i)\theta)\ge0\). Positive real parts ensure that the
adjusted upper argument is at most the lower one. The upper
bound \(D_H\) in this box can therefore use just the positive
midpoint change. Also, the derivative of the argument of each
individual half with respect to its own log radius is
nonnegative in the first box. Its numerator is the sum of
\((j-i)a_i a_j r^{i+j}\sin((j-i)\theta)\) for \(i<j\).
Since the two half radii are reciprocal, the derivative of
their phase sum is bounded in modulus by the maximum of the
two individual bounds, instead of their sum. We use the same
cap for both orientations, so \(C^*\) is halved in this box.

On all candidate intervals the replay proves positive real
parts and that each signed imaginary part is at most \(2h\).
Each half argument has magnitude at most
\(\arctan 2<3\pi/8\); hence the full phase difference is
strictly above \(-3\pi/2\). The replay gives its upper bound
strictly below \(\theta\), including the indicated first-box
improvement. The root equations require the difference to be
\(\theta\) modulo \(2\pi\), impossible in this interval.
Positive coefficients and even real positivity exclude the
real rays.

## Scope and uniform tail of the scalar checks

For \(m=13,\ldots,29\) use the exact scalar envelopes (1).
For \(m\ge30\) keep the first fourteen unhalved caps from
\(m=30\) at the outer radius. For the inner-radius calculation
use these same coefficients at radius one, bounding the
actual inner radius from above. The common geometric rate
on the tail is \(q=2/3\).
For even \(m=2h\) the logarithm of the unhalved outer cap
is convex in the index and vanishes at zero. The midpoint
bound \(e^{C/2-s h/4}<q^h\) follows at \(h=15\) and thereafter
since \(e^{-s/4}<q\). For odd \(m\) compare the previous even
degree at the same index. This works because the coefficient
cap and the positive radial factor both decrease at fixed
index with \(m\). The same bounds apply to the inner side.
The tail contribution to a coefficient difference uses one
cap because both arrays are nonnegative.

For each small geometry check, and the one uniform large check,
use 512 equal subintervals of \(\theta/\pi\in[0,1]\).
For the small phase comparisons use sixteen equal subintervals
of each candidate box (3). The sine and cosine bounds are
rational intervals obtained from Taylor's formula with
Lipschitz errors and a Machin bound for \(\pi\).
The uniform large phase estimate uses
\(h=1-\sum_{i\ge2}b_i\), \(B=\sum_i i b_i\),
\(C_2=\sum_i i^2 b_i\). Then
\(P_H=B/h\), \(E_H=E_0/h\),
\(D_H\le(E_1/h+BE_0/h^2)\theta\), and
\(C_H\le(C_2/h+B^2/h^2)\theta\).
Here \(E_0\) adds the signed change bounds. The bound
\(E_1\) weights positive changes by \(i\) and negative ones
at \(i\ge3\) by \(i/4\), since the first two sines are
nonnegative on these rays and \(-\sin x\le x/4\).
All these large bounds use the decreasing scalar estimate
(2) from \(n=30\). There are no finite moment tests.

## Exact scalar replay

~~~python
from fractions import Fraction as F
from functools import cache
from math import factorial

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

@cache
def qcap(m,k,right):
    return 1/exp_lo(S*k*(1-F(k,m))/2)

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
            norm_loss=shared
        else:
            norm_loss=F(0)
        coeff.append((max(cur,nxt),cur*min(1,norm_loss),center))
    if not right:
        assert sum(x[0] for x in coeff)<2
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
    for k,((b,dn,p),(lo,hi)) in enumerate(vals,1):
        if mode=='value':
            d,z=F(0),b
        elif mode=='euler':
            d,z=F(0),k*b
        elif mode=='pos':
            d,z=F(0),p
        else:
            d,z=dn,p
        cp,cm=max(0,hi),max(0,-lo)
        pos+=max(z*cp,d*cm)
        neg+=max(z*cm,d*cp)
    return pos,neg

def bounds(c,box_left,box_right,first):
    c_cos=list(zip(c,[trig_interval(k,box_left,box_right,True) for k in range(1,len(c)+1)]))
    c_sin=list(zip(c,[trig_interval(k,box_left,box_right,False) for k in range(1,len(c)+1)]))
    h=1-parts(c_cos,'value')[1]
    assert h>0,(box_left,box_right,float(h))
    ip,im=parts(c_sin,'value')
    ia=max(ip,im)
    assert ia < 2*h
    def real_log(mode):
        return max(parts(c_cos,mode))/h+max(parts(c_sin,mode))*ia/h**2
    P=real_log('euler')
    E=real_log('normal')
    mode='pos' if first else 'normal'
    rp,rm=parts(c_cos,mode)
    D=parts(c_sin,mode)[0]/h+max(rp*im,rm*ip)/h**2
    Q=max(parts(c_sin,'euler'))/h+max(parts(c_cos,'euler'))*ia/h**2
    return P,E,D,Q

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
            P,E,D,Q=map(sum,zip(*(bounds(c,a,b,k==1) for c in values)))
            assert P<n,(n,k,j,'slope',float(P))
            radial=(max(CL,CU)/n+E)/(1-P/n)
            angular=(D+Q*radial/(2*n if k==1 else n))/(a*pi_low)
            mx=[max(mx[0],radial),max(mx[1],angular)]
            assert radial>0
            assert angular<1,(n,k,j,'angle',float(angular))
    return mx


M,K=30,15
Q=F(2,3)
tail0=Q**K/(1-Q)
tail1=tail0*(K+Q/(1-Q))

@cache
def cosine_grid(k,t):
    return trig_interval(k,F(t,512),F(t+1,512),True)
@cache
def sine_grid(k,t):
    return trig_interval(k,F(t,512),F(t+1,512),False)

def geometry(m,large):
    out=exp_hi(CL/m)
    inn=F(1) if large else 1/exp_lo(CL/m)
    if large:
        o=[qcap(m,k,False)*out**k for k in range(1,K)]
        v=[qcap(m,k,False)*inn**k for k in range(1,K)]
        tail=tail0
    else:
        o=[half(m,k,False,out) for k in range(1,m//2+1)]
        v=[half(m,k,False,inn) for k in range(1,m//2+1)]
        tail=F(0)
    assert sum(o)+tail<2
    e=exp_lo(CL)
    worst=F(0)
    for t in range(512):
        cs=[cosine_grid(k,t) for k in range(1,len(o)+1)]
        ss=[sine_grid(k,t) for k in range(1,len(o)+1)]
        h=1-sum(b*max(0,-lo) for b,(lo,hi) in zip(v,cs))-tail
        ar=1+sum(b*max(0,hi) for b,(lo,hi) in zip(o,cs))+tail
        ai=max(sum(b*max(0,hi) for b,(lo,hi) in zip(o,ss)),
               sum(b*max(0,-lo) for b,(lo,hi) in zip(o,ss)))+tail
        assert h>0
        q=(ar*ar+ai*ai)/(e*h)**2
        assert q<1,(m,large,t,float(q))
        worst=max(worst,q)
    return worst


def global_bounds(n,right,large):
    if large:
        R=exp_hi((CU if right else CL)/n)
        t0,t1=tail0,tail1
        vals=[]
        for k in range(1,K):
            b=qcap(n,k,right)*R**k
            normalized=V*F(2*k*k+k,4*n**2)
            vals.append((b,b*normalized,F(0)))
    else:
        t0,t1=F(0),F(0)
        vals=half_envelopes(n,right)
    h=1-sum(b[0] for b in vals[1:])-t0
    A=sum(b[0] for b in vals)+t0
    B=sum(k*b[0] for k,b in enumerate(vals,1))+t1
    N0=sum(b[1]+b[2] for b in vals)+t0
    assert h>0
    if not right:
        assert A<2*h
    N1=sum(k*b[2] for k,b in enumerate(vals,1))+sum(F(k,4)*b[1] for k,b in enumerate(vals,1) if k>=3)+t1
    t2=tail0*(K*K+2*K*Q/(1-Q)+Q*(1+Q)/(1-Q)**2) if large else F(0)
    C=sum(k*k*b[0] for k,b in enumerate(vals,1))+t2
    return B/h,N0/h,N1/h+B*N0/h**2,C/h+B**2/h**2

def global_check(n,large):
    P,E,A,C=map(sum,zip(*(global_bounds(n,side,large) for side in (False,True))))
    assert P<n
    r=(max(CL,CU)/n+E)/(1-P/n)
    phase=A+C*r/n
    assert phase<1,(n,large,float(phase))
    return r,phase


for S,V,C in [
    (F(77,40),F(16,3),F(21,20)),
    (F(2),F(35,6),F(1)),
    (F(9,4),F(28,3),F(4,5)),
    (F(5,2),F(35,2),F(23,20)),
]:
    CL=CU=C
    qcap.cache_clear()
    assert 1/exp_lo(S/4)<Q
    assert qcap(M,M//2,True)*exp_hi(C/2)<Q**(M//2)
    worst=max(geometry(m,m==M) for m in range(13,M+1))
    alpha=max(check(n)[1] for n in range(13,M))
    large=global_check(M,True)[1]
    print(float(S),float(worst),float(alpha),float(large),flush=True)

~~~
