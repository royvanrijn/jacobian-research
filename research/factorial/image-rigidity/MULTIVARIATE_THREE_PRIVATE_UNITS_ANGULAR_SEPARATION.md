# All-degree angular separation with three private unit rows

Suppose an ordinary two-monomial profile has one private unit row
on one side and two private unit rows on the other. Its shared
exponents may have either orientation. Then its consecutive moment
polynomials have disjoint root arguments in every degree, including
under independent positive rescalings.

The
[one-sided shared reserve theorem](MULTIVARIATE_PRIVATE_ONE_TWO_ONE_SIDED_PHASE.md)
covers shared absolute curvature at least \(79/100\).
The argument below covers lower shared curvature analytically at
every terminal \(N\ge14\). The [terminal-thirteen result](MULTIVARIATE_TERMINAL_THIRTEEN_ALL_PROFILES.md)
and the preceding angular results supply terminals five through
thirteen. At terminal four, independent positive rescalings turn a
hypothetical common argument into a common unit-modulus root, and
the [terminal-four proof](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#terminal-4-for-arbitrary-supports)
uses only scale-invariant inequalities and the resulting common unit
quadratic factor. The both-private angle boxes separate terminal three,
and even real positivity separates terminal two. The bounded-degree checks
below are scalar inequalities on angle intervals. They neither
enumerate exponent profiles nor construct moment polynomials.

## Coefficient bounds and their adjacent change

Orient the side with two private units toward index \(k\). Reserving
one \(k!\) and one \((m-k)!\) cancels the binomial coefficient and
leaves the private factor \(k!\). After positive endpoint
normalization the left and right coefficients are thus bounded by
\[
                     f_{m,k}=\frac{k!}{(m!)^{k/m}},\qquad
                     f_{m,m-k},
\]
respectively. All normalized shared factors are at most one by
log convexity. Let \(S\) be the shared absolute curvature, and
\(V=\sum_{\rm shared}d_i^2/\min(a_i,b_i)\), \(d_i=b_i-a_i\).
Here \(S\le79/100\) gives \(|d_i|\le\min(a_i,b_i)\) and
\(V\le2S\le79/50\). Indeed an unequal row with larger gap
than its smaller exponent has curvature at least \(4/3\).

Let \(J(m)=\int_0^m[1-x\psi'(1+x)]\,dx\).
The integral and trapezoidal bounds on the positive trigamma series
give \(0\le1-x\psi'(1+x)\le\min(1,1/(2x))\). Hence
\[
         J(m)\le(1+\log(2m))/2\qquad(m\ge1/2).
\]
In particular \(J(m)\le12/5\) for \(m\le20\), and
\(J(m)/m^2\le(12/5)/20^2\) for \(m\ge20\) by monotonicity
of the bound divided by \(m^2\).

In either orientation the normalized coefficients at fixed \(k\)
are nonincreasing in real degree \(m\ge k\). Their negative
log derivatives are bounded by
\[
 \begin{split}
 D_L(m,k)&=\frac{k}{m}
                     +\frac{V(k^2/2+k/4)}{m^2},\\
 D_U(m,k)&=\frac{k(k-1)}{2m(m-k+1)}
                     +\frac{J(m)k+V(k^2/2+k/4)}{m^2}.          \tag{1}
 \end{split}
\]
For the private factor, integration by parts gives
\(m\psi(1+m)-\log\Gamma(1+m)=m-J(m)\).
Thus the left private derivative magnitude is bounded by
\(k/m\). The right magnitude is
\(\sum_{h=0}^{k-1}1/(m-h)-k/m+kJ(m)/m^2\), which is
nonnegative and bounded as in (1) by subtracting \(1/m\)
from each summand. For a shared row with smaller exponent
\(\ell\), write \(L_m(w)\) for its endpoint-normalized
logarithm per degree at \(w=k/m\). Its second derivative
is
\[
 L_m''(w)=d^2\int_0^\infty e^{-(a+dw)y}h(y/m)\,dy,\qquad
                            h(t)=t/(e^t-1).
\]
Since \(0\le-h'(t)\le1/2\), this gives
\(0\le L_m''\le d^2/\ell\) and
\(0\le\partial_m L_m''\le d^2/(2m^2\ell^2)\).
In differentiating \(mL_m(k/m)\), the derivative
\(L_m-wL_m'=-\int_0^w t L_m''(t)\,dt\) is nonpositive
and has magnitude at most \(d^2k^2/(2m^2\ell)\).
The remaining term \(m\partial_m L_m(w)\) is nonpositive
by the zero-endpoint Green kernel. The kernel has integral
\(w(1-w)/2\), giving magnitude at most \(d^2k/(4m^2\ell)\).
This proves the shared contributions to (1).

## Root annulus

Put \(C_L=9/10\), \(C_U=27/25\). Split an endpoint-normalized
polynomial as \(P_m(z)=L_m(z)+z^m U_m(1/z)\), with half weight
on the even midpoint. Every root has log radius in
\((-C_U/m,C_L/m)\) for \(m\ge13\).

To see this, we first bound \(\Re U_m\) at radii at most one.
Write \(u_k=[z^k]U_m\).
We have \(k u_k\ge(k+1)u_{k+1}\). For the private factor
the weighted ratio is \((1+1/k)(m!)^{1/m}/(m-k)\)
before the possible even middle halving.
For an oriented shared log-factorial function \(g(x)\) of the
index, \(g''\ge d^2/(mb+1)\) and
\(|g'''|\le d^3/(m^2\ell^2)\), where \(b\) denotes its larger
endpoint, \(d\) its absolute gap, and \(\ell\) its smaller
endpoint. These follow by integral comparison of the polygamma
series. The derivative at \(m/2\) after endpoint normalization
has magnitude at most \(d^3/(24\ell^2)\le S_i/12\) by
averaging Taylor's formula for \(g'\) around the midpoint,
where \(S_i=d^2/b\).
For \(k=1\) the shared factor makes a nonpositive contribution
to the successive log ratio, since
\(S_i/12-(m/2-2)d^2/(mb+1)<0\).
For \(k\ge2\) its possible positive contribution is at most
\(S_i/12\). The maximum of \((1+1/k)m/(m-k)\) for
\(2\le k\le m/2-1\) is at an endpoint by log convexity,
giving at most \(338/165\) for \(m\ge13\).
The private factor \((m!)^{1/m}/m=f_{m,m-1}\) decreases in
\(m\) by (1) with no shared rows. The two exact scalar checks
\[
    \tfrac{26}{12} f_{13,12}<1,\qquad
    \tfrac{338}{165} f_{13,12} e^{79/1200}<1
\]
therefore prove the asserted weighted decrease.

The Young cosine inequality
\(\sum_{k=1}^j\cos(k\theta)/k\ge-1\), applied by Abel
summation to the nonincreasing weights \(k u_k r^k\),
gives \(\Re U_m(r e^{i\theta})\ge1-u_1r\) for \(r\le1\).
[Brown and Yin, equation (1.2)](https://acta.hu/download.phtml?id=2604)
Let \(A_H(r)\) denote the coefficient cap sum for \(H-1\).
The checks below give
\[
 \begin{split}
 1+A_L(e^{C_L/m})&<e^{C_L}(1-u_1 e^{-C_L/m}),\\
 1+A_U(e^{C_U/m})&<e^{C_U}(1-A_L(e^{-C_U/m})).
 \end{split}
\]
Both lower modulus bounds on the right are positive, and hold
throughout the relevant reciprocal disks. Rouché in both
orientations gives the annulus.

## The paired equations on a candidate angle

By the
[both-private angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md),
a common nonreal argument must lie in \(0<\theta<\pi/2\) after
conjugation and in the same-label overlap
\[
              (2k-1)\pi/n<\theta<2k\pi/(n+2),\qquad 4k<n+2 .
\]
These statements apply to the separately normalized polynomials.
Write \(u_m\) for the log modulus of a root on that common ray.

At a fixed radius use the caps of each half on \(e^{C_H/n}\)
for both orders. These can be taken as a common cap \(b_i\)
at each index \(i\). By (1), the negative part of an adjacent
coefficient difference is bounded by the lower cap times the
derivative bound at \(n\), capped at the lower cap. The only
positive change can be an increase of the middle weight by
one half at index \(\lfloor(n+1)/2\rfloor\).
This also accounts for a newly introduced middle term.

For either half let \(h>0\) be a lower bound for its real part
on the angle interval, and let \(I_+,I_-\) bound the positive
and negative parts of its imaginary part. For a complex
variation \(X\) write \(R_X^\pm,I_X^\pm\) for bounds on its
signed real and imaginary parts. Directly from the quotient
formula,
\[
 \begin{split}
 |\Re(X/H)|&\le
  \max(R_X^+,R_X^-)/h+\max(I_X^+,I_X^-)\max(I_+,I_-)/h^2,\\
 \Im(X/H)&\le I_X^+/h+\max(R_X^+I_-,R_X^-I_+)/h^2,\\
 |\Im(X/H)|&\le
  \max(I_X^+,I_X^-)/h+\max(R_X^+,R_X^-)\max(I_+,I_-)/h^2 .
 \end{split}
\]
Apply the first two lines to the adjacent half difference
along the coefficient segment, giving modulus difference
\(E_H\) and an upper argument difference \(D_H\).
Apply the first and third lines to \(X=zH'(z)\), giving
log-radius derivative bounds \(P_H\) on log modulus and
\(C_H^*\) on argument. Take the sums \(P,E,D,C^*\) over
the two orientations. The real root equations give
\[
       n|u_N-u_n|\le
       r_*=\frac{\max(C_L,C_U)/n+E}{1-P/n}.
\]
Comparing the two phases at one radius and then integrating
the argument radius derivative shows that the full phase
difference is at most \(D+C^*r_*/n\).
The exact replay proves this is below \(\theta\) on all
candidate rays. The left half has cap sum below \(1/2\);
the right half has positive real part on these rays.
Hence the full phase difference is above \(-3\pi/2\).
The root equations instead require \(\theta\) modulo
\(2\pi\), which is impossible for \(0<\theta<\pi/2\).
The real rays are excluded by positivity and even-degree
real positivity.

## Scope of the scalar estimates

For \(13\le n\le17\), subdivide each candidate box in
\(\theta/\pi\) into sixteen equal intervals. Machin's
formula and alternating arctangent series bound \(\pi\),
and Taylor with the sine and cosine Lipschitz bounds
gives all signed intervals used above. This is a partition
of scalar analytic inequalities, independent of the
factorial profile.

For \(n=18,19\) and for \(n\ge20\) use
\(h=1-\sum_{i\ge2}b_i\), \(B=\sum_i i b_i\),
\(C=\sum_i i^2b_i\). In these rows
\(P_H=B/h\), \(E_H=E_0/h\),
\(D_H\le (E_1/h+BE_0/h^2)\theta\),
\(C_H^*\le(C/h+B^2/h^2)\theta\). Here \(E_0\) is the
sum of both change bounds, and \(E_1\) weights positive
changes by \(i\) and negative changes at \(i\ge3\) by
\(i/4\). The sine is nonnegative at \(i=1,2\) on the
candidate rays, and \(-\sin x\le x/4\) on the positive ray.

For \(m\ge20\) use the unhalved caps at \(20\) for the
first nine indices, by their fixed-index monotonicity.
The geometric tail rate is \(7/12\). Indeed on an even
half of length \(h\ge10\) the cap at the maximum radius
has endpoint value
\(e^{\max(C_L,C_U)/2}/\sqrt{\binom{2h}{h}}<(7/12)^h\).
The base \(h=10\) is checked and the consecutive ratios
are at most \(1/\sqrt3<7/12\). Log convexity of the cap
extends the bound to the entire half. An odd degree uses
the previous even unhalved caps. The right cap dominates
the left at the same radius by reflection of the strictly
decreasing trigamma across the midpoint. For the adjacent
tail differences use one common cap, since both arrays
have nonnegative entries. The large-order scalar sums and
the bounds in (1) are dominated by the calculation at
\(n=20\), using the displayed bound on \(J(m)/m^2\).

## Exact scalar replay

~~~python
from fractions import Fraction as F
from functools import cache
from math import factorial

V,Loss=F(79,50),F(12,5)
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
assert exp_lo(F(19,5))>40

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

assert F(26,12)*qcap(13,1,True)<1
assert F(338,165)*qcap(13,1,True)*exp_hi(F(79,1200))<1

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
        else:
            norm_loss=F(0)
        coeff.append((max(cur,nxt),cur*min(1,norm_loss),center))
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
    for k,((b,dn,p),(lo,hi)) in enumerate(vals,1):
        if mode=='value':
            d,z=F(0),b
        elif mode=='euler':
            d,z=F(0),k*b
        else:
            d,z=dn,p
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
    rp,rm=parts(c_cos,'normal')
    D=parts(c_sin,'normal')[0]/h+max(rp*im,rm*ip)/h**2
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
            P,E,D,Q=map(sum,zip(*(bounds(c,a,b) for c in values)))
            assert P<n,(n,k,j,'slope',float(P))
            radial=(max(CL,CU)/n+E)/(1-P/n)
            angular=(D+Q*radial/n)/(a*pi_low)
            mx=[max(mx[0],radial),max(mx[1],angular)]
            assert radial<F(1,2),(n,k,j,'radial',float(radial))
            assert angular<1,(n,k,j,'angle',float(angular))
    return mx

for n in range(13,18):
    print(n, *map(float,check(n)))

M,K,Q=20,10,F(7,12)
assert F(7,12)**2 > F(1,3)
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
        assert A<F(1,2)
    N1=sum(k*b[2] for k,b in enumerate(vals,1))+sum(F(k,4)*b[1] for k,b in enumerate(vals,1) if k>=3)+t1
    t2=tail0*(K*K+2*K*Q/(1-Q)+Q*(1+Q)/(1-Q)**2) if large else F(0)
    C=sum(k*k*b[0] for k,b in enumerate(vals,1))+t2
    return B/h,N0/h,N1/h+B*N0/h**2,C/h+B**2/h**2

def global_check(n,large):
    P,E,A,C=map(sum,zip(*(global_bounds(n,side,large) for side in (False,True))))
    assert P<n
    r=(max(CL,CU)/n+E)/(1-P/n)
    assert r<F(1,2),(n,large,float(r))
    phase=A+C*r/n
    assert phase<1,(n,large,float(phase))
    return r,phase

for m in range(13,M):
    geometry(m,False)
geometry(M,True)
for n in range(18,M):
    print('global small',n,*map(float,global_check(n,False)))
print('global large',*map(float,global_check(M,True)))

~~~
