#!/usr/bin/env python3
"""Exact Laguerre model and universal leading resultant coefficients.

No all-order coefficient positivity or unrestricted trinomial theorem is certified.
"""

from collections import defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, factorial

from flint import fmpq, fmpq_mpoly_ctx, fmpz_poly


@lru_cache(None)
def coupling(i, j, k):
    answer = 0
    for t in range(min(i, j, k) + 1):
        z = (i+j-k-t, i+k-j-t, j+k-i-t)
        if min(z) < 0 or any(v % 2 for v in z):
            continue
        p, q, r = [v//2 for v in z]
        n = p+q+r+t
        answer += factorial(n)*2**t // (factorial(p)*factorial(q)*factorial(r)*factorial(t))
    return answer


def poly(z):
    return fmpz_poly([z.get(i, 0) for i in range(max(z, default=0)+1)])


def filtered_moments(a, d, left, right):
    vec = {}
    for i in range(left+d+1):
        for j in range(right+d+1):
            def bn(n, k):
                return comb(n, k) if k <= n else 0
            p = bn(left, i)*bn(right+d, j)
            q = bn(left+a, i)*bn(right+d-a, j)
            r = bn(left+d, i)*bn(right, j)
            u, v = p-r, q-r
            if u or v:
                vec[i, j] = u, v
    moments = [defaultdict(int) for _ in range(7)]
    aa, bb, cc, dd, ee, ff, gg = moments
    for k, (u, v) in vec.items():
        aa[2*sum(k)] += u*u
        bb[2*sum(k)] += 2*u*v
        cc[2*sum(k)] += v*v
    for x, y, z in combinations_with_replacement(vec, 3):
        w = coupling(x[0], y[0], z[0])*coupling(x[1], y[1], z[1])
        if not w:
            continue
        w *= 1 if x == z else 3 if x == y or y == z else 6
        ui, vi = vec[x]
        uj, vj = vec[y]
        uk, vk = vec[z]
        n = sum(x)+sum(y)+sum(z)
        dd[n] += w*ui*uj*uk
        gg[n] += w*vi*vj*vk
        ff[n] += w*(ui*vj*vk+vi*uj*vk+vi*vj*uk)
        ee[n] += w*(ui*uj*vk+ui*vj*uk+vi*uj*uk)
    return tuple(map(poly, moments))


def filtered_resultant(a, d, left, right):
    aa, bb, cc, dd, ee, ff, gg = filtered_moments(a, d, left, right)
    ss = cc*cc*dd+aa*bb*gg-aa*cc*ff
    tt = cc*cc*ee-bb*cc*ff+(bb*bb-aa*cc)*gg
    result, remainder = divmod(cc*ss*ss-bb*ss*tt+aa*tt*tt, cc*cc)
    if remainder:
        raise ValueError('Incorrect resultant division')
    return result


RING = fmpq_mpoly_ctx.get(('L', 'H', 'a', 'b'))
L, H, A, B = RING.gens()
ZERO = RING.constant(0)
ONE = RING.constant(1)
EXPONENTS = ((L,H+A+B), (L+A,H+B), (L+A+B,H))


@lru_cache(None)
def binomial(column, row, n):
    out=ONE
    for j in range(n):
        out=out*(EXPONENTS[column][row]-j)/(j+1)
    return out


def symbolic_moments(limit):
    moments=[[ZERO for _ in range(limit+1)] for _ in range(7)]
    vec={}
    for n in range(1,limit+1):
        for i in range(n+1):
            j=n-i
            vals=[binomial(c,0,i)*binomial(c,1,j) for c in range(3)]
            vec[i,j]=(vals[0]-vals[2],vals[1]-vals[2])
    aa,bb,cc,dd,ee,ff,gg=moments
    for ix,(u,v) in vec.items():
        n=2*sum(ix)
        if n>limit:continue
        aa[n]+=u*u;bb[n]+=2*u*v;cc[n]+=v*v
    keys=list(vec)
    for x,y,z in combinations_with_replacement(keys,3):
        n=sum(x)+sum(y)+sum(z)
        if n>limit:continue
        w=coupling(x[0],y[0],z[0])*coupling(x[1],y[1],z[1])
        if not w:continue
        w*=1 if x==z else 3 if x==y or y==z else 6
        ui,vi=vec[x];uj,vj=vec[y];uk,vk=vec[z]
        dd[n]+=w*ui*uj*uk;gg[n]+=w*vi*vj*vk
        ff[n]+=w*(ui*vj*vk+vi*uj*vk+vi*vj*uk)
        ee[n]+=w*(ui*uj*vk+ui*vj*uk+vi*uj*uk)
    return moments


RESULTANT_TERMS = (
 ((3,0,0,0,0,0,2),1),((2,1,0,0,0,1,1),-1),((2,0,1,0,1,0,1),-2),
 ((2,0,1,0,0,2,0),1),((1,2,0,0,1,0,1),1),((1,1,1,1,0,0,1),3),
 ((1,1,1,0,1,1,0),-1),((1,0,2,1,0,1,0),-2),((1,0,2,0,2,0,0),1),
 ((0,3,0,1,0,0,1),-1),((0,2,1,1,0,1,0),1),((0,1,2,1,1,0,0),-1),
 ((0,0,3,2,0,0,0),1))


def multiply(left,right,limit):
    out={}
    for i,u in left.items():
        for j,v in right.items():
            if i+j<=limit:out[i+j]=out.get(i+j,ZERO)+u*v
    return {i:v for i,v in out.items() if v}


def resultant_jet(limit=20):
    ms=symbolic_moments(limit-10)
    if any(any(m[:2]) for m in ms[:3]) or any(any(m[:4]) for m in ms[3:]):
        raise ValueError('Unexpected low moment order')
    ms=[{i:v for i,v in enumerate(m) if v and i <= (limit-12 if n<3 else limit-10)}
        for n,m in enumerate(ms)]
    out={}
    for powers,scalar in RESULTANT_TERMS:
        term={0:ONE}
        for i,power in enumerate(powers):
            for _ in range(power):term=multiply(term,ms[i],limit)
        for n,v in term.items():out[n]=out.get(n,ZERO)+scalar*v
    return {n:v for n,v in out.items() if v}


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def verify_couplings():
    from fractions import Fraction

    def laguerre(n):
        return [Fraction((-1)**(n+k)*comb(n,k), factorial(k)) for k in range(n+1)]

    def multiply(left, right):
        out=[Fraction(0)]*(len(left)+len(right)-1)
        for i,u in enumerate(left):
            for j,v in enumerate(right):
                out[i+j]+=u*v
        return out

    checked=0
    for i,j,k in combinations_with_replacement(range(6),3):
        p=multiply(multiply(laguerre(i),laguerre(j)),laguerre(k))
        integral=sum(c*factorial(n) for n,c in enumerate(p))
        require(integral==coupling(i,j,k), 'The Laguerre coupling formula disagrees with direct integration.')
        rhs=sum(coupling(i-di,j-dj,k-dk) for di,dj,dk in ((1,1,0),(1,0,1),(0,1,1))
                if i>=di and j>=dj and k>=dk)
        if min(i,j,k)>=1:
            rhs+=2*coupling(i-1,j-1,k-1)
        if i+j+k:
            require(integral==rhs, 'The positive coupling recurrence failed.')
        checked+=1
    require(coupling(1,1,1)==2 and coupling(2,2,2)==10,
            'The small hand coupling controls failed.')
    return checked


def direct_factorial_moments(a,d,left,right):
    from fractions import Fraction

    offsets=(0,a,d)
    norms=[factorial(left+e)*factorial(right+d-e) for e in offsets]

    def moment(i,j,k):
        order=i+j+k
        off=a*j+d*k
        return Fraction(factorial(order*left+off)*factorial(order*(right+d)-off),
                        norms[0]**i*norms[1]**j*norms[2]**k)

    result=[]
    for order in (2,3):
        for y_degree in range(order+1):
            value=sum(comb(order,y_degree)*comb(order-y_degree,i)*comb(y_degree,j)
                      * (-1)**(order-i-j)*moment(i,j,order-i-j)
                      for i in range(order-y_degree+1) for j in range(y_degree+1))
            require(value.denominator==1, 'A normalized factorial moment is not integral.')
            result.append(int(value))
    return result


def verify_symbolic_jet():
    import sympy as sp

    aa,bb,cc,dd,ee,ff,gg,t=sp.symbols('A B C D E F G t')
    direct=sp.Poly(sp.resultant(aa+bb*t+cc*t*t,dd+ee*t+ff*t*t+gg*t**3,t),
                   aa,bb,cc,dd,ee,ff,gg)
    require(dict(direct.terms())==dict(RESULTANT_TERMS),
            'The retained Sylvester resultant coefficients changed.')
    jet=resultant_jet()
    span=A+B
    total=L+H+span
    vandermonde=(A*B*span)**6
    expected={
        18:162*vandermonde,
        19:324*total*vandermonde,
        20:vandermonde*(fmpq(27,2)*(47*total*total-74*total-1)
                         +81*(A*A+B*B+span*span)
                         +240*(3*(L-H)+A-B)**2),
    }
    require(set(jet)==set(expected), 'Unexpected order in the complete resultant jet through 20.')
    for order, value in expected.items():
        require(jet[order]==value, f'The universal resultant coefficient at order {order} changed.')
    # Keep the polynomial degree assertion in the computation itself.
    for order, value in jet.items():
        require(value.total_degree()<=order, 'A resultant coefficient exceeds its parameter degree bound.')


def verify_full_controls():
    profiles=((1,2,0,0),(1,2,0,10),(3,4,3,0),
              (1,4,0,3),(2,5,1,2),(6,12,0,0))
    coefficients=0
    for a,d,left,right in profiles:
        moments=filtered_moments(a,d,left,right)
        require([int(p(1)) for p in moments]==direct_factorial_moments(a,d,left,right),
                'The filtered moments at one disagree with direct factorial integration.')
        result=filtered_resultant(a,d,left,right)
        require(result.degree()<=12*(left+right+d), 'The full resultant exceeds its degree bound.')
        require(all(v>=0 for v in result.coeffs()), 'Full noise coefficient positivity failed in a retained control.')
        require(next((i for i,v in enumerate(result.coeffs()) if v),None)==18,
                'The initial noise order changed.')
        require(result[18]==162*(a*d*(d-a))**6, 'The universal leading coefficient failed in a full control.')
        require(result(1)>0, 'The full moment resultant vanishes in a control.')
        coefficients+=sum(bool(v) for v in result.coeffs())

    # The signed remainder has a negative coefficient in an admissible profile;
    # passing from S to the full resultant is essential for this proposed route.
    aa,bb,cc,dd,ee,ff,gg=filtered_moments(3,4,3,0)
    signed=cc*cc*dd+aa*bb*gg-aa*cc*ff
    require(signed[10]==-72576, 'The signed-remainder negative control changed.')

    base=filtered_resultant(1,2,0,0)
    z=fmpz_poly([0,1])
    require(base==2592*z**18*(z+1)**4*(z*z+4),
            'The explicit smallest resultant changed.')
    return len(profiles), coefficients


def main():
    couplings=verify_couplings()
    verify_symbolic_jet()
    profiles,coefficients=verify_full_controls()
    print(f'Passed {couplings} exact Laguerre coupling/integral controls and the positive recurrence.')
    print('Verified the universal resultant coefficients at orders 18, 19 and 20 over Q[L,H,a,b].')
    print(f'Passed {profiles} full filtered profiles ({coefficients} nonzero resultant coefficients), compared with direct factorial moments.')
    print('The signed-remainder negative control and the explicit span-two factorization passed.')
    print('All higher-order uniform noise coefficient positivity and the complete family theorem remain open.')


if __name__=='__main__':
    main()


