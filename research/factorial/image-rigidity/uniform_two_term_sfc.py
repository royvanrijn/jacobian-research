"""Exact finite certificates for the ordinary constant-plus-power SFC family.

Universal claims use ALL_UNIVARIATE_TWO_TERM_SFC.tex. This arithmetic module
constructs/checks the finite prefix and never extrapolates a finite scan.
"""
from math import comb,factorial,isqrt
from bisect import bisect_right
from functools import lru_cache
INDEX_CUTOFF=10000


def integer(x,minimum=0):
    if isinstance(x,bool) or not isinstance(x,int) or x<minimum:
        raise ValueError('integer outside the declared range')
    return x


def trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a


def add(a,b,p):
    out=[0]*max(len(a),len(b))
    for i,c in enumerate(a):out[i]=c
    for i,c in enumerate(b):out[i]=(out[i]+c)%p
    return trim(out)


def mul(a,b,p):
    out=[0]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b):out[i+j]=(out[i+j]+u*v)%p
    return trim(out)


def divrem(a,b,p):
    a=a[:];out=[0]*max(1,len(a)-len(b)+1);inv=pow(b[-1],-1,p)
    while len(a)>=len(b) and a!=[0]:
        shift=len(a)-len(b);c=a[-1]*inv%p;out[shift]=c
        for j,v in enumerate(b):a[shift+j]=(a[shift+j]-c*v)%p
        trim(a)
    return trim(out),a


def bezout(f,g,p):
    r0,r1=f,g;u0,u1=[1],[0];v0,v1=[0],[1]
    while r1!=[0]:
        q,r=divrem(r0,r1,p);r0,r1=r1,r
        u0,u1=u1,add(u0,[(-x)%p for x in mul(q,u1,p)],p)
        v0,v1=v1,add(v0,[(-x)%p for x in mul(q,v1,p)],p)
    scale=pow(r0[-1],-1,p)
    return [c*scale%p for c in r0],[c*scale%p for c in u0],[c*scale%p for c in v0]


@lru_cache(None)
def is_prime(p):
    return p>=2 and all(p%d for d in range(2,isqrt(p)+1))


def next_prime(k):
    k=integer(k,2)
    while not is_prime(k):k+=1
    return k


def primes_through(n):
    flags=bytearray([1])*(n+1);flags[:2]=b'\0\0'
    for p in range(2,isqrt(n)+1):
        if flags[p]:flags[p*p:n+1:p]=b'\0'*len(range(p*p,n+1,p))
    return [p for p in range(2,n+1) if flags[p]]


def is_prime_power_or_one(m):
    m=integer(m,1)
    if m==1:return True
    for p in range(2,isqrt(m)+1):
        if m%p==0:
            while m%p==0:m//=p
            return m==1
    return True


def finite_exponents():
    # m=6 already has its complete all-index proof, retained as an input.
    return [m for m in range(1,30) if m!=6 and not is_prime_power_or_one(m)]


def moment_coefficients(d,m,p):
    """Recurrence modulo p; requires p>m*d, so all denominators are units."""
    integer(d);integer(m,1);integer(p,2)
    if p<=m*d:raise ValueError('full-degree reduction requires p>m*d')
    out=[1];binomial=1;fact=1
    for j in range(1,d+1):
        binomial=binomial*(d-j+1)*pow(j,-1,p)%p
        for v in range(m*(j-1)+1,m*j+1):fact=fact*v%p
        out.append(binomial*fact%p)
    return out


@lru_cache(None)
def exact_factorial(k):return factorial(k)


def verify_record(row):
    n,m,p,d,mode,g,u,v=row
    assert is_prime(p)
    if mode=='preceding_prime_kernel':
        assert n==p+d and p>m*max(1,d)
    elif mode=='full_degree_modular':
        assert d==n and p>m*n
    else:raise AssertionError('unknown reduction mode')
    # Independent literal coefficient construction and multiplication check;
    # no trust in the coefficient recurrence or gcd output is required.
    f=[comb(d,j)*(exact_factorial(m*j)%p)%p for j in range(d+1)]
    derivative=trim([j*f[j]%p for j in range(1,len(f))] or [0])
    assert f[-1]!=0 and 1<=len(g)<=2 and g[-1]==1
    assert add(mul(u,f,p),mul(v,derivative,p),p)==g
    return True


def prefix_records():
    primes=primes_through(INDEX_CUTOFF-1)
    for m in finite_exponents():
        for n in range(2,INDEX_CUTOFF):
            if (1<<m)>=4*n*n:continue
            q=primes[bisect_right(primes,n)-1];d=n-q
            if q>m*max(1,d):p=q;mode='preceding_prime_kernel'
            else:p=next_prime(m*n+1000);d=n;mode='full_degree_modular'
            f=moment_coefficients(d,m,p)
            derivative=trim([j*f[j]%p for j in range(1,len(f))] or [0])
            g,u,v=bezout(f,derivative,p)
            if len(g)>2:
                # Additional good reduction primes are finite witnesses,
                # not an assumption that the first local test always works.
                p=next_prime(m*n+2000);d=n;mode='full_degree_modular'
                f=moment_coefficients(d,m,p)
                derivative=trim([j*f[j]%p for j in range(1,len(f))] or [0])
                g,u,v=bezout(f,derivative,p)
            row=[n,m,p,d,mode,g,u,v]
            verify_record(row)
            yield row
