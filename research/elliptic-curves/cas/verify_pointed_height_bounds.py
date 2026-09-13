"""Replay pointed-chart bounds with Fraction arithmetic and SymPy Sturm isolation.

No Sage imports, finite-field factorisation, minimisation or point search.
Previously proved prime inputs are bound to their retained certificate.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb, gcd, lcm
from pathlib import Path
import time


def require(value,message):
    if not value:
        raise ArithmeticError(message)


def trim(a):
    a = list(a)
    while len(a)>1 and not a[-1]: a.pop()
    return a or [0]


def add(a,b):
    out = [F(0)]*max(len(a),len(b))
    for i,c in enumerate(a):out[i]+=c
    for i,c in enumerate(b):out[i]+=c
    return trim(out)


def mul(a,b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i,c in enumerate(a):
        for j,d in enumerate(b):out[i+j]+=c*d
    return trim(out)


def homogeneous(a,M):
    u,v,w,z = map(F,M)
    require(u*z!=v*w,'singular map')
    out = [F(0)]*5
    for i,c in enumerate(a):
        for j in range(i+1):
            for k in range(5-i):
                out[j+k] += c*comb(i,j)*u**j*v**(i-j)*comb(4-i,k)*w**k*z**(4-i-k)
    return out


def affine(a,r,p):
    out=[F(0)]*len(a)
    for i,c in enumerate(a):
        for j in range(i+1):out[j]+=c*comb(i,j)*r**(i-j)*p**j
    return trim(out)


def valuation(a,p):
    a=F(a);require(a.denominator==1,'nonintegral local polynomial')
    if not a:return 10**9
    a=abs(a.numerator);v=0
    while a%p==0:a//=p;v+=1
    return v


def remainder(a,b,p):
    a=trim([int(c)%p for c in a]);b=trim([int(c)%p for c in b])
    require(b!=[0],'zero divisor in polynomial replay')
    while a!=[0] and len(a)>=len(b):
        q=a[-1]*pow(b[-1],-1,p)%p;shift=len(a)-len(b)
        for i,c in enumerate(b):a[i+shift]=(a[i+shift]-q*c)%p
        a=trim(a)
    return a


def pgcd(a,b,p):
    a,b=trim(a),trim(b)
    while b!=[0]:a,b=b,remainder(a,b,p)
    return trim([c*pow(int(a[-1]),-1,p)%p for c in a])


def power_x(exponent,modulus,p):
    value,base=[1],remainder([0,1],modulus,p)
    while exponent:
        if exponent&1:value=remainder(mul(value,base),modulus,p)
        base=remainder(mul(base,base),modulus,p);exponent//=2
    return value


def all_roots(f,g,roots,p):
    f,g=([int(c)%p for c in a] for a in (f,g))
    common=pgcd(f,g,p)
    if len(common)==1:
        require(not roots,'roots of a constant polynomial');return
    expected=pgcd(common,[int(c)%p for c in add(power_x(p,common,p),[0,-1])],p)
    require(roots==sorted(set(roots)) and all(0<=r<p for r in roots),'invalid residue list')
    product=[1]
    for r in roots:product=[int(c)%p for c in mul(product,[-r,1])]
    require(trim(product)==expected,'incomplete residue-root cover')


def local_node(node,f,g,q,p,cap,at=0):
    kind=node['kind']
    if kind=='budget-fallback':expected=cap
    elif kind=='nonsquare':
        v=min(valuation(c,p) for c in q)
        digits=3 if p==2 else 1
        require(valuation(q[0],p)==v and all(valuation(c,p)>=v+digits for c in q[1:]),'nonconstant local unit')
        unit=int(q[0])//p**v
        require(v%2 or (unit%8!=1 if p==2 else pow(unit,(p-1)//2,p)!=1),'unproved local nonsquare')
        expected=-1
    else:
        e=min(valuation(c,p) for c in f+g)
        if kind=='bezout-cap':
            require(at+e>=cap,'premature Bezout cap');expected=cap
        else:
            require(kind=='split' and node['content_valuation']==e and at+e<cap,'local content differs')
            f,g=([c/p**e for c in a] for a in (f,g))
            roots=node['roots'];all_roots(f,g,roots,p)
            require([c['residue'] for c in node['children']]==roots,'child cover differs')
            values=[at+e]
            for child in node['children']:
                r=child['residue']
                values.append(local_node(child['node'],affine(f,r,p),affine(g,r,p),affine(q,r,p),p,cap,at+e))
            expected=max(values)
    require(node['upper']==expected,'local bound differs')
    return expected


def real_check(n,d,q,L,U):
    import sympy as sp
    x=sp.Symbol('x')
    poly=lambda a:sp.Poly.from_list([sp.Rational(c.numerator,c.denominator) for c in a[::-1]],x,domain=sp.QQ)
    f,g,s=map(poly,(n,d,q))
    require(sp.gcd(s,s.diff()).degree()==0,'singular real domain')
    polynomials=[s,f-sp.Rational(L),f+sp.Rational(L),g-sp.Rational(L),g+sp.Rational(L),
                 f-sp.Rational(U),f+sp.Rational(U),g-sp.Rational(U),g+sp.Rational(U)]
    polynomials=[p for p in polynomials if not p.is_zero]
    intervals=sp.polys.polytools.intervals(polynomials,eps=sp.Rational(1,2**40),inf=-1,sup=1)
    samples={F(-1),F(1)};previous=F(-1)
    for (a,b),multiplicities in intervals:
        a,b=F(a),F(b)
        require(-1<=a<=b<=1 and previous<=a,'overlapping real root isolation')
        samples.update((a,b,(previous+a)/2));previous=b
    samples.add((previous+1)/2)
    def ev(a,x):
        y=F(0)
        for c in a[::-1]:y=y*x+c
        return y
    # These are sign cells cut out by *all* roots, not a numerical grid.
    # A simple root of s has a feasible side; isolated domain endpoints
    # are the rational +/-1 points explicitly checked here.
    for z in samples:
        if ev(q,z)>=0:
            value=max(abs(ev(n,z)),abs(ev(d,z)))
            require(L<=value<=U,'false real bound')
    return len(intervals)


def verify(packet,root):
    data=packet['input'];mapping=data['mapping']
    primes=Path(root)/data['prime_certificate']['path']
    require(hashlib.sha256(primes.read_bytes()).hexdigest()==data['prime_certificate']['sha256'],'prime certificate changed')
    proven=set(map(int,json.loads(primes.read_text())['proved_primes']))
    supplied=list(map(int,data['certified_primes']))
    require(supplied==sorted(set(supplied)) and set(supplied)<=proven,'uncertified prime inputs')
    A,B=map(F,data['curve'][3:]);raw=list(map(F,mapping['raw_coefficients']))
    require(list(map(F,data['curve'][:3]))==[0,0,0] and 4*A**3+27*B**2!=0,'bad curve')
    a,b=-raw[2]/6,-raw[1]/8
    require(raw==[-3*a*a-4*A,-8*b,-6*a,0,1] and b*b==a**3+A*a+B,'bad pointed anchor')
    require(list(map(F,packet['anchor']))==[a,b],'anchor changed')
    n=homogeneous([a**3+4*B,4*a*b,6*a*a+4*A,4*b,a],mapping['matrix'])
    d=homogeneous(raw,mapping['matrix']);q=list(map(F,mapping['discriminant_quartic']))
    ratio=F(mapping['square_ratio'])
    from math import isqrt
    require(ratio>0 and isqrt(ratio.numerator)**2==ratio.numerator and
            isqrt(ratio.denominator)**2==ratio.denominator and d==[ratio*c for c in q],'bad map identity')
    P=list(map(F,mapping['reduced_P']));Q=list(map(F,mapping['reduced_Q']))
    require(q==add([4*c for c in P],mul(Q,Q))+[F(0)]*(5-len(add([4*c for c in P],mul(Q,Q)))),'generalized quartic differs')
    scale=lcm(*(c.denominator for c in n+d));content=gcd(*(int(c*scale) for c in n+d))
    n,d=([c*scale/content for c in v] for v in (n,d))
    require(n==list(map(F,packet['numerator'])) and d==list(map(F,packet['denominator'])),'primitive x-map differs')
    for row,f,g in zip(packet['bezout'],(n,n[::-1]),(d,d[::-1])):
        u,v=map(lambda a:list(map(F,a)),(row['left'],row['right']))
        require(len(u)<=4 and len(v)<=4 and all(c.denominator==1 for c in u+v),'invalid homogeneous Bezout degrees')
        require(add(mul(u,f),mul(v,g))==[F(row['constant'])],'false Bezout identity')
    require(len(packet['bezout'])==2,'both projective charts required')
    C=lcm(*(int(r['constant']) for r in packet['bezout']))
    require(C>0 and C==int(packet['uniform_bezout_divisor']),'uniform divisor differs')
    remaining=C;refined=1;used=[]
    for row in packet['local_trees']:
        p=int(row['prime']);require(p in supplied and p not in used,'bad local prime')
        used.append(p);cap=valuation(remaining,p);require(cap>0 and cap==row['bezout_cap'],'local cap differs')
        remaining//=p**cap
        aa=local_node(row['affine'],n,d,q,p,cap)
        bb=local_node(row['infinity'],affine(n[::-1],0,p),affine(d[::-1],0,p),affine(q[::-1],0,p),p,cap)
        k=max(aa,bb);require(0<=k<=cap and k==row['upper'],'local upper differs')
        refined*=p**k
    refined*=remaining
    require(remaining==int(packet['unfactored_residual']) and refined==int(packet['finite_gcd_divisor']),'refined divisor differs')
    L,U=map(F,(packet['real']['lower'],packet['real']['upper']))
    require(0<L<=U,'invalid real enclosure')
    count=real_check(n,d,q,L,U)+real_check(n[::-1],d[::-1],q[::-1],L,U)
    radius=F(data['search_height'])**4*L/refined
    require(radius==F(packet['guaranteed_multiplicative_x_height_radius']),'coverage radius differs')
    return {'status':'PASS_FRACTION_AND_STURM_REPLAY','local_primes':len(used),
        'real_root_intervals':count,'finite_gcd_divisor':str(refined),
        'guaranteed_multiplicative_x_height_radius':str(radius),
        'boundary':'Exact rational/logarithmic bound definitions verified. Decimal log displays '
            'and claimed optimality are not proof inputs. Retained certified primality is reused.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--research-root',type=Path,default=Path(__file__).resolve().parents[2])
    args=parser.parse_args()
    require(not args.output.exists(),'preserve previous replay')
    start=time.process_time()
    receipt=verify(json.loads(args.packet.read_text()),args.research_root)
    receipt.update(cpu_seconds=time.process_time()-start,
        packet_sha256=hashlib.sha256(args.packet.read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(receipt['status'],receipt['cpu_seconds'],flush=True)
