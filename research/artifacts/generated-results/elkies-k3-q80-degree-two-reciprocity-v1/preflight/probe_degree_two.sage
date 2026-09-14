#!/usr/bin/env sage
"""Bounded degree-two norm-character preflight; no coefficient existence claim."""
from sage.all import GF, PolynomialRing, QQ
from collections import Counter, defaultdict
from pathlib import Path
import json
import resource
import time

resource.setrlimit(resource.RLIMIT_CPU,(60,65))
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
start=time.process_time()
ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parents[1]
source=json.loads((ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json').read_text())
p=131
F=GF(p)
R=PolynomialRing(F,'t'); t=R.gen()
F2=GF(p*p,'a',modulus=PolynomialRing(F,'a').gen()**2+1)
Q=PolynomialRing(F2,'x'); x=Q.gen()
K=R.fraction_field()
def poly(values): return R([F(QQ(v)) for v in values])
A=poly(source['weierstrass_model']['A_coefficients_low_to_high'])
B=poly(source['weierstrass_model']['B_coefficients_low_to_high'])
sections=[]
for row in source['sections']['records']:
    coords=[]
    for key in ['X','Y']:
        v=row[key]
        coords.append(K(poly(v['numerator_coefficients_low_to_high']))/poly(v['denominator_coefficients_low_to_high']))
    assert coords[1]**2==coords[0]**3+A*coords[0]+B
    sections.append(coords)
def encode(z):
    c=z.polynomial().list()
    return int(c[0])+(p*int(c[1]) if len(c)>1 else 0) if c else 0
def rational_at(f,b,w):
    n,d=f.numerator(),f.denominator()
    if b is not None:
        if not d(b): return None
        return n(b)/d(b)
    gap=d.degree()+w-n.degree()
    if gap<0:return None
    return F2(0) if gap else F2(n.leading_coefficient()/d.leading_coefficient())
def characters(b,roots,a,bb,degree):
    points=[]
    for X,Y in sections:
        u,v=rational_at(X,b,4),rational_at(Y,b,6)
        assert (u is None)==(v is None)
        assert u is None or v*v==u**3+a*u+bb
        points.append((u,v))
    result=[]
    for e in roots:
        code=0
        for i,(u,v) in enumerate(points):
            value=F2(1) if u is None else (3*e*e+a if u==e else u-e)
            assert value
            if degree==1:value=F(value)
            if not value.is_square():code|=1<<i
        result.append(code)
    return result
rational=[]; quadratic=[]; bad=[]
hist=Counter(); zeros=[]
for index in range(p):
    for j in range(0,(p+1)//2):
        if j==0:
            b=F2(index);degree=1
        else:
            b=F2(index)+j*F2.gen();degree=2
        a,bb=F2(A(b)),F2(B(b))
        r=(x**3+a*x+bb).roots(multiplicities=True)
        if degree==1:
            r=[(v,m) for v,m in r if v**p==v]
        roots=sorted([v for v,m in r for _ in range(int(m))],key=encode)
        hist[(degree,len(roots),bool(4*a**3+27*bb**2))]+=1
        if len(roots)!=3:continue
        codes=characters(b,roots,a,bb,degree)
        row={'t':encode(b),'roots':[encode(v) for v in roots],'codes':codes}
        assert codes[0]^codes[1]^codes[2]==0
        (rational if degree==1 else quadratic).append(row)
        if not (4*a**3+27*bb**2):bad.append(row)
        if degree==2 and codes.count(0)>=2:zeros.append(row)
a,bb=F2(A[8]),F2(B[12])
roots=sorted([v for v,m in (x**3+a*x+bb).roots() if v**p==v],key=encode)
if len(roots)==3:rational.append({'t':None,'roots':[encode(v) for v in roots],'codes':characters(None,roots,a,bb,1)})
groups=defaultdict(list)
for row in rational:
    for i,e in enumerate(row['roots']):
        for j,f in enumerate(row['roots']):
            if i!=j:groups[(row['codes'][i],row['codes'][j])].append([row['t'],e,f])
collisions=[v for v in groups.values() if len(v)>1]
result={'status':'PREFLIGHT','prime':p,'field_modulus':'a^2+1','histogram':[[*k,v] for k,v in sorted(hist.items())],
 'rational_split_fibres':rational,'quadratic_split_fibres':len(quadratic),'quadratic_zero_code_triples':zeros,
 'bad_quadratic_fibres':bad,'distinct_rational_ordered_pairs':sum(map(len,groups.values())),
 'rational_pair_collisions':collisions,'cpu_seconds':time.process_time()-start,
 'positive_endpoint_complete':False}
path=OUT/'preflight.json'
with path.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
print(json.dumps(result,sort_keys=True),flush=True)
