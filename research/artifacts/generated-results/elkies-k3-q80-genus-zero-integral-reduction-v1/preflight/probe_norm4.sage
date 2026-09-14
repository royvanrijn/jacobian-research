#!/usr/bin/env sage
"""Bounded equation-side scan of the known1313 unoriented norm4 classes."""
from sage.all import GF, QQ, ZZ, PolynomialRing, matrix, vector, pari
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
import json
import resource
import time

resource.setrlimit(resource.RLIMIT_CPU,(60,65))
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
start=time.process_time()
ROOT=Path(__file__).resolve().parents[4];OUT=Path(__file__).resolve().parents[1]
source=json.loads((ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json').read_text())
F=GF(131);R=PolynomialRing(F,'t');t=R.gen();K=R.fraction_field()
def poly(cs):return R([F(QQ(v)) for v in cs])
A=poly(source['weierstrass_model']['A_coefficients_low_to_high']);B=poly(source['weierstrass_model']['B_coefficients_low_to_high'])
G=matrix(ZZ,source['sections']['height_gram'])
basis=[]
for r in source['sections']['records']:
    pair=[]
    for key in ['X','Y']:
        f=r[key];pair.append(K(poly(f['numerator_coefficients_low_to_high']))/poly(f['denominator_coefficients_low_to_high']))
    assert pair[1]**2==pair[0]**3+A*pair[0]+B
    basis.append(tuple(pair))
def plus(P,Q):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u:
        if y==-v:return None
        s=(3*x*x+A)/(2*y)
    else:s=(v-y)/(u-x)
    z=s*s-x-u
    return (z,s*(x-z)-y)
@lru_cache(None)
def point(word):
    if not any(word):return None
    v=vector(ZZ,word);gv=G*v
    choices=[(2*(1 if c>0 else -1)*gv[i]-G[i,i],i) for i,c in enumerate(word) if c]
    i=max(choices)[1];sign=1 if word[i]>0 else -1
    rest=list(word);rest[i]-=sign
    x,y=basis[i]
    return plus(point(tuple(rest)),(x,sign*y))
enum=pari(G).qfminim(4)
assert int(enum[0])==2626
words=[]
for col in matrix(ZZ,enum[2]).columns():
    w=tuple(map(int,col));v=vector(ZZ,w);assert v*G*v==4
    if next(c for c in w if c)<0:w=tuple(-c for c in w)
    words.append(w)
words=sorted(set(words));assert len(words)==1313
sites=[2,9,12,27,32,35,38,46,48,63,74,75,103,105,107,115]
rows=[];repeated=[];all_roots=Counter()
for index,word in enumerate(words):
    P=point(word);assert P is not None
    x,y=P;assert x.denominator()==1 and y.denominator()==1
    x,y=R(x),R(y);assert x.degree()<=4 and y.degree()<=6
    assert y*y==x*x*x+A*x+B and (3*x*x+A).gcd(y)==1
    roots=[]
    for b in sites:
        if not y(b):
            n=1;d=y.derivative()
            while not d(b):n+=1;d=d.derivative();assert d
            row={'index':index,'word':list(word),'b':b,'root':int(x(b)),'order':n}
            roots.append([b,int(x(b)),n]);all_roots[b]+=1
            if n>1:repeated.append(row)
    rows.append({'index':index,'word':list(word),'x':[int(v) for v in x.list()],'y':[int(v) for v in y.list()],'roots':roots})
    if index%128==127:print(json.dumps({'compiled':index+1,'repeated_so_far':len(repeated),'cpu_seconds':time.process_time()-start}),flush=True)
with (OUT/'preflight-norm4.json').open('x') as f:json.dump({'status':'PREFLIGHT','rows':rows,'repeated':repeated,'root_counts':dict(all_roots),'cpu_seconds':time.process_time()-start},f,indent=2,sort_keys=True);f.write('\n')
print(json.dumps({'status':'PREFLIGHT','sections':len(rows),'repeated':repeated,'root_counts':dict(all_roots),'cpu_seconds':time.process_time()-start,'point_cache':str(point.cache_info())},sort_keys=True),flush=True)
