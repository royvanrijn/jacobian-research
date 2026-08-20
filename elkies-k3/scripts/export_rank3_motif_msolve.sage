from sage.all import *
from pathlib import Path
import argparse, random

ap=argparse.ArgumentParser()
ap.add_argument('--p',type=int,required=True)
ap.add_argument('--slices',type=int,default=5)
ap.add_argument('--seed',type=int,default=1)
ap.add_argument('--out',required=True)
ap.add_argument('--slice-width',type=int,default=6)
a=ap.parse_args()
if not is_prime(a.p) or a.p in (2,3,79): raise SystemExit('bad prime')
F=GF(a.p); random.seed(a.seed); set_random_seed(a.seed)

# Motif: P=1337,Q=736,R=2474; -P,-Q,-R; +/- (P+Q); -(P+R); -(P+Q+R).
names=[]
for z in ('P','Q','R'):
    names += [f'x{z}{i}' for i in range(5)]
    names += [f'y{z}{i}' for i in range(7)]
for z in ('PQ','PR','SR','TQ'):
    names += [f'm{z}{i}' for i in range(3)]
Rng=PolynomialRing(F,names,order='degrevlex')
g=Rng.gens_dict(); Rt=PolynomialRing(Rng,'t'); t=Rt.gen()
def pol(pre,d): return sum(g[f'{pre}{i}']*t^i for i in range(d+1))
def ce(P,d): return [Rng(P[i]) for i in range(d+1)]
def line(X1,Y1,X2,Y2,M): return Y2-Y1-M*(X2-X1)
def Avec(X1,Y1,X2,Y2,M): return M*(Y1+Y2)-(X1^2+X1*X2+X2^2)
def add(X1,Y1,X2,Y2,M):
    X3=M^2-X1-X2
    Y3=M*(X1-X3)-Y1
    return X3,Y3
x={z:pol('x'+z,4) for z in ('P','Q','R')}
y={z:pol('y'+z,6) for z in ('P','Q','R')}
m={z:pol('m'+z,2) for z in ('PQ','PR','SR','TQ')}
E=[]
E += ce(line(x['P'],y['P'],x['Q'],y['Q'],m['PQ']),6)
A=Avec(x['P'],y['P'],x['Q'],y['Q'],m['PQ'])
Sx,Sy=add(x['P'],y['P'],x['Q'],y['Q'],m['PQ'])
E += ce(line(x['P'],y['P'],x['R'],y['R'],m['PR']),6)
E += ce(Avec(x['P'],y['P'],x['R'],y['R'],m['PR'])-A,8)
Tx,Ty=add(x['P'],y['P'],x['R'],y['R'],m['PR'])
E += ce(line(Sx,Sy,x['R'],y['R'],m['SR']),6)
E += ce(Avec(Sx,Sy,x['R'],y['R'],m['SR'])-A,8)
U1x,U1y=add(Sx,Sy,x['R'],y['R'],m['SR'])
E += ce(line(Tx,Ty,x['Q'],y['Q'],m['TQ']),6)
E += ce(Avec(Tx,Ty,x['Q'],y['Q'],m['TQ'])-A,8)
U2x,U2y=add(Tx,Ty,x['Q'],y['Q'],m['TQ'])
E += ce(U1x-U2x,4)
E += ce(U1y-U2y,6)
E=[e for e in E if e!=0]; base=len(E)
vs=list(Rng.gens())
for _ in range(a.slices):
    chosen=random.sample(vs,min(a.slice_width,len(vs)))
    f=Rng(F.random_element())
    for v in chosen:
        c=F.random_element()
        while c==0: c=F.random_element()
        f += Rng(c)*v
    E.append(f)
out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
with out.open('w') as h:
    h.write(','.join(names)+'\n'+str(a.p)+'\n')
    for i,e in enumerate(E):
        h.write(str(Rng(e)).replace('**','^'))
        h.write(',\n' if i+1<len(E) else '\n')
print(f'MOTIFMSOLVE|p={a.p}|seed={a.seed}|slices={a.slices}|vars={len(names)}|base_eqs={base}|eqs={len(E)}|out={out}')
