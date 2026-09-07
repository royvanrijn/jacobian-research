#!/usr/bin/env python3
"""Scan the S9 orbit of the rank-10 lattice fingerprint on fixed Coxeter roots.

The raw Coxeter-9 clique has Gram 2(I+J), so the numerical reconstruction does
not intrinsically identify the nine lattice labels.  A rank-10 extension vector
therefore has an S9 orbit of possible pairing fingerprints against the 45
explicit Coxeter minimal sections V_i and D_ij=V_i-V_j.

For each distinct fingerprint and each supplied Coxeter root, this script keeps
that surface fixed and solves only for the added section Q.  All quadratic line
slopes except one are eliminated linearly at every residual evaluation.  This
makes the scan cheap enough to rank the symmetry orbit before another expensive
full continuation.

This is a *proximity/ranking* scan, not a proof that Q exists on the fixed root:
the true rank-10 locus is codimension one inside the rank-9 Coxeter moduli.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path
import argparse, csv, json
import numpy as np
from scipy.optimize import least_squares

BASE = Path(__file__).resolve().parents[1]
PAIRS = [(i,j) for i in range(9) for j in range(i+1,9)]
PI = {p:n for n,p in enumerate(PAIRS)}
TRIPLES = [(i,j,k) for i in range(9) for j in range(i+1,9) for k in range(j+1,9)]
TINC=np.zeros((84,9)); EINC=np.zeros((36,9))
for r,t in enumerate(TRIPLES): TINC[r,list(t)]=1
for r,(i,j) in enumerate(PAIRS): EINC[r,[i,j]]=1
TP=np.linalg.inv(TINC.T@TINC)@TINC.T
EP=np.linalg.inv(EINC.T@EINC)@EINC.T

def pmul(a,b): return np.convolve(a,b)
def p3(a): return pmul(pmul(a,a),a)
def pslope(s,i,j): return -s[PI[(j,i)]] if i>j else s[PI[(i,j)]]

def reconstruct(s):
    f=np.empty((84,5))
    for r,(i,j,k) in enumerate(TRIPLES):
        mij,mik,mjk=pslope(s,i,j),pslope(s,i,k),pslope(s,j,k)
        f[r]=pmul(mik,mij)+pmul(mik,mjk)-pmul(mij,mjk)
    x=TP@f
    q=np.empty((36,7))
    for r,(i,j) in enumerate(PAIRS): q[r]=pmul(s[r],x[i]-x[j])
    y=EP@q
    ae=np.empty((36,9))
    for r,(i,j) in enumerate(PAIRS):
        ae[r]=pmul(s[r],y[i]-y[j])-pmul(x[i],x[i])-pmul(x[i],x[j])-pmul(x[j],x[j])
    A=ae.mean(0)
    bn=np.empty((9,13))
    for i in range(9): bn[i]=pmul(y[i],y[i])-p3(x[i])-pmul(A,x[i])
    B=bn.mean(0)
    return x,y,A,B

def derived45(s):
    x,y,A,B=reconstruct(s)
    xs=[x[i].copy() for i in range(9)]; ys=[y[i].copy() for i in range(9)]
    labels=[f'V{i}' for i in range(9)]
    for r,(i,j) in enumerate(PAIRS):
        m=s[r]
        xd=pmul(m,m)-x[i]-x[j]
        yd=-y[i]-pmul(m,xd-x[i])
        xs.append(xd); ys.append(yd); labels.append(f'D{i}{j}')
    return np.asarray(xs),np.asarray(ys),labels,A,B

def conv_matrix(dx):
    M=np.zeros((7,3))
    for j in range(3): M[j:j+5,j]=dx
    return M

def unique_permutations(values):
    counts=Counter(values); keys=sorted(counts); n=len(values); cur=[0]*n
    def rec(pos):
        if pos==n:
            yield tuple(cur); return
        for k in keys:
            if counts[k]:
                counts[k]-=1; cur[pos]=k
                yield from rec(pos+1)
                counts[k]+=1
    yield from rec(0)

def fingerprint45(pv):
    fp=list(map(int,pv))
    fp += [int(pv[i]-pv[j]) for i,j in PAIRS]
    return tuple(fp)

def anchors_from_fp(fp):
    return [(i,int(p)) for i,p in enumerate(fp) if abs(int(p))==2]

def choose_primary(anchors):
    verts=[a for a in anchors if a[0]<9]
    return verts[0] if verts else anchors[0]

def line_projection_error(lhs,dx):
    M=conv_matrix(dx)
    m,*_=np.linalg.lstsq(M,lhs,rcond=None)
    return lhs-M@m

def solve_fixed_surface(xs,ys,A,B,anchors,rng,max_nfev,restarts):
    primary=choose_primary(anchors)
    ordered=[primary]+[a for a in anchors if a!=primary]
    xscale=max(1e-6,float(np.median(np.linalg.norm(xs[:9],axis=1))))
    mscale=max(.1,float(np.median(np.linalg.norm(xs[:9],axis=1))**0.5))
    best=None
    for rr in range(restarts):
        coeff=rng.normal(size=9); coeff/=max(np.linalg.norm(coeff),1e-12)
        x0=coeff@xs[:9] + rng.normal(0,.05*xscale,size=5)
        m0=rng.normal(0,mscale,size=3)
        z0=np.concatenate([x0,m0])
        def residual(z,limit=None):
            xq=z[:5]; m=z[5:8]
            ai,p=ordered[0]
            yq=pmul(m,xq-xs[ai])-(p/2.)*ys[ai]
            pieces=[.5*(pmul(yq,yq)-p3(xq)-pmul(A,xq)-B)]
            use=len(ordered) if limit is None else min(len(ordered),limit)
            for aj,pj in ordered[1:use]:
                lhs=yq+(pj/2.)*ys[aj]
                pieces.append(line_projection_error(lhs,xq-xs[aj]))
            # keep Q distinct from the primary anchor
            d=np.linalg.norm(xq-xs[ai])/xscale
            pieces.append(np.asarray([1e-4*max(0.,1e-4-d)]))
            return np.concatenate(pieces)
        z=z0
        for lim in (min(4,len(ordered)),None):
            res=least_squares(lambda zz,L=lim:residual(zz,L),z,method='trf',jac='2-point',
                max_nfev=max_nfev,ftol=1e-11,xtol=1e-11,gtol=1e-11,x_scale='jac')
            z=res.x
        xq=z[:5]; m=z[5:8]; ai,p=ordered[0]
        yq=pmul(m,xq-xs[ai])-(p/2.)*ys[ai]
        curve=float(np.max(np.abs(pmul(yq,yq)-p3(xq)-pmul(A,xq)-B)))
        line=0.
        for aj,pj in ordered[1:]:
            e=line_projection_error(yq+(pj/2.)*ys[aj],xq-xs[aj])
            line=max(line,float(np.max(np.abs(e))))
        score=max(curve,line)
        if best is None or score<best[0]: best=(score,curve,line,xq,yq)
    return best

ap=argparse.ArgumentParser()
ap.add_argument('--root',type=Path,action='append',default=None,
    help='Coxeter root directory; repeat option for multiple roots')
ap.add_argument('--chain',type=Path,default=BASE/'results/rank17-extension-chain-v1')
ap.add_argument('--restarts',type=int,default=1)
ap.add_argument('--max-nfev',type=int,default=180)
ap.add_argument('--top',type=int,default=30)
ap.add_argument('--seed',type=int,default=20260820)
ap.add_argument('--out',type=Path,default=BASE/'results/rank10-fingerprint-scan-v1')
args=ap.parse_args()
roots=args.root or [
    BASE/'results/coxeter9-slope-numeric-v1/root-000001',
    BASE/'results/coxeter9-slope-numeric-v1/root-000029',
]
chain=args.chain.resolve(); out=args.out.resolve(); out.mkdir(parents=True,exist_ok=True)
H=np.loadtxt(BASE/'data/lattice/short_vector_basis_gram.txt',dtype=np.int64)
Bsel=np.loadtxt(chain/'selected-vectors.txt',dtype=np.int64)
Vlat=Bsel[:9]; qvec=Bsel[9]
pv_orig=tuple(map(int,((qvec@H)@Vlat.T).tolist()))
print('original V-pairing signature =',' '.join(f'{p:+d}' for p in pv_orig))
# The Q -> -Q orientation is geometrically equivalent for x; include both and dedupe.
fp_seen={};
for pv in unique_permutations(pv_orig):
    for sgn in (1,-1):
        qpv=tuple(sgn*p for p in pv); fp=fingerprint45(qpv)
        fp_seen.setdefault(fp,qpv)
fingerprints=list(fp_seen.items())
print('distinct S9/Q-sign fingerprints =',len(fingerprints))
identity_fp=fingerprint45(pv_orig)
rows=[]; rng=np.random.default_rng(args.seed)
for root0 in roots:
    root=root0.resolve()
    if not (root/'slopes.txt').exists():
        print('SKIP root missing:',root); continue
    slopes=np.loadtxt(root/'slopes.txt').reshape(36,3)
    xs,ys,labels,A,B=derived45(slopes)
    print('\nROOT',root.name,'fingerprints=',len(fingerprints),flush=True)
    for idx,(fp,pv) in enumerate(fingerprints,1):
        anchors=anchors_from_fp(fp)
        if len(anchors)<2: continue
        score,curve,line,xq,yq=solve_fixed_surface(xs,ys,A,B,anchors,rng,args.max_nfev,args.restarts)
        rows.append({
            'root':root.name,'fingerprint_index':idx,'identity':int(fp==identity_fp),
            'score':score,'curve_max':curve,'line_max':line,'anchor_count':len(anchors),
            'v_pairings':' '.join(map(str,pv)),
            'anchors':' '.join(f'{labels[i]}:{p:+d}' for i,p in anchors),
        })
        if idx%50==0:
            best=min(r['score'] for r in rows if r['root']==root.name)
            print(f'SCAN|root={root.name}|done={idx}/{len(fingerprints)}|best={best:.3e}',flush=True)
rows.sort(key=lambda r:r['score'])
with (out/'fingerprints.tsv').open('w',newline='') as f:
    fields=list(rows[0]); w=csv.DictWriter(f,delimiter='\t',fieldnames=fields); w.writeheader(); w.writerows(rows)
print('\nTOP FINGERPRINTS')
for r in rows[:args.top]:
    print(f"FP|root={r['root']}|idx={r['fingerprint_index']}|identity={r['identity']}|score={r['score']:.3e}|curve={r['curve_max']:.3e}|line={r['line_max']:.3e}|anchors={r['anchor_count']}|V={r['v_pairings']}")
print('saved =',out/'fingerprints.tsv')
