#!/usr/bin/env python3
"""Diagnose the tangent dimension of the refined rank-10 reconstruction.

The first rank-10 refinement reported Jacobian nullities by applying relative
thresholds directly to an unscaled finite-difference Jacobian.  In the current
Weierstrass gauge the polynomial coefficients span many orders of magnitude,
so those threshold counts are not reliable dimension estimates.

This diagnostic rebuilds the rank-10 equations with *explicit* quadratic slope
variables for all nine Coxeter anchors (rather than eliminating eight slopes by
least squares), then:

* reconstructs the secondary slopes at the saved solution;
* forms central finite-difference Jacobians at several step sizes;
* applies iterative row/column equilibration before SVD;
* reports the smallest singular values and largest spectral gaps;
* compares raw and equilibrated nullity counts.

Exact rank is invariant under nonsingular row/column scaling; the purpose of
equilibration is only to make floating-point rank decisions less dependent on
coefficient scale.  Stable spectral gaps across step sizes are more meaningful
than any single hard threshold.
"""
from __future__ import annotations

from pathlib import Path
import argparse
import json
import numpy as np

BASE = Path(__file__).resolve().parents[1]
PAIRS=[(i,j) for i in range(9) for j in range(i+1,9)]
PI={p:n for n,p in enumerate(PAIRS)}
TRIPLES=[(i,j,k) for i in range(9) for j in range(i+1,9) for k in range(j+1,9)]
TINC=np.zeros((84,9)); EINC=np.zeros((36,9))
for r,t in enumerate(TRIPLES): TINC[r,list(t)]=1.0
for r,(i,j) in enumerate(PAIRS): EINC[r,[i,j]]=1.0
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
    x=TP@f; rx=f-TINC@x
    q=np.empty((36,7))
    for r,(i,j) in enumerate(PAIRS): q[r]=pmul(s[r],x[i]-x[j])
    y=EP@q; ry=q-EINC@y
    ae=np.empty((36,9))
    for r,(i,j) in enumerate(PAIRS):
        ae[r]=pmul(s[r],y[i]-y[j])-pmul(x[i],x[i])-pmul(x[i],x[j])-pmul(x[j],x[j])
    A=ae.mean(0); rA=ae-A
    bn=np.empty((9,13))
    for i in range(9): bn[i]=pmul(y[i],y[i])-p3(x[i])-pmul(A,x[i])
    B=bn.mean(0); rB=bn-B
    return x,y,A,B,(rx,ry,rA,rB)


def derived45(s):
    x,y,A,B,blocks=reconstruct(s)
    xs=[x[i].copy() for i in range(9)]; ys=[y[i].copy() for i in range(9)]
    labels=[f'V{i}' for i in range(9)]
    for r,(i,j) in enumerate(PAIRS):
        m=s[r]
        xd=pmul(m,m)-x[i]-x[j]
        yd=-y[i]-pmul(m,xd-x[i])
        xs.append(xd); ys.append(yd); labels.append(f'D{i}{j}')
    return np.asarray(xs),np.asarray(ys),labels,A,B,blocks


def fixed_from_meta(meta):
    sign=int(meta['gauge_sign']); sp=tuple(sorted(map(int,meta.get('second_pair',[0,2]))))
    return {3*PI[(0,1)]:float(sign),3*PI[(0,1)]+1:0.,3*PI[(0,1)]+2:1.,3*PI[sp]+1:1.}
def freepos(fixed): return [i for i in range(108) if i not in fixed]
def pack(s,fixed):
    f=s.reshape(-1); return np.asarray([f[i] for i in freepos(fixed)])
def unpack(z,fixed):
    f=np.empty(108); p=freepos(fixed); f[p]=z
    for i,v in fixed.items(): f[i]=v
    return f.reshape(36,3)


def fingerprint45(pv):
    fp=list(map(int,pv)); fp += [int(pv[i]-pv[j]) for i,j in PAIRS]; return tuple(fp)
def anchors_from_pv(pv): return [(i,int(p)) for i,p in enumerate(fingerprint45(pv)) if abs(int(p))==2]
def choose_primary(anchors):
    vertices=[a for a in anchors if a[0]<9]; return vertices[0] if vertices else anchors[0]


def conv_matrix(dx):
    M=np.zeros((7,3))
    for j in range(3): M[j:j+5,j]=dx
    return M


def explicit_state(saved_state,fixed,anchors):
    """Expand saved [104 base, 5 xQ, 3 primary slope] to all 9 slopes."""
    base=saved_state[:104]
    xq=saved_state[104:109]
    m0=saved_state[109:112]
    s=unpack(base,fixed)
    xs,ys,_,A,B,_=derived45(s)
    primary=choose_primary(anchors)
    ordered=[primary]+[a for a in anchors if a!=primary]
    ai,p=ordered[0]
    yq=pmul(m0,xq-xs[ai])-(p/2.0)*ys[ai]
    ms=[m0]
    for aj,pj in ordered[1:]:
        lhs=yq+(pj/2.0)*ys[aj]
        M=conv_matrix(xq-xs[aj])
        m,*_=np.linalg.lstsq(M,lhs,rcond=None)
        ms.append(m)
    return np.concatenate([base,xq,np.asarray(ms).reshape(-1)])


def raw_explicit(z,fixed,anchors):
    base=z[:104]; xq=z[104:109]; ms=z[109:].reshape(len(anchors),3)
    s=unpack(base,fixed)
    xs,ys,_,A,B,blocks=derived45(s)
    primary=choose_primary(anchors)
    ordered=[primary]+[a for a in anchors if a!=primary]
    ai,p=ordered[0]
    yq=pmul(ms[0],xq-xs[ai])-(p/2.0)*ys[ai]
    pieces=[blocks[0].ravel(),blocks[1].ravel(),0.5*blocks[2].ravel(),0.25*blocks[3].ravel()]
    pieces.append(0.5*(pmul(yq,yq)-p3(xq)-pmul(A,xq)-B))
    for k,(aj,pj) in enumerate(ordered[1:],start=1):
        pieces.append(yq+(pj/2.0)*ys[aj]-pmul(ms[k],xq-xs[aj]))
    return np.concatenate(pieces)


def jacobian(z,fixed,anchors,step):
    f0=raw_explicit(z,fixed,anchors)
    J=np.empty((len(f0),len(z)))
    for c in range(len(z)):
        h=step*max(1.0,abs(float(z[c])))
        a=z.copy(); b=z.copy(); a[c]+=h; b[c]-=h
        J[:,c]=(raw_explicit(a,fixed,anchors)-raw_explicit(b,fixed,anchors))/(2.0*h)
    return J


def equilibrate(J,iterations=8,floor=1e-300):
    """Symmetric-ish Ruiz equilibration; returns scaled J and scales."""
    A=np.asarray(J,dtype=float).copy()
    rs=np.ones(A.shape[0]); cs=np.ones(A.shape[1])
    for _ in range(iterations):
        rn=np.linalg.norm(A,axis=1)
        rf=1.0/np.sqrt(np.maximum(rn,floor))
        A*=rf[:,None]; rs*=rf
        cn=np.linalg.norm(A,axis=0)
        cf=1.0/np.sqrt(np.maximum(cn,floor))
        A*=cf[None,:]; cs*=cf
    return A,rs,cs


def report_spectrum(tag,sv,nvars):
    smax=float(sv[0]); rel=sv/max(smax,1e-300)
    print(tag+'|smax=%.6e|smin=%.6e|cond=%.6e'%(smax,float(sv[-1]),smax/max(float(sv[-1]),1e-300)))
    for e in (4,5,6,7,8,9,10,11,12):
        null=nvars-int(np.count_nonzero(rel>10.0**(-e)))
        print('%s|null1e-%d=%d'%(tag,e,null))
    tail=min(20,len(sv))
    print(tag+'|tail='+' '.join('%.3e'%x for x in sv[-tail:]))
    # A tangent dimension d corresponds to a gap between sv[-d-1] and sv[-d].
    gaps=[]
    for d in range(1,min(18,len(sv)-1)+1):
        hi=float(sv[-d-1]); lo=float(sv[-d])
        gaps.append((hi/max(lo,1e-300),d,hi,lo))
    gaps.sort(reverse=True)
    print(tag+'|largest_bottom_gaps='+' '.join('d%d:%.3e'%(d,g) for g,d,_,_ in gaps[:8]))


ap=argparse.ArgumentParser()
ap.add_argument('--root',type=Path,default=BASE/'results/coxeter9-slope-numeric-v1/root-000029')
ap.add_argument('--refined',type=Path,default=BASE/'results/rank10-winning-refine-v1')
ap.add_argument('--v-pairings',default='0,0,2,2,2,1,1,1,1')
ap.add_argument('--steps',default='1e-4,3e-5,1e-5,3e-6,1e-6,3e-7,1e-7')
ap.add_argument('--out',type=Path,default=BASE/'results/rank10-tangent-diagnostics-v1')
args=ap.parse_args()
root=args.root.resolve(); refined=args.refined.resolve(); out=args.out.resolve(); out.mkdir(parents=True,exist_ok=True)
meta=json.loads((root/'candidate.json').read_text())
fixed=fixed_from_meta(meta)
pv=tuple(map(int,args.v_pairings.split(','))); anchors=anchors_from_pv(pv)
state=np.load(refined/'state.npy')
if len(state)!=112: raise SystemExit('expected saved rank10 state of length 112, got %d'%len(state))
z=explicit_state(state,fixed,anchors)
raw=raw_explicit(z,fixed,anchors)
print('explicit_variables =',len(z))
print('explicit_residuals =',len(raw))
print('anchors =',len(anchors))
print('explicit_raw_max = %.6e'%float(np.max(np.abs(raw))))
print('expected_generic_rank10_moduli_dimension = 8')
print()
steps=[float(x) for x in args.steps.split(',') if x.strip()]
summary=[]
for step in steps:
    print('STEP %.1e'%step,flush=True)
    J=jacobian(z,fixed,anchors,step)
    sv=np.linalg.svd(J,compute_uv=False)
    Jeq,rs,cs=equilibrate(J)
    sve=np.linalg.svd(Jeq,compute_uv=False)
    report_spectrum('RAW',sv,len(z))
    report_spectrum('EQUIL',sve,len(z))
    np.savetxt(out/('sv-raw-%.0e.txt'%step),sv.reshape(1,-1),fmt='%.17g')
    np.savetxt(out/('sv-equilibrated-%.0e.txt'%step),sve.reshape(1,-1),fmt='%.17g')
    # Keep the strongest bottom-gap candidate for a compact cross-step table.
    gaps=[]
    for d in range(1,min(18,len(sve)-1)+1):
        gaps.append((float(sve[-d-1])/max(float(sve[-d]),1e-300),d))
    g,d=max(gaps)
    summary.append((step,d,g,float(sve[-1]),float(sve[0])))
    print()
print('CROSS-STEP EQUILIBRATED GAP SUMMARY')
for step,d,g,smin,smax in summary:
    print('GAP|step=%.1e|best_d=%d|gap=%.6e|smax=%.6e|smin=%.6e'%(step,d,g,smax,smin))
print('saved =',out)
