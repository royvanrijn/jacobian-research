#!/usr/bin/env python3
"""Refine the rank-10 fingerprint while explicitly excluding Coxeter collisions.

The first healthy-looking rank-10 refinement was later found to sit extremely
close to an existing Coxeter section (V3).  In the current gauge many y/A/B
coefficients are tiny, so near-collision can make several unrelated quadratic
line fits look artificially excellent in absolute error.

This script keeps the same target fingerprint

    (0,0,2,2,2,1,1,1,1)

but adds a scale-free distinctness guard.  Let d45 be the minimum coefficient
norm distance from x_Q to the 45 explicit Coxeter x-polynomials and let s45 be
the median pairwise distance among those 45 polynomials.  We require

    d45 / s45 >= --min-distinct-ratio.

The final diagnostic also avoids an arbitrary line-error threshold.  It reports

    edge_gap = (best non-target quadratic-line error)
               / (worst required-edge error).

A genuine fingerprint should have a large edge_gap.  A value near 1 means the
supposed target edges are numerically indistinguishable from accidental line
fits and the section is not certified.
"""
from __future__ import annotations

from pathlib import Path
import argparse
import json
import numpy as np
from scipy.optimize import least_squares

BASE = Path(__file__).resolve().parents[1]
PAIRS = [(i,j) for i in range(9) for j in range(i+1,9)]
PI = {p:n for n,p in enumerate(PAIRS)}
TRIPLES = [(i,j,k) for i in range(9) for j in range(i+1,9) for k in range(j+1,9)]
TINC=np.zeros((84,9)); EINC=np.zeros((36,9))
for r,t in enumerate(TRIPLES): TINC[r,list(t)]=1.0
for r,(i,j) in enumerate(PAIRS): EINC[r,[i,j]]=1.0
TP=np.linalg.inv(TINC.T@TINC)@TINC.T
EP=np.linalg.inv(EINC.T@EINC)@EINC.T


def pmul(a,b): return np.convolve(a,b)
def p3(a): return pmul(pmul(a,a),a)
def pder(a): return np.asarray([k*a[k] for k in range(1,len(a))],dtype=float) if len(a)>1 else np.zeros(1)
def padd(a,b,sa=1.,sb=1.):
    z=np.zeros(max(len(a),len(b))); z[:len(a)]+=sa*a; z[:len(b)]+=sb*b; return z

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
    verts=[a for a in anchors if a[0]<9]; return verts[0] if verts else anchors[0]


def conv_matrix(dx):
    M=np.zeros((7,3))
    for j in range(3): M[j:j+5,j]=dx
    return M

def projection_error(lhs,dx):
    M=conv_matrix(dx); m,*_=np.linalg.lstsq(M,lhs,rcond=None); return lhs-M@m

def line_abs_error(xq,yq,xp,yp,sign):
    return float(np.max(np.abs(projection_error(yq+sign*yp,xq-xp))))


def branch_norms(A,B):
    disc=padd(p3(A),pmul(B,B),4.,27.)
    jv=padd(pmul(pder(A),B),pmul(A,pder(B)),3.,-2.)
    return float(np.linalg.norm(disc)),float(np.linalg.norm(jv))


def pairwise_scale(xs):
    ds=[]
    for i in range(len(xs)):
        for j in range(i+1,len(xs)):
            d=float(np.linalg.norm(xs[i]-xs[j]))
            if d>0: ds.append(d)
    return max(float(np.median(ds)),1e-300)


def decode(z,fixed,anchors):
    s=unpack(z[:104],fixed)
    xs,ys,labels,A,B,blocks=derived45(s)
    xq=z[104:109]; m0=z[109:112]
    primary=choose_primary(anchors); ordered=[primary]+[a for a in anchors if a!=primary]
    ai,p=ordered[0]
    yq=pmul(m0,xq-xs[ai])-(p/2.)*ys[ai]
    return s,xs,ys,labels,A,B,blocks,xq,m0,yq,ordered


def raw_blocks(z,fixed,anchors):
    s,xs,ys,labels,A,B,blocks,xq,m0,yq,ordered=decode(z,fixed,anchors)
    out=[blocks[0].ravel(),blocks[1].ravel(),.5*blocks[2].ravel(),.25*blocks[3].ravel()]
    out.append(.5*(pmul(yq,yq)-p3(xq)-pmul(A,xq)-B))
    for aj,pj in ordered[1:]:
        out.append(projection_error(yq+(pj/2.)*ys[aj],xq-xs[aj]))
    return out


def diagnostics(z,fixed,anchors,root_disc,root_jv):
    s,xs,ys,labels,A,B,blocks,xq,m0,yq,ordered=decode(z,fixed,anchors)
    base=max(float(np.max(np.abs(b))) for b in blocks)
    curve=float(np.max(np.abs(pmul(yq,yq)-p3(xq)-pmul(A,xq)-B)))
    expected=dict(anchors)
    target_errors=[]; nontarget_errors=[]; wrong_errors=[]
    for i in range(45):
        ep=line_abs_error(xq,yq,xs[i],ys[i],+1.)
        em=line_abs_error(xq,yq,xs[i],ys[i],-1.)
        if i in expected:
            want=expected[i]
            target_errors.append(ep if want==2 else em)
            wrong_errors.append(em if want==2 else ep)
        else:
            nontarget_errors.append(min(ep,em))
    target_max=max(target_errors)
    nontarget_min=min(nontarget_errors)
    wrong_min=min(wrong_errors)
    edge_gap=nontarget_min/max(target_max,1e-300)
    wrong_gap=wrong_min/max(target_max,1e-300)
    scale=pairwise_scale(xs)
    dists=np.linalg.norm(xs-xq[None,:],axis=1)
    nearest=int(np.argmin(dists)); distinct_ratio=float(dists[nearest]/scale)
    dn,jn=branch_norms(A,B)
    return {
        'raw_max':max(base,curve,target_max),
        'base_max':base,'curve_max':curve,'target_line_max':target_max,
        'nontarget_line_min':nontarget_min,'wrong_sign_min':wrong_min,
        'edge_gap':edge_gap,'wrong_sign_gap':wrong_gap,
        'distinct_ratio':distinct_ratio,'nearest_label':labels[nearest],
        'nearest_x_distance':float(dists[nearest]),'coxeter_x_scale':scale,
        'delta_ratio':dn/max(root_disc,1e-300),'jvar_ratio':jn/max(root_jv,1e-300),
    },(s,xs,ys,A,B,xq,yq)


ap=argparse.ArgumentParser()
ap.add_argument('--root',type=Path,default=BASE/'results/coxeter9-slope-numeric-v1/root-000029')
ap.add_argument('--previous',type=Path,default=BASE/'results/rank10-winning-refine-v1')
ap.add_argument('--v-pairings',default='0,0,2,2,2,1,1,1,1')
ap.add_argument('--restarts',type=int,default=16)
ap.add_argument('--max-nfev',type=int,default=3500)
ap.add_argument('--seed',type=int,default=20260820)
ap.add_argument('--min-distinct-ratio',type=float,default=2e-2)
ap.add_argument('--distinct-weight',type=float,default=1e-5)
ap.add_argument('--min-edge-gap',type=float,default=1e2)
ap.add_argument('--min-wrong-sign-gap',type=float,default=1e2)
ap.add_argument('--min-branch-ratio',type=float,default=1e-3)
ap.add_argument('--branch-weight',type=float,default=1e-7)
ap.add_argument('--out',type=Path,default=BASE/'results/rank10-nondegenerate-v1')
args=ap.parse_args()
root=args.root.resolve(); previous=args.previous.resolve(); out=args.out.resolve(); out.mkdir(parents=True,exist_ok=True)
pv=tuple(map(int,args.v_pairings.split(',')))
if len(pv)!=9: raise SystemExit('--v-pairings must contain nine integers')
anchors=anchors_from_pv(pv)
meta=json.loads((root/'candidate.json').read_text())
fixed=fixed_from_meta(meta)
slopes0=np.loadtxt(root/'slopes.txt').reshape(36,3)
free0=pack(slopes0,fixed)
_,_,A0,B0,_=reconstruct(slopes0); root_disc,root_jv=branch_norms(A0,B0)
if (previous/'state.npy').exists():
    zprev=np.load(previous/'state.npy')
    if len(zprev)!=112: raise SystemExit('previous state must have length 112')
else:
    raise SystemExit(f'missing previous state: {previous}/state.npy')

# Report the collision in the previous solution explicitly.
prevstats,_=diagnostics(zprev,fixed,anchors,root_disc,root_jv)
print('PREVIOUS')
for k in ('raw_max','distinct_ratio','nearest_label','edge_gap','wrong_sign_gap','target_line_max','nontarget_line_min','delta_ratio','jvar_ratio'):
    print(f'{k} = {prevstats[k]}')
print()

rng=np.random.default_rng(args.seed)
best=None
for restart in range(args.restarts):
    z0=zprev.copy()
    # First run lets the guard push the old solution away from collision.  Later
    # starts perturb Q on the natural Coxeter x scale and the base only weakly.
    _,xs0,_,_,_,_,_,_,_,_,_=decode(z0,fixed,anchors)
    xscale=pairwise_scale(xs0)
    if restart:
        z0[:104]+=rng.normal(0,2e-5,size=104)
        z0[104:109]+=rng.normal(0,.05*xscale,size=5)
        z0[109:112]+=rng.normal(0,.02*max(1.,float(np.linalg.norm(z0[109:112]))),size=3)

    def fun(z):
        pieces=list(raw_blocks(z,fixed,anchors))
        _,xs,_,_,A,B,_,xq,_,_,_=decode(z,fixed,anchors)
        scale=pairwise_scale(xs)
        dmin=float(np.min(np.linalg.norm(xs-xq[None,:],axis=1)))/scale
        dn,jn=branch_norms(A,B)
        dr=dn/max(root_disc,1e-300); jr=jn/max(root_jv,1e-300)
        pieces.append(np.asarray([
            args.distinct_weight*max(0.,args.min_distinct_ratio-dmin),
            args.branch_weight*max(0.,args.min_branch_ratio-dr),
            args.branch_weight*max(0.,args.min_branch_ratio-jr),
        ]))
        return np.concatenate([np.asarray(p).reshape(-1) for p in pieces])

    res=least_squares(fun,z0,method='trf',jac='2-point',max_nfev=args.max_nfev,
        ftol=1e-14,xtol=1e-14,gtol=1e-14,x_scale='jac',verbose=0)
    stats,data=diagnostics(res.x,fixed,anchors,root_disc,root_jv)
    distinct=stats['distinct_ratio']>=.95*args.min_distinct_ratio
    healthy=stats['delta_ratio']>=args.min_branch_ratio and stats['jvar_ratio']>=args.min_branch_ratio
    separated=stats['edge_gap']>=args.min_edge_gap and stats['wrong_sign_gap']>=args.min_wrong_sign_gap
    usable=distinct and healthy and separated
    # Algebraic residual first after all structural validation gates pass.
    key=(stats['raw_max'],-min(stats['edge_gap'],1e12),-min(stats['wrong_sign_gap'],1e12))
    isbest=usable and (best is None or key<best[0])
    print('NONDEG|restart=%d|raw=%.3e|base=%.3e|curve=%.3e|line=%.3e|distinct=%.3e|nearest=%s|edge_gap=%.3e|wrong_gap=%.3e|nontarget=%.3e|delta=%.3e|jvar=%.3e|usable=%d%s'%(
        restart,stats['raw_max'],stats['base_max'],stats['curve_max'],stats['target_line_max'],
        stats['distinct_ratio'],stats['nearest_label'],stats['edge_gap'],stats['wrong_sign_gap'],stats['nontarget_line_min'],
        stats['delta_ratio'],stats['jvar_ratio'],int(usable),'|BEST' if isbest else ''),flush=True)
    if isbest: best=(key,res.x.copy(),stats,data)

if best is None: raise SystemExit('no healthy, non-colliding, incidence-separated rank-10 candidate found')
_,z,stats,data=best; s,xs,ys,A,B,xq,yq=data
np.save(out/'state.npy',z)
np.savetxt(out/'slopes.txt',s,fmt='%.17g')
np.savetxt(out/'x-coxeter45.txt',xs,fmt='%.17g')
np.savetxt(out/'y-coxeter45.txt',ys,fmt='%.17g')
np.savetxt(out/'x-q.txt',xq.reshape(1,-1),fmt='%.17g')
np.savetxt(out/'y-q.txt',yq.reshape(1,-1),fmt='%.17g')
np.savetxt(out/'A.txt',A.reshape(1,-1),fmt='%.17g')
np.savetxt(out/'B.txt',B.reshape(1,-1),fmt='%.17g')
(out/'candidate.json').write_text(json.dumps({'version':2,'kind':'rank10-nondegenerate-refinement','v_pairings':list(pv),**stats},indent=2,sort_keys=True)+'\n')
print('\nBEST NONDEGENERATE RANK10')
for k,v in stats.items(): print(k,'=',v)
print('saved =',out)
