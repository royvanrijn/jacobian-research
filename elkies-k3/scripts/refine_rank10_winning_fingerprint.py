#!/usr/bin/env python3
"""Jointly refine the winning rank-10 Coxeter fingerprint on root-000029.

The fixed-surface S9 scan found an essentially machine-precision extra section
for the fingerprint

    (0,0,2,2,2,1,1,1,1)

against the nine numerical Coxeter generators.  This script first recovers that
section robustly on the fixed root, then allows the 104 gauge-fixed Coxeter
slope coefficients to move together with the extra section.

Secondary chord slopes are eliminated by linear least squares, so the joint
nonlinear state has only

    104 Coxeter variables + 5 x_Q coefficients + 3 primary-slope coefficients.

The script reports raw algebraic residuals, non-isotrivial branch ratios, and
the finite-difference Jacobian nullity of the combined rank-10 system.  It also
exports all full-lattice Coxeter permutations compatible with the winning
fingerprint; these are the candidates to distinguish at rank 11.
"""
from __future__ import annotations

from collections import Counter
from itertools import permutations
from pathlib import Path
import argparse
import csv
import json
import numpy as np
from scipy.optimize import least_squares

BASE = Path(__file__).resolve().parents[1]
PAIRS = [(i,j) for i in range(9) for j in range(i+1,9)]
PI = {p:n for n,p in enumerate(PAIRS)}
TRIPLES = [(i,j,k) for i in range(9) for j in range(i+1,9) for k in range(j+1,9)]
TINC=np.zeros((84,9)); EINC=np.zeros((36,9))
for r,t in enumerate(TRIPLES): TINC[r,list(t)] = 1.0
for r,(i,j) in enumerate(PAIRS): EINC[r,[i,j]] = 1.0
TP=np.linalg.inv(TINC.T@TINC)@TINC.T
EP=np.linalg.inv(EINC.T@EINC)@EINC.T


def pmul(a,b): return np.convolve(a,b)
def p3(a): return pmul(pmul(a,a),a)
def pder(a):
    return np.asarray([k*a[k] for k in range(1,len(a))],dtype=float) if len(a)>1 else np.zeros(1)
def padd(a,b,sa=1.,sb=1.):
    z=np.zeros(max(len(a),len(b))); z[:len(a)] += sa*a; z[:len(b)] += sb*b; return z

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


def base_raw(free,fixed):
    *_,blocks=reconstruct(unpack(free,fixed))
    return np.concatenate([blocks[0].ravel(),blocks[1].ravel(),.5*blocks[2].ravel(),.25*blocks[3].ravel()])


def conv_matrix(dx):
    M=np.zeros((7,3))
    for j in range(3): M[j:j+5,j]=dx
    return M

def projection_error(lhs,dx):
    M=conv_matrix(dx); m,*_=np.linalg.lstsq(M,lhs,rcond=None); return lhs-M@m


def branch_polys(A,B):
    disc=padd(p3(A),pmul(B,B),4.,27.)
    jv=padd(pmul(pder(A),B),pmul(A,pder(B)),3.,-2.)
    return disc,jv

def branch_norms(A,B):
    d,j=branch_polys(A,B); return float(np.linalg.norm(d)),float(np.linalg.norm(j))


def fingerprint45(pv):
    fp=list(map(int,pv)); fp += [int(pv[i]-pv[j]) for i,j in PAIRS]; return tuple(fp)
def anchors_from_pv(pv):
    fp=fingerprint45(pv)
    return [(i,int(p)) for i,p in enumerate(fp) if abs(int(p))==2]
def choose_primary(anchors):
    verts=[a for a in anchors if a[0]<9]; return verts[0] if verts else anchors[0]


def decode_joint(z,fixed):
    s=unpack(z[:104],fixed)
    xs,ys,labels,A,B,blocks=derived45(s)
    xq=z[104:109]; m0=z[109:112]
    return s,xs,ys,labels,A,B,blocks,xq,m0


def joint_raw(z,fixed,anchors):
    s,xs,ys,labels,A,B,blocks,xq,m0=decode_joint(z,fixed)
    primary=choose_primary(anchors)
    ordered=[primary]+[a for a in anchors if a!=primary]
    ai,p=ordered[0]
    yq=pmul(m0,xq-xs[ai])-(p/2.)*ys[ai]
    pieces=[blocks[0].ravel(),blocks[1].ravel(),.5*blocks[2].ravel(),.25*blocks[3].ravel()]
    pieces.append(.5*(pmul(yq,yq)-p3(xq)-pmul(A,xq)-B))
    for aj,pj in ordered[1:]:
        pieces.append(projection_error(yq+(pj/2.)*ys[aj],xq-xs[aj]))
    return np.concatenate(pieces)


def diagnostics(z,fixed,anchors,root_disc,root_jv):
    s,xs,ys,labels,A,B,blocks,xq,m0=decode_joint(z,fixed)
    primary=choose_primary(anchors); ordered=[primary]+[a for a in anchors if a!=primary]
    ai,p=ordered[0]; yq=pmul(m0,xq-xs[ai])-(p/2.)*ys[ai]
    base=max(float(np.max(np.abs(b))) for b in blocks)
    curve=float(np.max(np.abs(pmul(yq,yq)-p3(xq)-pmul(A,xq)-B)))
    line=0.
    for aj,pj in ordered[1:]:
        line=max(line,float(np.max(np.abs(projection_error(yq+(pj/2.)*ys[aj],xq-xs[aj])))))
    dn,jn=branch_norms(A,B)
    return {
        'raw_max':max(base,curve,line),'base_max':base,'curve_max':curve,'line_max':line,
        'delta_ratio':dn/max(root_disc,1e-300),'jvar_ratio':jn/max(root_jv,1e-300),
        'slope_rank':int(np.linalg.matrix_rank(s,tol=1e-9)),
    },(s,xs,ys,A,B,xq,yq)


def fixed_surface_start(slopes,pv,rng,restarts,max_nfev):
    xs,ys,labels,A,B,_=derived45(slopes); anchors=anchors_from_pv(pv)
    primary=choose_primary(anchors); ordered=[primary]+[a for a in anchors if a!=primary]
    xscale=max(1e-6,float(np.median(np.linalg.norm(xs[:9],axis=1))))
    mscale=max(.1,float(np.median(np.abs(slopes))))
    best=None
    for rr in range(restarts):
        coeff=rng.normal(size=9); coeff/=max(np.linalg.norm(coeff),1e-12)
        z=np.concatenate([coeff@xs[:9]+rng.normal(0,.05*xscale,size=5),rng.normal(0,mscale,size=3)])
        def f(zz):
            xq,m=zz[:5],zz[5:8]; ai,p=ordered[0]
            yq=pmul(m,xq-xs[ai])-(p/2.)*ys[ai]
            pieces=[.5*(pmul(yq,yq)-p3(xq)-pmul(A,xq)-B)]
            for aj,pj in ordered[1:]: pieces.append(projection_error(yq+(pj/2.)*ys[aj],xq-xs[aj]))
            return np.concatenate(pieces)
        res=least_squares(f,z,method='trf',jac='2-point',max_nfev=max_nfev,
            ftol=1e-13,xtol=1e-13,gtol=1e-13,x_scale='jac')
        score=float(np.max(np.abs(f(res.x))))
        if best is None or score<best[0]: best=(score,res.x.copy())
    return best


def jacobian_singular_values(z,fixed,anchors,step):
    f0=joint_raw(z,fixed,anchors); J=np.empty((len(f0),len(z)))
    for c in range(len(z)):
        h=step*max(1.,abs(float(z[c]))); a=z.copy(); b=z.copy(); a[c]+=h; b[c]-=h
        J[:,c]=(joint_raw(a,fixed,anchors)-joint_raw(b,fixed,anchors))/(2*h)
    return np.linalg.svd(J,compute_uv=False)


def compatible_mappings(original_pv,target_pv):
    # mapping[j] is the original lattice Coxeter index represented by numerical V_j.
    out=[]
    for qsign in (1,-1):
        src=[qsign*p for p in original_pv]
        groups_src={v:[i for i,x in enumerate(src) if x==v] for v in set(src)}
        groups_dst={v:[j for j,x in enumerate(target_pv) if x==v] for v in set(target_pv)}
        if {k:len(v) for k,v in groups_src.items()} != {k:len(v) for k,v in groups_dst.items()}:
            continue
        keys=sorted(groups_src)
        partial=[None]*9
        def rec(k):
            if k==len(keys):
                out.append((qsign,tuple(partial))); return
            val=keys[k]; srcidx=groups_src[val]; dstidx=groups_dst[val]
            for perm in permutations(srcidx):
                for j,i in zip(dstidx,perm): partial[j]=i
                rec(k+1)
        rec(0)
    return out

ap=argparse.ArgumentParser()
ap.add_argument('--root',type=Path,default=BASE/'results/coxeter9-slope-numeric-v1/root-000029')
ap.add_argument('--chain',type=Path,default=BASE/'results/rank17-extension-chain-v1')
ap.add_argument('--v-pairings',default='0,0,2,2,2,1,1,1,1')
ap.add_argument('--fixed-restarts',type=int,default=16)
ap.add_argument('--fixed-max-nfev',type=int,default=1200)
ap.add_argument('--joint-max-nfev',type=int,default=3000)
ap.add_argument('--joint-restarts',type=int,default=4)
ap.add_argument('--seed',type=int,default=20260820)
ap.add_argument('--jacobian-step',type=float,default=1e-6)
ap.add_argument('--out',type=Path,default=BASE/'results/rank10-winning-refine-v1')
args=ap.parse_args()
root=args.root.resolve(); chain=args.chain.resolve(); out=args.out.resolve(); out.mkdir(parents=True,exist_ok=True)
pv=tuple(map(int,args.v_pairings.split(',')))
if len(pv)!=9: raise SystemExit('--v-pairings must contain nine integers')
meta=json.loads((root/'candidate.json').read_text()); slopes0=np.loadtxt(root/'slopes.txt').reshape(36,3)
fixed=fixed_from_meta(meta); free0=pack(slopes0,fixed)
_,_,A0,B0,_=reconstruct(slopes0); root_disc,root_jv=branch_norms(A0,B0)
anchors=anchors_from_pv(pv); labels=[f'V{i}' for i in range(9)]+[f'D{i}{j}' for i,j in PAIRS]
print('root =',root)
print('V fingerprint =',' '.join(f'{p:+d}' for p in pv))
print('anchors =',' '.join(f'{labels[i]}:{p:+d}' for i,p in anchors))
print('anchor count =',len(anchors))
print('root raw =',meta.get('raw_max'))
print('root branch norms: delta=%.6e jvar=%.6e'%(root_disc,root_jv))

# Record the remaining S9 ambiguity for rank 11.
H=np.loadtxt(BASE/'data/lattice/short_vector_basis_gram.txt',dtype=np.int64)
Bsel=np.loadtxt(chain/'selected-vectors.txt',dtype=np.int64)
orig=tuple(map(int,((Bsel[9]@H)@Bsel[:9].T).tolist()))
mappings=compatible_mappings(orig,pv)
print('compatible full-lattice Coxeter mappings =',len(mappings))
with (out/'compatible-coxeter-mappings.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['q_sign','numeric_to_lattice'])
    for qsign,m in mappings: w.writerow([qsign,' '.join(map(str,m))])

rng=np.random.default_rng(args.seed)
fixed_score,fixed_z=fixed_surface_start(slopes0,pv,rng,args.fixed_restarts,args.fixed_max_nfev)
print('fixed-surface recovered score = %.6e'%fixed_score)

best=None
for restart in range(args.joint_restarts):
    z0=np.concatenate([free0,fixed_z])
    if restart:
        z0[:104] += rng.normal(0,1e-5,size=104)
        z0[104:] += rng.normal(0,1e-5,size=8)
    # A tiny proximity term only selects the nearest point on the positive-dimensional rank-10 locus.
    def fun(z):
        raw=joint_raw(z,fixed,anchors)
        prox=1e-14*(z[:104]-free0)
        return np.concatenate([raw,prox])
    res=least_squares(fun,z0,method='trf',jac='2-point',max_nfev=args.joint_max_nfev,
        ftol=1e-14,xtol=1e-14,gtol=1e-14,x_scale='jac',verbose=0)
    stats,data=diagnostics(res.x,fixed,anchors,root_disc,root_jv)
    score=stats['raw_max']
    healthy=stats['delta_ratio']>1e-3 and stats['jvar_ratio']>1e-3
    print('REFINE|restart=%d|score=%.3e|base=%.3e|curve=%.3e|line=%.3e|delta_ratio=%.3e|jvar_ratio=%.3e|healthy=%d'%(
        restart,score,stats['base_max'],stats['curve_max'],stats['line_max'],stats['delta_ratio'],stats['jvar_ratio'],int(healthy)),flush=True)
    if healthy and (best is None or score<best[0]): best=(score,res.x.copy(),stats,data)

if best is None: raise SystemExit('no healthy joint rank-10 refinement found')
score,z,stats,data=best; s,xs,ys,A,B,xq,yq=data
sv=jacobian_singular_values(z,fixed,anchors,args.jacobian_step)
smax=float(sv[0]); nullities={}
for e in (6,8,10,12):
    tol=(10.**(-e))*max(smax,1e-300); nullities[f'nullity_1e-{e}']=len(z)-int(np.count_nonzero(sv>tol))
print('JACOBIAN|null6=%d|null8=%d|null10=%d|null12=%d|smax=%.3e|smin=%.3e'%(
    nullities['nullity_1e-6'],nullities['nullity_1e-8'],nullities['nullity_1e-10'],nullities['nullity_1e-12'],smax,float(sv[-1])))

np.save(out/'state.npy',z); np.savetxt(out/'slopes.txt',s,fmt='%.17g'); np.savetxt(out/'x-coxeter45.txt',xs,fmt='%.17g'); np.savetxt(out/'y-coxeter45.txt',ys,fmt='%.17g'); np.savetxt(out/'x-q.txt',xq.reshape(1,-1),fmt='%.17g'); np.savetxt(out/'y-q.txt',yq.reshape(1,-1),fmt='%.17g'); np.savetxt(out/'A.txt',A.reshape(1,-1),fmt='%.17g'); np.savetxt(out/'B.txt',B.reshape(1,-1),fmt='%.17g'); np.savetxt(out/'jacobian-singular-values.txt',sv.reshape(1,-1),fmt='%.17g')
payload={'version':1,'kind':'rank10-winning-fingerprint-refinement','source_root':str(root),'v_pairings':list(pv),'anchors':[(labels[i],p) for i,p in anchors],'compatible_mappings':len(mappings),**stats,**nullities}
(out/'candidate.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
print('\nBEST RANK10 REFINEMENT')
for k,v in stats.items(): print(k,'=',v)
for k,v in nullities.items(): print(k,'=',v)
print('compatible_mappings =',len(mappings))
print('saved =',out)
