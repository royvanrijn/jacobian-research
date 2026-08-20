#!/usr/bin/env python3
"""Rank-9 -> rank-10 continuation using the full 45-section Coxeter scaffold.

The first continuation experiment only anchored the new lattice section against
the nine Coxeter generators.  But the slope reconstruction gives explicit
coordinates for 36 further minimal sections D_ij = V_i - V_j.  The target
rank-10 lattice vector has a known height pairing with all 45 sections, so all
|pairing|=2 links can be imposed as quadratic chord-line identities.

This solver also initializes the Coxeter surface only along the numerical
9-dimensional tangent space of the non-isotrivial root, instead of perturbing
all 104 gauge-fixed variables.  Weak scale-relative branch guards prevent the
least-squares solve from escaping to the ubiquitous isotrivial/cuspidal locus.
"""
from __future__ import annotations
from pathlib import Path
import argparse, json
import numpy as np
from scipy.optimize import least_squares

BASE=Path(__file__).resolve().parents[1]
PAIRS=[(i,j) for i in range(9) for j in range(i+1,9)]
PI={p:n for n,p in enumerate(PAIRS)}
TRIPLES=[(i,j,k) for i in range(9) for j in range(i+1,9) for k in range(j+1,9)]
TINC=np.zeros((84,9)); EINC=np.zeros((36,9))
for r,t in enumerate(TRIPLES): TINC[r,list(t)]=1
for r,(i,j) in enumerate(PAIRS): EINC[r,[i,j]]=1
TP=np.linalg.inv(TINC.T@TINC)@TINC.T
EP=np.linalg.inv(EINC.T@EINC)@EINC.T

def pmul(a,b): return np.convolve(a,b)
def p3(a): return pmul(pmul(a,a),a)
def pder(a): return np.asarray([k*a[k] for k in range(1,len(a))]) if len(a)>1 else np.zeros(1)
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
def base_raw(z,fixed):
    *_,blocks=reconstruct(unpack(z,fixed))
    return np.concatenate([blocks[0].ravel(),blocks[1].ravel(),.5*blocks[2].ravel(),.25*blocks[3].ravel()])
def branch_norms(A,B):
    disc=padd(p3(A),pmul(B,B),4.,27.)
    jv=padd(pmul(pder(A),B),pmul(A,pder(B)),3.,-2.)
    return float(np.linalg.norm(disc)),float(np.linalg.norm(jv))
def derived45(s):
    x,y,A,B,blocks=reconstruct(s)
    xs=[x[i].copy() for i in range(9)]; ys=[y[i].copy() for i in range(9)]
    labels=[f'V{i}' for i in range(9)]
    for r,(i,j) in enumerate(PAIRS):
        m=s[r]; xd=pmul(m,m)-x[i]-x[j]
        yd=-y[i]-pmul(m,xd-x[i])
        xs.append(xd); ys.append(yd); labels.append(f'D{i}{j}')
    return np.asarray(xs),np.asarray(ys),labels,A,B,blocks
def lattice45(chain_vectors):
    V=chain_vectors[:9]; nodes=[V[i].copy() for i in range(9)]
    for i,j in PAIRS: nodes.append(V[i]-V[j])
    return np.asarray(nodes,dtype=np.int64)
def tangent_basis(rootfree,fixed,step=2e-6,nullity=9):
    m=len(base_raw(rootfree,fixed)); J=np.empty((m,len(rootfree)))
    for c in range(len(rootfree)):
        h=step*max(1.,abs(rootfree[c])); a=rootfree.copy(); b=rootfree.copy(); a[c]+=h; b[c]-=h
        J[:,c]=(base_raw(a,fixed)-base_raw(b,fixed))/(2*h)
    _,sv,Vt=np.linalg.svd(J,full_matrices=False)
    return Vt[-nullity:],sv
def decode(z,fixed,nanchors):
    s=unpack(z[:104],fixed); xs,ys,labels,A,B,blocks=derived45(s)
    xq=z[104:109]; slopes=z[109:].reshape(nanchors,3)
    return s,xs,ys,labels,A,B,blocks,xq,slopes
def make_residual(z,fixed,anchors,root_disc,root_jv,args,limit=None,proximity=None):
    n=len(anchors); s,xs,ys,labels,A,B,blocks,xq,ms=decode(z,fixed,n)
    use=n if limit is None else min(n,limit)
    ai,pairing,_=anchors[0]; yq=pmul(ms[0],xq-xs[ai])-(pairing/2.)*ys[ai]
    pieces=[blocks[0].ravel(),blocks[1].ravel(),.5*blocks[2].ravel(),.25*blocks[3].ravel()]
    pieces.append(args.curve_weight*(pmul(yq,yq)-p3(xq)-pmul(A,xq)-B))
    for k in range(1,use):
        ai,pairing,_=anchors[k]
        pieces.append(yq+(pairing/2.)*ys[ai]-pmul(ms[k],xq-xs[ai]))
    dn,jn=branch_norms(A,B); dr=dn/max(root_disc,1e-300); jr=jn/max(root_jv,1e-300)
    pieces.append(np.asarray([args.branch_weight*max(0.,args.min_branch_ratio-dr),args.branch_weight*max(0.,args.min_branch_ratio-jr)]))
    if proximity is not None and args.proximity_weight>0: pieces.append(args.proximity_weight*(z[:104]-proximity))
    return np.concatenate(pieces)
def diagnostics(z,fixed,anchors,root_disc,root_jv):
    n=len(anchors); s,xs,ys,labels,A,B,blocks,xq,ms=decode(z,fixed,n)
    ai,pairing,_=anchors[0]; yq=pmul(ms[0],xq-xs[ai])-(pairing/2.)*ys[ai]
    base=max(float(np.max(np.abs(b))) for b in blocks)
    curve=float(np.max(np.abs(pmul(yq,yq)-p3(xq)-pmul(A,xq)-B)))
    line=0.
    for k,(ai,pairing,_) in enumerate(anchors[1:],start=1):
        e=yq+(pairing/2.)*ys[ai]-pmul(ms[k],xq-xs[ai]); line=max(line,float(np.max(np.abs(e))))
    dn,jn=branch_norms(A,B); dr=dn/max(root_disc,1e-300); jr=jn/max(root_jv,1e-300)
    return base,curve,line,dr,jr,(s,xs,ys,A,B,xq,yq)

ap=argparse.ArgumentParser()
ap.add_argument('--root',type=Path,default=BASE/'results/coxeter9-slope-numeric-v1/root-000001')
ap.add_argument('--chain',type=Path,default=BASE/'results/rank17-extension-chain-v1')
ap.add_argument('--restarts',type=int,default=24)
ap.add_argument('--max-nfev',type=int,default=2200)
ap.add_argument('--seed',type=int,default=20260820)
ap.add_argument('--tangent-scale',type=float,default=.25)
ap.add_argument('--curve-weight',type=float,default=.5)
ap.add_argument('--branch-weight',type=float,default=1e-6)
ap.add_argument('--min-branch-ratio',type=float,default=1e-2)
ap.add_argument('--proximity-weight',type=float,default=1e-10)
ap.add_argument('--out',type=Path,default=BASE/'results/rank10-full-coxeter-v2')
args=ap.parse_args()
root=args.root.resolve(); chain=args.chain.resolve(); out=args.out.resolve(); out.mkdir(parents=True,exist_ok=True)
meta=json.loads((root/'candidate.json').read_text()); roots=np.loadtxt(root/'slopes.txt').reshape(36,3)
fixed=fixed_from_meta(meta); rootfree=pack(roots,fixed)
_,_,A0,B0,_=reconstruct(roots); root_disc,root_jv=branch_norms(A0,B0)
H=np.loadtxt(BASE/'data/lattice/short_vector_basis_gram.txt',dtype=np.int64)
Bsel=np.loadtxt(chain/'selected-vectors.txt',dtype=np.int64); qvec=Bsel[9]
L45=lattice45(Bsel); pairings=np.asarray([(qvec@H)@v for v in L45],dtype=int)
labels=[f'V{i}' for i in range(9)]+[f'D{i}{j}' for i,j in PAIRS]
anchors=[(i,int(pairings[i]),labels[i]) for i in range(45) if abs(int(pairings[i]))==2]
anchors.sort(key=lambda a:(0 if a[0]<9 else 1,a[0]))
if len(anchors)<2: raise SystemExit('too few full-Coxeter anchors')
print('target signed index =',np.loadtxt(chain/'selected-signed-indices.txt',dtype=int).reshape(-1)[9])
print('full Coxeter abs2 anchors =',len(anchors))
print('anchors =',' '.join(f'{lab}:{p:+d}' for _,p,lab in anchors))
N,sv=tangent_basis(rootfree,fixed)
print('base tangent nullity used =',len(N))
print('jacobian singular tail =',' '.join(f'{x:.3e}' for x in sv[-12:]))
print('root branch norms: delta=%.6e jvar=%.6e'%(root_disc,root_jv)); print()
rng=np.random.default_rng(args.seed); best=None; n=len(anchors)
for restart in range(args.restarts):
    b0=rootfree.copy() if restart==0 else rootfree+rng.normal(0,args.tangent_scale,size=len(N))@N
    proj=least_squares(lambda z:base_raw(z,fixed),b0,method='trf',jac='2-point',max_nfev=500,ftol=1e-12,xtol=1e-12,gtol=1e-12,x_scale='jac')
    b0=proj.x; s0=unpack(b0,fixed); xs0,ys0,_,_,_,_=derived45(s0)
    z=np.zeros(104+5+3*n); z[:104]=b0
    coeff=rng.normal(size=9); coeff/=max(np.linalg.norm(coeff),1e-12); z[104:109]=coeff@xs0[:9]+rng.normal(0,.2,size=5)
    slope_scale=max(.25,float(np.median(np.abs(s0)))); z[109:]=rng.normal(0,slope_scale,size=3*n)
    stages=[]
    for lim in (2,5,None):
        fun=lambda zz,L=lim:make_residual(zz,fixed,anchors,root_disc,root_jv,args,L,rootfree)
        res=least_squares(fun,z,method='trf',jac='2-point',max_nfev=args.max_nfev,ftol=1e-12,xtol=1e-12,gtol=1e-12,x_scale='jac')
        z=res.x; stages.append((lim,int(res.nfev),float(np.linalg.norm(res.fun))))
    base,curve,line,dr,jr,data=diagnostics(z,fixed,anchors,root_disc,root_jv)
    score=max(base,curve,line); usable=dr>=args.min_branch_ratio and jr>=args.min_branch_ratio
    isbest=usable and (best is None or score<best[0])
    print(f'FULL45|restart={restart}|score={score:.3e}|base={base:.3e}|curve={curve:.3e}|line={line:.3e}|delta_ratio={dr:.3e}|jvar_ratio={jr:.3e}|usable={int(usable)}'+('|BEST' if isbest else ''),flush=True)
    if isbest: best=(score,z.copy(),base,curve,line,dr,jr,stages,data)
if best is None: raise SystemExit('no branch-preserving full-Coxeter rank-10 root found')
score,z,base,curve,line,dr,jr,stages,data=best; s,xs,ys,A,B,xq,yq=data
np.save(out/'state.npy',z); np.savetxt(out/'slopes.txt',s,fmt='%.17g'); np.savetxt(out/'x-coxeter45.txt',xs,fmt='%.17g'); np.savetxt(out/'y-coxeter45.txt',ys,fmt='%.17g'); np.savetxt(out/'x-q.txt',xq.reshape(1,-1),fmt='%.17g'); np.savetxt(out/'y-q.txt',yq.reshape(1,-1),fmt='%.17g'); np.savetxt(out/'A.txt',A.reshape(1,-1),fmt='%.17g'); np.savetxt(out/'B.txt',B.reshape(1,-1),fmt='%.17g')
(out/'candidate.json').write_text(json.dumps({'score':score,'base_max':base,'curve_max':curve,'line_max':line,'delta_ratio':dr,'jvar_ratio':jr,'anchors':[(lab,p) for _,p,lab in anchors],'stages':stages},indent=2)+'\n')
print('\nBEST FULL-COXETER RANK10'); print('score =',score); print('delta_ratio =',dr); print('jvar_ratio =',jr); print('saved =',out)
