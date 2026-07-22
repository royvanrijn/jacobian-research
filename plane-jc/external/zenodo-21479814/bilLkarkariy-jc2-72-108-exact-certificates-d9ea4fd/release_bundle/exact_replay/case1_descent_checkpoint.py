#!/usr/bin/env python3
from pathlib import Path
import pickle,time,gc,sys
import exact_core as ec
ROOT=Path(__file__).resolve().parent
CK=ROOT/'case1_checkpoint.pkl'

def bracket(Pi,i,Qj,j):
    return ec.tadd(ec.tscale(ec.tmul(Pi,ec.tder(Qj)),ec.K(i)),
                   ec.tscale(ec.tmul(ec.tder(Pi),Qj),ec.K(-j)))
def support_P(r): return list(range(0,8+r+1))
def support_Q(s): return list(range(1,13)) if s==0 else list(range(0,12+s+1))
def ppterms(v): return sum(len(q.d) for q in v)

def encK(k):
    return [(int(k.p[i].p),int(k.p[i].q)) for i in range(len(k.p))]
def decK(v):
    from flint import fmpq,fmpq_poly
    return ec.K(fmpq_poly([fmpq(n,d) for n,d in v]))
def encPP(q): return [(m,encK(c)) for m,c in q.d.items()]
def decPP(v): return ec.PP({tuple(m):decK(c) for m,c in v})
def encBands(X): return {i:[[encPP(q) for q in band] for band in [v]][0] for i,v in X.items()}
def decBands(X): return {int(i):[decPP(q) for q in v] for i,v in X.items()}
def save_state():
    raw=dict(P=encBands(P),Q=encBands(Q),np=np,all_comp=[encPP(q) for q in all_comp],labels=labels,done=done)
    CK.write_bytes(pickle.dumps(raw,protocol=5))
if CK.exists():
    st=pickle.loads(CK.read_bytes())
    P,Q,np,all_comp,labels,done=decBands(st['P']),decBands(st['Q']),st['np'],[decPP(q) for q in st['all_comp']],st['labels'],st['done']
    print('resume',done,'np',np,'comps',len(all_comp),flush=True)
else:
    B,E,C,F,np=ec.high_layers(); P={2:[ec.PP.const(x) for x in ec.Avec],1:B,0:C}; Q={3:[ec.PP.const(x) for x in ec.Dvec],2:E,1:F}
    all_comp=[];labels=[];done=[]
    save_state()

for k in [1,0,-1,-2,-3]:
    if k in done: continue
    T=time.time(); r=k-2;s=k-1;ps=support_P(r);qs=support_Q(s);ncols=len(ps)+len(qs)
    known=[]
    pairs=[]
    for i,Pi in list(P.items()):
        for j,Qj in list(Q.items()):
            if i+j-1==k:
                q0=time.time(); b=bracket(Pi,i,Qj,j); known=ec.tadd(known,b)
                pairs.append((i,j,len(b),ppterms(b),time.time()-q0))
                print(' k',k,'pair',i,j,'len',len(b),'terms',ppterms(b),'sec',round(time.time()-q0,2),flush=True)
    maxdeg=max(len(known)-1,max(ps)+11,7+max(qs))
    rows=[[ec.K(0) for _ in range(ncols)] for _ in range(maxdeg+1)]
    for col,m in enumerate(ps):
        for j,Dj in enumerate(ec.Dvec):
            deg=m+j-1
            if 0<=deg<=maxdeg: rows[deg][col]+=Dj*(r*j-3*m)
    off=len(ps)
    for jj,n in enumerate(qs):
        for i,Ai in enumerate(ec.Avec):
            deg=i+n-1
            if 0<=deg<=maxdeg: rows[deg][off+jj]+=Ai*(2*n-s*i)
    keep=[];rhs=[]
    for deg,row in enumerate(rows):
        kn=known[deg] if deg<len(known) else ec.PP()
        if any(row) or kn:keep.append(row);rhs.append(-kn)
    print(' k',k,'known total terms',ppterms(known),'solve start',flush=True)
    sol,free,comp,np2=ec.linear_solve(keep,rhs,np)
    Pv=[ec.PP() for _ in range(max(ps)+1)]
    for col,m in enumerate(ps):Pv[m]=sol[col]
    Qv=[ec.PP() for _ in range(max(qs)+1)]
    for jj,n in enumerate(qs):Qv[n]=sol[off+jj]
    P[r]=Pv;Q[s]=Qv
    all_comp.extend(comp);labels.extend([k]*len(comp));np=np2;done.append(k)
    print(f'DONE k={k} rows={len(keep)} cols={ncols} rank={ncols-len(free)} nullity={len(free)} compat={len(comp)} params={np} elapsed={time.time()-T:.1f}s',flush=True)
    save_state()
    gc.collect()
print('ALL DONE labels',labels,'np',np,flush=True)
# serialize equations and Singular inputs
assert np==6 and [labels.count(k) for k in [1,0,-1,-2,-3]]==[0,2,2,4,5]
def mod_sing(poly):
    terms=[]
    for i in range(len(poly)):
        c=poly[i]
        if not c:continue
        num=int(c.p);den=int(c.q);cs=str(num) if den==1 else f'({num}/{den})';mon='' if i==0 else ('u' if i==1 else f'u^{i}')
        term=cs if not mon else mon if c==1 else '-'+mon if c==-1 else f'{cs}*{mon}'
        terms.append(term)
    return '+'.join(terms).replace('+-','-')
minpoly=mod_sing(ec.MOD)
for tag,select in [('case1_all13_exact',list(range(len(all_comp)))),('case1_tail9_exact',[i for i,k in enumerate(labels) if k<=-2])]:
    eqs=[all_comp[i].sing(6) for i in select]
    p=ROOT/(tag+'.sing')
    p.write_text('LIB "resources.lib";\nResources::setcores(1);\nLIB "nfmodstd.lib";\nring R=(0,u),(r,s,h,u1,u2,u3),dp;\nminpoly='+minpoly+';\noption(redSB);\nideal I=\n'+',\n'.join(eqs)+';\nprint("INPUT_SIZE="+string(size(I)));\nideal J=nfmodStd(I);\nprint("SIZE="+string(size(J)));\nJ;\nif(size(J)==1 && J[1]==1){print("UNIT");}else{print("NONUNIT");}\nquit;\n')
    print(tag,p.stat().st_size,flush=True)
(ROOT/'case1_residuals_exact.txt').write_text('\n'.join(f'k={k} R{i}={q.sing(6)}' for i,(k,q) in enumerate(zip(labels,all_comp))))
