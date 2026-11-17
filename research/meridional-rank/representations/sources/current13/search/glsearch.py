"""Representations of a knot group into GL(4,q), q prime, with meridians -> pseudo-reflections
I + v f (f(v) = lam - 1), seed vectors e_1..e_4, seed functionals the rows of A (det A != 0)."""
import ast,time,sys,json,itertools,numpy as np
sys.path.insert(0,'.')
from quandle import crossings, plan; from wirt import load
OFF=[(i,j) for i in range(4) for j in range(4) if i!=j]
PATH=[(0,1),(1,2),(2,3)]
FREE=[p for p in OFF if p not in PATH]
def det4(A,q):
    # A: (K,4,4) int64 -> det mod q via Leibniz (24 terms)
    s=np.zeros(A.shape[0],np.int64)
    for p in itertools.permutations(range(4)):
        inv=sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])
        t=np.ones(A.shape[0],np.int64)
        for i in range(4): t=(t*A[:,i,p[i]])%q
        s=(s+(-1)**inv*t)%q
    return s
def search_gl(n,cr,q,lam,chunk=500_000,want=1):
    s,steps=plan(n,cr,4)
    linv=pow(lam,-1,q)
    found=[]
    for pathbits in itertools.product((0,1),repeat=3):
        total=q**len(FREE)
        for start in range(0,total,chunk):
            idx=np.arange(start,min(total,start+chunk),dtype=np.int64)
            K=len(idx)
            A=np.zeros((K,4,4),np.int64)
            for i in range(4): A[:,i,i]=(lam-1)%q
            for k,(i,j) in enumerate(FREE): A[:,i,j]=(idx//q**k)%q
            for b,(i,j) in zip(pathbits,PATH): A[:,i,j]=b
            keep=det4(A,q)!=0; A=A[keep]; K=len(A)
            if K==0: continue
            V=[None]*n; F=[None]*n
            for i,si in enumerate(s):
                v=np.zeros((K,4),np.int64); v[:,i]=1; V[si]=v; F[si]=A[:,i,:].copy()
            def conj(o,a,e):
                fv=(F[o]*V[a]).sum(1)%q          # f_o(v_a)
                af=(F[a]*V[o]).sum(1)%q          # f_a(v_o)
                if e==1:
                    v=(V[a]+fv[:,None]*V[o])%q; f=(F[a]-(af*linv%q)[:,None]*F[o])%q
                else:
                    v=(V[a]-(fv*linv%q)[:,None]*V[o])%q; f=(F[a]+af[:,None]*F[o])%q
                return v,f
            for t,o,src,e in steps:
                V[t],F[t]=conj(o,src,e)
            ok=np.ones(K,bool)
            for o,u,v,e in cr:
                vv,ff=conj(o,u,e)
                # compare rank-one matrices vv ff^T and V[v] F[v]^T
                d=(vv[:,:,None]*ff[:,None,:]-V[v][:,:,None]*F[v][:,None,:])%q
                ok&=(d.reshape(K,16)==0).all(1)
                if not ok.any(): break
            for k in np.nonzero(ok)[0][:want]:
                found.append({'A':A[k].tolist(),'V':[V[a][k].tolist() for a in range(n)],'F':[F[a][k].tolist() for a in range(n)]})
            if len(found)>=want: return found
    return found
if __name__=='__main__':
    q=int(sys.argv[1]); lam=int(sys.argv[2]); todo=json.load(open(sys.argv[3])); out=sys.argv[4]
    rows,I=load(); S=set(todo); sel=[x for x in rows if x[I['name']] in S]
    res={}; t=time.time()
    for x in sel:
        pd=ast.literal_eval(x[I['pd_notation']]); n,cr=crossings(pd)
        f=search_gl(n,cr,q,lam)
        if f: res[x[I['name']]]={'q':q,'lam':lam,'kind':'gl',**f[0]}
        json.dump(res,open(out,'w'))
    print('gl',q,lam,len(sel),len(res),round(time.time()-t,1),flush=True)
