from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="Second-generation low-q neighbor search from q=90 MW7 frame.")
ap.add_argument("--qmin",type=int,default=2)
ap.add_argument("--qmax",type=int,default=16)
ap.add_argument("--enum-cap",type=int,default=4000)
ap.add_argument("--report",type=int,default=20)
a=ap.parse_args()

FRAME=Path("elkies-k3/data/fibrations/q90_mw7_frame.txt")
M=matrix(ZZ,[[ZZ(x) for x in ln.split()] for ln in FRAME.read_text().splitlines() if ln.strip()])
assert M.nrows()==17 and M.det()==948
U=matrix(ZZ,[[0,1],[1,0]])
NS=block_diagonal_matrix(U,-M)
print(f"Q90NBR|stage=start|detM={M.det()}|q={a.qmin}..{a.qmax}",flush=True)

def qform_from_hessian(H):
    co=[]
    for i in range(H.nrows()):
        for j in range(i,H.ncols()):
            co.append(H[i,i]//2 if i==j else H[i,j])
    return QuadraticForm(ZZ,H.nrows(),co)

# Exact-norm Cholesky DFS sampler.
def vectors_of_qnorm(n):
    target=ZZ(2*n)
    MR=M.change_ring(RR)
    Rchol=MR.cholesky().transpose()
    dim=M.nrows(); x=[ZZ(0)]*dim; found=[]; seen=set()
    eps=RR("1e-10")
    def rec(k,partial):
        if len(found)>=a.enum_cap: return
        if k<0:
            v=vector(ZZ,x)
            if v!=0 and ZZ(v*M*v)==target:
                tup=tuple(v); neg=tuple(-v); canon=min(tup,neg)
                if canon not in seen:
                    seen.add(canon); found.append(vector(ZZ,canon))
            return
        off=sum(Rchol[k,j]*RR(x[j]) for j in range(k+1,dim))
        rem=RR(target)-partial
        if rem < -eps: return
        rad=sqrt(max(RR(0),rem)); diag=Rchol[k,k]
        lo=ceil((-rad-off)/diag-eps); hi=floor((rad-off)/diag+eps)
        center=ZZ(round(-off/diag)); vals=[]
        if lo<=center<=hi: vals.append(center)
        d=1
        while center-d>=lo or center+d<=hi:
            if center-d>=lo: vals.append(center-d)
            if center+d<=hi: vals.append(center+d)
            d+=1
        for z in vals:
            x[k]=ZZ(z)
            y=diag*RR(z)+off
            rec(k-1,partial+y*y)
            if len(found)>=a.enum_cap: break
        x[k]=ZZ(0)
    rec(dim-1,RR(0))
    print(f"Q90NBR|stage=enumerate|q={n}|pairs={len(found)}",flush=True)
    return found

def primitive_div(v):
    return gcd([abs(ZZ(x)) for x in NS*v])

def bezout_vector_for_pairing(f):
    p=list(NS*f); cur=ZZ(0); vec=[ZZ(0)]*19
    for i,pi in enumerate(p):
        if pi==0: continue
        gg,s,t=xgcd(cur,ZZ(pi))
        vec=[s*x for x in vec]; vec[i]+=t; cur=gg
    if abs(cur)!=1: return None
    if cur==-1: vec=[-x for x in vec]
    return vector(ZZ,vec)

def build_frame(f):
    g=bezout_vector_for_pairing(f)
    if g is None: return None
    gsq=ZZ(g*NS*g)
    if gsq%2: return None
    g0=g-(gsq//2)*f
    if g0*NS*g0!=0 or f*NS*g0!=1: return None
    K=matrix(ZZ,[list(f*NS),list(g0*NS)]).right_kernel_matrix()
    P=-(K*NS*K.transpose())
    if not P.is_positive_definite(): return None
    return P

def root_data(P):
    Q=qform_from_hessian(P)
    arr=Q.short_vector_list_up_to_length(2,True)
    half=arr[1] if len(arr)>1 else []
    if not half: return (0,0,1)
    R=matrix(ZZ,[list(r) for r in half]); rank=R.rank()
    L=matrix(ZZ,[list(r) for r in list(half)+[tuple(-vector(ZZ,r)) for r in half]]).row_module()
    B=L.basis_matrix(); G=B*P*B.transpose()
    return rank,2*len(half),abs(G.det())

results=[]; seen=set()
for q in range(a.qmin,a.qmax+1):
    vv=vectors_of_qnorm(q)
    for v in vv:
        for aa in divisors(q):
            bb=q//aa
            f=vector(ZZ,[ZZ(aa),ZZ(bb)]+list(v))
            if f*NS*f!=0 or primitive_div(f)!=1: continue
            P=build_frame(f)
            if P is None: continue
            rr,rc,rd=root_data(P)
            try: Pred=P.LLL_gram()
            except Exception: Pred=P
            key=(rr,rc,rd,tuple(Pred.list()))
            if key in seen: continue
            seen.add(key)
            mw=17-rr
            results.append((mw,rr,rc,rd,q,aa,bb,tuple(v),P))
            print(f"Q90NBR|cand={len(results)}|q={q}|ab={aa},{bb}|root_rank={rr}|roots={rc}|rootdet={rd}|MW={mw}",flush=True)

results.sort(key=lambda x:(x[0],-x[2],x[3],x[4]))
print(f"Q90NBR|stage=summary|unique={len(results)}",flush=True)
for i,row in enumerate(results[:a.report],1):
    mw,rr,rc,rd,q,aa,bb,v,P=row
    print(f"Q90NBR_BEST|rank={i}|MW={mw}|root_rank={rr}|roots={rc}|rootdet={rd}|q={q}|ab={aa},{bb}|v={v}",flush=True)
    print(f"Q90NBR_BEST|frame_gram={P}",flush=True)

out=Path("artifacts/local/elkies-k3/q90-second-neighbors.txt")
out.parent.mkdir(parents=True,exist_ok=True)
with out.open("w") as h:
    for i,row in enumerate(results[:a.report],1):
        mw,rr,rc,rd,q,aa,bb,v,P=row
        h.write(f"rank={i} MW={mw} root_rank={rr} roots={rc} rootdet={rd} q={q} ab={aa},{bb} v={v}\nframe=\n{P}\n\n")
print(f"Q90NBR|stage=done|out={out}",flush=True)
