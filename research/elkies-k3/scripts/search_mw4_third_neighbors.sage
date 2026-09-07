from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="Third-generation low-q search from the MW4 frame.")
ap.add_argument("--qmin",type=int,default=2)
ap.add_argument("--qmax",type=int,default=12)
ap.add_argument("--enum-cap",type=int,default=5000)
ap.add_argument("--report",type=int,default=20)
a=ap.parse_args()

FRAME=Path("elkies-k3/data/fibrations/q90_mw4_frame.txt")
M=matrix(ZZ,[[ZZ(x) for x in ln.split()] for ln in FRAME.read_text().splitlines() if ln.strip()])
assert M.det()==948
U=matrix(ZZ,[[0,1],[1,0]])
NS=block_diagonal_matrix(U,-M)
print(f"MW4NBR|stage=start|q={a.qmin}..{a.qmax}|det={M.det()}",flush=True)

def qform(H):
    co=[]
    for i in range(H.nrows()):
        for j in range(i,H.ncols()):
            co.append(H[i,i]//2 if i==j else H[i,j])
    return QuadraticForm(ZZ,H.nrows(),co)

def vectors_of_qnorm(n):
    target=ZZ(2*n)
    MR=M.change_ring(RR); Rchol=MR.cholesky().transpose()
    dim=17; x=[ZZ(0)]*dim; found=[]; seen=set(); eps=RR("1e-10")
    def rec(k,partial):
        if len(found)>=a.enum_cap: return
        if k<0:
            v=vector(ZZ,x)
            if v!=0 and ZZ(v*M*v)==target:
                can=min(tuple(v),tuple(-v))
                if can not in seen:
                    seen.add(can); found.append(vector(ZZ,can))
            return
        off=sum(Rchol[k,j]*RR(x[j]) for j in range(k+1,dim))
        rem=RR(target)-partial
        if rem < -eps:return
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
            x[k]=ZZ(z); y=diag*RR(z)+off
            rec(k-1,partial+y*y)
            if len(found)>=a.enum_cap: break
        x[k]=ZZ(0)
    rec(dim-1,RR(0))
    print(f"MW4NBR|stage=enumerate|q={n}|pairs={len(found)}",flush=True)
    return found

def primitive_div(v): return gcd([abs(ZZ(x)) for x in NS*v])

def bezout(f):
    p=list(NS*f); cur=ZZ(0); vec=[ZZ(0)]*19
    for i,pi in enumerate(p):
        if pi==0:continue
        gg,s,t=xgcd(cur,ZZ(pi)); vec=[s*x for x in vec]; vec[i]+=t; cur=gg
    if abs(cur)!=1:return None
    if cur==-1:vec=[-x for x in vec]
    return vector(ZZ,vec)

def frame(f):
    g=bezout(f)
    if g is None:return None
    gsq=ZZ(g*NS*g)
    if gsq%2:return None
    g0=g-(gsq//2)*f
    if g0*NS*g0!=0 or f*NS*g0!=1:return None
    K=matrix(ZZ,[list(f*NS),list(g0*NS)]).right_kernel_matrix()
    P=-(K*NS*K.transpose())
    return P if P.is_positive_definite() else None

def roots(P):
    A=qform(P).short_vector_list_up_to_length(2,True)
    half=A[1] if len(A)>1 else []
    if not half:return 0,0,1
    R=matrix(ZZ,[list(r) for r in half]); rr=R.rank()
    signed=list(half)+[tuple(-vector(ZZ,r)) for r in half]
    L=matrix(ZZ,[list(r) for r in signed]).row_module()
    B=L.basis_matrix(); G=B*P*B.transpose()
    return rr,2*len(half),abs(G.det())

res=[]; seen=set()
for q in range(a.qmin,a.qmax+1):
    for v in vectors_of_qnorm(q):
        for aa in divisors(q):
            bb=q//aa
            f=vector(ZZ,[ZZ(aa),ZZ(bb)]+list(v))
            if f*NS*f!=0 or primitive_div(f)!=1:continue
            P=frame(f)
            if P is None:continue
            rr,rc,rd=roots(P)
            try:key=(rr,rc,rd,tuple(P.LLL_gram().list()))
            except:key=(rr,rc,rd,tuple(P.list()))
            if key in seen:continue
            seen.add(key)
            mw=17-rr
            res.append((mw,rr,rc,rd,q,aa,bb,tuple(v),P))
            if mw<=4:
                print(f"MW4NBR|cand={len(res)}|q={q}|ab={aa},{bb}|root_rank={rr}|roots={rc}|rootdet={rd}|MW={mw}",flush=True)

res.sort(key=lambda z:(z[0],-z[2],z[3],z[4]))
print(f"MW4NBR|stage=summary|unique={len(res)}",flush=True)
for i,z in enumerate(res[:a.report],1):
    mw,rr,rc,rd,q,aa,bb,v,P=z
    print(f"MW4NBR_BEST|rank={i}|MW={mw}|root_rank={rr}|roots={rc}|rootdet={rd}|q={q}|ab={aa},{bb}|v={v}",flush=True)
    print(f"MW4NBR_BEST|frame_gram={P}",flush=True)
