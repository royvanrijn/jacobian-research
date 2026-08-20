from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser()
ap.add_argument("--frame",required=True)
ap.add_argument("--name",required=True)
args=ap.parse_args()

F=matrix(ZZ,[[ZZ(x) for x in ln.split()] for ln in Path(args.frame).read_text().splitlines() if ln.strip()])
assert F.nrows()==17 and F.det()==948 and F.is_positive_definite()

def qform(H):
    co=[]
    for i in range(H.nrows()):
        for j in range(i,H.ncols()):
            co.append(H[i,i]//2 if i==j else H[i,j])
    return QuadraticForm(ZZ,H.nrows(),co)

half=[vector(ZZ,r) for r in qform(F).short_vector_list_up_to_length(2,True)[1]]
signed=half+[-r for r in half]
Rmod=matrix(ZZ,[list(r) for r in signed]).row_module()
RB=Rmod.basis_matrix()
GR=RB*F*RB.transpose()
rr=RB.rank()
print(f"BRANCH|name={args.name}|stage=roots|rank={rr}|roots={2*len(half)}|rootdet={abs(GR.det())}",flush=True)

# exact essential intersection
C=(RB*F).right_kernel_matrix()
GC=C*F*C.transpose()
print(f"BRANCH|name={args.name}|stage=intersection|rank={C.rank()}|det={GC.det()}|gram={GC}",flush=True)

A=block_matrix([[RB],[C]],subdivide=False)
idx=abs(A.det())
D,U,V=A.smith_form()
smith=[abs(ZZ(D[i,i])) for i in range(17)]
print(f"BRANCH|name={args.name}|stage=glue|index={idx}|smith={smith}",flush=True)

Ainv=A.inverse()
def mod1(q):
    q=QQ(q); return q-floor(q)
def key(z):
    c=vector(QQ,z)*Ainv
    return tuple(mod1(x) for x in c)

zero=vector(ZZ,[0]*17)
cosets={key(zero):zero}; queue=[zero]; head=0
while head<len(queue) and len(cosets)<idx:
    z=queue[head]; head+=1
    for i in range(17):
        e=vector(ZZ,[0]*17); e[i]=1
        for s in (1,-1):
            y=z+s*e; k=key(y)
            if k not in cosets:
                cosets[k]=y; queue.append(y)
                if len(cosets)==idx: break
        if len(cosets)==idx: break
print(f"BRANCH|name={args.name}|stage=cosets|count={len(cosets)}",flush=True)

GRinv=GR.inverse(); GCinv=GC.inverse()
def project(x):
    x=vector(QQ,x)
    coeff=(x*F*RB.transpose())*GRinv
    return x-coeff*RB
def ccoords(p):
    return (vector(QQ,p)*F*C.transpose())*GCinv

proj=[ccoords(project(z)) for z in cosets.values()]
rows=[vector(QQ,[1 if i==j else 0 for i in range(C.rank())]) for j in range(C.rank())] + proj
den=lcm([QQ(x).denominator() for r in rows for x in r])
Zrows=matrix(ZZ,[[ZZ(den*x) for x in r] for r in rows])
Bz=Zrows.row_module().basis_matrix()
B=Bz.change_ring(QQ)/den
H=B*GC*B.transpose()
print(f"BRANCH|name={args.name}|stage=MW|rank={H.rank()}|height_gram={H}|det={H.det()}",flush=True)

hden=lcm([QQ(x).denominator() for x in H.list()])
ZH=(hden*H).change_ring(ZZ)
print(f"BRANCH|name={args.name}|stage=scaled|den={hden}|gram={ZH}",flush=True)
try:
    T=matrix(ZZ,pari(ZH).qflllgram())
    cand=[]
    for X in (T,T.transpose()):
        if X.nrows()==H.nrows():
            Red=X*ZH*X.transpose()
            if Red.is_symmetric() and Red.det()==ZH.det():
                cand.append((max(abs(v) for v in Red.list()),Red))
    if cand:
        _,Red=min(cand,key=lambda z:z[0])
        print(f"BRANCH|name={args.name}|stage=reduced|scaled_gram={Red}|scale={hden}",flush=True)
except Exception as e:
    print(f"BRANCH|name={args.name}|reduce_warning={e}",flush=True)

out=Path(f"artifacts/local/elkies-k3/{args.name}-analysis.txt")
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(f"root_rank={rr}\nroot_count={2*len(half)}\nroot_det={abs(GR.det())}\nMW={H}\ndet={H.det()}\n")
print(f"BRANCH|name={args.name}|stage=done|out={out}",flush=True)
