from sage.all import *
from pathlib import Path

OUT=Path("artifacts/local/elkies-k3")
OUT.mkdir(parents=True,exist_ok=True)

# Hyperbolic plane, positive Cartan matrices.
U=matrix(ZZ,[[0,1],[1,0]])
E8=matrix(ZZ,[
[2,-1,0,0,0,0,0,0],
[-1,2,-1,0,0,0,0,0],
[0,-1,2,-1,0,0,0,-1],
[0,0,-1,2,-1,0,0,0],
[0,0,0,-1,2,-1,0,0],
[0,0,0,0,-1,2,-1,0],
[0,0,0,0,0,-1,2,0],
[0,0,-1,0,0,0,0,2],
])
A2=matrix(ZZ,[[2,-1],[-1,2]])
assert E8.det()==1 and A2.det()==3

# CM NS lattice for Utsumi No.1: U + E8(-1)^2 + A2(-1).
NS=block_diagonal_matrix(U,-E8,-E8,-A2)
assert NS.nrows()==20 and NS.det()==-3
print(f"CMEMBED|stage=cm_ns|rank={NS.nrows()}|det={NS.det()}",flush=True)

# Primitive vector of E8 norm 316. Coordinates are in the simple-root basis above.
# Found deterministically/verified below; its pairing gcd is 1, hence primitive/divisibility 1
# in the unimodular E8 summand.
ve8=vector(ZZ,[9,-2,-10,-7,-4,4,3,-7])
assert ve8*E8*ve8==316
assert gcd([abs(x) for x in E8*ve8])==1

w=vector(ZZ,[0]*20)
for i,x in enumerate(ve8):
    w[2+i]=x
assert w*NS*w==-316
div=gcd([abs(x) for x in NS*w])
assert div==1
print(f"CMEMBED|stage=extra_class|w={tuple(w)}|square={w*NS*w}|div={div}",flush=True)

# Integral orthogonal complement S = w^\perp in NS.
row=matrix(ZZ,1,20,list(w*NS))
K=row.right_kernel_matrix()
assert K.nrows()==19
S=K*NS*K.transpose()
print(f"CMEMBED|stage=generic_ns|rank={S.nrows()}|det={S.det()}",flush=True)
assert abs(S.det())==948

# Compare discriminant form with the target S19=U+(-M).
Mpath=Path("elkies-k3/data/lattice/rank17_gram.txt")
M=matrix(ZZ,[[ZZ(x) for x in ln.split()] for ln in Mpath.read_text().splitlines() if ln.strip() and not ln.startswith("#")])
TARGET=block_diagonal_matrix(U,-M)
assert TARGET.det()==S.det()

def class_order(Hinv,z):
    x=Hinv*vector(ZZ,z)
    return lcm([QQ(q).denominator() for q in x])

def frac2(q):
    q=QQ(q)
    return q-2*floor(q/2)

def cyclic_generator_q(H):
    Hi=H.inverse(); n=abs(H.det())
    # Search tiny sparse vectors for full-order discriminant class.
    for i in range(H.nrows()):
        z=[0]*H.nrows(); z[i]=1
        if class_order(Hi,z)==n:
            zz=vector(ZZ,z)
            return frac2(zz*Hi*zz)
    for i in range(H.nrows()):
      for j in range(i+1,H.nrows()):
       for si in (-1,1):
        for sj in (-1,1):
         z=[0]*H.nrows(); z[i]=si; z[j]=sj
         if class_order(Hi,z)==n:
            zz=vector(ZZ,z)
            return frac2(zz*Hi*zz)
    raise RuntimeError("no cyclic generator found")

qS=cyclic_generator_q(S)
qT=cyclic_generator_q(TARGET)
n=948
hits=[u for u in range(1,n) if gcd(u,n)==1 and frac2(QQ(u*u)*qS-qT)==0]
print(f"CMEMBED|stage=discform|qS={qS}|qTarget={qT}|unit_hits={hits}",flush=True)
assert hits

# Roots inherited from the CM No.1 trivial lattice.
# The entire second E8 and A2 are orthogonal to w, hence survive in S.
# Test whether any roots from the first E8 survive as well: vectors r with r^2=2 and <r,ve8>=0.
Q8=QuadraticForm(ZZ,8,[E8[i,i]//2 if i==j else E8[i,j] for i in range(8) for j in range(i,8)])
sv=Q8.short_vector_list_up_to_length(2,True)
roots=sv[1]  # Q=1 means E8 norm 2
orth=[]
for r in roots:
    rv=vector(ZZ,r)
    if rv*E8*ve8==0:
        orth.append(tuple(rv))
print(f"CMEMBED|stage=first_E8_roots|total={len(roots)}|orthogonal_to_w={len(orth)}",flush=True)

# Root rank from inherited configuration:
# second E8 contributes 8, A2 contributes 2, plus rank of first-E8 roots orthogonal to w.
if orth:
    Rmat=matrix(ZZ,[list(r) for r in orth])
    rrank=Rmat.rank()
else:
    rrank=0
inherited_root_rank=8+2+rrank
mw_rank=19-2-inherited_root_rank
print(f"CMEMBED|stage=fibration_prediction|firstE8_root_rank={rrank}|inherited_root_rank={inherited_root_rank}|predicted_MW_rank={mw_rank}",flush=True)

# Save matrices and a concise result.
(OUT/"generic-ns-from-disc3-cm.txt").write_text(
    "w="+" ".join(map(str,w))+"\n"
    +"S=\n"+"\n".join(" ".join(map(str,row)) for row in S.rows())+"\n"
    +f"qS={qS}\nqTarget={qT}\nunit_hits={hits}\n"
    +f"firstE8_root_rank={rrank}\ninherited_root_rank={inherited_root_rank}\npredicted_MW_rank={mw_rank}\n"
)
print("CMEMBED|stage=done|status=OK",flush=True)
