from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser()
ap.add_argument("--Tgram",default="artifacts/local/elkies-k3/cm-t2-candidates.txt")
ap.add_argument("--Tid",type=int,default=0)
ap.add_argument("--search-bound",type=int,default=20)
ap.add_argument("--out",default="artifacts/local/elkies-k3/clifford-to-k3-map.txt")
a=ap.parse_args()

D=ZZ(6); LEVEL=ZZ(79); N=D*LEVEL

def load_T(path,tid):
    rows=[]
    for ln in Path(path).read_text().splitlines():
        if ln.startswith(f"TGRAM|{tid}|"):
            rows.append([ZZ(x) for x in ln.split("|",2)[2].split()])
    if len(rows)!=3: raise SystemExit(f"cannot load T={tid}")
    return matrix(ZZ,rows)

HT=load_T(a.Tgram,a.Tid)
QT=HT.change_ring(QQ)/2
print(f"CLMAP|stage=T|detH={HT.det()}|detQ={QT.det()}|four_detQ={4*QT.det()}",flush=True)

A=QuaternionAlgebra(D)
M=A.maximal_order(); mb=list(M.basis())

def nr(x):
    try:return QQ(x.reduced_norm())
    except:return QQ(x.norm())
def tr(x):
    try:return QQ(x.reduced_trace())
    except:return QQ(x.trace())
def qcoords(x):
    try:return vector(QQ,x.coefficient_tuple())
    except:
        try:return vector(QQ,x.list())
        except:return vector(QQ,tuple(x))
def prim4(c): return gcd([abs(ZZ(x)) for x in c])==1

def order_from_alpha(alpha):
    ai=alpha.inverse()
    M2=A.quaternion_order([ai*b*alpha for b in mb])
    O=M.intersection(M2)
    try: idx=O.free_module().index_in(M.free_module())
    except: return None
    return O if idx==LEVEL else None

O=None
for B in range(1,7):
  for c0 in range(-B,B+1):
   for c1 in range(-B,B+1):
    for c2 in range(-B,B+1):
     for c3 in range(-B,B+1):
      c=(c0,c1,c2,c3)
      if c==(0,0,0,0) or not prim4(c): continue
      x=sum((ZZ(c[i])*mb[i] for i in range(4)),A(0))
      if abs(nr(x))!=LEVEL: continue
      OO=order_from_alpha(x)
      if OO is not None: O=OO; break
     if O is not None: break
    if O is not None: break
   if O is not None: break
  if O is not None: break
if O is None: raise RuntimeError("failed Eichler order")
ob=list(O.basis())
print(f"CLMAP|stage=order|basis={' ; '.join(map(str,ob))}",flush=True)

# Dual basis for Trd(x*y)
P=matrix(QQ,4,4,lambda i,j:tr(ob[i]*ob[j]))
Dcoef=P.inverse().transpose()
dual=[sum((Dcoef[r,j]*ob[j] for j in range(4)),A(0)) for r in range(4)]
for i in range(4):
 for j in range(4):
  assert tr(ob[i]*dual[j])==(1 if i==j else 0)
print("CLMAP|stage=dual|status=OK",flush=True)

# Integral kernel of trace on O^sharp
tv=vector(QQ,[tr(x) for x in dual])
den=lcm([x.denominator() for x in tv])
ti=vector(ZZ,[ZZ(den*x) for x in tv])
K=matrix(ZZ,1,4,list(ti)).right_kernel_matrix()
sharp0=[sum((ZZ(K[r,j])*dual[j] for j in range(4)),A(0)) for r in range(3)]
print(f"CLMAP|sharp0_basis={' ; '.join(map(str,sharp0))}",flush=True)

# q = -N*nrd to match signature of K3 T
G=matrix(QQ,3,3)
for i in range(3): G[i,i]=-N*nr(sharp0[i])
for i in range(3):
 for j in range(i+1,3):
  G[i,j]=G[j,i]=-N*(nr(sharp0[i]+sharp0[j])-nr(sharp0[i])-nr(sharp0[j]))/2

print(f"CLMAP|stage=form|G={G}|detG={G.det()}|four_det={4*G.det()}",flush=True)
# G is the symmetric matrix of the quadratic form q(x)=x*G*x.
# For an integral quadratic form G may be half-integral off diagonal.
# The integral lattice Hessian is HG = 2*G.
HGq = 2*G
if not all(QQ(x).denominator() == 1 for x in HGq.list()):
    raise RuntimeError("inverse-Clifford Hessian 2G is not integral")
HG = HGq.change_ring(ZZ)
print(f"CLMAP|HG_det={HG.det()}|HT_det={HT.det()}",flush=True)

# Exact rank-3 isometry search U*HG*U^T=HT.
B=a.search_bound
targets=sorted(set(int(HT[i,i]) for i in range(3)))
vecs={n:[] for n in targets}
for x in range(-B,B+1):
 for y in range(-B,B+1):
  for z in range(-B,B+1):
   v=vector(ZZ,[x,y,z])
   if v==0: continue
   q=int(v*HG*v)
   if q in vecs: vecs[q].append(v)
print("CLMAP|candidate_counts|"+"|".join(f"{n}:{len(vs)}" for n,vs in vecs.items()),flush=True)

rows=[None]*3; found=[None]
def bt(i):
    if i==3:
        U=matrix(ZZ,rows)
        if abs(U.det())==1 and U*HG*U.transpose()==HT:
            found[0]=U; return True
        return False
    for v in vecs[int(HT[i,i])]:
        if all(v*HG*rows[j]==HT[i,j] for j in range(i)):
            rows[i]=v
            if bt(i+1): return True
    return False
bt(0)
U=found[0]
if U is None:
    print(f"CLMAP|NO_ISOMETRY_IN_BOX|bound={B}",flush=True)
else:
    print(f"CLMAP|JACKPOT|U={U}|det={U.det()}",flush=True)
    print(f"CLMAP|Uinv={U.inverse()}",flush=True)

Path(a.out).parent.mkdir(parents=True,exist_ok=True)
Path(a.out).write_text(f"HG={HG}\nHT={HT}\nU={U}\n")
print(f"CLMAP|stage=done|out={a.out}",flush=True)
