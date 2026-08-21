from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="Transport a Gross CM embedding into recovered K3 transcendental coordinates.")
ap.add_argument("--Tgram",default="artifacts/local/elkies-k3/cm-t2-candidates.txt")
ap.add_argument("--Tid",type=int,default=0)
ap.add_argument("--Ufile",default="artifacts/local/elkies-k3/clifford-class-match.txt")
ap.add_argument("--target",type=int,choices=(3,24,43),default=3)
ap.add_argument("--out",default="artifacts/local/elkies-k3/cm-delta3-k3-vector.txt")
a=ap.parse_args()

D=ZZ(6); LEVEL=ZZ(79); N=D*LEVEL

def load_T(path,tid):
    rows=[]
    for ln in Path(path).read_text().splitlines():
        if ln.startswith(f"TGRAM|{tid}|"):
            rows.append([ZZ(x) for x in ln.split("|",2)[2].split()])
    if len(rows)!=3: raise SystemExit(f"cannot load T={tid}")
    return matrix(ZZ,rows)

def parse_U(path,tid):
    txt=Path(path).read_text()
    marker=f"JACKPOT T={tid}"
    pos=txt.find(marker)
    if pos<0: raise SystemExit(f"cannot find T={tid} jackpot in {path}")
    sub=txt[pos:]
    upos=sub.find("U=")
    if upos<0: raise SystemExit("U= missing")
    chunk=sub[upos+2:]
    rows=[]
    for ln in chunk.splitlines():
        cleaned = ln.replace("["," ").replace("]"," ").strip()
        if not cleaned:
            continue

        vals=[]
        for tok in cleaned.split():
            try:
                vals.append(ZZ(tok))
            except (TypeError, ValueError):
                pass

        if len(vals) == 3:
            rows.append(vals)
            if len(rows)==3:
                break
    if len(rows)!=3: raise SystemExit(f"could not parse U rows: {rows}")
    return matrix(ZZ,rows)

HT=load_T(a.Tgram,a.Tid)
U=parse_U(a.Ufile,a.Tid)
print(f"CMK3|stage=input|Tid={a.Tid}|detHT={HT.det()}|detU={U.det()}",flush=True)

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

# Dual O^sharp for Trd(xy).
P=matrix(QQ,4,4,lambda i,j:tr(ob[i]*ob[j]))
Dcoef=P.inverse().transpose()
dual=[sum((Dcoef[r,j]*ob[j] for j in range(4)),A(0)) for r in range(4)]

tv=vector(QQ,[tr(x) for x in dual])
den=lcm([x.denominator() for x in tv])
ti=vector(ZZ,[ZZ(den*x) for x in tv])
K=matrix(ZZ,1,4,list(ti)).right_kernel_matrix()
sharp0=[sum((ZZ(K[r,j])*dual[j] for j in range(4)),A(0)) for r in range(3)]

# Inverse-Clifford quadratic form qC = -N*nrd on sharp0.
GC=matrix(QQ,3,3)
for i in range(3): GC[i,i]=-N*nr(sharp0[i])
for i in range(3):
 for j in range(i+1,3):
  GC[i,j]=GC[j,i]=-N*(nr(sharp0[i]+sharp0[j])-nr(sharp0[i])-nr(sharp0[j]))/2
HC=(2*GC).change_ring(ZZ)

assert U*HC*U.transpose()==HT
print(f"CMK3|stage=isometry|status=OK|HCdet={HC.det()}",flush=True)

# Exact primitive Gross vectors previously found.
gens=A.gens(); ii,jj,kk=gens[0],gens[1],gens[2]
beta_by_target = {
    3: -3*ii-jj+kk,
    24: -78*ii+38*jj-24*kk,
    # Gross coordinates (-11,5,-1) in the exact basis
    # (i+j+37*k, 2*j+112*k, 158*k).
    43: -11*ii-jj-5*kk,
}
beta=beta_by_target[a.target]
assert nr(beta)==a.target
print(f"CMK3|stage=beta|beta={beta}|norm={nr(beta)}",flush=True)

# Functional lambda(y)=Trd(beta*y) on (O^sharp)^0.
lam=vector(QQ,[tr(beta*y) for y in sharp0])
print(f"CMK3|functional={tuple(lam)}",flush=True)

# HC is the bilinear Hessian for qC. The representing vector r satisfies
#     r * HC * z^T = lambda(z)
# for all integral coordinate vectors z.
# Thus r = lambda * HC^{-1}.
rC=lam*HC.inverse()
print(f"CMK3|clifford_rational_vector={tuple(rC)}",flush=True)

# Coordinate convention from class matching:
# rows of U are the HT basis vectors expressed in HC coordinates:
#     x_K3 -> z_Clifford = x_K3 * U
# therefore x_K3 = z_Clifford * U^{-1}.
rT=rC*U.inverse()
print(f"CMK3|k3_rational_vector={tuple(rT)}",flush=True)

# Primitive integral vector spanning the same rational line.
dens=[QQ(x).denominator() for x in rT]
scale=lcm(dens)
vi=vector(ZZ,[ZZ(scale*x) for x in rT])
g=gcd([abs(x) for x in vi if x])
if g>1: vi=vector(ZZ,[x//g for x in vi])
# canonical sign
for x in vi:
    if x!=0:
        if x<0: vi=-vi
        break

normH=ZZ(vi*HT*vi)
div=gcd([abs(ZZ(x)) for x in HT*vi])
print(f"CMK3|Delta=-{a.target}|k3_vector={tuple(vi)}|normH={normH}|qT={QQ(normH)/2}|div={div}",flush=True)

# Orthogonal complement in T.
row=matrix(ZZ,1,3,list(vi*HT))
KB=row.right_kernel_matrix()
if KB.nrows()!=2: raise RuntimeError("orthogonal complement rank !=2")
G2=KB*HT*KB.transpose()
print(f"CMK3|complement_det={G2.det()}|gram={G2}",flush=True)

# Reduce binary form associated to even complement Hessian.
if G2[0,0]%2==0 and G2[1,1]%2==0:
    aa=ZZ(G2[0,0]//2); bb=ZZ(G2[0,1]); cc=ZZ(G2[1,1]//2)
    cont=gcd([abs(aa),abs(bb),abs(cc)])
    ap,bp,cp=aa//cont,bb//cont,cc//cont
    try:
        red=BinaryQF([ap,bp,cp]).reduced_form()
        redt=tuple(map(int,red))
    except Exception:
        redt=(int(ap),int(bp),int(cp))
    disc=bp*bp-4*ap*cp
    print(f"CMK3|binary={redt}|content={cont}|binary_disc={disc}",flush=True)

# Determinant identity for primitive vector.
rhs = QQ(abs(HT.det()) * abs(normH)) / QQ(div*div)
print(f"CMK3|det_identity|lhs={G2.det()}|rhs={rhs}|ok={QQ(G2.det())==rhs}",flush=True)

Path(a.out).parent.mkdir(parents=True,exist_ok=True)
Path(a.out).write_text(
    f"beta={beta}\n"
    f"functional={tuple(lam)}\n"
    f"clifford_rational_vector={tuple(rC)}\n"
    f"k3_rational_vector={tuple(rT)}\n"
    f"k3_vector={tuple(vi)}\n"
    f"normH={normH}\n"
    f"div={div}\n"
    f"complement={G2}\n"
)
print(f"CMK3|stage=done|out={a.out}",flush=True)
