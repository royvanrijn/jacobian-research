from sage.all import *
from pathlib import Path
import argparse
from collections import deque

ap=argparse.ArgumentParser(description="Compute Atkin-Lehner action on norm-3 and norm-24 Gross vectors for Eichler order (6,79).")
ap.add_argument("--search-bound",type=int,default=12,help="coefficient box for normalizer search")
ap.add_argument("--hit-bound",type=int,default=80,help="Gross-coordinate box to find seed vectors")
ap.add_argument("--targets",default="3,24")
ap.add_argument("--out",default="artifacts/local/elkies-k3/atkin-lehner-cm-orbits.txt")
a=ap.parse_args()

D=ZZ(6); N=ZZ(79)
targets=[ZZ(x) for x in a.targets.split(",") if x.strip()]
A=QuaternionAlgebra(D)
M=A.maximal_order()
mb=list(M.basis())

def qcoords(x):
    try: return vector(QQ,x.coefficient_tuple())
    except Exception:
        try: return vector(QQ,x.list())
        except Exception: return vector(QQ,tuple(x))

def nr(x):
    try: return QQ(x.reduced_norm())
    except Exception: return QQ(x.norm())

def tr(x):
    try: return QQ(x.reduced_trace())
    except Exception:
        try: return QQ(x.trace())
        except Exception: return QQ(2*qcoords(x)[0])

def primitive4(c):
    return gcd([abs(ZZ(x)) for x in c])==1

# Reconstruct exactly the Eichler order used by the successful previous experiment.
# Search a primitive maximal-order element of norm +-79; first one giving index 79 wins.
def order_from_alpha(alpha):
    ainv=alpha.inverse()
    M2=A.quaternion_order([ainv*b*alpha for b in mb])
    try:
        if not M2.is_maximal(): return None
    except Exception:
        pass
    O=M.intersection(M2)
    try:
        idx=O.free_module().index_in(M.free_module())
    except Exception:
        return None
    return O if idx==N else None

O=None; alpha_seed=None
for B in range(1,7):
  for c0 in range(-B,B+1):
   for c1 in range(-B,B+1):
    for c2 in range(-B,B+1):
     for c3 in range(-B,B+1):
      c=(c0,c1,c2,c3)
      if c==(0,0,0,0) or not primitive4(c): continue
      x=sum((ZZ(c[i])*mb[i] for i in range(4)),A(0))
      if abs(nr(x))!=N: continue
      OO=order_from_alpha(x)
      if OO is not None:
          O=OO; alpha_seed=x
          break
     if O is not None: break
    if O is not None: break
   if O is not None: break
  if O is not None: break

if O is None: raise RuntimeError("failed to reconstruct level-79 Eichler order")
ob=list(O.basis())
print(f"AL|stage=order|alpha79_seed={alpha_seed}|order={O}",flush=True)

# Coordinate matrices for membership and normalizer tests.
OB=matrix(QQ,[qcoords(b) for b in ob])
OBinv=OB.inverse()
MB=matrix(QQ,[qcoords(b) for b in mb])
MBinv=MB.inverse()

def coords_in_basis(x,Binv):
    return qcoords(x)*Binv

def integral_vec(v):
    return all(QQ(x).denominator()==1 for x in v)

def normalizer_matrix(alpha):
    # Matrix of conjugation x -> alpha^-1 x alpha on O basis.
    ai=alpha.inverse()
    rows=[]
    for b in ob:
        c=coords_in_basis(ai*b*alpha,OBinv)
        if not integral_vec(c): return None
        rows.append([ZZ(x) for x in c])
    C=matrix(ZZ,rows)
    if abs(C.det())!=1: return None
    return C

# Search elements in O first, then maximal order, with reduced norm +-p.
def find_normalizer(p,bound):
    seen=set()
    for source_name,basis in (("O",ob),("M",mb)):
      for B in range(1,bound+1):
        rng=range(-B,B+1)
        for c0 in rng:
         for c1 in rng:
          for c2 in rng:
           for c3 in rng:
            c=(c0,c1,c2,c3)
            if c==(0,0,0,0) or not primitive4(c): continue
            # Shell only, avoids repeating previous boxes.
            if max(abs(x) for x in c)!=B: continue
            x=sum((ZZ(c[i])*basis[i] for i in range(4)),A(0))
            nx=nr(x)
            if abs(nx)!=p: continue
            key=tuple(qcoords(x))
            if key in seen: continue
            seen.add(key)
            C=normalizer_matrix(x)
            if C is None: continue
            return source_name,c,x,nx,C
    return None

norms={}
for p in (2,3,79):
    print(f"AL|stage=search_normalizer|p={p}|bound={a.search_bound}",flush=True)
    r=find_normalizer(ZZ(p),a.search_bound)
    if r is None:
        print(f"AL|NORMALIZER_MISSING|p={p}",flush=True)
        continue
    source,c,x,nx,C=r
    norms[p]=(x,C)
    print(f"AL|normalizer|p={p}|source={source}|coeffs={c}|norm={nx}|alpha={x}|det_action={C.det()}",flush=True)

# Gross lattice O^T.
tau=[2*b-tr(b)*A.one() for b in ob]
V=[]
for x in tau:
    c=qcoords(x)
    assert c[0]==0
    V.append(vector(QQ,c[1:]))
den=lcm([QQ(x).denominator() for v in V for x in v])
Z=matrix(ZZ,[[ZZ(den*x) for x in v] for v in V])
Bz=Z.row_module().basis_matrix()
B=(Bz.change_ring(QQ))/den
Binv=B.inverse()
gens=A.gens(); i,j,k=gens[0],gens[1],gens[2]
def from3(v): return QQ(v[0])*i+QQ(v[1])*j+QQ(v[2])*k
gb=[from3(B.row(r)) for r in range(3)]

G=matrix(QQ,3,3)
for r in range(3): G[r,r]=nr(gb[r])
for r in range(3):
 for s in range(r+1,3):
    G[r,s]=G[s,r]=(nr(gb[r]+gb[s])-nr(gb[r])-nr(gb[s]))/2
print(f"AL|gross_gram={G}|det={G.det()}",flush=True)

def gross_coords(beta):
    c=qcoords(beta)
    if c[0]!=0: raise RuntimeError("not trace zero")
    z=vector(QQ,c[1:])*Binv
    if not integral_vec(z):
        raise RuntimeError(f"normalizer image escaped Gross lattice: {z}")
    return vector(ZZ,[ZZ(x) for x in z])

def qg(v):
    vv=vector(ZZ,v); return QQ(vv*G*vv)

def prim3(v): return gcd([abs(ZZ(x)) for x in v])==1

def seed_for_norm(n,bound):
    for x in range(-bound,bound+1):
     for y in range(-bound,bound+1):
      for z in range(-bound,bound+1):
       if x==y==z==0 or not prim3((x,y,z)): continue
       if qg((x,y,z))==n:
           return vector(ZZ,[x,y,z])
    return None

# Exact Gross-lattice action matrices.
gross_actions={}
for p,(alpha,C4) in norms.items():
    ai=alpha.inverse()
    rows=[]
    for b in gb:
        rows.append(list(gross_coords(ai*b*alpha)))
    W=matrix(ZZ,rows)
    # row coordinates: v -> v*W
    assert W*G*W.transpose()==G
    gross_actions[p]=W
    print(f"AL|gross_action|p={p}|det={W.det()}|W={W}",flush=True)

def canon_pm(v):
    t=tuple(map(int,v))
    nt=tuple(-x for x in t)
    return min(t,nt)

def inspect_actions(seed):
    """
    Raw conjugation on Gross vectors is NOT a finite Atkin-Lehner orbit:
    w_p is involutive only modulo O^x-conjugacy.  Do not BFS these matrices.
    Instead record first images and squares; embedding-class reduction comes next.
    """
    out = []

    for p,W in gross_actions.items():
        v1 = seed * W
        v2 = v1 * W

        rel1 = (
            "same" if tuple(v1) == tuple(seed)
            else "minus" if tuple(v1) == tuple(-seed)
            else "other"
        )

        rel2 = (
            "same" if tuple(v2) == tuple(seed)
            else "minus" if tuple(v2) == tuple(-seed)
            else "other"
        )

        out.append((p, v1, v2, rel1, rel2))

    return out

lines=[]
for n in targets:
    seed=seed_for_norm(n,a.hit_bound)
    if seed is None:
        print(f"AL|target={n}|NO_SEED",flush=True); continue
    beta=sum((ZZ(seed[r])*gb[r] for r in range(3)),A(0))
    print(f"AL|target={n}|seed={tuple(seed)}|beta={beta}|norm={nr(beta)}",flush=True)

    actions = inspect_actions(seed)

    for p,v1,v2,rel1,rel2 in actions:
        print(
            f"ALACTION|target={n}|p={p}"
            f"|relation={rel1}|image={tuple(v1)}"
            f"|square_relation={rel2}|square_image={tuple(v2)}",
            flush=True,
        )
        lines.append(
            f"ALACTION|target={n}|p={p}"
            f"|relation={rel1}|image={tuple(v1)}"
            f"|square_relation={rel2}|square_image={tuple(v2)}"
        )


Path(a.out).parent.mkdir(parents=True,exist_ok=True)
Path(a.out).write_text("\n".join(lines)+"\n")
print(f"AL|stage=done|out={a.out}",flush=True)
