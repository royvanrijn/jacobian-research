from sage.all import *
from pathlib import Path
import argparse, math, collections, ast

ap = argparse.ArgumentParser(description="Enumerate primitive negative-vector / Heegner candidates in the recovered generic transcendental lattice.")
ap.add_argument("--candidates", default="artifacts/local/elkies-k3/cm-t2-candidates.txt")
ap.add_argument("--T", type=int, default=0, help="ternary representative id; all current survivors are in the same genus/class")
ap.add_argument("--bound", type=int, default=60, help="coordinate box for primitive v")
ap.add_argument("--norm-max", type=int, default=5000, help="only keep -v^2 <= this")
ap.add_argument("--top", type=int, default=100)
ap.add_argument("--out", default="artifacts/local/elkies-k3/heegner-orbits.txt")
args = ap.parse_args()

def load_T(path, tid):
    rows=[]
    for ln in Path(path).read_text().splitlines():
        if ln.startswith(f"TGRAM|{tid}|"):
            rows.append([ZZ(x) for x in ln.split("|",2)[2].split()])
    if len(rows)!=3:
        raise SystemExit(f"could not find 3 rows for T={tid}")
    return matrix(ZZ,rows)

def primitive(v):
    return gcd([abs(ZZ(x)) for x in v]) == 1

def divisibility(H,v):
    w=H*v
    return gcd([abs(ZZ(x)) for x in w])

def complement_gram(H,v):
    # Integer solutions x with <x,v>=0.
    row=matrix(ZZ,1,3,list(v * H))
    K=row.right_kernel_matrix()
    if K.nrows()!=2:
        return None
    G=K*H*K.transpose()
    if G.det() <= 0:
        return None
    # Ensure positive definite (rank 2).
    if G[0,0] <= 0:
        K[0] = -K[0]
        G=K*H*K.transpose()
    if G[0,0] <= 0 or G.det() <= 0:
        return None
    return K,G

def bqf_data(G):
    # Even binary lattice G = [[2a,b],[b,2c]].
    assert G[0,0] % 2 == 0 and G[1,1] % 2 == 0
    a=ZZ(G[0,0]//2); b=ZZ(G[0,1]); c=ZZ(G[1,1]//2)
    content=gcd([abs(a),abs(b),abs(c)])
    ap,bp,cp=a//content,b//content,c//content
    D=bp*bp-4*ap*cp  # primitive quadratic-order discriminant
    try:
        red=BinaryQF([ap,bp,cp]).reduced_form()
        key=tuple(map(int,red))
    except Exception:
        key=(int(ap),int(bp),int(cp))
    return dict(a=a,b=b,c=c,content=content,primitive=(ap,bp,cp),
                order_disc=ZZ(D),reduced=key)

H=load_T(args.candidates,args.T)
assert H.det()==-948
print(f"HEEGNER|stage=start|T={args.T}|det={H.det()}|bound={args.bound}|norm_max={args.norm_max}", flush=True)

# Deduplicate by robust orbit surrogate:
# (norm, divisibility, primitive reduced binary complement).
# For a rank-(2,1) lattice this is much sharper than complement discriminant alone.
groups={}
count=0
B=args.bound
for x in range(-B,B+1):
  for y in range(-B,B+1):
    for z in range(-B,B+1):
      v=vector(ZZ,[x,y,z])
      if not primitive(v):
        continue
      n=ZZ(v*H*v)
      if n>=0 or -n>args.norm_max:
        continue
      d=ZZ(divisibility(H,v))
      cg=complement_gram(H,v)
      if cg is None:
        continue
      K,G=cg
      # determinant identity is a useful exact consistency check.
      lhs=QQ(G.det())
      rhs = QQ(abs(H.det()) * (-n)) / QQ(d*d)
      if lhs != rhs:
        raise RuntimeError(f"det identity failed: {lhs} != {rhs}, v={v}, d={d}")
      bd=bqf_data(G)
      key=(int(n),int(d),bd["reduced"],int(bd["content"]),int(bd["order_disc"]))
      score=(abs(int(n)),int(d),int(G.det()),max(abs(int(q)) for q in v))
      old=groups.get(key)
      rec=dict(v=tuple(map(int,v)),norm=int(n),div=int(d),G=G,K=K,bd=bd,score=score)
      if old is None or score < old["score"]:
          groups[key]=rec
      count += 1

rows=sorted(groups.values(), key=lambda r:(r["score"][0], r["G"].det(), -r["div"], r["bd"]["reduced"]))
print(f"HEEGNER|stage=enumerated|vectors={count}|orbit_surrogates={len(rows)}",flush=True)

lines=[]
for i,r in enumerate(rows[:args.top],1):
    G=r["G"]; bd=r["bd"]
    # field fundamental discriminant derived from the primitive form discriminant
    D=ZZ(bd["order_disc"])
    try:
        Kfld=QuadraticField(D)
        Dfund=ZZ(Kfld.discriminant())
        conductor=ZZ(sqrt(ZZ(D//Dfund))) if Dfund != 0 and D % Dfund == 0 and D//Dfund >= 0 and is_square(D//Dfund) else None
    except Exception:
        Dfund=None; conductor=None
    s=(f"HEEGNER|rank={i}|v={r['v']}|norm={r['norm']}|div={r['div']}"
       f"|comp_det={G.det()}|binary={bd['reduced']}|content={bd['content']}"
       f"|order_disc={D}|field_disc={Dfund}|conductor={conductor}")
    print(s,flush=True); lines.append(s)
    lines.append("HEEGNERGRAM|rank=%d|%s"%(i,";".join(",".join(map(str,row)) for row in G.rows())))

Path(args.out).parent.mkdir(parents=True,exist_ok=True)
Path(args.out).write_text("\n".join(lines)+"\n")
print(f"HEEGNER|stage=done|out={args.out}",flush=True)
print("HEEGNER|note=These are O(T)-orbit SURROGATES, not a proof of orbit equivalence. Same (norm,div,reduced complement) is a very strong invariant; exact O(T) orbit certification can follow for a short list.",flush=True)
