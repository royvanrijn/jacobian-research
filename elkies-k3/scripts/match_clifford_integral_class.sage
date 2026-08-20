from sage.all import *
from pathlib import Path
import argparse, math

ap=argparse.ArgumentParser()
ap.add_argument("--candidates",default="artifacts/local/elkies-k3/cm-t2-candidates.txt")
ap.add_argument("--HG-file",default="artifacts/local/elkies-k3/clifford-to-k3-map.txt")
ap.add_argument("--bound",type=int,default=200)
ap.add_argument("--out",default="artifacts/local/elkies-k3/clifford-class-match.txt")
a=ap.parse_args()

def load_candidates(path):
    Ts={}
    for ln in Path(path).read_text().splitlines():
        if ln.startswith("TGRAM|"):
            _,sid,row=ln.split("|",2)
            Ts.setdefault(int(sid),[]).append([ZZ(x) for x in row.split()])
    return {k:matrix(ZZ,v) for k,v in Ts.items() if len(v)==3}

# Parse HG from the previous result file, whose first line starts HG=[...].
# Sage matrix repr spans lines, so parse between HG= and HT=.
txt=Path(a.HG_file).read_text()
chunk=txt.split("HG=",1)[1].split("\nHT=",1)[0].strip()
# Matrix textual repr is not eval-safe; parse integer rows from brackets.
rows=[]
for ln in chunk.splitlines():
    vals=[ZZ(x) for x in ln.replace("["," ").replace("]"," ").split()]
    if vals: rows.append(vals)
if len(rows)!=3:
    raise SystemExit(f"failed to parse HG from {a.HG_file}: got {rows}")
HG=matrix(ZZ,rows)
print(f"CLASSMATCH|HG={HG}|det={HG.det()}|bound={a.bound}",flush=True)

Ts=load_candidates(a.candidates)
print(f"CLASSMATCH|candidate_count={len(Ts)}",flush=True)

# Fast representations q(v)=target using quadratic formula in z.
def reps_of_norm(H,n,B):
    h00,h01,h02=map(ZZ,(H[0,0],H[0,1],H[0,2]))
    h11,h12,h22=map(ZZ,(H[1,1],H[1,2],H[2,2]))
    out=[]
    for x in range(-B,B+1):
      for y in range(-B,B+1):
        aa=h22
        bb=2*(h02*x+h12*y)
        cc=h00*x*x+2*h01*x*y+h11*y*y-n
        disc=bb*bb-4*aa*cc
        if disc<0: continue
        sd=isqrt(int(disc))
        if sd*sd!=disc: continue
        den=2*aa
        for num in (-bb+sd,-bb-sd):
            if den and num%den==0:
                z=ZZ(num//den)
                if abs(z)<=B:
                    v=vector(ZZ,[x,y,z])
                    if v!=0 and v*H*v==n:
                        out.append(v)
    # dedupe
    seen=set(); ans=[]
    for v in out:
        t=tuple(v)
        if t not in seen:
            seen.add(t); ans.append(v)
    return ans

def find_U(Hsrc,Hdst,B):
    # rows u_i in source coordinates must realize Hdst pairings.
    norms=[int(Hdst[i,i]) for i in range(3)]
    pools={}
    for n in set(norms):
        pools[n]=reps_of_norm(Hsrc,n,B)
    counts={n:len(pools[n]) for n in pools}
    rows=[None]*3
    def bt(i):
        if i==3:
            U=matrix(ZZ,rows)
            if abs(U.det())==1 and U*Hsrc*U.transpose()==Hdst:
                return U
            return None
        for v in pools[norms[i]]:
            ok=True
            for j in range(i):
                if v*Hsrc*rows[j] != Hdst[i,j]:
                    ok=False; break
            if not ok: continue
            rows[i]=v
            r=bt(i+1)
            if r is not None: return r
        return None
    return bt(0),counts

lines=[]
for tid in sorted(Ts):
    H=Ts[tid]
    if H.det()!=HG.det():
        continue
    U,counts=find_U(HG,H,a.bound)
    print(f"CLASSMATCH|T={tid}|counts={counts}|match={U is not None}",flush=True)
    if U is not None:
        print(f"CLASSMATCH|JACKPOT|T={tid}|detU={U.det()}|U={U}",flush=True)
        lines.append(f"JACKPOT T={tid} detU={U.det()}\nU={U}\n")
        # Don't stop: see whether several stored candidates are same class.
Path(a.out).parent.mkdir(parents=True,exist_ok=True)
Path(a.out).write_text("\n".join(lines))
print(f"CLASSMATCH|stage=done|jackpots={len(lines)}|out={a.out}",flush=True)
