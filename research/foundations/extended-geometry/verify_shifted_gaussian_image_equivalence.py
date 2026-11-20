#!/usr/bin/env python3
"""Exact finite controls for the shifted Gaussian Image equivalence.

Only the Python standard library is used. Sparse substitution and Gaussian
Wick evaluation are separate from original falling-factorial contraction.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import factorial
from pathlib import Path
import json
import sys

if not __debug__:
    raise SystemExit("Exact certificate replay requires assertions; refuse -O/-OO.")

HERE = Path(__file__).resolve().parent
PROOF = "SHIFTED_GAUSSIAN_IMAGE_EQUIVALENCE.md"
INPUTS = "shifted_gaussian_image_equivalence_inputs.json"
RECEIPT = "shifted_gaussian_image_equivalence_verification.json"

# Exact Q(i), deliberately represented independently of the original contraction.
class C:
    __slots__ = ("r", "i")
    def __init__(self, r=0, i=0):
        if isinstance(r, C):
            self.r, self.i = r.r, r.i
        else:
            self.r, self.i = F(r), F(i)
    def __add__(self, b):
        b=C(b); return C(self.r+b.r, self.i+b.i)
    __radd__=__add__
    def __neg__(self): return C(-self.r,-self.i)
    def __sub__(self,b): return self+-C(b)
    def __rsub__(self,b): return C(b)+-self
    def __mul__(self,b):
        b=C(b); return C(self.r*b.r-self.i*b.i,self.r*b.i+self.i*b.r)
    __rmul__=__mul__
    def __eq__(self,b):
        b=C(b); return self.r==b.r and self.i==b.i
    def __bool__(self): return bool(self.r or self.i)
    def __repr__(self): return f"C({self.r!s},{self.i!s})"

ZERO=C(0)
ONE=C(1)
I=C(0,1)


def clean(p): return {e:C(c) for e,c in p.items() if c}
def add(*ps):
    d=defaultdict(C)
    for p in ps:
        for e,c in p.items(): d[e]=d[e]+c
    return clean(d)
def scale(p,c): return clean({e:v*c for e,v in p.items()})
def mul(p,q):
    d=defaultdict(C)
    for a,c in p.items():
        for b,v in q.items():
            e=tuple(x+y for x,y in zip(a,b)); d[e]=d[e]+c*v
    return clean(d)
def one(n): return {(0,)*n:ONE}
def var(n,j):
    e=[0]*n; e[j]=1; return {tuple(e):ONE}
def pw(p,k,n):
    r=one(n)
    while k:
        if k&1: r=mul(r,p)
        k//=2
        if k: p=mul(p,p)
    return r

def substitute(p,vs,nout):
    maxs=[max((e[j] for e in p),default=0) for j in range(len(vs))]
    powers=[]
    for v,d in zip(vs,maxs):
        vv=[one(nout)]
        for _ in range(d): vv.append(mul(vv[-1],v))
        powers.append(vv)
    result={}
    for e,c in p.items():
        term={(0,)*nout:c}
        for j,k in enumerate(e):
            if k: term=mul(term,powers[j][k])
        result=add(result,term)
    return result

def image(p,n):
    """Original definition E(xi^alpha s^beta)=partial^alpha s^beta."""
    d=defaultdict(C)
    for e,c in p.items():
        a,b=e[:n],e[n:]
        if any(x>y for x,y in zip(a,b)): continue
        factor=1
        for x,y in zip(a,b): factor*=factorial(y)//factorial(y-x)
        k=tuple(y-x for x,y in zip(a,b)); d[k]=d[k]+c*factor
    return clean(d)

def gaussian(p,n):
    """Independent standard real Gaussian Wick moments, leaving z parameters."""
    d=defaultdict(C)
    for e,c in p.items():
        if any(k%2 for k in e[:2*n]): continue
        factor=1
        for k in e[:2*n]:
            for odd in range(1,k,2): factor*=odd
        k=e[2*n:]; d[k]=d[k]+c*factor
    return clean(d)

def shifted_vs(n,mode="correct"):
    nv=3*n
    xs=[var(nv,j) for j in range(n)]
    ys=[var(nv,n+j) for j in range(n)]
    zs=[var(nv,2*n+j) for j in range(n)]
    vs=[scale(add(xs[j],scale(ys[j],-I)),F(1,2)) for j in range(n)]
    if mode=="wrong_half": vs=[scale(v,2) for v in vs]
    ws=[add(xs[j],scale(ys[j],I if mode!="wrong_sign" else -I)) for j in range(n)]
    ss=[add(ws[j],zs[j]) if mode!="omitted_shift" else ws[j] for j in range(n)]
    return vs+ss

def inverse_vs(n):
    # Output variables (xi,s,z), input variables (X,Y,z).
    nv=3*n
    xis=[var(nv,j) for j in range(n)]
    ss=[var(nv,n+j) for j in range(n)]
    zs=[var(nv,2*n+j) for j in range(n)]
    ws=[add(ss[j],scale(zs[j],-1)) for j in range(n)]
    xs=[add(xis[j],scale(ws[j],F(1,2))) for j in range(n)]
    ys=[scale(add(ws[j],scale(xis[j],-2)),-I*F(1,2)) for j in range(n)]
    return xs+ys+zs

def load_poly(terms,n):
    p={}
    for e,c in terms:
        assert len(e)==2*n and all(isinstance(k,int) and k>=0 for k in e)
        e=tuple(e); p=add(p,{e:C(F(c))})
    return p

def rational_payload(p):
    result=[]
    for e,c in sorted(p.items()):
        assert c.i==0, (e,c)
        result.append([list(e),str(c.r)])
    return result

def expect_equal(p,q,label):
    assert p==q, (label,rational_payload(p),rational_payload(q))

def sha(path): return sha256(path.read_bytes()).hexdigest()

def find_root(arg):
    if arg:
        root=Path(arg).resolve()
        assert (root/"research/MATH_STATUS.json").is_file()
        return root
    for candidate in [HERE,*HERE.parents]:
        if (candidate/"research/MATH_STATUS.json").is_file(): return candidate
    raise SystemExit("Use --repo-root when replaying this scratch copy outside the repository.")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path)
    ap.add_argument("--write-receipt",action="store_true")
    ap.add_argument("--repo-root",type=Path)
    args=ap.parse_args()
    root=find_root(args.repo_root)
    ledgerpath=HERE/"shifted_gaussian_primary_source_ledger.json"
    sourceledger=json.loads(ledgerpath.read_text())
    validated=[]
    for entry in sourceledger["files"]:
        base=HERE if entry["location"]=="packet" else root
        source=base/entry["path"]
        assert sha(source)==entry["sha256"],"Source pin mismatch: "+entry["path"]
        validated.append({"path":entry["path"],"sha256":entry["sha256"],"role":entry["role"]})
    inputs=json.loads((HERE/INPUTS).read_text())
    assert inputs["schema"]==1
    monomial_controls=0
    for spec in inputs["monomial_grids"]:
        n,bound=spec["n"],spec["max_exponent"]
        vs=shifted_vs(n)
        for exps in product(range(bound+1),repeat=2*n):
            h={exps:ONE}
            expect_equal(image(h,n),gaussian(substitute(h,vs,3*n),n),f"monomial-{n}-{exps}")
            monomial_controls+=1
    mixed_controls=0
    mixed_results=[]
    for spec in inputs["phase_marker_fixtures"]:
        n=spec["n"]; f=load_poly(spec["phase"],n); q=load_poly(spec["marker"],n)
        sf=substitute(f,shifted_vs(n),3*n)
        sq=substitute(q,shifted_vs(n),3*n)
        fo=one(2*n); fg=one(3*n)
        for m in range(spec["max_m"]+1):
            orig=image(mul(q,fo),n)
            shifted=gaussian(mul(sq,fg),n)
            expect_equal(orig,shifted,f"mixed-{spec['name']}-{m}")
            mixed_results.append({"name":spec["name"],"m":m,"output":rational_payload(orig)})
            mixed_controls+=1
            fo=mul(fo,f); fg=mul(fg,sf)
    inverse_controls=0
    for n in inputs["inverse_dimensions"]:
        back=inverse_vs(n)
        # Compose S with its inverse on every coordinate; parameters are fixed.
        sv=shifted_vs(n)+[var(3*n,2*n+j) for j in range(n)]
        for j,expr in enumerate(back):
            expect_equal(substitute(expr,sv,3*n),var(3*n,j),f"inverse-X-Y-z-{n}-{j}")
            inverse_controls+=1
        for j,expr in enumerate(sv):
            expect_equal(substitute(expr,back,3*n),var(3*n,j),f"inverse-xi-s-z-{n}-{j}")
            inverse_controls+=1
        # General Gaussian markers really pull back to original polynomial markers.
        for j in range(2*n):
            a=add(pw(var(3*n,j),3,3*n),scale(var(3*n,2*n),I))
            expect_equal(substitute(substitute(a,back,3*n),sv,3*n),a,f"inverse-marker-{n}-{j}")
            inverse_controls+=1
    finite_marker_degree_controls=0
    for spec in inputs["finite_marker_degree_grids"]:
        n,B=spec["n"],spec["B"]
        back=inverse_vs(n)
        sv=shifted_vs(n)+[var(3*n,2*n+j) for j in range(n)]
        for exps in product(range(B+1),repeat=2*n):
            if sum(exps)>B: continue
            a={tuple(exps)+(0,)*n:ONE}
            original=substitute(a,back,3*n)
            assert all(sum(e[:2*n])<=B for e in original)
            expect_equal(substitute(original,sv,3*n),a,f"finite-marker-degree-{n}-{B}-{exps}")
            finite_marker_degree_controls+=1
    # Controls that the affine shift and covariance normalization are essential.
    s=var(2,1)
    xis=mul(var(2,0),s)
    corrupt=[]
    for mode,h in [("omitted_shift",s),("wrong_half",xis),("wrong_sign",xis)]:
        assert image(h,1)!=gaussian(substitute(h,shifted_vs(1,mode),3),1)
        corrupt.append(mode)
    # A concrete known pure unsafe phase: literal original bilinears, not a CT proxy.
    n=2; nv=4
    xi1,xi2,s1,s2=[var(nv,j) for j in range(nv)]
    R=add(mul(xi1,s1),mul(xi2,s2)); Z=mul(xi1,s2)
    W=scale(mul(xi2,s1),2); T=add(mul(xi1,s1),scale(mul(xi2,s2),-1))
    f=mul(add(R,Z),add(mul(pw(R,2,nv),W),scale(mul(add(scale(R,2),Z),pw(T,2,nv)),F(-1,2))))
    assert len(f)==16
    known=[]; fo=one(nv)
    for m in range(1,inputs["known_unsafe_max_m"]+1):
        fo=mul(fo,f)
        assert image(fo,n)=={}
        marked=image(mul(Z,fo),n)
        denom=1
        for j in range(1,2*m+2,2): denom*=j
        value=F(factorial(4*m+2)*factorial(m),denom)
        expect_equal(marked,{(0,0):C(value)},f"known-unsafe-original-{m}")
        known.append({"m":m,"pure":[],"marked":str(value),"original_support":len(fo)})
    # Actual shifted-Gaussian expansion for both q_F and its bad marker, m=1.
    sf=substitute(f,shifted_vs(n),3*n)
    sz=substitute(Z,shifted_vs(n),3*n)
    expect_equal(gaussian(sf,n),{},"known-pure-shift")
    expect_equal(gaussian(mul(sz,sf),n),{(0,0):C(240)},"known-bad-marker-shift")
    # Simple contrasting phases; finite controls only, with proof-level classes explicit.
    assert image({(0,1):ONE},1)=={(1,):ONE}  # an unbalanced output is not discarded
    assert image({(1,0):ONE},1)=={}
    assert image({(0,0):ONE},1)=={(0,):ONE}
    pin_files={
        "proof":HERE/PROOF,"checker":Path(__file__).resolve(),"inputs":HERE/INPUTS,
        **{x:root/x for x in inputs["retained_proof_inputs"]},
    }
    receipt={
        "schema":1,"result":"PASS","assurance":"finite exact internal controls; unrestricted statements use written proofs",
        "source_hashes":{k:sha(v) for k,v in pin_files.items()},
        "primary_source_ledger_sha256":sha(ledgerpath),"validated_source_pins":validated,
        "monomial_controls":monomial_controls,"mixed_phase_marker_controls":mixed_controls,
        "affine_inverse_controls":inverse_controls,"finite_marker_degree_controls":finite_marker_degree_controls,"negative_controls":corrupt,
        "known_pure_unsafe_original_powers":known,
        "known_pure_unsafe_shifted_controls":2,"known_phase_shifted_support":len(sf),
        "mixed_outputs":mixed_results,
    }
    dest=args.output or HERE/RECEIPT
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:receipt[k] for k in ("result","monomial_controls","mixed_phase_marker_controls","affine_inverse_controls","finite_marker_degree_controls","negative_controls","known_pure_unsafe_shifted_controls")}))
    print(str(dest))

if __name__=="__main__": main()
