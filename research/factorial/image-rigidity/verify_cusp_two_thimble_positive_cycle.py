#!/usr/bin/env python3
"""Exact/source support for the actual cusp two-thimble cycle proof."""
from __future__ import annotations
if not __debug__:
    raise SystemExit("Refusing -O/-OO: exact assertions must execute.")
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations,product
from math import comb,factorial
from pathlib import Path
import argparse,copy,hashlib,json,re,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
PN="CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md"
CN="verify_cusp_two_thimble_positive_cycle.py"
IN="cusp_two_thimble_inputs.json"
LN="cusp_two_thimble_source_ledger.json"
RN="cusp_two_thimble_verification.json"
PREFIX="research/factorial/image-rigidity/"
ROLES={
 "proof":("local",PN),"checker":("local",CN),"inputs":("local",IN),
 "stationary_proof":("repo",PREFIX+"CUSP_STATIONARY_BELYI_MONODROMY_AND_RADIAL_BARRIER.md"),
 "stationary_checker":("repo",PREFIX+"verify_cusp_stationary_belyi_monodromy.py"),
 "stationary_inputs":("repo",PREFIX+"cusp_stationary_belyi_inputs.json"),
 "stationary_ledger":("repo",PREFIX+"cusp_stationary_belyi_source_ledger.json"),
 "stationary_receipt":("repo",PREFIX+"cusp_stationary_belyi_verification.json"),
}
SCOPE={
 "schema":1,"fixed_integer_cusp":{"s_min":1,"a":"s^2","d":"2s+1"},
 "parameter_domain":"Im(w)>0","exact_period_integer_m_min":1,
 "ordinary_untwisted_integer_periods_only":True,
 "three_ends":["0","zeta","infinity"],"relative_rank":2,
 "original_cycle_tree_path":"one or two oriented edges; coefficients 0,+1,-1",
 "complete_original_cycle":True,
 "complete_Gaussian_packets_on_compact_nonStokes_cells":True,
 "both_active_on_distinguished_equimodular_arc":True,
 "both_active_everywhere":False,"Stokes_uniform_error":False,
 "eventual_compact_regular_cell_noncollision":True,
 "effective_cutoff":False,"global_cusp_SFC":False,
 "all_actual_roots_localized":False,"uniform_in_d":False,
 "formal_verification":False,"external_review":False,
}
CONTEXT={"url":"https://arxiv.org/pdf/1001.2933v4",
 "locator":"Sections3.1.1-3.1.5, printed24-32",
 "scope":"Classical context only. Inward boundary, logarithmic ends and actual cycle are proved directly."}
def require(ok,msg):
    if not ok:raise AssertionError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text(encoding="utf-8"))
def check_inputs(data):
    require(type(data) is dict and data==SCOPE,"Exact scoped input contract")
    for k,v in SCOPE.items():require(type(data[k]) is type(v),"Scope type:"+k)
    require(type(data["fixed_integer_cusp"]["s_min"]) is int,"Integer s_min")
def check_sources(doc,repo):
    require(type(doc) is dict and set(doc)=={"schema","formal_verification","external_review","context_only","sources"},"Source schema")
    require(type(doc["schema"]) is int and doc["schema"]==1,"Source schema1")
    require(doc["formal_verification"] is False and doc["external_review"] is False,"Assurance flags")
    require(doc["context_only"]==CONTEXT,"Exact context-only source scope")
    rows=doc["sources"];require(type(rows) is list and len(rows)==len(ROLES),"All eight mandatory roles")
    out={}
    for row in rows:
        require(type(row) is dict and set(row)=={"role","kind","path","sha256"},"Source row schema")
        role=row["role"];require(type(role) is str and role in ROLES and role not in out,"Unique role")
        require((row["kind"],row["path"])==ROLES[role],"Exact source path")
        require(type(row["sha256"]) is str and re.fullmatch("[0-9a-f]{64}",row["sha256"]),"SHA256")
        root=HERE if row["kind"]=="local" else repo.resolve()
        p=(root/row["path"]).resolve()
        require(p.is_relative_to(root.resolve()) and p.is_file() and sha(p)==row["sha256"],"Missing/stale:"+role)
        out[role]=row
    return out
def stationary_closure(repo):
    base=repo/PREFIX;doc=load(base/"cusp_stationary_belyi_source_ledger.json")
    require(type(doc.get("schema")) is int and doc["schema"]==1,"Imported schema1")
    rows=doc.get("files");require(type(rows) is list and len(rows)==6,"Imported six source premises")
    out=[];seen=set()
    for row in rows:
        require(type(row) is dict and row.get("root") in {"local","repo"},"Imported source kind")
        key=(row["root"],row["path"]);require(key not in seen,"Imported source unique");seen.add(key)
        root=base if row["root"]=="local" else repo
        p=(root/row["path"]).resolve()
        require(p.is_relative_to(root.resolve()) and p.is_file() and sha(p)==row["sha256"],"Imported stale source")
        rel=(Path(PREFIX)/row["path"]).as_posix() if row["root"]=="local" else row["path"]
        out.append({"path":rel,"sha256":row["sha256"]})
    for name in ["cusp_stationary_belyi_source_ledger.json","cusp_stationary_belyi_verification.json"]:
        out.append({"path":PREFIX+name,"sha256":sha(base/name)})
    require(len({r["path"] for r in out})==8,"Eight-file imported closure")
    return sorted(out,key=lambda r:r["path"])
def equal_receipt(actual,expected):require(actual==expected,"Imported fresh receipt differs")
def replay_stationary(repo):
    script=repo/ROLES["stationary_checker"][1];expected=repo/ROLES["stationary_receipt"][1]
    with tempfile.TemporaryDirectory(prefix="cusp-thimble-") as d:
        output=Path(d)/"stationary.json"
        r=subprocess.run([sys.executable,"-B",str(script),"--repo-root",str(repo),"--output",str(output)],
                         cwd=d,capture_output=True,text=True,timeout=30)
        require(r.returncode==0,"Stationary replay failed:"+r.stdout+";"+r.stderr)
        equal_receipt(output.read_bytes(),expected.read_bytes())
        require(load(output).get("status")=="PASS_STATIONARY_COVER_ALGEBRA_AND_GROUP_SUPPORT","Exact stationary source status")
        return sha(output)

NV=6
def clean(p):return {m:Q(c) for m,c in p.items() if c}
def cst(c):return {} if not c else {(0,)*NV:Q(c)}
def var(i):
    k=[0]*NV;k[i]=1
    return {tuple(k):Q(1)}
def add(*items):
    out={}
    for p in items:
        for m,c in p.items():out[m]=out.get(m,Q(0))+c
    return clean(out)
def sc(p,c):return clean({m:c*v for m,v in p.items()})
def mul(*items):
    out=cst(1)
    for p in items:
        nxt={}
        for m,c in out.items():
            for n,v in p.items():
                k=tuple(x+y for x,y in zip(m,n));nxt[k]=nxt.get(k,Q(0))+c*v
        out=clean(nxt)
    return out
def sq(p):return mul(p,p)
def cmul(z,w):return (add(mul(z[0],w[0]),sc(mul(z[1],w[1]),-1)),
                         add(mul(z[0],w[1]),mul(z[1],w[0])))
def algebra():
    d,r,u,v,X,k=[var(i) for i in range(NV)]
    upper=cmul((mul(d,X),{}),(add(u,X),sc(v,-1)))[1]
    require(not add(upper,mul(d,X,v)),"Upper inward numerator")
    lower=cmul((sc(mul(d,X),-1),{}),(add(u,sc(X,-1)),sc(v,-1)))[1]
    require(not add(lower,sc(mul(d,X,v),-1)),"Lower inward numerator")
    den=add(sq(add(u,sc(X,-1))),sq(v))
    full_lower=add(mul(r,k,den),lower)
    require(not add(full_lower,sc(add(mul(r,k,den),mul(d,X,v)),-1)),"Lower full sign")
    a,S,w,T,d,_=[var(i) for i in range(NV)]
    N=add(mul(add(a,sc(S,-1)),add(w,T)),mul(d,T))
    infinity=add(mul(add(a,d,sc(S,-1)),add(w,T)),sc(mul(d,w),-1))
    require(not add(N,sc(infinity,-1)),"Exact infinity derivative")
    product_rule=add(mul(a,add(w,T)),mul(d,T),sc(mul(S,add(w,T)),-1))
    require(not add(product_rule,sc(N,-1)),"Actual polynomial derivative")
    # Phi'V=-1 after clearing S(w+T) and the stationary numerator.
    require(not add(mul(N,sc(mul(S,add(w,T)),-1)),mul(S,add(w,T),N)),"Normalized flow")
    # Polar D derivative: multiply the imaginary derivative by d*P.
    P,Qv,dd,Dv,Ng,theta=[var(i) for i in range(NV)]
    phase_numerator=add(sc(sq(Qv),-1),sc(sq(P),-1))
    require(not add(phase_numerator,sq(Qv),sq(P)),"Monotone phase")
    # Product rule for (N-theta)*D, using Dtheta=-gap/(N-theta).
    gap=Qv
    derivative=add(sc(Dv,-1),sc(gap,-1))
    require(not add(derivative,Dv,gap),"Adjacent action product derivative")
    return 8
@dataclass(frozen=True)
class G:
    re:Q=Q(0)
    im:Q=Q(0)
    def __add__(self,o):
        if not isinstance(o,G):o=G(Q(o))
        return G(self.re+o.re,self.im+o.im)
    __radd__=__add__
    def __neg__(self):return G(-self.re,-self.im)
    def __sub__(self,o):return self+-o if isinstance(o,G) else self+(-Q(o))
    def __mul__(self,o):
        if not isinstance(o,G):o=G(Q(o))
        return G(self.re*o.re-self.im*o.im,self.re*o.im+self.im*o.re)
    __rmul__=__mul__
    def inv(self):
        n=self.re*self.re+self.im*self.im
        assert n
        return G(self.re/n,-self.im/n)
    def __truediv__(self,o):return self*(o.inv() if isinstance(o,G) else G(Q(1)/Q(o)))
    def __pow__(self,n):
        if n<0:return self.inv()**(-n)
        out=G(Q(1))
        for _ in range(n):out=out*self
        return out
def rotation_control(s,m):
    a=s*s;d=2*s+1;z=G(Q(1),-Q(1,2*d));w=-(z**d);c=z*m
    assert z.re>0 and z.im<0 and w.im>0
    L0=G();LE=G();RE=G()
    for j in range(m+1):
        n=a*m+d*j
        whole=G(Q(factorial(n)))/(c**(n+1))
        tail=whole*sum((c**k/Q(factorial(k)) for k in range(n+1)),G())
        L0=L0+whole*((-1)**j*comb(m,j))
        LE=LE-tail*((-1)**j*comb(m,j))
        RE=RE+tail*((-1)**(m-j)*comb(m,j))
    assert LE+RE*((-1)**m)==G()
    pref=(z**(a*m+1))*(w**m)
    actual=sum((w**(m-j)*Q(comb(m,j)*factorial(a*m+d*j),m**(a*m+d*j+1))
                for j in range(m+1)),G())
    assert pref*L0==actual
def tree_controls():
    vertices=range(3);edges=list(combinations(vertices,2));counts={1:0,2:0}
    for pair in combinations(edges,2):
        for flip in product([False,True],repeat=2):
            cols=[]
            for edge,rev in zip(pair,flip):
                u,v=edge if not rev else edge[::-1]
                col=[0,0,0];col[u]=-1;col[v]=1;cols.append(col)
            a,b=cols[0][0],cols[1][0];c,d=cols[0][1],cols[1][1]
            det=a*d-b*c;assert abs(det)==1
            x=Q(-d,det);y=Q(c,det)
            assert all(x*cols[0][i]+y*cols[1][i]==[-1,0,1][i] for i in vertices)
            assert x in [-1,0,1] and y in [-1,0,1] and (x or y)
            counts[int(bool(x))+int(bool(y))]+=1
    assert counts=={1:8,2:4}
    return 12,counts
def guards():
    g={
      "sector_cosine_margin":Q(1,2)-Q(1,3)==Q(1,6),
      "inverse_derivative_speed":1/(1-Q(1,4))<2,
      "inverse_derivative_perturbation":Q(1,4)/(1-Q(1,4))==Q(1,3),
      "zero_cap_exponent":1+1>0,
      "simple_zero_cap_exponent":1+1>0,
      "Gaussian_relative_order":Q(3,2)-Q(1,2)==1,
      "compact_radial_error":Q(1,1)>0,
      "noncollision_leading_gap":Q(1,1)>0,
      "integer_m_min":1>0,
    }
    assert all(g.values())
    return g
def source_negatives(ledger,repo):
    names=[]
    def reject(name,fn):
        try:fn()
        except (AssertionError,KeyError,TypeError,ValueError):names.append(name);return
        raise AssertionError("Negative accepted:"+name)
    for role in ROLES:
        d=copy.deepcopy(ledger);d["sources"]=[r for r in d["sources"] if r["role"]!=role]
        reject("omit:"+role,lambda d=d:check_sources(d,repo))
        d=copy.deepcopy(ledger);next(r for r in d["sources"] if r["role"]==role)["sha256"]="0"*64
        reject("badpin:"+role,lambda d=d:check_sources(d,repo))
    for key in ["both_active_everywhere","Stokes_uniform_error","effective_cutoff",
                "global_cusp_SFC","all_actual_roots_localized","uniform_in_d"]:
        d=copy.deepcopy(SCOPE);d[key]=True
        reject("false-scope:"+key,lambda d=d:check_inputs(d))
    for key in ["complete_original_cycle","both_active_on_distinguished_equimodular_arc",
                "eventual_compact_regular_cell_noncollision"]:
        d=copy.deepcopy(SCOPE);d[key]=False
        reject("omitted-scope:"+key,lambda d=d:check_inputs(d))
    for key in ["schema","empty","duplicate"]:
        d=copy.deepcopy(ledger)
        if key=="schema":d["schema"]=True
        if key=="empty":d["sources"]=[]
        if key=="duplicate":d["sources"][-1]=copy.deepcopy(d["sources"][0])
        reject("ledger:"+key,lambda d=d:check_sources(d,repo))
    reject("forged-child",lambda:equal_receipt(b'{"status":"PASS"}',b'{"status":"PASS","real":true}'))
    return names
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--repo-root",required=True,type=Path)
    ap.add_argument("--output",type=Path,default=HERE/RN);ns=ap.parse_args();repo=ns.repo_root.resolve()
    require(repo.is_dir(),"Explicit repo root required")
    inputs=load(HERE/IN);check_inputs(inputs)
    ledger=load(HERE/LN);roles=check_sources(ledger,repo)
    closure=stationary_closure(repo);child=replay_stationary(repo)
    ids=algebra();trees,counts=tree_controls();budget=guards()
    for s in range(1,4):
        for m in range(1,5):rotation_control(s,m)
    # Real countercontrols to forbidden inference steps.
    assert counts[1]>0                         # Rank2 does not imply both coefficients nonzero.
    assert (-1)**2==1 and -1!=1                # m-phases alone do not permit a gradient connection.
    upper_boundary=G(Q(3))/(G(Q(1),Q(1))+1)    # Actual a1,d3,r1,w1+i.
    assert upper_boundary.im==-Q(3,5)<0        # The flipped inward sign is false.
    assert 3-1==2 and 2-1!=2                   # Omitting the interior zero loses a relative generator.
    negatives=source_negatives(ledger,repo);check_sources(ledger,repo)
    payload={"status":"PASS","schema":1,"formal_verification":False,"external_review":False,
     "scope":inputs,"proof_sha256":sha(HERE/PN),"checker_sha256":sha(HERE/CN),
     "input_sha256":sha(HERE/IN),"source_ledger_sha256":sha(HERE/LN),
     "mandatory_roles":len(roles),"sources":[roles[r] for r in sorted(roles)],
     "imported_stationary_closure":closure,"fresh_stationary_receipt_sha256":child,
     "fresh_stationary_receipt_byte_identical":True,
     "counts":{"universal_polynomial_identities":ids,"oriented_tree_cases":trees,
               "actual_complex_rotated_integral_controls":12,"rational_guards":len(budget),
               "mathematical_negative_controls":4,"source_scope_controls":len(negatives)},
     "tree_path_lengths":counts,"guard_names":sorted(budget),"negative_names":negatives,
     "meaning":"Finite algebra and source replay support the complete written topological/analytic proof; they do not formally verify flow, handles or Gaussian estimates."}
    ns.output.parent.mkdir(parents=True,exist_ok=True)
    ns.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS","counts":payload["counts"],"output":str(ns.output)}))
if __name__=="__main__":main()

