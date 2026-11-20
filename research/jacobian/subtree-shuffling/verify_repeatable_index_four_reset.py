#!/usr/bin/env python3
"""Exact symbolic and literal checks for a repeatable cubic index-four reset."""
if not __debug__:
    raise SystemExit("Refusing optimized Python: exact assertion checks are required.")
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse, copy, hashlib, json
HERE=Path(__file__).resolve().parent

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def clean(p):return {m:c for m,c in p.items() if c}
def add(*ps):
    r={}
    for p in ps:
        for m,c in p.items():r[m]=r.get(m,F(0))+c
    return clean(r)
def scale(p,c):return clean({m:v*c for m,v in p.items()})
def term(*ids,c=1):return {tuple(sorted(ids)):F(c)}
def mul(p,q):
    r={}
    for m,c in p.items():
        for n,d in q.items():
            k=tuple(sorted(m+n));r[k]=r.get(k,F(0))+c*d
    return clean(r)
def derivative(p,i):
    r={}
    for m,c in p.items():
        if i in m:
            n=list(m);n.remove(i)
            r[tuple(n)]=c*m.count(i)
    return r
def product(*ps):
    r=term()
    for p in ps:r=mul(r,p)
    return r
def jacobian(terms):
    r={}
    for t in terms:
        m=t["inputs"]
        for i,k in Counter(m).items():
            n=list(m);n.remove(i);n=tuple(n)
            q=r.setdefault((t["output"],i),{})
            q[n]=q.get(n,F(0))+F(t["coefficient"])*k
    return {ij:clean(q) for ij,q in r.items() if clean(q)}
def matmul(a,b):
    byrow={}
    for (i,j),q in b.items():byrow.setdefault(i,[]).append((j,q))
    r={}
    for (i,k),p in a.items():
        for j,q in byrow.get(k,[]):
            r[(i,j)]=add(r.get((i,j),{}),mul(p,q))
    return {ij:q for ij,q in r.items() if q}
def rowtimes(v,j):
    r={}
    for (i,k),p in j.items():
        if i in v:r[k]=add(r.get(k,{}),mul(v[i],p))
    return {i:p for i,p in r.items() if p}
def drow(p):return {i:derivative(p,i) for i in sorted({i for m in p for i in m})}
def monomial_power(p,n):
    r=term()
    for _ in range(n):r=mul(r,p)
    return r
def substitute(p,values):
    r={}
    for m,c in p.items():
        q=term(c=c)
        for i in m:q=mul(q,values.get(i,term(i)))
        r=add(r,q)
    return r
def encode(p):
    return [{"inputs":list(m),"coefficient":str(c)} for m,c in sorted(p.items())]
def vector_add(*rows):
    r={}
    for row in rows:
        for i,p in row.items():r[i]=add(r.get(i,{}),p)
    return {i:p for i,p in r.items() if p}
def vector_scale(row,p):return {i:mul(q,p) for i,q in row.items() if mul(q,p)}
def module_action(row,transitions):
    return vector_add(*(vector_scale(transitions.get(i,{}),p) for i,p in row.items()))
def module_power(row,transitions,k):
    for _ in range(k):row=module_action(row,transitions)
    return row
def module_substitute(row,values):
    return {i:substitute(p,values) for i,p in row.items() if substitute(p,values)}

def universal_controls():
    # Coefficient variables u,w,p,V,z,y,t,alpha,beta,gamma,delta,epsilon.
    u,w,p,V,z,y,t,al,be,ga,de,ep=[term(i) for i in range(12)]
    tr={"du":{},"dw":{"s":term()},"s":{"r":term()},"r":{},
        "dp":{"du":scale(mul(u,w),2),"dw":mul(u,u)},
        "v":{"a":term()},"a":{"b":term()},"b":{"c":term()},"c":{},
        "dz":{"du":mul(al,V),"v":mul(al,u)},
        "dy":{"dw":mul(be,V),"v":mul(be,w)},
        "dt":{"du":add(product(ga,w,z),scale(product(de,u,y),2)),
              "dw":product(ga,u,z),"dz":product(ga,u,w),"dy":product(de,u,u),
              "dp":mul(ep,V),"v":mul(ep,p)}}
    fourth=module_power({"dt":term()},tr,4)
    A=add(mul(de,be),ep);B=add(mul(ga,al),mul(de,be))
    assert fourth=={"r":product(u,u,V,A),"b":product(u,u,w,B),"c":mul(ep,p)}
    third=module_power({"dt":term()},tr,3)
    assert third=={"r":product(ga,u,z),"s":product(u,u,V,A),
                   "a":product(u,u,w,B),"b":mul(ep,p)}
    # Standard valid gadget: alpha=beta=gamma=epsilon=1, delta=-1.
    good={7:term(),8:term(),9:term(),10:term(c=-1),11:term()}
    g3=module_substitute(third,good)
    assert g3=={"r":mul(u,z),"b":p}
    g4=module_substitute(fourth,good)
    assert g4=={"c":p}
    # d(p*z-u*t)*J^3 is zero when c=0; in fact c is not needed here.
    dnext={"dp":z,"dz":p,"du":scale(t,-1),"dt":scale(u,-1)}
    assert not module_substitute(module_power(dnext,tr,3),good)
    assert module_substitute(module_power(dnext,tr,2),good)=={"du":scale(product(u,u,w,V),-3)}
    # epsilon=lambda while other four parameters are standard.
    lam={7:term(),8:term(),9:term(),10:term(c=-1)}
    f=module_substitute(fourth,lam)
    assert f=={"r":product(u,u,V,add(ep,term(c=-1))),"c":mul(ep,p)}
    # Exact inverse reset, independent of v and W.
    # u,w,p_input,Z,Y,T,v use variables 0,...,6.
    u,W,p0,Z,Y,T,v=[term(i) for i in range(7)]
    P=add(p0,product(u,u,W))
    zi=add(Z,mul(u,v));yi=add(Y,mul(W,v))
    ti=add(T,product(u,W,zi),scale(product(u,u,yi),-1),mul(P,v))
    assert ti==add(T,product(u,W,Z),scale(product(u,u,Y),-1),mul(P,v))
    reset=add(mul(P,zi),scale(mul(u,ti),-1))
    expected=add(mul(p0,Z),scale(mul(u,T),-1),product(u,u,u,Y))
    assert reset==expected
    # Changed coefficient leaves exactly the inherited term.
    lam=term(7)
    ti_lam=add(T,product(u,W,zi),scale(product(u,u,yi),-1),product(lam,P,v))
    assert add(mul(P,zi),scale(mul(u,ti_lam),-1))==add(expected,product(add(term(),scale(lam,-1)),u,P,v))
    return {"universal_fourth_row":{k:encode(q) for k,q in fourth.items()},
            "inverse_reset":encode(reset),"symbolic_identities":11}

def validate_input(d):
    assert set(d)=={"schema_version","base_relative_path","u","w","p","seed_quadratic",
                   "iterations","base_dimension","base_inverse_degree","seed_quadratic_weight",
                   "base_w_weight","lambda_negative_controls"}
    assert d["schema_version"]==1 and d["base_dimension"]==58
    assert [d[k] for k in ("u","w","p")]==[0,47,48]
    assert d["base_inverse_degree"]==207 and d["seed_quadratic_weight"]==90
    assert d["base_w_weight"]==29 and d["iterations"]==5
    assert d["lambda_negative_controls"]==["0","2"]
    q={};seen=set()
    for a in d["seed_quadratic"]:
        assert set(a)=={"inputs","coefficient"}
        m=a["inputs"];c=F(a["coefficient"])
        assert len(m)==2 and m==sorted(m) and all(type(i) is int and 0<=i<58 for i in m)
        assert tuple(m) not in seen and c
        seen.add(tuple(m));q[tuple(m)]=c
    assert q=={(0,51):F(1),(0,53):F(-4),(48,49):F(9)}
    return q
def validate_terms(terms,weights,n):
    assert len(weights)==n and all(type(x) is int and x>0 for x in weights)
    seen=set()
    for a in terms:
        assert set(a)=={"output","inputs","coefficient"}
        i=a["output"];m=a["inputs"];c=F(a["coefficient"])
        assert type(i) is int and 0<=i<n and len(m)==3
        assert all(type(k) is int and 0<=k<i for k in m) and m==sorted(m) and c
        assert (i,tuple(m)) not in seen;seen.add((i,tuple(m)))
        assert sum(weights[k] for k in m)==weights[i]
def append_gadget(terms,n,V,u,w,p,lam=F(1)):
    z,y,t=n,n+1,n+2
    zs=mul(term(u),V);ys=mul(term(w),V)
    ts=add(term(u,w,z),term(u,u,y,c=-1),scale(mul(term(p),V),lam))
    new=copy.deepcopy(terms)
    for i,q in [(z,zs),(y,ys),(t,ts)]:
        for a in encode(q):new.append({"output":i,**a})
    return new,add(term(p,z),term(u,t,c=-1))
def evaluate(p,values):
    out=F(0)
    for m,c in p.items():
        for i in m:c*=values[i]
        out+=c
    return out
def inverse_at(terms,n,targets):
    vals=[F(targets.get(i,0)) for i in range(n)]
    rows={}
    for a in terms:rows.setdefault(a["output"],[]).append(a)
    for i in range(n):
        for a in rows.get(i,[]):
            vals[i]+=evaluate({tuple(a["inputs"]):F(a["coefficient"])},vals)
    return vals
def literal_controls(d,base):
    terms=base["H_terms"];weights=base["weights"]
    assert base["schema_version"]==1 and base["dimension"]==58
    assert len(terms)==1248 and weights[0]==1 and weights[47]==29 and max(weights)==207
    assert weights.count(1)==1 and weights.count(207)==1 and weights[57]==207
    validate_terms(terms,weights,58)
    u,w,p=d["u"],d["w"],d["p"]
    assert not any(a["output"]==u for a in terms)
    assert [a for a in terms if a["output"]==p]==[
        {"output":48,"inputs":[0,0,47],"coefficient":"1"}]
    V=validate_input(d)
    assert {sum(weights[i] for i in m) for m in V}=={90}
    j=jacobian(terms);j2=matmul(j,j);j3=matmul(j2,j);j4=matmul(j3,j)
    assert j3 and not j4 and not any(i==w for i,k in j3)
    assert not rowtimes(drow(V),j3)
    r={k:q for (i,k),q in j2.items() if i==w}
    assert r and V
    # Literal inverse degree207 follows weights and the preserved old top coefficient.
    vals=inverse_at(terms,58,{0:1})
    TOP=F(8803690708443294969,95329784685296)
    assert vals[57]==-2*TOP**7 and vals[47]==TOP
    profiles=[]
    allterms=terms;allweights=list(weights);n=58
    originalV=V
    for k in range(1,d["iterations"]+1):
        vweight=90+32*(k-1)
        allterms,V=append_gadget(allterms,n,V,u,w,p)
        allweights += [1+vweight,29+vweight,31+vweight]
        n+=3
        validate_terms(allterms,allweights,n)
        jj=jacobian(allterms);jj2=matmul(jj,jj);jj3=matmul(jj2,jj);jj4=matmul(jj3,jj)
        assert not jj4 and jj3 and not rowtimes(drow(V),jj3)
        assert not any(i==w for i,col in jj3)
        profiles.append({"iterations":k,"dimension":n,"cubic_terms":len(allterms),
                         "jacobian_cube_entries":len(jj3),"jacobian_fourth_entries":len(jj4),
                         "largest_added_weight":allweights[-1],
                         "returned_quadratic_weight":vweight+32})
    # All newly appended polynomials really satisfy the inverse reset on full input assignments.
    inverse_checks=0
    for fixture in range(3):
        targets={0:1,48:fixture,49:2*fixture,51:fixture,53:-fixture}
        for i in range(58,n):targets[i]=((i+fixture)%5)-2
        values=inverse_at(allterms,n,targets)
        V=originalV
        for k in range(d["iterations"]):
            z,y,t=58+3*k,59+3*k,60+3*k
            V=add(term(p,z),term(u,t,c=-1))
            expected=F(targets.get(p,0))*targets[z]-targets[0]*targets[t]+targets[0]**3*targets[y]
            assert evaluate(V,values)==expected
            inverse_checks+=1
    negatives=[]
    for lam in map(F,d["lambda_negative_controls"]):
        bad,_=append_gadget(terms,58,originalV,u,w,p,lam)
        b=jacobian(bad);b2=matmul(b,b);b4=matmul(matmul(b2,b),b)
        expected={k:scale(mul(term(u,u),mul(originalV,q)),lam-1) for k,q in r.items()}
        actual={col:q for (i,col),q in b4.items() if i==60}
        assert actual==expected and actual
        assert all(i==60 for i,col in b4)
        negatives.append({"lambda":str(lam),"fourth_row_nonzero_columns":len(actual)})
    rejected=0
    for kind in ["duplicate","nonquadratic","zero","coefficient"]:
        bad=copy.deepcopy(d)
        if kind=="duplicate":bad["seed_quadratic"].append(copy.deepcopy(bad["seed_quadratic"][0]))
        elif kind=="nonquadratic":bad["seed_quadratic"][0]["inputs"].append(1)
        elif kind=="zero":bad["seed_quadratic"][0]["coefficient"]="0"
        else:bad["seed_quadratic"][0]["coefficient"]="2"
        try:validate_input(bad)
        except (AssertionError,ValueError):rejected+=1
    assert rejected==4
    return {"base_dimension":58,"base_inverse_degree":207,"first_block_degree_upper_bound":121,
            "later_block_degree_upper_bound":35,"all_iterated_inverse_degrees":207,
            "literal_profiles":profiles,"literal_inverse_reset_checks":inverse_checks,
            "lambda_negative_controls":negatives,"malformed_inputs_rejected":rejected}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root",type=Path)
    ap.add_argument("--input",type=Path,default=HERE/"repeatable_index_four_reset_inputs.json")
    ap.add_argument("--ledger",type=Path,default=HERE/"repeatable_index_four_reset_source_ledger.json")
    ap.add_argument("--output",type=Path,default=HERE/"repeatable_index_four_reset_verification.json")
    a=ap.parse_args()
    if a.repo_root is None:
        a.repo_root=next((p for p in [HERE,*HERE.parents,Path.cwd(),*Path.cwd().parents]
                         if (p/"research/MATH_STATUS.json").is_file()),None)
        if a.repo_root is None:ap.error("--repo-root is required outside the repository")
    repo=a.repo_root.resolve()
    led=json.loads(a.ledger.read_text());assert led["schema_version"]==1
    required={"proof","checker","input","base_proof","base_input","primary_context"}
    assert {s["role"] for s in led["sources"]}==required and len(led["sources"])==len(required)
    paths={}
    for s in led["sources"]:
        root=a.ledger.parent if s["scope"]=="packet" else repo
        p=(root/s["path"]).resolve()
        assert p.is_file() and digest(p)==s["sha256"],"Source pin mismatch: "+str(p)
        paths[s["role"]]=p
    assert paths["checker"]==Path(__file__).resolve() and paths["input"]==a.input.resolve()
    d=json.loads(a.input.read_text());validate_input(d)
    assert paths["base_input"]==(repo/d["base_relative_path"]).resolve()
    base=json.loads(paths["base_input"].read_text())
    result={"schema_version":1,"status":"PASS","source_ledger_sha256":digest(a.ledger),
            "scope":"Universal reset identities plus five literal iterates; no inverse-degree amplification",
            **universal_controls(),**literal_controls(d,base)}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result))
if __name__=="__main__":main()
