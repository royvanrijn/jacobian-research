#!/usr/bin/env python3
"""Exact literal verification of a 58-variable cubic index-four map of inverse degree207."""
if not __debug__:
    raise SystemExit("Refusing optimized Python: exact assertion checks are required.")
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse,copy,hashlib,json
HERE=Path(__file__).resolve().parent
TOP=F(8803690708443294969,95329784685296)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def jacobian(terms):
    out={}
    for t in terms:
        ids=t["inputs"]
        for i,m in Counter(ids).items():
            mon=list(ids);mon.remove(i)
            p=out.setdefault((t["output"],i),{})
            key=tuple(mon);p[key]=p.get(key,F(0))+F(t["coefficient"])*m
    return {ij:{m:c for m,c in p.items() if c} for ij,p in out.items() if any(p.values())}

def multiply(a,b):
    byrow={}
    for (r,c),p in b.items():byrow.setdefault(r,[]).append((c,p))
    out={}
    for (i,k),p in a.items():
        for j,q in byrow.get(k,[]):
            r=out.setdefault((i,j),{})
            for u,c in p.items():
                for v,d in q.items():
                    m=tuple(sorted(u+v));r[m]=r.get(m,F(0))+c*d
    return {ij:{m:c for m,c in p.items() if c} for ij,p in out.items() if any(p.values())}

def poly_mul(p,q):
    r={}
    for u,c in p.items():
        for v,d in q.items():
            mon=tuple(sorted(u+v));r[mon]=r.get(mon,F(0))+c*d
    return {m:c for m,c in r.items() if c}

def poly_derivative(p,i):
    r={}
    for mon,c in p.items():
        k=mon.count(i)
        if k:
            ids=list(mon);ids.remove(i);r[tuple(ids)]=c*k
    return r


def validate(d):
    assert set(d)=={"schema_version","dimension","weights","H_terms"}
    assert d["schema_version"]==1 and d["dimension"]==58
    w=d["weights"];assert len(w)==58
    assert all(type(x) is int and x>0 for x in w)
    assert w.count(1)==1 and w[0]==1 and w.count(207)==1 and w[-1]==207
    assert max(w[:-1])<207 and max(w[:48])==29 and w[47]==29
    assert w[-10:]==[31,59,87,89,117,89,119,147,177,207]
    terms=d["H_terms"];assert len(terms)==1248
    seen=set()
    for t in terms:
        assert set(t)=={"output","inputs","coefficient"}
        o=t["output"];ids=t["inputs"]
        assert type(o) is int and 0<=o<58
        assert isinstance(ids,list) and len(ids)==3 and ids==sorted(ids)
        assert all(type(i) is int and 0<=i<o for i in ids)
        assert (o,tuple(ids)) not in seen;seen.add((o,tuple(ids)))
        assert isinstance(t["coefficient"],str) and F(t["coefficient"])!=0
        assert w[o]==sum(w[i] for i in ids)
    assert len([t for t in terms if t["output"]<48])==1213
    return terms,w

def polynomial_derivation(p,rows):
    out={}
    for i,q in rows.items():
        for mon,c in poly_mul(poly_derivative(p,i),q).items():out[mon]=out.get(mon,F(0))+c
    return {m:c for m,c in out.items() if c}

def universal_controls(terms):
    # u,w,p1,p2,p3,A,B,E,D,C,R,S = 0,...,11.
    mp={0:0,47:1,**{i:i-46 for i in range(48,58)}}
    ps={}
    for t in terms:
        if t["output"]>=48:
            ps.setdefault(mp[t["output"]],{})[tuple(sorted(mp[i] for i in t["inputs"]))]=F(t["coefficient"])
    dr={2:{(0,0):F(1)},3:{(0,1):F(2)},4:{(1,1):F(3)},
        5:{(0,3):F(3),(1,2):F(-18)},6:{(1,3):F(6),(0,4):F(-2)},
        7:{(0,3):F(3)},8:{(2,3):F(-9),(0,5):F(-1)},
        9:{(0,6):F(12),(1,5):F(8),(2,4):F(24),(3,3):F(-21)},
        10:{(3,5):F(6),(2,6):F(24)}}
    du={8:{(0,0,1,1,1,1):F(15)},9:{(0,1,1,1,1,1):F(-24)},
        10:{(0,1,1,1,1,3):F(-126),(0,0,1,1,1,4):F(60),(1,1,1,1,1,2):F(72)}}
    for i in range(5,12):assert not polynomial_derivation(ps[i],dr)
    assert not polynomial_derivation(ps[11],du)
    return 8

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input",type=Path,default=HERE/"cubic_index_four_syzygy58.json")
    ap.add_argument("--ledger",type=Path,default=HERE/"cubic_index_four_syzygy58_source_ledger.json")
    ap.add_argument("--output",type=Path,default=HERE/"cubic_index_four_syzygy58_verification.json")
    a=ap.parse_args()
    ledger=json.loads(a.ledger.read_text())
    assert ledger["schema_version"]==1 and len(ledger["sources"])==4
    assert {s["role"] for s in ledger["sources"]}=={"proof","literal_input","checker","primary_context"}
    for item in ledger["sources"]:
        p=(a.ledger.parent/item["path"]).resolve()
        assert p.is_file() and digest(p)==item["sha256"],"Source pin mismatch: "+str(p)
        if item["role"]=="literal_input":assert p==a.input.resolve()
        if item["role"]=="checker":assert p==Path(__file__).resolve()
    d=json.loads(a.input.read_text());terms,w=validate(d)
    symbolic=universal_controls(terms)
    base=[t for t in terms if t["output"]<48]
    c=jacobian(base);c2=multiply(c,c);c3=multiply(c2,c)
    assert not c3 and c2[(47,0)][(2,2,2,10)]==36
    assert c2[(47,1)][(0,1,1,30)]==F(393016350032033013,3404635167332)
    assert all(t["output"]!=0 for t in base)
    j=jacobian(terms);j2=multiply(j,j);j3=multiply(j2,j);j4=multiply(j3,j)
    assert j3 and not j4 and j3[(48,0)][(0,0,2,2,2,10)]==36
    vals=[F(1)]+[F(0)]*57
    for i in range(1,58):
        for t in terms:
            if t["output"]==i:
                x,y,z=t["inputs"];vals[i]+=F(t["coefficient"])*vals[x]*vals[y]*vals[z]
    assert vals[47]==TOP
    assert vals[-10:]==[TOP,TOP**2,TOP**3,-5*TOP**3,TOP**4,
                         TOP**3,-TOP**4,-5*TOP**5,-TOP**6,-2*TOP**7]
    bad=copy.deepcopy(terms)
    k=next(k for k,t in enumerate(bad) if t["output"]==57)
    bad[k]["coefficient"]=str(F(bad[k]["coefficient"])+1)
    b=jacobian(bad);b2=multiply(b,b);b4=multiply(multiply(b2,b),b)
    assert b4
    rejected=0
    for mode in ("duplicate","zero","index","weight"):
        bad=copy.deepcopy(d)
        if mode=="duplicate":bad["H_terms"][-1]=copy.deepcopy(bad["H_terms"][-2])
        elif mode=="zero":bad["H_terms"][-1]["coefficient"]="0"
        elif mode=="index":bad["H_terms"][-1]["inputs"][0]=58
        else:bad["weights"][-1]=209
        try:validate(bad)
        except (AssertionError,ValueError):rejected+=1
    assert rejected==4
    receipt={"schema_version":1,"status":"PASS","source_ledger_sha256":digest(a.ledger),
      "dimension":58,"cubic_terms":len(terms),"universal_polynomial_cancellations":symbolic,
      "base_jacobian_cube_entries":len(c3),"full_jacobian_cube_nonzero_entries":len(j3),
      "full_jacobian_fourth_power_entries":len(j4),"nilpotency_index":4,
      "inverse_degree":207,"inverse_line_top_coefficient":str(vals[-1]),
      "coefficient_perturbation_fourth_power_nonzero_entries":len(b4),"malformed_inputs_rejected":rejected}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt))
if __name__=="__main__":main()
