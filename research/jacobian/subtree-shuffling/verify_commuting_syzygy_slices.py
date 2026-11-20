#!/usr/bin/env python3
"""Exact localized syzygy slices and the post-S scalar obstruction."""
if not __debug__:raise SystemExit("Refusing optimized Python: exact assertion checks are required.")
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from math import factorial
import argparse,hashlib,json,copy
HERE=Path(__file__).resolve().parent
NV=11
ZERO=(0,)*NV
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mon(spec,c=1):
    e=[0]*NV
    for i,k in spec.items():e[i]=k
    return {tuple(e):F(c)} if c else {}
def plus(*ps):
    out={}
    for p in ps:
        for e,c in p.items():out[e]=out.get(e,F(0))+c
    return {e:c for e,c in out.items() if c}
def scale(p,c):return {e:F(c)*d for e,d in p.items() if c*d}
def times(p,q):
    out={}
    for e,c in p.items():
        for f,d in q.items():
            k=tuple(x+y for x,y in zip(e,f));out[k]=out.get(k,F(0))+c*d
    return {e:c for e,c in out.items() if c}
def diff(p,i):
    out={}
    for e,c in p.items():
        if e[i]:
            k=list(e);k[i]-=1;out[tuple(k)]=c*e[i]
    return out
def apply(p,rows):return plus(*(times(diff(p,i),q) for i,q in rows.items()))
def poly_power(p,n):
    assert n>=0
    out=mon({})
    for _ in range(n):out=times(out,p)
    return out
def rows():
    dr={2:mon({0:2}),3:mon({0:1,1:1},2),4:mon({1:2},3),
      5:plus(mon({0:1,3:1},3),mon({1:1,2:1},-18)),
      6:plus(mon({1:1,3:1},6),mon({0:1,4:1},-2)),
      7:mon({0:1,3:1},3),
      8:plus(mon({2:1,3:1},-9),mon({0:1,5:1},-1)),
      9:plus(mon({0:1,6:1},12),mon({1:1,5:1},8),mon({2:1,4:1},24),mon({3:2},-21)),
      10:plus(mon({3:1,5:1},6),mon({2:1,6:1},24))}
    du={8:mon({0:2,1:4},15),9:mon({0:1,1:5},-24),
      10:plus(mon({0:1,1:4,3:1},-126),mon({0:2,1:3,4:1},60),mon({1:5,2:1},72))}
    return dr,du
def projection(p,t,dr):
    out={};current=p
    for k in range(8):
        if not current:return out
        out=plus(out,scale(times(poly_power(t,k),current),F((-1)**k,factorial(k))))
        current=apply(current,dr)
    raise AssertionError("Unexpected coordinate nilpotency length")
def inverse_line(p):
    # Values of u,w,p1,p2,p3,A,B,E,D,C,R, with formal w=W.
    values=[(1,0,1),(0,1,1),(2,1,1),(1,2,1),(0,3,1),
            (2,3,-5),(1,4,1),(2,3,1),(3,4,-1),(2,5,-5),(3,6,-1)]
    out={}
    for e,c in p.items():
        a=b=0;q=c
        for k,(r,s,d) in zip(e,values):a+=k*r;b+=k*s;q*=F(d)**k
        out[(a,b)]=out.get((a,b),F(0))+q
    return {e:c for e,c in out.items() if c}
def algebra_checks():
    dr,du=rows()
    for i in range(NV):
        z=mon({i:1})
        assert apply(apply(z,dr),du)==apply(apply(z,du),dr)
    t=mon({0:-2,2:1})
    assert apply(t,dr)==mon({}) and not apply(t,du)
    normalized={i:projection(mon({i:1}),t,dr) for i in range(3,11)}
    explicit={
      3:plus(mon({3:1}),mon({0:-1,1:1,2:1},-2)),
      4:plus(mon({4:1}),mon({0:-2,1:2,2:1},-3)),
      5:plus(mon({5:1}),mon({0:-1,2:1,3:1},-3),mon({0:-2,1:1,2:2},12)),
      6:plus(mon({6:1}),mon({0:-2,1:1,2:1,3:1},-6),mon({0:-1,2:1,4:1},2),mon({0:-3,1:2,2:2},3)),
      7:plus(mon({7:1}),mon({0:-1,2:1,3:1},-3),mon({0:-2,1:1,2:2},3)),
      8:plus(mon({8:1}),mon({0:-1,2:1,5:1}),mon({0:-2,2:2,3:1},3),mon({0:-3,1:1,2:3},4))}
    for i,p in normalized.items():
        assert not apply(p,dr)
        assert p.get(tuple(int(j==i) for j in range(NV)))==1
        for e in p:
            assert not any(e[j] for j in range(i+1,NV))
            assert e[i] in (0,1)
            if e[i]:assert sum(abs(x) for x in e)==1
    for i,p in explicit.items():assert normalized[i]==p
    a=mon({0:2,1:4},15)
    assert apply(normalized[8],du)==a
    assert apply(normalized[9],du)==mon({0:1,1:5},-24)
    assert apply(normalized[10],du)==du[10]
    v=times(normalized[8],mon({0:-2,1:-4},F(1,15)))
    assert not apply(v,dr) and apply(v,du)==mon({})
    joint={i:normalized[i] for i in (3,4,5,6,7)}
    for i in (9,10):joint[i]=plus(normalized[i],scale(times(v,apply(normalized[i],du)),-1))
    for i,p in joint.items():
        assert not apply(p,dr) and not apply(p,du)
        assert p.get(tuple(int(j==i) for j in range(NV)))==1
        assert all(not any(e[j] for j in range(i+1,NV)) for e in p)
    assert [inverse_line(normalized[i]) for i in (3,4,5,6)]==[
      {(1,2):F(-1)},{(0,3):F(-2)},{(2,3):F(4)},{}]
    assert inverse_line(normalized[8])=={(3,4):F(1)}
    # Actual algebraic corruption controls, rather than sampled equality tests.
    wrong=plus(normalized[5],mon({0:-2,1:1,2:2},-1))
    assert apply(wrong,dr)
    broken=copy.deepcopy(dr);broken[2]=scale(broken[2],2)
    assert apply(apply(mon({10:1}),broken),du)!=apply(apply(mon({10:1}),du),broken)
    return {"commutator_generator_identities":NV,"slice_identities":4,
      "normalized_coordinates":8,"joint_kernel_coordinates":7,
      "algebraic_negative_controls":2,"line_values_checked":5}

def jacobian(terms):
    out={}
    for t in terms:
        for i,k in Counter(t["inputs"]).items():
            ids=list(t["inputs"]);ids.remove(i)
            p=out.setdefault((t["output"],i),{});key=tuple(ids)
            p[key]=p.get(key,F(0))+F(t["coefficient"])*k
    return {ij:{m:c for m,c in p.items() if c} for ij,p in out.items() if any(p.values())}
def multiply(a,b):
    byrow={}
    for (i,j),p in b.items():byrow.setdefault(i,[]).append((j,p))
    out={}
    for (i,k),p in a.items():
        for j,q in byrow.get(k,[]):
            r=out.setdefault((i,j),{})
            for u,c in p.items():
                for v,d in q.items():
                    key=tuple(sorted(u+v));r[key]=r.get(key,F(0))+c*d
    return {ij:{m:c for m,c in p.items() if c} for ij,p in out.items() if any(p.values())}
def matrix_checks(data):
    assert data["dimension"]==58 and len(data["H_terms"])==1248
    j=jacobian(data["H_terms"]);j2=multiply(j,j);j3=multiply(j2,j);j4=multiply(j3,j)
    assert j3 and not j4
    base=jacobian([t for t in data["H_terms"] if t["output"]<48])
    b2=multiply(base,base)
    assert not multiply(b2,base)
    assert b2[(47,1)][(0,1,1,30)]==F(393016350032033013,3404635167332)
    # Verify EVERY table entry against the actual full58-variable cube.
    dr,du=rows()
    local_to_global={0:0,1:47,**{i:i+46 for i in range(2,NV)}}
    def ordinary(p):
        out={}
        for e,c in p.items():
            assert all(k>=0 for k in e)
            ids=tuple(sorted(local_to_global[i] for i,k in enumerate(e) for _ in range(k)))
            out[ids]=out.get(ids,F(0))+c
        return {m:c for m,c in out.items() if c}
    def product(p,q):
        x=multiply({(0,0):p},{(0,0):q})
        return x.get((0,0),{})
    table_checks=0
    for local_row in range(2,NV):
        actual_row=local_to_global[local_row]
        d=ordinary(dr.get(local_row,{}));e=ordinary(du.get(local_row,{}))
        for col in range(58):
            want=product(d,b2.get((47,col),{}))
            if col==0:
                for mon,c in e.items():want[mon]=want.get(mon,F(0))+c
            want={mon:c for mon,c in want.items() if c}
            assert j3.get((actual_row,col),{})==want
            table_checks+=1
    assert table_checks==522
    assert all(not p for (i,col),p in j3.items() if 48<=i<57 and col>=48)
    # Necessary-only negative control: Q_S=0 alone does not preserve Jfourth.
    q_gradient={(0,48):{(48,48):F(3)}} # Q=p1^3, so Q_S=0.
    bad_fourth_row=multiply(q_gradient,j3)
    bad_coefficient=F(1179049050096099039,3404635167332)
    assert bad_fourth_row[(0,1)][(0,0,0,1,1,30,48,48)]==bad_coefficient
    expected={
      48:{(0,47,47,47,47,49):F(-108),(47,47,47,47,47,48):F(216)},
      49:{(0,0,47,47,47,49):F(216),(0,47,47,47,47,48):F(-432)},
      50:{(0,0,0,47,47,49):F(-108),(0,0,47,47,47,48):F(216)}}
    for col,p in expected.items():
        assert j3[(57,col)]==p
        assert not any(q for (i,c),q in j3.items() if c==col and i!=57)
    omega={(0,48):{(47,47):F(1)},(0,49):{(0,47):F(-2)},(0,50):{(0,0):F(1)}}
    assert not multiply(omega,j)
    wrong=copy.deepcopy(omega);wrong[(0,49)][(0,47)]=F(-1)
    assert multiply(wrong,j)
    return {"jacobian_cube_nonzero_entries":len(j3),"jacobian_fourth_power_entries":len(j4),
      "pre_S_row_column_identities":table_checks,"pre_S_new_fibre_columns_zero":True,
      "post_S_columns_checked":3,"omega_J_entries":0,"omega_coefficient_negative_control":True,
      "negative_converse_Q":"p1^3","negative_converse_Q_S":0,
      "negative_converse_fourth_row_nonzero_entries":len(bad_fourth_row),
      "negative_converse_coefficient":str(bad_coefficient)}
def detect_repo(arg):
    if arg:return arg.resolve()
    for p in list(HERE.parents)+[Path.cwd(),*Path.cwd().parents]:
        if (p/"research/MATH_STATUS.json").is_file():return p
    raise SystemExit("Pass --repo-root for the canonical literal input.")
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root",type=Path)
    ap.add_argument("--input",type=Path,default=HERE/"commuting_syzygy_slice_inputs.json")
    ap.add_argument("--ledger",type=Path,default=HERE/"commuting_syzygy_slice_source_ledger.json")
    ap.add_argument("--output",type=Path,default=HERE/"commuting_syzygy_slice_verification.json")
    a=ap.parse_args();repo=detect_repo(a.repo_root)
    ledger=json.loads(a.ledger.read_text());assert ledger["schema_version"]==1
    assert {s["role"] for s in ledger["sources"]}=={"proof","checker","input","canonical_proof","literal_data","classical_linear_context","classical_covariant_context"}
    for item in ledger["sources"]:
        base=repo/"research" if item["scope"]=="repo" else a.ledger.parent
        p=(base/item["path"]).resolve()
        assert p.is_file() and digest(p)==item["sha256"],"Source pin mismatch: "+str(p)
        if item["role"]=="checker":assert p==Path(__file__).resolve()
        if item["role"]=="input":assert p==a.input.resolve()
    inp=json.loads(a.input.read_text())
    assert inp["schema_version"]==1 and inp["dimension"]==58 and inp["pre_S_dimension"]==57
    assert inp["u"]==0 and inp["w"]==47 and inp["auxiliary_coordinates"]==list(range(48,58))
    data_path=repo/"research"/inp["canonical_literal_data"]
    assert data_path.resolve().is_relative_to(repo.resolve())
    out={"status":"PASS","source_ledger_sha256":digest(a.ledger),"algebra":algebra_checks(),
         "literal_matrix":matrix_checks(json.loads(data_path.read_text()))}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2)+"\n");print(json.dumps(out))
if __name__=="__main__":main()
