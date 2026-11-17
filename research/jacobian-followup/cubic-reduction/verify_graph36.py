#!/usr/bin/env python3
"""Independent exact dual-certificate audit. No SymPy, FLINT or rank computations."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
from math import comb
import hashlib,json
from verify import add,mul,variable,compose,decode
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def downset(exponents):
    return set(e for top in exponents for e in product(*(range(q+1) for q in top)))
def sliced(F,j):
    n=len(F);g=add(variable(j,n),F[j],Q(-1))
    assert not any(e[j] for e in g)
    sub=[variable(k,n) for k in range(n)];sub[j]=g
    kept={i:compose(f,sub,n) for i,f in enumerate(F) if i!=j}
    zero=(0,)*n
    for i,f in kept.items():
        assert not f.get(zero,0) and not any(e[j] for e in f)
        assert [f.get(tuple(int(a==b) for a in range(n)),Q(0)) for b in range(n)]==[Q(b==i) for b in range(n)]
    return kept
def verify_functional(F,data,degree_bound):
    j=data["pivot"]-1;i=data["component"]-1;n=len(F)
    assert i!=j and 3<=j<12
    kept=sliced(F,j);allowed=[a-1 for a in data["allowed_outputs"]]
    assert allowed==[a for a in kept if a!=i] and len(allowed)==10
    weights={tuple(row["exponent"]):Q(row["weight"]) for row in data["functional"]}
    assert len(weights)==len(data["functional"]) and all(weights.values())
    assert all(len(e)==n and min(e)>=0 and sum(e)>3 and not e[j] for e in weights)
    D=max(map(sum,weights));T=sorted(downset(weights));index={e:k for k,e in enumerate(T)}
    assert all(tuple(e[a]-int(a==b) for a in range(n)) in index for e in T for b in range(n) if e[b])
    assert sum((w*kept[i].get(e,Q(0)) for e,w in weights.items()),Q(0))==1
    # Multiplication by each raw output is a sparse linear action on this quotient.
    actions=[]
    for a in allowed:
        action={}
        assert not kept[a].get((0,)*n,0)
        for e,k in index.items():
            entries=[]
            for f,c in kept[a].items():
                g=tuple(q+r for q,r in zip(e,f))
                if g in index:entries.append((index[g],c))
            if entries:action[k]=entries
        actions.append(action)
    lambda_by_index={index[e]:w for e,w in weights.items()}
    checked=0;zero_words=0;m=len(actions)
    def apply(action,vector):
        out={}
        for a,c in vector.items():
            for b,d in action.get(a,[]):out[b]=out.get(b,Q(0))+c*d
        return {a:c for a,c in out.items() if c}
    def visit(vector,depth,last):
        nonlocal checked,zero_words
        if depth==degree_bound:return
        for a in range(last,m):
            child=apply(actions[a],vector);next_depth=depth+1
            if not child:
                zero_words+=sum(comb(m-a+h-1,h) for h in range(degree_bound-next_depth+1))
                continue
            checked+=1
            assert sum((w*child.get(k,Q(0)) for k,w in lambda_by_index.items()),Q(0))==0,(j+1,i+1,next_depth,a)
            visit(child,next_depth,a)
    visit({index[(0,)*n]:Q(1)},0,0)
    expected=sum(comb(m+d-1,d) for d in range(1,degree_bound+1))
    assert checked+zero_words==expected
    return dict(pivot=j+1,component=i+1,functional_terms=len(weights),maximum_jet=D,
       quotient_dimension=len(T),target_degree_checked=degree_bound,total_target_monomials=expected,
       nonzero_words_checked=checked,zero_words_covered=zero_words)
def positive_control():
    # Exact BCR2 square completion: the degree-four graph defect is removable.
    n=12;x=variable(0,n);y=variable(1,n);a=variable(3,n);b=variable(7,n)
    x2=mul(x,x);y2=mul(y,y);raw=add(a,mul(x2,y2));generator=add(b,mul(x,y))
    corrected=add(raw,mul(generator,generator),Q(-1))
    assert max(map(sum,raw))==4 and max(map(sum,corrected))==3
    expected=add(add(a,mul(mul(x,y),b),Q(-2)),mul(b,b),Q(-1))
    assert corrected==expected
    # The tempting functional on the single defect coefficient fails to kill P=Y8^2.
    e=(2,2)+(0,)*10
    assert raw[e]==mul(generator,generator)[e]==1
def main():
    source=HERE/"cubic12-integer-witness.json";d=json.loads(source.read_text())
    F=[decode(row,12) for row in d["coefficients"]]
    assert (d["dimension"],d["degree"],d["nonlinear_total"])==(12,3,36)
    raw_manifest=json.loads((HERE/"graph36/checkpoint.json").read_text())
    wide_manifest=json.loads((HERE/"adaptive36-wide/checkpoint.json").read_text())
    assert raw_manifest["status"]==wide_manifest["status"]=="complete"
    assert {c["pivot"] for c in raw_manifest["cases"]}==set(range(4,13))
    assert {c["pivot"] for c in wide_manifest["cases"]}==set(range(4,13))
    finite=[];unbounded=[];hashes={}
    for j in range(4,13):
        raw_path=HERE/"graph36"/("pivot-%02d.json"%j);raw=json.loads(raw_path.read_text())
        assert raw["input_sha256"]==sha(source) and raw["source_sha256"]==sha(HERE/"graph36.py")
        sep=next(row for row in raw["tests"] if row["status"]=="separated")
        raw_input=dict(pivot=j,**sep);finite.append(verify_functional(F,raw_input,3))
        wide_path=HERE/"adaptive36-wide"/("pivot-%02d.json"%j);wide=json.loads(wide_path.read_text())
        assert wide["status"]=="arbitrary-degree-obstructed" and wide["input_sha256"]==sha(source)
        assert wide["source_sha256"]==sha(HERE/"graph36_adaptive_wide.py")
        assert wide["algebra_sha256"]==sha(HERE/"graph36_all_degree.py")
        assert wide["linear_solver_sha256"]==sha(HERE/"graph36.py")
        D=max(sum(row["exponent"]) for row in wide["functional"])
        unbounded.append(verify_functional(F,wide,D))
        hashes[str(j)]=dict(raw_sha256=sha(raw_path),arbitrary_degree_sha256=sha(wide_path))
        print("verified",j,"jet",D,"words",unbounded[-1]["total_target_monomials"],flush=True)
    positive_control()
    receipt=dict(status="pass",source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/"verify.py"),
      input_sha256=sha(source),finite_target_degree_cases=finite,arbitrary_target_degree_cases=unbounded,
      case_hashes=hashes,all_nine_arbitrary_degree_certificates=True,
      normalized_raw_source_jets_checked=True,dual_annihilation_and_nonzero_defect=True,
      commutative_target_monomial_coverage_checked=True,zero_branch_pruning_coverage_checked=True,
      higher_degree_terms_excluded_by_positive_order=True,
      known_BCR2_square_completion_control=True,
      scope="For each of nine fixed literal graph deletions in KCUB3, one raw retained component cannot be made cubic by any polynomial in the other ten raw retained outputs, of any degree. No sequential-stage, general source-change or minimum-dimension exclusion.")
    (HERE/"graph36-validation.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in receipt.items() if k not in ("finite_target_degree_cases","arbitrary_target_degree_cases","case_hashes")}))
if __name__=="__main__":main()

