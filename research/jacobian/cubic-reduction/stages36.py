#!/usr/bin/env python3
"""Find cubic-pullback target relations that mix an obstructed raw coordinate."""
from pathlib import Path
from itertools import combinations_with_replacement
import hashlib,json,time
import sympy as s
from flint import fmpq_mat
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    source=HERE/"cubic12-integer-witness.json";d=json.loads(source.read_text());X=s.symbols(d["variables"])
    F=[s.expand(sum(s.Rational(t["coefficient"])*s.prod(v**n for v,n in zip(X,t["exponent"])) for t in row)) for row in d["coefficients"]]
    folder=HERE/"stages36";folder.mkdir(exist_ok=True);cases=[]
    for j in range(3,12):
        start=time.monotonic();g=s.expand(X[j]-F[j]);kept={a:s.expand(f.subs(X[j],g)) for a,f in enumerate(F) if a!=j}
        case=json.loads((HERE/"adaptive36-wide"/("pivot-%02d.json"%(j+1))).read_text());i=case["component"]-1
        labels=list(kept);words=[ix for degree in range(1,4) for ix in combinations_with_replacement(labels,degree)]
        products=[s.Poly(s.prod(kept[a] for a in ix),*X) for ix in words]
        columns=[{e:c for e,c in p.terms() if sum(e)>3} for p in products]
        monomials=sorted({e for p in columns for e in p})
        matrix=fmpq_mat([[str(p.get(e,0)) for p in columns] for e in monomials])
        rr,rank=matrix.rref();pivots=[next(k for k in range(len(words)) if rr[a,k]) for a in range(rank)]
        free=[k for k in range(len(words)) if k not in pivots];bad_slots=[k for k,w in enumerate(words) if i in w]
        candidates=[];mixed_basis=0
        for f in free:
            vector={f:s.Rational(1)}
            for a,p in enumerate(pivots):
                c=-s.Rational(str(rr[a,f]))
                if c:vector[p]=c
            if not any(k in vector for k in bad_slots):continue
            mixed_basis+=1
            # Drop linear terms in already cubic raw outputs; they are trivial kernel directions.
            vector={k:c for k,c in vector.items() if len(words[k])!=1 or products[k].total_degree()>3}
            used={a for k in vector for a in words[k]}
            eligible=[a for a in labels if a!=i and a not in used and s.Poly(kept[a],*X).total_degree()<=3]
            correction=s.Poly(sum(c*products[k].as_expr() for k,c in vector.items()),*X)
            assert correction.total_degree()<=3
            polynomial=[dict(indices=[a+1 for a in words[k]],coefficient=str(c)) for k,c in vector.items()]
            candidates.append(dict(eligible_targets=[a+1 for a in eligible],target_polynomial=polynomial,
                pullback=[dict(exponent=list(e),coefficient=str(c)) for e,c in correction.terms()]))
        result=dict(status="complete",pivot=j+1,obstructed_component=i+1,target_degree_bound=3,
            matrix_rows=len(monomials),matrix_columns=len(words),rank=rank,nullity=len(free),
            mixed_kernel_basis_vectors=mixed_basis,candidates=candidates,
            source_sha256=sha(Path(__file__)),input_sha256=sha(source),seconds=round(time.monotonic()-start,3),
            boundary="Complete raw polynomial target space through degree3. Zero mixed projection proves no such raw cubic-pullback relation; eligible candidates are elementary first steps, not completed dimension reductions.")
        (folder/("pivot-%02d.json"%(j+1))).write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
        cases.append(dict(pivot=j+1,file="pivot-%02d.json"%(j+1),mixed_basis=mixed_basis,eligible=sum(bool(c["eligible_targets"]) for c in candidates)))
        (folder/"checkpoint.json").write_text(json.dumps(dict(status="running",source_sha256=sha(Path(__file__)),input_sha256=sha(source),cases=cases),indent=2,sort_keys=True)+"\n")
        print(json.dumps({k:result[k] for k in ["pivot","obstructed_component","rank","nullity","mixed_kernel_basis_vectors","seconds"]}),flush=True)
    (folder/"checkpoint.json").write_text(json.dumps(dict(status="complete",source_sha256=sha(Path(__file__)),input_sha256=sha(source),cases=cases),indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()

