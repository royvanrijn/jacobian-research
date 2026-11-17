#!/usr/bin/env python3
"""Standard-library exact derivative bound for the fixed padded companion tensor."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
from verify import decode
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    source=HERE/"de-bondt12-witness.json";base=json.loads(source.read_text())
    initial=[decode(row,12) for row in base["coefficients"]]
    diagonal=[Q(1)]*12;diagonal[10]=Q(1,2);K=[]
    for i,p in enumerate(initial):
        row={}
        for e,c in p.items():
            factor=Q(1)/diagonal[i]
            for a,power in zip(diagonal,e):factor*=a**power
            row[e]=c*factor
        K.append(row)
    H=[]
    for i in range(12):
        p={e+(0,)*5+(1,):c for e,c in K[i].items() if sum(e)==2}
        if i<5:
            e=[0]*18;e[12+i]=1;e[17]=2;p[tuple(e)]=Q(-1)
        H.append(p)
    H += [{e+(0,)*6:c for e,c in K[i].items() if sum(e)==3} for i in range(5)]
    H += [{(0,)*17+(3,):Q(1)}]
    derivatives=[];labels=[]
    for i,p in enumerate(H):
        for j in range(18):
            d={}
            for e,c in p.items():
                if e[j]:
                    ee=list(e);ee[j]-=1;d[tuple(ee)]=c*e[j]
            derivatives.append(d);labels.append([i,j])
    monomials=sorted(set().union(*(set(d) for d in derivatives)))
    rows=[[d.get(e,Q(0)) for e in monomials] for d in derivatives]
    # Incremental Gaussian elimination with independent derivative-row and column witnesses.
    echelon={};independent=[]
    for label,row0 in zip(labels,rows):
        row=list(row0)
        for pivot,basis in sorted(echelon.items()):
            factor=row[pivot]
            if factor:row=[a-factor*b for a,b in zip(row,basis)]
        nonzero=next((j for j,c in enumerate(row) if c),None)
        if nonzero is None:continue
        value=row[nonzero];echelon[nonzero]=[c/value for c in row];independent.append(label)
    rank=len(echelon);assert rank==34
    primary=json.loads((HERE/"fixed-tensor-bound-primary.json").read_text())
    assert primary["source_sha256"]==sha(HERE/"compute_fixed_tensor_bound.py")
    assert primary["base_source_sha256"]==base["source_sha256"]
    tensor_hash=hashlib.sha256(json.dumps([[[list(e),str(c)] for e,c in sorted(p.items())] for p in H],separators=(",",":")).encode()).hexdigest()
    assert primary["tensor_sha256"]==tensor_hash
    assert primary["ambient_quadratic_monomials"]==len(monomials)
    assert primary["nonzero_derivatives"]==sum(any(r) for r in rows)
    assert json.loads((HERE/"fixed-tensor-bound-primary.json").read_text())["rank"]==rank
    receipt=dict(status="pass",source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/"verify.py"),
      base_input_sha256=sha(source),primary_diagnostic_sha256=sha(HERE/"fixed-tensor-bound-primary.json"),
      tensor_sha256=tensor_hash,rank=rank,derivatives=len(rows),nonzero_derivatives=sum(any(r) for r in rows),
      ambient_quadratic_monomials=len(monomials),independent_derivative_rows=independent,
      pivot_monomials=[list(monomials[j]) for j in sorted(echelon)],
      bound="Every simultaneous cube representation of this fixed padded18-output tensor requires at least34 cubes, over any characteristic-zero extension.",
      scope="Fixed tensor only. No lower bound on arbitrary Keller counterexamples, alternative companions, changed padding, special-form dimensions or the global research goal.")
    (HERE/"fixed-tensor-bound-validation.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in receipt.items() if k not in ["independent_derivative_rows","pivot_monomials"]}))
if __name__=="__main__":main()
