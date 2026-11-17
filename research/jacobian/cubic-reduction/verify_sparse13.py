#!/usr/bin/env python3
"""Independent end-to-end replay of the 42-term cubic witness. Standard library only."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
from verify13 import replay
from verify import add,variable,compose,derivative,identity,madd,mmul,determinant3,decode,evaluate
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    raw=json.loads((HERE/"reuse-witness.json").read_text())
    assert raw["source_sha256"]==sha(HERE/"search_reuse.py")
    for name,digest in raw["helper_sha256"].items():assert sha(HERE/name)==digest
    F,points,_,stats=replay(raw);n=len(F);assert n==13
    base=json.loads((HERE/"cubic13-witness.json").read_text())
    assert base["source_sha256"]==sha(HERE/"construct13.py")
    assert base["discovery_sha256"]==sha(HERE/"reuse-witness.json")
    assert F==[decode(row,n) for row in base["coefficients"]]
    data=HERE/"cubic13-sparse-witness.json";d=json.loads(data.read_text())
    assert d["source_sha256"]==sha(HERE/"construct_sparse13.py")
    assert d["base_sha256"]==sha(HERE/"cubic13-witness.json")
    candidate=HERE/"cubic13-sparse-candidate.json"
    assert d["candidate_sha256"]==sha(candidate)
    c=json.loads(candidate.read_text())
    assert c["source_sha256"]==sha(HERE/"optimize13.py")
    assert c["input_sha256"]==sha(HERE/"cubic13-witness.json") and c["moves"]==d["conjugations"]
    for move in d["conjugations"]:
        i,j,t=move["i"],move["j"],Q(move["t"])
        assert 0<=i<n and 0<=j<n and i!=j
        substitution=[variable(k,n) for k in range(n)]
        substitution[i]=add(substitution[i],substitution[j],t)
        transformed=[compose(f,substitution,n) for f in F]
        transformed[i]=add(transformed[i],transformed[j],-t);F=transformed
        for point in points:point[i]-=t*point[j]
    assert F==[decode(row,n) for row in d["coefficients"]]
    assert F==[decode(row,n) for row in c["coefficients"]]
    counts=[sum(sum(e)>=2 for e in f) for f in F]
    assert counts==d["nonlinear_counts"] and sum(counts)==d["nonlinear_total"]==42
    assert max(sum(e) for f in F for e in f)==d["degree"]==3
    assert points==[[Q(v) for v in p] for p in d["points"]]
    images=[[evaluate(f,p) for f in F] for p in points]
    assert len({tuple(p) for p in points})==3
    assert images[0]==images[1]==images[2]==[Q(v) for v in d["common_image"]]
    J=[[derivative(f,j) for j in range(n)] for f in F];zero=(0,)*n
    assert all(not f.get(zero,0) for f in F)
    assert [[p.get(zero,Q(0)) for p in row] for row in J]==[[Q(i==j) for j in range(n)] for i in range(n)]
    A=[r[:3] for r in J[:3]];B=[r[3:] for r in J[:3]]
    C=[r[:3] for r in J[3:]];D=[r[3:] for r in J[3:]];I=identity(n-3,n)
    N=madd(D,I,Q(-1));N2=mmul(N,N);N3=mmul(N2,N)
    assert not any(p for row in N3 for p in row)
    inv=madd(madd(I,N,Q(-1)),N2)
    assert mmul(D,inv)==I and mmul(inv,D)==I
    S=madd(A,mmul(B,mmul(inv,C)),Q(-1))
    assert determinant3(S,n)=={zero:Q(1)}
    receipt=dict(status="pass",source_sha256=sha(Path(__file__)),
        arithmetic_sha256=sha(HERE/"verify.py"),replay_sha256=sha(HERE/"verify13.py"),
        primary_source_sha256=sha(HERE/"construct_sparse13.py"),witness_sha256=sha(data),
        base_sha256=sha(HERE/"cubic13-witness.json"),discovery_sha256=sha(HERE/"reuse-witness.json"),
        dimension=n,degree=3,nonlinear_counts=counts,nonlinear_total=42,
        determinant_one=True,identity_derivative=True,collision_points=3,
        elementary_base_moves=stats,conjugations=d["conjugations"],
        base_determinant_verified=True,entire_transformations_replayed=True,
        all_variables_symbolic=True,full_polynomial_schur_determinant=True,
        carrier_nilpotence_exponent=3,
        scope="Cubic dimension13 upper bound with42 nonlinear occurrences in displayed tangent-normalized coordinates; collision fibre is not the origin fibre. No minimum dimension or support claim. Three-variable bounds [10,48] unchanged.")
    (HERE/"sparse-validation13.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt))
if __name__=="__main__":main()

