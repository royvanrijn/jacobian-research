#!/usr/bin/env python3
"""Independent rational audit of the 36-term twelve-variable cubic conjugate."""
from pathlib import Path

import tempfile as _check_tempfile

_CHECK_OUTPUT_DIRECTORY = None

def _check_output(original):
    """Write new runs to a temporary directory; preserve deposited receipts."""
    global _CHECK_OUTPUT_DIRECTORY
    if _CHECK_OUTPUT_DIRECTORY is None:
        _CHECK_OUTPUT_DIRECTORY = Path(_check_tempfile.mkdtemp(prefix="math-check-"))
        print("Verification output:", _CHECK_OUTPUT_DIRECTORY)
    return _CHECK_OUTPUT_DIRECTORY / Path(original).name


from fractions import Fraction as Q
import hashlib,json
from verify import add,variable,compose,derivative,identity,madd,mmul,determinant3,decode,evaluate
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def audit(F):
    n=len(F);zero=(0,)*n
    assert all(not f.get(zero,0) for f in F)
    J=[[derivative(f,j) for j in range(n)] for f in F]
    assert [[p.get(zero,Q(0)) for p in row] for row in J]==[[Q(i==j) for j in range(n)] for i in range(n)]
    A=[r[:3] for r in J[:3]];B=[r[3:] for r in J[:3]]
    C=[r[:3] for r in J[3:]];D=[r[3:] for r in J[3:]]
    I=identity(n-3,n);N=madd(D,I,Q(-1));N2=mmul(N,N);N3=mmul(N2,N)
    assert not any(p for row in N3 for p in row)
    inv=madd(madd(I,N,Q(-1)),N2)
    assert mmul(D,inv)==I and mmul(inv,D)==I
    S=madd(A,mmul(B,mmul(inv,C)),Q(-1))
    assert determinant3(S,n)=={zero:Q(1)}
def main():
    path=HERE/"cubic12-sparse-witness.json";d=json.loads(path.read_text())
    base_path=HERE/"de-bondt12-witness.json";base=json.loads(base_path.read_text())
    cand_path=HERE/"de-bondt12-sparse-candidate.json";candidate=json.loads(cand_path.read_text())
    F=[decode(row,12) for row in base["coefficients"]];audit(F)
    assert sum(sum(e)>=2 for f in F for e in f)==40
    points=[[Q(v) for v in p] for p in base["points"]]
    assert [evaluate(f,points[0]) for f in F]==[evaluate(f,points[1]) for f in F]
    assert d["conjugations"]==candidate["moves"]
    for move in d["conjugations"]:
        i,j,t=move["i"],move["j"],Q(move["t"]);assert i!=j and 0<=i<12 and 0<=j<12
        sub=[variable(k,12) for k in range(12)];sub[i]=add(sub[i],sub[j],t)
        F=[compose(f,sub,12) for f in F];F[i]=add(F[i],F[j],-t)
        for point in points:point[i]-=t*point[j]
    assert F==[decode(row,12) for row in d["coefficients"]]
    assert F==[decode(row,12) for row in candidate["coefficients"]]
    audit(F)
    counts=[sum(sum(e)>=2 for e in f) for f in F]
    assert counts==d["nonlinear_counts"] and sum(counts)==d["nonlinear_total"]==36
    assert max(sum(e) for f in F for e in f)==d["degree"]==3
    assert d["dimension"]==12 and points==[[Q(v) for v in p] for p in d["points"]]
    images=[[evaluate(f,p) for f in F] for p in points]
    assert points[0]!=points[1] and images[0]==images[1]==[Q(v) for v in d["common_image"]]
    receipt=dict(status="pass",source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/"verify.py"),
      primary_source_sha256=sha(HERE/"construct_sparse12.py"),witness_sha256=sha(path),
      base_sha256=sha(base_path),base_source_sha256=base["source_sha256"],
      dimension=12,degree=3,baseline_nonlinear_total=40,nonlinear_total=36,nonlinear_counts=counts,
      determinant_one=True,identity_derivative=True,collision_points=2,
      baseline_and_final_full_polynomial_determinants=True,all_variables_symbolic=True,
      both_carrier_inverse_identities=True,carrier_nilpotence_exponent=3,
      conjugations=d["conjugations"],
      scope="Sparsity improvement40 to36 in the already verified BCRDB1 twelve-variable cubic map. No dimension improvement, global optimum or priority claim. Three-variable normalized bounds [10,48] unchanged.")
    (_check_output(HERE/"sparse-validation12.json")).write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt))
if __name__=="__main__":main()

