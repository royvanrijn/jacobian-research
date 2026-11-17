#!/usr/bin/env python3
"""Independent integer-coefficient cubic verification, including its rational base replay."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
import verify_sparse12 as rational
from verify import decode,evaluate
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    rational.main()
    path=HERE/"cubic12-integer-witness.json";d=json.loads(path.read_text())
    assert d["source_sha256"]==sha(HERE/"construct_integer12.py")
    base_path=HERE/"cubic12-sparse-witness.json";base=json.loads(base_path.read_text())
    assert d["rational_witness_sha256"]==sha(base_path)
    diagonal=[Q(v) for v in d["diagonal"]]
    assert len(diagonal)==12 and all(diagonal)
    assert diagonal==[Q(1)]*8+[Q(1,2),Q(1),Q(1,2),Q(1)]
    old=[decode(row,12) for row in base["coefficients"]];F=[]
    for i,f in enumerate(old):
        row={}
        for e,c in f.items():
            scale=Q(1)/diagonal[i]
            for v,n in zip(diagonal,e):scale*=v**n
            row[e]=c*scale
        F.append(row)
    assert F==[decode(row,12) for row in d["coefficients"]]
    assert all(c.denominator==1 for f in F for c in f.values())
    rational.audit(F)
    points=[[Q(v)/scale for v,scale in zip(p,diagonal)] for p in base["points"]]
    assert points==[[Q(v) for v in p] for p in d["points"]]
    assert all(v.denominator==1 for p in points for v in p) and points[0]!=points[1]
    images=[[evaluate(f,p) for f in F] for p in points]
    assert images[0]==images[1]==[Q(v) for v in d["common_image"]]
    counts=[sum(sum(e)>=2 for e in f) for f in F]
    assert counts==d["nonlinear_counts"] and sum(counts)==36
    assert max(sum(e) for f in F for e in f)==3 and d["dimension"]==12
    receipt=dict(status="pass",source_sha256=sha(Path(__file__)),
        arithmetic_sha256=sha(HERE/"verify.py"),rational_checker_sha256=sha(HERE/"verify_sparse12.py"),
        primary_source_sha256=sha(HERE/"construct_integer12.py"),witness_sha256=sha(path),
        rational_validation_sha256=sha(HERE/"sparse-validation12.json"),
        dimension=12,degree=3,baseline_nonlinear_total=40,nonlinear_total=36,nonlinear_counts=counts,
        integral_coefficients=True,integral_collision_points=True,determinant_one=True,identity_derivative=True,
        collision_points=2,full_polynomial_determinant=True,all_variables_symbolic=True,
        source_and_target_conjugation_replayed=True,
        scope="36-term integral conjugate of BCRDB1's12-variable cubic. Known dimension retained; no dimension record or optimality. Earlier17/13 experiments are dominated.")
    (HERE/"integer-validation12.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt))
if __name__=="__main__":main()

