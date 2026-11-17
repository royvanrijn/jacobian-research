#!/usr/bin/env python3
"""Build the89-variable integer cubic-linear Keller map from the exact cube identity."""
from pathlib import Path
from math import lcm
import hashlib,importlib.util,json
import sympy as s
HERE=Path(__file__).resolve().parent;DB=HERE.parent/"de-bondt-2610.24889"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cube(v):return v.applyfunc(lambda a:a**3)
def main():
    folder=HERE/"power89";folder.mkdir(exist_ok=True)
    cp=HERE/"power-basis-anneal/candidate.json";candidate=json.loads(cp.read_text())
    checkpoint=json.loads((HERE/"power-basis-anneal/checkpoint.json").read_text())
    assert checkpoint["status"]=="candidate" and candidate["source_sha256"]==sha(HERE/"power_basis_anneal.py")
    for name,digest in checkpoint["source_inputs"].items():assert sha(DB/name)==digest
    L=s.Matrix(candidate["forms"]);rawB=s.Matrix([[s.Rational(c) for c in row] for row in candidate["B"]])
    n=L.rows;assert n==candidate["size"]==89 and L.cols==rawB.rows==18
    spec=importlib.util.spec_from_file_location("db89_source",DB/"verify.py");module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    _,K,base_points=module.cubic12(2);H,C=module.homogeneous18(K);Z=module.X
    padded=list(H);padded[-1]=Z[-1]**3
    for i in range(18):
        identity=sum(rawB[i,j]*(sum(L[j,k]*Z[k] for k in range(18)))**3 for j in range(n))
        assert s.Poly(s.expand(identity-padded[i]),*Z).is_zero
    scale=lcm(*(int(c.q) for c in rawB[:17,:]));assert scale
    B=rawB.copy()
    for i in range(17):
        for j in range(n):B[i,j]*=scale
    assert all(c.q==1 for c in B)
    padding=candidate["padding_index"];e=s.zeros(18,1);e[17]=1
    assert L.row(padding)==e.T and B.row(17)==s.eye(n).row(padding)
    L0=L.copy();L0[padding,:]=s.zeros(1,18);A=L0*B
    assert all(c.q==1 for c in A)
    _,pivots=B.rref();assert len(pivots)==18
    inverse=B[:,list(pivots)].inv();C0=s.zeros(n,18)
    for i,j in enumerate(pivots):C0[j,:]=inverse[i,:]
    assert B*C0==s.eye(18)
    lifted=[s.Matrix(p+[-scale*a for a in module.evaluate(C,p)]+[s.Rational(1,scale)]) for p in base_points]
    base=lambda z:z+scale*s.Matrix(module.evaluate(H,list(z)))
    assert base(lifted[0])==base(lifted[1])
    shear=B[:,padding]-e;assert shear[17]==0
    u=C0*lifted[0];trial=C0*lifted[1];image=u+cube(A*u);v=image-cube(A*trial)
    assert B*u==lifted[0] and B*v==lifted[1] and u!=v
    assert v+cube(A*v)==image
    assert B*image==base(lifted[0])-shear*lifted[0][17]**3
    payload=dict(status="exact-primary-pass",source_sha256=sha(Path(__file__)),candidate_sha256=sha(cp),
       base_source_sha256=sha(DB/"verify.py"),dimension=n,coefficient_ring="Z",scale=scale,
       forms=[[int(c) for c in row] for row in L.tolist()],
       raw_B=[[str(c) for c in row] for row in rawB.tolist()],
       B=[[int(c) for c in row] for row in B.tolist()],
       A=[[int(c) for c in row] for row in A.tolist()],
       padding_index=padding,right_inverse_columns=list(pivots),
       right_inverse_minor=[[str(c) for c in row] for row in inverse.tolist()],
       base_points=[[str(c) for c in p] for p in base_points],
       lifted_points=[[str(c) for c in p] for p in lifted],
       u=[str(c) for c in u],v=[str(c) for c in v],image=[str(c) for c in image],
       target_shear=[str(c) for c in shear],
       map="F(Y)=Y+(A*Y)^{circ3}",
       determinant_proof="Exact shared-cube identity, companion Schur determinant and Sylvester identity; no89-variable determinant expansion.",
       scope="Explicit integer cubic-linear witness in89 variables, improving the previous90-variable construction. Minimum dimension and priority are not claimed.")
    (folder/"witness.json").write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({key:payload[key] for key in ["status","dimension","scale","padding_index","scope"]}))
if __name__=="__main__":main()

