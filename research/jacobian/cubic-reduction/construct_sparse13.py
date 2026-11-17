"""Independently construct and check the 42-term cubic thirteen-variable witness."""
from pathlib import Path
import hashlib,json
import sympy as s
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    base=json.loads((HERE/"cubic13-witness.json").read_text())
    candidate=json.loads((HERE/"cubic13-sparse-candidate.json").read_text())
    assert candidate["source_sha256"]==sha(HERE/"optimize13.py")
    assert candidate["input_sha256"]==sha(HERE/"cubic13-witness.json")
    X=s.symbols(base["variables"]);n=len(X)
    F=s.Matrix([sum(s.Rational(t["coefficient"])*s.prod(v**e for v,e in zip(X,t["exponent"])) for t in row) for row in base["coefficients"]])
    points=[s.Matrix([s.Rational(c) for c in p]) for p in base["points"]]
    for move in candidate["moves"]:
        i,j=move["i"],move["j"];t=s.Rational(move["t"]);assert i!=j
        S=s.eye(n);S[i,j]=t
        F=(S.inv()*F.subs(dict(zip(X,S*s.Matrix(X))),simultaneous=True)).applyfunc(s.expand)
        points=[S.inv()*p for p in points]
    expected=s.Matrix([sum(s.Rational(t["coefficient"])*s.prod(v**e for v,e in zip(X,t["exponent"])) for t in row) for row in candidate["coefficients"]])
    assert (F-expected).applyfunc(s.expand)==s.zeros(n,1)
    J=F.jacobian(X);assert J.subs(dict.fromkeys(X,0))==s.eye(n)
    A=J[:3,:3];B=J[:3,3:];C=J[3:,:3];D=J[3:,3:];N=D-s.eye(n-3)
    assert (N**3).applyfunc(s.expand)==s.zeros(n-3)
    inv=(s.eye(n-3)-N+N*N).applyfunc(s.expand)
    assert (D*inv).applyfunc(s.expand)==s.eye(n-3)
    Schur=(A-B*inv*C).applyfunc(s.expand)
    assert s.expand(Schur.det(method="berkowitz"))==1
    images=[[s.expand(q.subs(dict(zip(X,p)))) for q in F] for p in points]
    assert images[0]==images[1]==images[2]
    coefficients=[[dict(exponent=list(e),coefficient=str(c)) for e,c in s.Poly(q,*X).terms()] for q in F]
    counts=[sum(sum(t["exponent"])>=2 for t in row) for row in coefficients]
    assert sum(counts)==42
    d=dict(status="exact-primary-pass",source_sha256=sha(Path(__file__)),
      base_sha256=sha(HERE/"cubic13-witness.json"),candidate_sha256=sha(HERE/"cubic13-sparse-candidate.json"),
      variables=list(map(str,X)),dimension=n,degree=max(s.Poly(q,*X).total_degree() for q in F),
      nonlinear_counts=counts,nonlinear_total=sum(counts),coefficients=coefficients,
      points=[[str(c) for c in p] for p in points],common_image=list(map(str,images[0])),
      determinant="1",identity_derivative=True,carrier_nilpotence_exponent=3,
      conjugations=candidate["moves"],scope="Explicit degree-three dimension13 existence witness with42 displayed nonlinear occurrences; no optimality claim.")
    (HERE/"cubic13-sparse-witness.json").write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in d.items() if k not in ("coefficients","points")}))
if __name__=="__main__":main()

