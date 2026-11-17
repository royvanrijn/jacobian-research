"""Construct the 36-term conjugate of the verified de Bondt twelve-variable map."""
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


import hashlib,importlib.util,json
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    base_path=HERE/"de-bondt12-witness.json";base=json.loads(base_path.read_text())
    spec=importlib.util.spec_from_file_location("de_bondt_primary",source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    _,F,base_points=module.cubic12();X=module.X[:12];n=12;F=s.Matrix(F)
    expected=s.Matrix([sum(s.Rational(t["coefficient"])*s.prod(v**e for v,e in zip(X,t["exponent"])) for t in row) for row in base["coefficients"]])
    assert (F-expected).applyfunc(s.expand)==s.zeros(n,1)
    candidate_path=HERE/"de-bondt12-sparse-candidate.json";candidate=json.loads(candidate_path.read_text())
    points=[s.Matrix(point) for point in base_points]
    for move in candidate["moves"]:
        i,j=move["i"],move["j"];t=s.Rational(move["t"]);assert i!=j
        S=s.eye(n);S[i,j]=t
        F=(S.inv()*F.subs(dict(zip(X,S*s.Matrix(X))),simultaneous=True)).applyfunc(s.expand)
        points=[S.inv()*p for p in points]
    expected=s.Matrix([sum(s.Rational(t["coefficient"])*s.prod(v**e for v,e in zip(X,t["exponent"])) for t in row) for row in candidate["coefficients"]])
    assert (F-expected).applyfunc(s.expand)==s.zeros(n,1)
    J=F.jacobian(X);assert J.subs(dict.fromkeys(X,0))==s.eye(n)
    A=J[:3,:3];B=J[:3,3:];C=J[3:,:3];D=J[3:,3:];N=D-s.eye(9)
    power=s.eye(9);inverse=power;nilpotence=None
    for exponent in range(1,10):
        power=(power*(-N)).applyfunc(s.expand)
        if power==s.zeros(9):nilpotence=exponent;break
        inverse=(inverse+power).applyfunc(s.expand)
    assert nilpotence and (D*inverse).applyfunc(s.expand)==s.eye(9)
    Schur=(A-B*inverse*C).applyfunc(s.expand)
    assert s.expand(Schur.det(method="berkowitz"))==1
    images=[[s.expand(f.subs(dict(zip(X,p)))) for f in F] for p in points]
    assert images[0]==images[1]
    coefficients=[[dict(exponent=list(e),coefficient=str(c)) for e,c in s.Poly(f,*X).terms()] for f in F]
    counts=[sum(sum(t["exponent"])>=2 for t in row) for row in coefficients];assert sum(counts)==36
    d=dict(status="exact-primary-pass",source_sha256=sha(Path(__file__)),base_sha256=sha(base_path),
      base_source_sha256=sha(source),candidate_sha256=sha(candidate_path),
      variables=list(map(str,X)),dimension=n,degree=max(s.Poly(f,*X).total_degree() for f in F),
      nonlinear_counts=counts,nonlinear_total=sum(counts),coefficients=coefficients,
      points=[[str(v) for v in p] for p in points],common_image=list(map(str,images[0])),
      determinant="1",identity_derivative=True,carrier_nilpotence_exponent=nilpotence,
      conjugations=candidate["moves"],scope="36-term conjugate of the already verified12-variable cubic BCRDB1 map. Dimension unchanged; no optimality or novelty claim.")
    (_check_output(HERE/"cubic12-sparse-witness.json")).write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in d.items() if k not in ("coefficients","points")}))
if __name__=="__main__":main()

