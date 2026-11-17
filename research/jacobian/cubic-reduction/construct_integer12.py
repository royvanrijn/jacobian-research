"""Clear the two residual halves by a diagonal conjugation, preserving36 terms."""
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


import hashlib,json
import sympy as s
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    path=HERE/"cubic12-sparse-witness.json";base=json.loads(path.read_text());X=s.symbols(base["variables"]);n=12
    F=s.Matrix([sum(s.Rational(t["coefficient"])*s.prod(v**e for v,e in zip(X,t["exponent"])) for t in row) for row in base["coefficients"]])
    diagonal=[s.Rational(1)]*n;diagonal[8]=diagonal[10]=s.Rational(1,2);D=s.diag(*diagonal)
    F=(D.inv()*F.subs(dict(zip(X,D*s.Matrix(X))),simultaneous=True)).applyfunc(s.expand)
    points=[D.inv()*s.Matrix([s.Rational(v) for v in p]) for p in base["points"]]
    coefficients=[[dict(exponent=list(e),coefficient=str(c)) for e,c in s.Poly(f,*X).terms()] for f in F]
    assert all(s.Rational(t["coefficient"]).q==1 for row in coefficients for t in row)
    assert all(v.q==1 for p in points for v in p)
    counts=[sum(sum(t["exponent"])>=2 for t in row) for row in coefficients];assert sum(counts)==36
    assert F.jacobian(X).subs(dict.fromkeys(X,0))==s.eye(n)
    images=[[s.expand(f.subs(dict(zip(X,p)))) for f in F] for p in points];assert images[0]==images[1]
    d=dict(status="exact-primary-pass",source_sha256=sha(Path(__file__)),rational_witness_sha256=sha(path),
      dimension=n,degree=max(s.Poly(f,*X).total_degree() for f in F),variables=list(map(str,X)),coefficients=coefficients,
      nonlinear_counts=counts,nonlinear_total=36,diagonal=list(map(str,diagonal)),
      points=[[str(v) for v in p] for p in points],common_image=list(map(str,images[0])),
      integral_coefficients=True,integral_collision_points=True,determinant="1",identity_derivative=True,
      determinant_justification="Invertible constant conjugation of the independently checked rational36-term map.",
      scope="Integer36-term conjugate of BCRDB1's already verified12-variable cubic. Dimension unchanged, no optimality claim.")
    (_check_output(HERE/"cubic12-integer-witness.json")).write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in d.items() if k not in ("coefficients","points")}))
if __name__=="__main__":main()

