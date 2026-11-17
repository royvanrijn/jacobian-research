"""Exact finite controls for PROOF.md; not a proof of its universal statements."""
import hashlib
import json
import platform
import time
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parent
started = time.monotonic()
checks = []

def check(name, condition):
    if not condition:
        raise AssertionError(name)
    checks.append(name)

u, v, X, Y, a, b, t = s.symbols("u v X Y a b t")
pair = s.Matrix([[a, b], [b, a*u]])
check("quadratic arbitrary functional determinant",
      s.expand(pair.det() - (a*a*u-b*b)) == 0)
check("adjugate common-denominator identity",
      pair.adjugate()*pair == pair.det()*s.eye(2))
coefficient_pair = s.Matrix([[0, 1], [1, 0]])
check("coefficient functional nondegenerate even in characteristic two",
      coefficient_pair.det() == -1)
# In characteristic two, t^2=u has trace(1)=trace(t)=0,
# while coefficient-of-t sends t to 1.
mul_t = s.Matrix([[0, u], [1, 0]])
check("purely inseparable degree-two field trace vanishes",
      int(s.eye(2).trace()) % 2 == 0 and mul_t.trace() == 0)
cubic_pair = s.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
check("degree-three coefficient pairing nondegenerate",
      cubic_pair.det() == -1)
as_pair = s.Matrix([[0, 1], [1, 1]])
check("Artin-Schreier degree-two coefficient pairing",
      int(as_pair.det()) % 2 == 1)
as_map = t*t+t
check("Artin-Schreier distinct preimages in characteristic two",
      s.Poly(s.expand(as_map.subs(t, t+1)-as_map), t, modulus=2).is_zero)
check("Artin-Schreier etale derivative",
      s.Poly(s.diff(as_map,t), t, modulus=2).as_expr() == 1)

chart = {X: 1/u, Y: u**3*v-u**2}
generators = [Y, X*Y, X**2*Y, X**3*Y+X]
expected = [u**3*v-u**2, u**2*v-u, u*v-1, v]
for index, (g, value) in enumerate(zip(generators, expected),1):
    check(f"Wright generator chart {index}",
          s.cancel(g.subs(chart)-value) == 0)
jacobian = s.det(s.Matrix([[s.diff(chart[X],u),s.diff(chart[X],v)],
                           [s.diff(chart[Y],u),s.diff(chart[Y],v)]]))
check("Wright volume-form chart", s.cancel(jacobian+u)==0)

def order_u(expr):
    numerator, denominator = s.fraction(s.cancel(expr))
    def minimum(poly):
        return min(monomial[0] for monomial, coefficient
                   in s.Poly(poly,u).terms() if coefficient != 0)
    return minimum(numerator)-minimum(denominator)

check("Wright missing source coordinate pole", order_u(chart[X]) == -1)
for i, g in enumerate(expected,1):
    check(f"Wright generator nonnegative valuation {i}", order_u(g) >= 0)
for power in range(1,7):
    # These are controls, not a substitute for the proof for arbitrary c.
    check(f"multiplier Y^{power} fails on X^{2*power+1}",
          order_u((Y**power*X**(2*power+1)).subs(chart)) == -1)
for degree in range(8):
    value = Y*s.diff(X**(degree+3),X,degree)
    check(f"finite derivative pole control order {degree}",
          order_u(value.subs(chart)) == -1)

E = lambda f: s.expand(f).subs(X,0)
bxy = X*Y
operator = E
for depth in range(9):
    check(f"nonlocal evaluation commutator depth {depth}",
          s.expand(operator(s.Integer(1))-(-bxy)**depth) == 0)
    previous = operator
    operator = lambda f, previous=previous: s.expand(
        previous(bxy*f)-bxy*previous(f))
check("nonlocal evaluation fails B-linearity",
      E(bxy) != bxy*E(s.Integer(1)))
sample = X**3+X*Y+Y**2+1
origin = lambda f: s.expand(f).subs({X:0,Y:0})
for i, g in enumerate(generators,1):
    check(f"torsion residue-module linearity control {i}",
          origin(g*sample) == origin(g)*origin(sample))
check("transcendental-extension evaluation functional",
      (t*(X**2+1)).subs(X,0) == t)

manifest = ROOT/"sources"/"manifest.json"
receipt = {
    "schema_version": 1,
    "result": "pass",
    "assurance": "Exact representative identities only. Universal finiteness, finite-target and all-order obstruction statements have written proofs in PROOF.md; this is not formal verification or independent whole-proof review.",
    "checks": checks,
    "check_count": len(checks),
    "python": platform.python_version(),
    "sympy": s.__version__,
    "runtime_seconds": round(time.monotonic()-started,3),
    "sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
               for p in [ROOT/"PROOF.md", Path(__file__), manifest]},
    "closes_JC2": False,
    "closes_HC4": False
}
(ROOT/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"result":receipt["result"],"checks":len(checks),
                  "runtime_seconds":receipt["runtime_seconds"]}))
