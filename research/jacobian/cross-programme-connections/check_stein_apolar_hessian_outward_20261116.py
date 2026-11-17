"""Finite controls for the outward-transfer note; the written arguments are universal."""
from pathlib import Path
import hashlib
import json
import sympy as s

root = Path(__file__).resolve().parent
x, y, t = s.symbols("x y t")
def gaussian_moment(n):
    return 0 if n % 2 else s.factorial2(n - 1)

stein_checks = 0
for n in range(15):
    mu = gaussian_moment(3*n)
    expected = (3*n+1)*(3*n+3)*(3*n+5)*mu
    assert gaussian_moment(3*n+6) == expected
    f = y**n
    Sf = 27*y**3*s.diff(f,y,3)+162*y**2*s.diff(f,y,2)+177*y*s.diff(f,y)+(15-y*y)*f
    value = sum(c*gaussian_moment(3*e[0]) for e,c in s.Poly(s.expand(Sf),y).terms())
    assert value == 0
    stein_checks += 2

apolar_checks = 0
P = x**2+y**3
nilpotence = []
for j in range(7):
    poly = s.expand(P**j)
    nonzero_powers = []
    for k in range(9):
        value = s.diff(poly, x, k, y, k)
        expected = any(2*a >= k and 3*(j-a) >= k for a in range(j+1))
        assert (value != 0) == expected
        if value != 0:
            nonzero_powers.append(k)
        apolar_checks += 1
    nu = 1+max(nonzero_powers)
    formula = 1+max(min(2*a,3*(j-a)) for a in range(j+1))
    assert nu == formula
    nilpotence.append(nu)
for j in range(7):
    inverted = 1+max([0]+[k for k in range(1,9)
                         if s.ceiling(s.Rational(k,2))+s.ceiling(s.Rational(k,3)) <= j])
    assert max(nilpotence[:j+1]) == inverted
    apolar_checks += 1
# Uniformly nilpotent branch: d_y annihilates every power of P=x.
assert all(s.diff(x**j,y) == 0 for j in range(7))

xs = s.symbols("a b c d")
ps = s.symbols("u v w z")
h = xs[0]*xs[2]+xs[1]*xs[3]+xs[0]**4
inverse = [ps[2],ps[3],ps[0]-4*ps[2]**3,ps[1]]
dual = s.expand(sum(u*v for u,v in zip(ps,inverse))-h.subs(dict(zip(xs,inverse))))
assert dual == ps[0]*ps[2]+ps[1]*ps[3]-ps[2]**4
assert [s.diff(dual,p) for p in ps] == inverse
assert s.hessian(h,xs).det() == s.hessian(dual,ps).det() == 1

note = root / "STEIN_APOLAR_AND_HESSIAN_OUTWARD_BRIDGES_2026_11_16.tex"
dependencies = [
    root.parent/"moment-vanishing-audit/BOUNDED_BERNSTEIN_STEPS_AND_RATIONAL_SURVIVAL_SLOPES.md",
    root.parent/"../hessian/hc4-triangular-graphs/RICCI_FLAT_HC4_AND_ONE_DIMENSIONAL_JACOBIAN_IMAGES.md",
]
receipt = {
    "date":"2026-11-16",
    "status":"finite_exact_controls_pass",
    "stein_checks":stein_checks,
    "apolar_checks":apolar_checks,
    "apolar_nilpotence_indices_example":nilpotence,
    "legendre_control":"nonlinear polynomial dual, inverse gradient and reciprocal Hessian determinant",
    "proof_sha256":hashlib.sha256(note.read_bytes()).hexdigest(),
    "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "dependency_sha256":{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies},
    "scope":"All-degree corollaries rest on written arguments and retained input theorems, not these finite examples.",
    "formal_verification":False,"external_review":False,"registry_promotion":False,
    "novelty_claim":False,"agent_messages_sent":0
}
out=root/"stein-apolar-hessian-outward-controls-2026-11-16.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(receipt,indent=2))
