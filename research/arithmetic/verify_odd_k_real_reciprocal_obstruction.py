#!/usr/bin/env python3
"""Exact finite controls for the uniform odd-k real reciprocal obstruction.

The companion note, not the bounded controls, proves the all-odd-k theorem.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import platform

import sympy as sp

HERE = Path(__file__).resolve().parent
v = sp.symbols("v", positive=True)
w = sp.symbols("w")


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def reconstruct(k, sigma):
    sigma = sp.Poly(sp.expand(sigma), v)
    require(sigma.degree() == 2*k and sigma.LC() == 1, "Wrong critical degree or leading coefficient")
    L = 4*k-1
    return sp.expand(sum(L*sigma.nth(j)*v**j/sp.Integer(2*j-1) for j in range(2*k+1)))


def symbolic_control(k):
    coeffs = sp.symbols("a0:"+str(k))
    s = w**k + sum(coeffs[j]*w**j for j in range(k))
    sigma = sp.expand(v**k*s.subs(w, v+1/v))
    F = reconstruct(k, sigma)
    L = 4*k-1
    require(sp.expand(2*v*sp.diff(F,v)-F-L*sigma) == 0, "Euler critical equation")
    require(F.subs(v,0) == -L, "Negative constant coefficient")
    require(sp.expand(v**(2*k)*sigma.subs(v,1/v)-sigma) == 0, "Critical reciprocity")
    A = (F+v*F.subs(v,1/v))/sp.sqrt(v)
    B = (F-v*F.subs(v,1/v))/sp.sqrt(v)
    scale = 2*v**(2*k)*sp.sqrt(v)
    require(sp.expand(scale*sp.diff(A,v)-L*sigma*(v**(2*k-1)-1)) == 0,
            "A derivative with positive weight")
    require(sp.expand(scale*sp.diff(B,v)-L*sigma*(v**(2*k-1)+1)) == 0,
            "B derivative with positive weight")
    G = (F+L)/sp.sqrt(v)
    require(sp.expand(2*v*sp.sqrt(v)*sp.diff(G,v)-L*(sigma-1)) == 0,
            "Integral antiderivative identity")
    require(sp.limit(G,v,0,dir="+") == 0, "Integral boundary constant")
    return {"k":k, "degree":4*k, "symbolic_coefficients":k,
            "Euler_reciprocity_derivative_and_endpoint_identities":True}


def endpoint_control(k):
    roots = [sp.Rational(1,j+2) for j in range(k)]
    sigma = sp.prod((v-a)*(v-1/a) for a in roots)
    F = reconstruct(k, sigma)
    alpha = min(roots)
    value = sp.cancel(F.subs(v,alpha))
    require(value < 0, "First critical value sign")
    require(sp.expand(sigma.subs(v,0)) == 1, "Positive reciprocal constant")
    require(sp.expand(sigma-sp.prod((1-v/a)*(1-v*a) for a in roots)) == 0,
            "Product normalization used by the strict integral inequality")
    return {"k":k, "smallest_critical_root":str(alpha),
            "F_at_smallest_critical_root":str(value), "strictly_negative":True}


def rolle_compatible(word):
    p = [i for i,x in enumerate(word) if x == "P"]
    q = [i for i,x in enumerate(word) if x == "Q"]
    # Consecutive P zeros require a Q derivative zero between them.
    if any(not any(a<j<b for j in q) for a,b in zip(p,p[1:])):
        return False
    # B(0)=0 requires a P before the first Q, and between successive Q.
    if q and not any(i<q[0] for i in p):
        return False
    if any(not any(a<i<b for i in p) for a,b in zip(q,q[1:])):
        return False
    return True


def order_controls():
    controls = []
    for k in range(1,14):
        compatible = ["".join(word) for word in product("PQ",repeat=k) if rolle_compatible(word)]
        expected = "".join("P" if i%2 == 0 else "Q" for i in range(k))
        require(compatible == [expected], "Rolle conditions did not force alternation")
        controls.append({"k":k, "unique_Rolle_compatible_word":expected,
                         "last_sign":"P" if k%2 else "Q"})
    return controls


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt",type=Path)
    args = parser.parse_args()
    require(sp.__version__ == "1.14.0", "Expected exact Sympy version")
    report = {
        "schema_version":1, "result":"PASS",
        "python":platform.python_version(), "sympy":sp.__version__,
        "checker_sha256":sha256(Path(__file__).read_bytes()).hexdigest(),
        "symbolic_controls":[symbolic_control(k) for k in (1,2,3,5,7)],
        "exact_endpoint_controls":[endpoint_control(k) for k in (1,2,3,5,7)],
        "Rolle_order_controls":order_controls(),
        "proof_boundary":(
            "The uniform proof uses the strict product inequality on the entire first-root interval, "
            "the convergent antiderivative, Rolle counts for arbitrary k, and the last-gap sign. "
            "Bounded symbolic/formula and sign-order checks supplement that proof. "
            "No existence is certified for even k or for nonreciprocal models."
        ),
        "formal_verification":False, "independent_replay":False, "external_review":False,
    }
    text = json.dumps(report,indent=2)+"\n"
    if args.write_receipt:
        args.write_receipt.write_text(text,encoding="utf-8")
    print(text,end="")


if __name__ == "__main__":
    main()
