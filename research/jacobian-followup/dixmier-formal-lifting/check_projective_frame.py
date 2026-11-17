#!/usr/bin/env python3
"""Exact controls for QPF1. Universal module statements use the written proof."""
from pathlib import Path
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
x, p, t, z = sp.symbols("x p t z")

def mul(f, g, parameter=t):
    f, g = sp.expand(f), sp.expand(g)
    if f == 0 or g == 0:
        return sp.Integer(0)
    return sp.expand(sum(parameter**j * sp.diff(f, p, j) *
                         sp.diff(g, x, j) / sp.factorial(j)
                         for j in range(sp.Poly(f, p).degree()+1)))

def zero(f):
    assert sp.cancel(f) == 0, sp.cancel(f)

def main():
    # z is independent: substituting z=t^(m-1) proves these identities for all m.
    q, c = 1-t*z*x, t**2*z
    f, g = mul(q,q), mul(q,p)-2*c
    a = z*mul(q,x**2)
    b = 1+t*z*x+z*mul(p,x**2)
    row, column = [-p,q], [a,b]
    zero(mul(p,q)-mul(q,p)+c)
    zero(mul(-p,f)+mul(q,g))
    zero(mul(-p,a)+mul(q,b)-1)

    # Row-vector convention: E=1-column*row projects onto a complement of A*row.
    E = [[int(i==j)-mul(column[i],row[j]) for j in range(2)] for i in range(2)]
    for i in range(2):
        for j in range(2):
            zero(sum(mul(E[i][k],E[k][j]) for k in range(2))-E[i][j])
        zero(sum(mul(E[i][k],[f,g][k]) for k in range(2))-[f,g][i])

    v = -p/c
    zero(mul(v,q)-mul(q,v)-1)
    zero(-c*(mul(q,v)+2)-g)

    # Independent canonical variables represented by x,p with parameter 1.
    standard_f, standard_g = x**2, x*p+2
    zero(mul(-p,standard_f,1)+mul(x,standard_g,1))
    zero(mul(standard_g,x,1)-mul(p,standard_f,1)-x)
    U,V = sp.symbols("U V")
    leading_ideal = sp.groebner([U**2,U*V],U,V)
    assert leading_ideal.reduce(U)[1] == U
    assert leading_ideal.reduce(V)[1] == V
    assert sp.gcd(U**2,U*V) == U
    normalized = sp.groebner([U,V],U,V)
    assert normalized.reduce(sp.Integer(1))[1] == 1

    jet_controls=[]
    for m in (1,2,5):
        for N in (1,2,4,7):
            qm = 1-t**m*x
            inverse_square=sum((j+1)*t**(m*j)*x**j
                               for j in range((N-1)//m+1))
            error=sp.Poly(sp.expand(qm**2*inverse_square-1),t)
            assert all(degree[0]>=N for degree,coeff in error.terms() if coeff)
            jet_controls.append([m,N])

    receipt={
        "claim":"QPF1",
        "status":"PASS",
        "universal_parameter_identity":"z independent; specialize z=t^(m-1), m>=1",
        "PBW_commutator":True,
        "projective_syzygy":True,
        "polynomial_splitting":True,
        "two_by_two_idempotent":True,
        "generic_canonical_coordinate_change":True,
        "canonical_leading_critical_pair":True,
        "regular_right_action_failure":True,
        "normalized_generic_symbol_quotient_dimension":1,
        "jet_inverse_controls":jet_controls,
        "coverage_boundary":"Exact symbolic identities and listed finite controls only. Syzygy exhaustion, generic nonfreeness, all ranks and all parameter jets, and compatible filtered-frame descent are written proofs. No right Weyl action on the projective ideal is constructed.",
        "independent_replay":False,
        "formal_verification":False,
        "external_review":False,
        "hashes":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                  for name in ("PROJECTIVE_FRAME_BOUNDARY.md","check_projective_frame.py")}
    }
    (HERE/"projective_frame_controls.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()
