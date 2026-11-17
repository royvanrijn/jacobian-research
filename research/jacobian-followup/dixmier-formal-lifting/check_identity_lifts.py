#!/usr/bin/env python3
"""QNI1 exact full-family control; the general implications use written proofs."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sympy as sp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
INPUT=ROOT/"archive/non-elliptic/scripts/verify_counterexample_independent.py"
spec=importlib.util.spec_from_file_location("retained_foundational_input",INPUT)
retained=importlib.util.module_from_spec(spec)
spec.loader.exec_module(retained)
x,y,z,t=sp.symbols("x y z t")
p1,p2,p3=sp.symbols("p1 p2 p3")
positions=(x,y,z)
momenta=(p1,p2,p3)

def zero(F):
    assert sp.expand(F)==0,sp.expand(F)

def rational_zero(F):
    assert sp.cancel(F)==0,sp.cancel(F)

def convert(poly):
    return sp.Add(*(sp.Integer(c)*sp.prod(v**e for v,e in zip(positions,m))
                    for m,c in poly.items()))

def degree(F):
    return int(sp.Poly(F,*positions,*momenta).total_degree())

def poisson(F,G):
    return sp.expand(sum(sp.diff(F,p)*sp.diff(G,q)-sp.diff(F,q)*sp.diff(G,p)
                         for q,p in zip(positions,momenta)))

def main():
    raw=[convert(P) for P in retained.mapping]
    f=sp.Matrix([raw[2]/2,raw[1],raw[0]])
    origin=dict.fromkeys(positions,0)
    J=f.jacobian(positions)
    assert f.subs(origin)==sp.zeros(3,1)
    assert J.subs(origin)==sp.eye(3)
    zero(J.det()-1)
    A=J.adjugate().applyfunc(sp.expand)
    for v in J*A-sp.eye(3):
        zero(v)
    for v in A*J-sp.eye(3):
        zero(v)
    scale={v:t*v for v in positions}
    ft=(f.subs(scale,simultaneous=True)/t).applyfunc(sp.expand)
    At=A.subs(scale,simultaneous=True).applyfunc(sp.expand)
    assert ft.subs(t,0)==sp.Matrix(positions)
    assert At.subs(t,0)==sp.eye(3)
    assert ft.jacobian(positions)==J.subs(scale,simultaneous=True)
    for F in list(ft)+list(At):
        sp.Poly(F,*positions,t) # no negative parameter powers
    for v in ft.jacobian(positions)*At-sp.eye(3):
        zero(v)
    P=At.T*sp.Matrix(momenta)
    relations=0
    # Each P is first-order. These Poisson identities are the exact
    # normalized normal-order Weyl commutators, without higher terms.
    for i in range(3):
        for j in range(3):
            zero(poisson(P[i],ft[j])-int(i==j))
            relations+=1
    for i in range(3):
        for j in range(i+1,3):
            zero(poisson(P[i],P[j]))
            zero(poisson(ft[i],ft[j]))
            relations+=2
    assert relations==15
    degrees=[degree(a) for a in list(ft)+list(P)]
    assert degrees==[4,6,7,12,10,9]
    for point in retained.points:
        scaled={v:sp.Rational(c.numerator,c.denominator)/t
                for v,c in zip(positions,point)}
        image=ft.subs(scaled,simultaneous=True)
        for a,b in zip(image,[0,0,-1/(4*t)]):
            rational_zero(a-b)
    leading=(ft-sp.Matrix(positions)).diff(t).subs(t,0)
    assert leading==sp.Matrix([0,3*x*z,4*y**2])
    H=(leading.T*sp.Matrix(momenta))[0]
    quantum_leading=(P-sp.Matrix(momenta)).diff(t).subs(t,0)
    for i in range(3):
        zero(poisson(H,positions[i])-leading[i])
        zero(poisson(H,momenta[i])-quantum_leading[i])

    receipt={
        "claim":"QNI1","status":"PASS",
        "normalized_tangent_identity_and_jacobian":True,
        "inverse_Jacobian_column_degrees":
            [max(int(sp.Poly(A[j,i],*positions).total_degree()) for j in range(3))
             for i in range(3)],
        "full_relations_as_polynomial_t_identities":relations,
        "forward_generator_spatial_degrees":degrees,
        "uniform_forward_bound":12,
        "scaled_rational_collision_points":len(retained.points),
        "first_Hamiltonian":str(H),
        "parameter_scope":"Relations hold for independent t and Weyl h. Specializing t=h^N gives the family for every N>=1; N>=2 also ensures h-adic unit innerness.",
        "coverage_boundary":"Exact polynomial-family identities and the retained collision are checked. Generic nonsurjectivity, all-degree inverse nontermination, formal graph equivalence, failure of universal bounded comparison and DC_n=>P_n=>JC_n use written proofs. The extension is not independently replayed or formally verified; no rank-two counterexample is produced.",
        "independent_replay":False,"formal_verification":False,"external_review":False,
        "hashes":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                  for name in ("IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md",
                               "check_identity_lifts.py")},
        "retained_input":{
            "claim":"F1","path":str(INPUT.relative_to(ROOT)),
            "sha256":hashlib.sha256(INPUT.read_bytes()).hexdigest(),
            "scope":"Imports the retained exact map and rational collision; importing replays its small sparse determinant control. This extension does not inherit F1's formal-verification status."}
    }
    (HERE/"identity_lift_controls.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()
