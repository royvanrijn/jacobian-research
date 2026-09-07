#!/usr/bin/env python3
"""Check the exact local identities used by the HC4 affine-leaf obstruction.

The universal polynomial-degree argument and the passage to affine leaves
are written proofs. This checker reuses the published frame constructor,
verifies literal curvature combinations, and checks the rational ODE identity.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "scripts/verify_hc4_motion_frame_transport.py"
OUTPUT = ROOT / "artifacts/generated-results/hc4-negative-motion-polynomial-obstruction-v1.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location("hc4_transport_audit", AUDIT)
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    ns = audit.load_frame_constructor()
    a, p, q = ns["a"], ns["p0_geometric"], ns["q"]
    parameters = ns["parameters"]
    derivative = ns["derivative_of"]
    h = parameters[12]+parameters[3]+parameters[7]
    relations = [p+a, q+a, h]
    equations = audit.branch_equations(ns, relations)
    substitution = {p: -a, q: -a, parameters[3]: -parameters[12]-parameters[7]}

    certificates = {
        "e1_a": {0: sp.Rational(1, 2), 5: -1, 10: sp.Rational(-1, 2)},
        "e2_a": {9: -1, 14: -1, 16: -1, 26: -1},
    }
    targets = {
        "e1_a": derivative(a, 0),
        "e2_a": derivative(a, 1)-sp.Rational(3, 2)*a**2,
    }
    for name, coefficients in certificates.items():
        combination = sum(c*equations[i] for i, c in coefficients.items())
        assert sp.expand(combination-targets[name]) == 0

    connection = [m.subs(substitution, simultaneous=True) for m in ns["connection"]]
    assert connection[0][0, :] == sp.zeros(1, 4)
    assert connection[0][1, :] == sp.Matrix([[-a, 0, 0, 0]])
    assert connection[1][0, :] == sp.Matrix([[sp.Rational(5, 2)*a, 0, 0, 0]])
    assert connection[1][1, 1] == a/2
    assert connection[1][1, 2:] == sp.zeros(1, 2)

    # On an affine leaf, e1=f*ds and e2=h*ds+v*dt, with f_s=v_s=0.
    # These are the exact t-derivatives forced by the checked identities.
    f, v, alpha = sp.symbols("f v alpha", nonzero=True)
    rules = {f: 5*alpha*f/(2*v), v: alpha/2, alpha: 3*alpha**2/(2*v)}

    def Dt(expression):
        return sp.factor(sum(sp.diff(expression, variable)*value
                             for variable, value in rules.items()))

    n = f/v
    dn = Dt(n)
    assert sp.factor(dn-2*alpha*f/v**2) == 0
    ode = sp.factor(2*n*Dt(dn)-3*dn**2)
    assert ode == 0
    degree = sp.symbols("d", integer=True)
    leading_factor = sp.factor(2*degree*(degree-1)-3*degree**2)
    assert leading_factor == -degree*(degree+2)

    # Normalized Jordan top coefficient matches the Hessian cofactor line.
    S, N = ns["S"], ns["N"]
    first = sp.eye(4)[:, 0]
    assert N**3*S.inv() == first*first.T
    assert ns["T"].adjugate() == -first*first.T

    paths = [
        Path(__file__).resolve(), AUDIT,
        ROOT / "scripts/verify_hc4_affine_plane_prolongation.py",
    ]
    data = {
        "format": "hc4-negative-motion-polynomial-obstruction-v1",
        "theorem": "HC4MRA2",
        "source_sha256": {
            path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths
        },
        "branch": "p=q=-a, a!=0, p12+p3+p7=0",
        "curvature_certificates": {
            name: {"target": str(targets[name]),
                   "equation_coefficients": {str(i): str(c) for i, c in row.items()}}
            for name, row in certificates.items()
        },
        "affine_leaf_connection": {
            "nabla_e1_e1": [str(x) for x in connection[0][0, :]],
            "nabla_e1_e2": [str(x) for x in connection[0][1, :]],
            "nabla_e2_e1": [str(x) for x in connection[1][0, :]],
            "nabla_e2_e2": [str(x) for x in connection[1][1, :]],
        },
        "affine_coordinate_derivatives": {str(k): str(val) for k, val in rules.items()},
        "polynomial_matrix_coefficient": str(n),
        "coefficient_derivative": str(dn),
        "ode_residual": str(ode),
        "positive_degree_leading_factor": str(leading_factor),
        "boundary": (
            "Written global proof: N is polynomial, so its coefficient n on an affine "
            "leaf is polynomial. The ODE forces n constant in characteristic zero, "
            "contradicting a!=0. No finite-jet inconsistency or HC4 counterexample is claimed."
        ),
        "normalization_boundary": (
            "The written Piola/quasi-translation argument makes the cofactor-normalized "
            "e1 affinely parallel along itself; the matrix identities check its algebraic scale."
        ),
    }
    serialized = json.dumps(data, indent=2)+"\n"
    if args.write:
        OUTPUT.write_text(serialized)
    else:
        assert OUTPUT.read_text() == serialized, "artifact differs; inspect before --write"
    print("PASS literal curvature certificates e1(a)=0 and e2(a)=3*a^2/2")
    print("PASS affine-leaf connection and cofactor-normalization matrices")
    print("PASS 2*n*n''-3*(n')^2=0 and leading factor -d*(d+2)")
    print("PASS " + ("wrote " if args.write else "byte-identical ") + str(OUTPUT))
    print("Global polynomiality and affine-leaf restriction remain explicit written proof steps")


if __name__ == "__main__":
    main()
