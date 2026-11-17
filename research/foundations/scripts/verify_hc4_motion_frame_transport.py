#!/usr/bin/env python3
"""Audit the HC4 determinant transport and prolong the surviving sign branches.

The transport is checked directly by exact matrix algebra. The local branch
calculation reuses only the published first-order frame/curvature constructor;
it does not import the extra d(pq)=0 assumption or run its output writer.
A compatible finite jet is not a polynomial or formal all-order HC4 solution.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/verify_hc4_affine_plane_prolongation.py"
OUTPUT = ROOT / "artifacts/generated-results/hc4-motion-frame-transport-v1.json"


def matrix_rows(matrix):
    return [[str(sp.factor(value)) for value in row] for row in matrix.tolist()]


def check_transport():
    a, b, p, q, h = sp.symbols("a b p q h", nonzero=True)
    metric = sp.zeros(4)
    for i in range(4):
        metric[i, 3-i] = 1
    jordan = sp.zeros(4)
    for i in range(3):
        jordan[i, i+1] = 1
    tensor = metric*jordan
    # Freeze the adapted frame and take z=e1, q=(e2,e3,e4).
    passive = tensor[1:, 1:]
    M = passive.inv()
    kernel_derivative = sp.Matrix([[0, a, b], [0, 0, a], [0, 0, 0]])
    third_z = -kernel_derivative.T*passive
    assert third_z == third_z.T
    B = sp.simplify(M*third_z*M)
    assert B == sp.Matrix([[-b, -a, 0], [-a, 0, 0], [0, 0, 0]])
    # Target rows theta2,theta3 modulo theta4; p-coordinates reverse q.
    U = sp.Matrix([[-q, 0, 0], [-h, -p, 0], [0, 0, 0]])
    assert U[:2, :2].det() == p*q

    # Source q=C*q', z=c*z'; determinant is one and actual u=e3 is fixed.
    C = sp.Matrix([[1, b/(2*a), 0], [0, 1, 0], [0, 0, -a]])
    c = -1/a
    D = sp.diag(c, C)
    assert sp.factor(D.det()) == 1
    u = sp.Matrix([0, 0, 1])
    assert sp.simplify(c*C.T*u) == u
    B_normal = sp.simplify(c*C.inv()*B*C.inv().T)
    assert B_normal == sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    U_normal = sp.simplify(c*C.T*U*C.inv().T)
    assert sp.factor(U_normal[:2, :2].det()) == p*q/a**2
    # Already b=0 gives a metric- and volume-preserving counterexample
    # to the assertion that the smaller determinant has a constant factor.
    D0 = D.subs(b, 0)
    assert sp.simplify(D0.T*metric*D0) == metric
    assert sp.simplify(D0.inv()*jordan*D0) != jordan
    return {
        "adapted_graph_hessian": matrix_rows(B),
        "adapted_projective_derivative": matrix_rows(U),
        "determinant_one_source_change": matrix_rows(D),
        "normalized_graph_hessian": matrix_rows(B_normal),
        "normalized_projective_derivative": matrix_rows(U_normal),
        "normalized_projective_determinant": str(p*q/a**2),
        "b_zero_preserves_metric": True,
        "change_does_not_preserve_normalized_jordan_matrix": True,
        "conclusion": "Frozen normalization controls pq/a^2, not pq; d(pq)=0 does not follow.",
    }


def load_frame_constructor():
    source = SOURCE.read_text()
    start = source.index("n = 4\n")
    stop = source.index(
        "derivative_matrix, compatibility = eliminate_derivatives(curvature_equations)"
    )
    namespace = {"sp": sp}
    exec(compile(source[start:stop], str(SOURCE), "exec"), namespace)
    return namespace


def branch_equations(ns, relations):
    parameters = ns["parameters"]
    # Keep the geometric a coordinate, eliminating p, q, and then p3.
    substitution = sp.solve(
        relations, [ns["p0_geometric"], ns["q"], parameters[3]], dict=True
    )[0]
    equations = [
        sp.expand(f.subs(substitution, simultaneous=True))
        for f in ns["curvature_equations"]
    ]
    equations += [
        ns["derivative_of"](relation, direction)
        for relation in relations for direction in range(4)
    ]
    return equations


def check_branches():
    ns = load_frame_constructor()
    a, p, q = ns["a"], ns["p0_geometric"], ns["q"]
    parameters = ns["parameters"]
    derivative = ns["derivative_of"]
    derivative_variables = ns["derivative_variables"]
    assert (str(a), str(p), str(q)) == ("p4", "p0", "p8")

    plus = branch_equations(ns, [p-a, q-a])
    matrix, rhs = sp.linear_eq_to_matrix(plus, derivative_variables)
    witnesses = [
        vector for vector in matrix.T.nullspace()
        if sp.factor((vector.T*rhs)[0]) == 4*a**2
    ]
    assert witnesses
    vector = witnesses[0]
    # Literal polynomial combination: all derivative terms cancel.
    assert sp.expand((vector.T*sp.Matrix(plus))[0]) == -4*a**2

    minus = branch_equations(ns, [p+a, q+a])
    _, compatibility = ns["eliminate_derivatives"](minus)
    h = parameters[12]+parameters[3]+parameters[7]
    assert compatibility and all(
        sp.cancel(f/(a*h)).is_number and sp.cancel(f/(a*h)) != 0
        for f in compatibility
    )
    minus_next = branch_equations(ns, [p+a, q+a, h])
    _, next_compatibility = ns["eliminate_derivatives"](minus_next)
    assert next_compatibility

    # A rational finite jet satisfying the full curvature equations,
    # both negative-sign identities, h=0, and their first derivatives.
    point = dict.fromkeys(parameters, sp.Integer(0))
    point.update({a: sp.Integer(1), p: sp.Integer(-1), q: sp.Integer(-1)})
    equations_at_point = [
        f.subs(point) for f in ns["curvature_equations"]
    ] + [
        derivative(relation, direction)
        for relation in (p+a, q+a, h) for direction in range(4)
    ]
    solution = next(iter(sp.linsolve(equations_at_point, derivative_variables)))
    free = set().union(*(value.free_symbols for value in solution))
    jet = dict(zip(derivative_variables, [
        value.subs(dict.fromkeys(free, 0)) for value in solution
    ]))
    assert all(sp.expand(f.subs(jet)) == 0 for f in equations_at_point)
    assert all(f.subs(point) == 0 for f in next_compatibility)
    d_pq = [derivative(p*q, i).subs(point).subs(jet) for i in range(4)]
    d_ratio = [derivative(p*q/a**2, i).subs(point).subs(jet) for i in range(4)]
    assert d_pq == [0, 3, 0, 0]
    assert d_ratio == [0, 0, 0, 0]
    return {
        "geometric_parameter_names": {"a": str(a), "p": str(p), "q": str(q)},
        "positive_sign_equation_count": len(plus),
        "positive_sign_left_kernel_witness": [str(v) for v in vector],
        "positive_sign_combination": str(-4*a**2),
        "positive_sign_conclusion": "a!=0 and p=q=a are inconsistent after differentiating the branch identities.",
        "negative_sign_compatibility": [str(sp.factor(f)) for f in compatibility],
        "negative_sign_extra_linear_relation": str(h),
        "negative_sign_next_compatibility": [str(sp.factor(f)) for f in next_compatibility],
        "finite_jet": {
            "parameters": {str(k): str(v) for k, v in point.items()},
            "nonzero_directional_derivatives": {str(k): str(v) for k, v in jet.items() if v},
            "checked_equations": len(equations_at_point),
            "d_pq": [str(v) for v in d_pq],
            "d_normalized_ratio": [str(v) for v in d_ratio],
            "boundary": "Compatible finite connection jet only; no all-order integrability or polynomial Hessian pair is claimed.",
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = {
        "format": "hc4-motion-frame-transport-v1",
        "theorem": "HC4MRA1",
        "source_sha256": {
            path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (SOURCE, Path(__file__).resolve())
        },
        "transport": check_transport(),
        "branches": check_branches(),
        "status_effect": "HC4MR1 and HC4MR2 are partial; the negative maximal-motion sign remains open.",
    }
    serialized = json.dumps(data, indent=2)+"\n"
    if args.write:
        OUTPUT.write_text(serialized)
    else:
        assert OUTPUT.read_text() == serialized, "artifact differs; inspect before --write"
    print("PASS determinant-one frame transport gives pq/a^2")
    print("PASS positive-sign closure using differentiated branch identities")
    print("PASS 108 exact negative-sign finite-jet equations with d(pq) nonzero")
    print("PASS " + ("wrote " if args.write else "byte-identical ") + str(OUTPUT))
    print("The surviving finite jet is not an HC4 counterexample or an all-order solution")


if __name__ == "__main__":
    main()
