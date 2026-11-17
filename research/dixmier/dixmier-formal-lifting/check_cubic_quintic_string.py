#!/usr/bin/env python3
"""QSO35 exact finite Ore, complete principal-balance and resonance controls."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import tempfile

import sympy as sp

HERE = Path(__file__).resolve().parent
PDF = HERE / "sources/hayford-2405.03260v2.pdf"
PDF_SHA256 = "3763cddb3e723b8685c7974ad54a0873323c7f076b7cc6cc96705d4d97824c5c"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert sha(PDF) == PDF_SHA256, "Primary-source PDF pin mismatch"
    hierarchy_audit_path = HERE / "cubic_hierarchy_source_audit.json"
    hierarchy_audit = json.loads(hierarchy_audit_path.read_text())
    for source in hierarchy_audit["sources"]:
        assert sha(HERE / source["cached_file"]) == source["sha256"], "Hierarchy source pin mismatch"
    x = sp.symbols("x")
    U, V = sp.Function("U")(x), sp.Function("V")(x)
    c4, c2, c1 = sp.symbols("c4 c2 c1")

    def mul(left, right):
        out = {}
        for i, ai in left.items():
            for j, bj in right.items():
                for k in range(i+1):
                    n = i+j-k
                    out[n] = out.get(n, 0)+ai*sp.binomial(i, k)*sp.diff(bj, x, k)
        return {n: sp.expand(coefficient) for n, coefficient in out.items()}

    def comm(left, right):
        lr, rl = mul(left, right), mul(right, left)
        out = {}
        for n in set(lr) | set(rl):
            coefficient = sp.simplify(lr.get(n, 0)-rl.get(n, 0))
            if coefficient != 0:
                out[n] = coefficient
        return out

    def same(left, right):
        return all(sp.simplify(left.get(n, 0)-right.get(n, 0)) == 0
                   for n in set(left) | set(right))

    def add(operators):
        result = {}
        for scalar, operator in operators:
            for n, coefficient in operator.items():
                result[n] = result.get(n, 0)+scalar*coefficient
        return result

    p = {3: sp.Integer(1), 1: -sp.Rational(3, 2)*U,
         0: -sp.Rational(3, 4)*sp.diff(U, x)+sp.Rational(3, 2)*V}
    b1 = {1: sp.Integer(1)}
    b2 = {2: sp.Integer(1), 0: -U}
    b4 = {4: sp.Integer(1), 2: -2*U, 1: 2*(V-sp.diff(U, x)),
          0: sp.diff(V, x)+U**2/2-sp.Rational(5, 6)*sp.diff(U, x, 2)}
    b5 = {5: sp.Integer(1), 3: -sp.Rational(5, 2)*U,
          2: sp.Rational(5, 2)*(V-sp.Rational(3, 2)*sp.diff(U, x)),
          1: sp.Rational(5, 4)*(U**2-sp.Rational(7, 3)*sp.diff(U, x, 2)+2*sp.diff(V, x)),
          0: sp.Rational(5, 4)*(sp.Rational(4, 3)*sp.diff(V, x, 2)+U*sp.diff(U, x)
                                -sp.Rational(2, 3)*sp.diff(U, x, 3)-2*U*V)}
    for operator in (b1, b2, b4, b5):
        assert set(comm(p, operator)) <= {0, 1}
    assert comm(p, p) == {}
    f = sp.Function("f")(x)
    for order in range(6):
        assert sp.simplify(comm(p, {order: f})[order+2]-3*sp.diff(f, x)) == 0

    A = sp.diff(V, x, 2)/2-sp.Rational(3, 2)*U*V
    B = (sp.diff(U, x, 4)/12-sp.Rational(3, 4)*U*sp.diff(U, x, 2)
         -sp.Rational(3, 8)*sp.diff(U, x)**2+sp.Rational(3, 2)*V**2+U**3/2)
    H1 = (-sp.diff(U, x, 4)/6+sp.Rational(5, 4)*U*sp.diff(U, x, 2)
          +sp.Rational(15, 16)*sp.diff(U, x)**2-sp.Rational(5, 8)*U**3
          -sp.Rational(15, 4)*V**2)
    H2 = (sp.diff(V, x, 4)/6-sp.Rational(5, 4)*U*sp.diff(V, x, 2)
          -sp.Rational(5, 8)*sp.diff(U, x)*sp.diff(V, x)
          -sp.Rational(5, 8)*sp.diff(U, x, 2)*V+sp.Rational(15, 8)*U**2*V)
    J1 = H1-2*c4*A-3*c2*V+sp.Rational(3, 2)*c1*U
    J2 = H2-c4*B+c2*(sp.Rational(3, 4)*U**2-sp.diff(U, x, 2)/4)-sp.Rational(3, 2)*c1*V
    q = add([(1, b5), (c4, b4), (c2, b2), (c1, b1)])
    expected = {1: sp.diff(J1, x), 0: sp.diff(J1, x, 2)/2+sp.diff(J2, x)}
    assert same(comm(p, q), expected)
    corrupted = dict(q)
    corrupted[0] += U*V/4
    assert not same(comm(p, corrupted), expected), "Corrupted B5 coefficient passed"

    top_matrix = sp.Matrix([[sp.diff(z, derivative) for derivative in
                             (sp.diff(U, x, 4), sp.diff(V, x, 4))]
                            for z in (J1, J2)])
    assert top_matrix == sp.Matrix([[-sp.Rational(1, 6), 0],
                                    [-c4/12, sp.Rational(1, 6)]])
    assert sp.simplify(top_matrix.det()) == -sp.Rational(1, 36)

    t = sp.symbols("t", positive=True)
    r, eps, du, dv, u, v = sp.symbols("r eps du dv u v")
    ut = u/t**2+eps*du*t**(r-2)
    vt = v/t**3+eps*dv*t**(r-3)
    h1t = (-sp.diff(ut, t, 4)/6+sp.Rational(5, 4)*ut*sp.diff(ut, t, 2)
           +sp.Rational(15, 16)*sp.diff(ut, t)**2-sp.Rational(5, 8)*ut**3
           -sp.Rational(15, 4)*vt**2)
    h2t = (sp.diff(vt, t, 4)/6-sp.Rational(5, 4)*ut*sp.diff(vt, t, 2)
           -sp.Rational(5, 8)*sp.diff(ut, t)*sp.diff(vt, t)
           -sp.Rational(5, 8)*sp.diff(ut, t, 2)*vt+sp.Rational(15, 8)*ut**2*vt)
    balances = [sp.factor(z.subs({eps: 0, t: 1})) for z in (h1t, h2t)]
    principal = sorted((a, b) for a, b in sp.solve(balances, (u, v)) if a != 0)
    assert principal == [(2, 0), (4, -4), (4, 4), (8, -8), (8, 8), (16, 0)]
    variation = [sp.expand(sp.diff(z, eps).subs({eps: 0, t: 1})) for z in (h1t, h2t)]
    matrix = sp.Matrix([[sp.diff(z, coefficient) for coefficient in (du, dv)]
                        for z in variation])
    W = r*r-8*r+18
    displayed = sp.Matrix([
        [-sp.prod(r-j for j in (2, 3, 4, 5))/6+sp.Rational(5, 4)*u*W-sp.Rational(15, 8)*u*u,
         -sp.Rational(15, 2)*v],
        [-sp.Rational(5, 8)*v*(r*r-8*r+36-6*u),
         sp.prod(r-j for j in (3, 4, 5, 6))/6-sp.Rational(5, 4)*u*W+sp.Rational(15, 8)*u*u]
    ])
    assert all(sp.simplify(z) == 0 for z in matrix-displayed)
    expected_roots = {
        (2, 0): [-1, 1, 3, 3, 5, 5, 7, 9],
        (4, -4): [-2, -1, 2, 3, 5, 6, 9, 10],
        (4, 4): [-2, -1, 2, 3, 5, 6, 9, 10],
        (8, -8): [-4, -2, -1, 3, 5, 9, 10, 12],
        (8, 8): [-4, -2, -1, 3, 5, 9, 10, 12],
        (16, 0): [-7, -4, -2, -1, 9, 10, 12, 15]
    }
    determinants = {}
    for pair in principal:
        determinant = sp.factor(matrix.det().subs({u: pair[0], v: pair[1]}))
        root_mult = sp.roots(determinant, r)
        assert all(root.is_Integer for root in root_mult)
        roots = sorted(int(root) for root, multiplicity in root_mult.items()
                       for _ in range(multiplicity))
        assert roots == expected_roots[pair] and len(roots) == 8
        determinants[str(pair)] = {"determinant": str(determinant),
                                    "integer_roots_with_multiplicity": roots}

    output = args.output or (Path(tempfile.mkdtemp(prefix="qso35-controls-")) / "cubic_quintic_string_controls.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    record = {
        "claim": "QSO35", "status": "PASS",
        "environment": {"python": platform.python_version(), "sympy": sp.__version__},
        "finite_normal_form_commutator_order_at_most_one": "PASS for B1,B2,B4,B5",
        "leading_coefficient_extraction_for_orders_0_through_5": "3*fprime at order n+2",
        "all_lower_flows": ["c4*B4", "c2*B2", "c1*B1"],
        "universal_commutator_identity": "[P,Q]=J1prime*E+J1second/2+J2prime",
        "universal_commutator_check": "PASS",
        "corrupted_B5_coefficient_rejected": True,
        "principal_balance_exhaustion": [[int(a), int(b)] for a,b in principal],
        "resonance_matrix_rederived_by_ODE_first_variation": "PASS",
        "resonance_determinants": determinants,
        "bounded_jet_leading_matrix": [["-1/6", "0"], ["-c4/12", "1/6"]],
        "bounded_jet_matrix_determinant": "-1/36",
        "primary_source": {
            "author": "Nathan Hayford",
            "arxiv": "2405.03260v2",
            "url": "https://arxiv.org/pdf/2405.03260v2",
            "cached_file": str(PDF.relative_to(HERE)), "sha256": PDF_SHA256,
            "inspected_location": "Appendix A.2: positive fractional powers Q^(4/3), Q^(5/3), PDF page 39.",
            "import_scope": "Comparison of displayed finite operator formulas only; no algebraic-solution exclusion theorem imported."
        },
        "hashes": {name: sha(HERE/name) for name in [
            "CUBIC_QUINTIC_STRING_OBSTRUCTION.md", "check_cubic_quintic_string.py",
            "CUBIC_QUARTIC_STRING_OBSTRUCTION.md", "QUADRATIC_OPERATOR_MATE_CLASSIFICATION.md",
            "cubic_hierarchy_source_audit.json"]},
        "general_cubic_source_boundary": {
            "source_audit": "cubic_hierarchy_source_audit.json",
            "cached_source_pins": "PASS",
            "Miwa_regular_nonresonance": "Section 2 equation (2.22) is required",
            "Malgrange_scalar_reconstruction": "Application obligation; meromorphic connection forms alone do not give meromorphic normalized scalar coefficients"
        },
        "written_proof_coverage": [
            "Arbitrary algebraic coefficient order-(3,5) normalization with every lower constant flow.",
            "Time primitive constructed inside the coefficient field.",
            "All finite-pole exponent balances exhausted by strict pole comparisons.",
            "No fractional Puiseux exponent at a finite pole or bounded finite point.",
            "Curve-cover rationality by Riemann-Hurwitz and rational infinity contradiction.",
            "Exact-R h-regular prescribed Long cubic S requires mate order at least seven using polynomial leading-order subtraction and QSO2/QSO34/QSO35."
        ],
        "coverage_boundary": "Exact finite Ore identities, exhaustive principal-coefficient balances, resonance factors and bounded-jet leading matrix. Universal field normalization, all algebraic places, fractional-index induction, rationality, infinity and valued Long arguments are written proofs. No order-(3,7) theorem, all-cubic-partner exclusion, finite Long lift or DC1/DC2 resolution is claimed.",
        "independent_replay": False, "formal_verification": False, "external_review": False
    }
    output.write_text(json.dumps(record, indent=2)+"\n", encoding="utf-8")
    print("QSO35 exact controls PASS")
    print("Verification output:", output)


if __name__ == "__main__":
    main()
