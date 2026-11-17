#!/usr/bin/env python3
"""Exact QSO34 Ore and resonance controls; universal geometry is a written proof."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
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

    x = sp.symbols("x")
    u, v, U, V = [sp.Function(q)(x) for q in ("u", "v", "U", "V")]
    c, c1, c0 = sp.symbols("c c1 c0")

    def mul(left, right):
        out = {}
        for i, ai in left.items():
            for j, bj in right.items():
                for k in range(i + 1):
                    degree = i + j - k
                    out[degree] = out.get(degree, 0) + ai * sp.binomial(i, k) * sp.diff(bj, x, k)
        return {i: sp.expand(t) for i, t in out.items()}

    def comm(left, right):
        lr, rl = mul(left, right), mul(right, left)
        result = {}
        for i in set(lr) | set(rl):
            coefficient = sp.simplify(lr.get(i, 0) - rl.get(i, 0))
            if coefficient != 0:
                result[i] = coefficient
        return result

    def same(left, right):
        return all(sp.simplify(left.get(i, 0) - right.get(i, 0)) == 0
                   for i in set(left) | set(right))

    # No pseudodifferential expansion is imported: use finite normal-order multiplication.
    p0 = {3: sp.Integer(1), 1: u, 0: v}
    q0 = {4: sp.Integer(1), 2: 4*u/3+c, 1: 2*sp.diff(u, x)/3+4*v/3+c1,
          0: 2*sp.diff(u, x, 2)/9+2*u*u/9+2*sp.diff(v, x)/3+2*c*u/3+c0}
    generic_comm = comm(p0, q0)
    assert set(generic_comm) <= {0, 1}

    def normal_pair(uu, vv):
        pp = {3: sp.Integer(1), 1: -sp.Rational(3, 2)*uu,
              0: -sp.Rational(3, 4)*sp.diff(uu, x)+sp.Rational(3, 2)*vv}
        qq = {4: sp.Integer(1), 2: c-2*uu, 1: 2*(vv-sp.diff(uu, x)),
              0: sp.diff(vv, x)+uu*uu/2-sp.Rational(5, 6)*sp.diff(uu, x, 2)-c*uu}
        return pp, qq

    pshift, qshift = normal_pair(-2*u/3, 2*v/3-sp.diff(u, x)/3+c1/2)
    expected_p = dict(p0)
    expected_p[0] += 3*c1/4
    expected_q = dict(q0)
    expected_q[0] -= c0
    assert same(pshift, expected_p)
    assert same(qshift, expected_q)

    p, q = normal_pair(U, V)
    I1 = sp.diff(V, x, 2)/2-sp.Rational(3, 2)*U*V+sp.Rational(3, 2)*c*V
    I2 = (sp.diff(U, x, 4)/12-sp.Rational(3, 4)*U*sp.diff(U, x, 2)
          -sp.Rational(3, 8)*sp.diff(U, x)**2+sp.Rational(3, 2)*V**2+U**3/2
          -c*(3*U**2-sp.diff(U, x, 2))/4)
    expected = {1: -2*sp.diff(I1, x), 0: -sp.diff(I1, x, 2)-sp.diff(I2, x)}
    assert same(comm(p, q), expected)
    corrupted = dict(q)
    corrupted[0] += sp.diff(U, x, 2)/6
    assert not same(comm(p, corrupted), expected), "Corrupted coefficient unexpectedly passed"

    # Re-derive principal balances and matrix by symbolic first variation of the ODE.
    t = sp.symbols("t", positive=True)
    r, eps, alpha, beta, uu, vv = sp.symbols("r eps alpha beta uu vv")
    utest = uu/t**2+eps*alpha*t**(r-2)
    vtest = vv/t**3+eps*beta*t**(r-3)
    f1 = sp.diff(vtest, t, 2)/2-sp.Rational(3, 2)*utest*vtest
    f2 = (sp.diff(utest, t, 4)/12-sp.Rational(3, 4)*utest*sp.diff(utest, t, 2)
          -sp.Rational(3, 8)*sp.diff(utest, t)**2
          +sp.Rational(3, 2)*vtest**2+utest**3/2)
    balances = [sp.expand(f.subs({t: 1, eps: 0})) for f in (f1, f2)]
    solutions = sorted((a, b) for a, b in sp.solve(balances, (uu, vv)) if a != 0)
    assert solutions == [(2, 0), (4, -4), (4, 4), (10, 0)]
    variation = [sp.expand(sp.diff(f, eps).subs({eps: 0, t: 1})) for f in (f1, f2)]
    matrix = sp.Matrix([[sp.diff(f, a) for a in (alpha, beta)] for f in variation])
    A = (sp.prod(r-j for j in (2, 3, 4, 5))/12
         -sp.Rational(3, 4)*uu*((r-2)*(r-3)+6)
         +sp.Rational(3, 2)*uu*(r-2)+sp.Rational(3, 2)*uu*uu)
    displayed = sp.Matrix([[-sp.Rational(3, 2)*vv, ((r-3)*(r-4)-3*uu)/2],
                           [A, 3*vv]])
    assert all(sp.simplify(z) == 0 for z in matrix-displayed)

    factors = {}
    expected_roots = {(2, 0): (-1, 1, 3, 4, 6, 8),
                      (4, -4): (-2, -1, 3, 4, 8, 9),
                      (4, 4): (-2, -1, 3, 4, 8, 9),
                      (10, 0): (-5, -2, -1, 8, 9, 12)}
    for principal in solutions:
        determinant = sp.factor(-matrix.det().subs({uu: principal[0], vv: principal[1]}))
        roots = sp.roots(determinant, r)
        assert len(roots) == 6 and all(root.is_Integer and multiplicity == 1
                                       for root, multiplicity in roots.items())
        assert tuple(sorted(roots)) == expected_roots[principal]
        factors[str(principal)] = {"minus_determinant": str(determinant),
                                  "all_roots": [int(root) for root in sorted(roots)]}

    output = args.output or (Path(tempfile.mkdtemp(prefix="qso34-controls-")) / "cubic_quartic_string_controls.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    record = {
        "claim": "QSO34", "status": "PASS",
        "normalization_high_order_equations": "PASS",
        "constant_lower_flow_removal": "PASS",
        "universal_commutator_identity": "[P,Q]=-2*I1prime*E-I1second-I2prime",
        "universal_commutator_identity_check": "PASS",
        "corrupted_quadratic_derivative_coefficient_rejected": True,
        "principal_balance_exhaustion": [[int(a), int(b)] for a,b in solutions],
        "resonance_matrices_rederived_by_ODE_linearization": "PASS",
        "resonance_determinants": factors,
        "primary_source": {
            "author": "Nathan Hayford",
            "title": "The Ising Model Coupled to 2D Gravity: Higher-order Painleve Equations/The (3,4) String Equation",
            "arxiv": "2405.03260v2", "version_date": "2025-06-16",
            "url": "https://arxiv.org/pdf/2405.03260v2",
            "cached_file": str(PDF.relative_to(HERE)),
            "sha256": PDF_SHA256,
            "inspected_locations": "Appendix A.2, equations (A.8), (A.12), Proposition A.1 and (A.13), PDF pages 38-39; introductory equation (1.1).",
            "import_scope": "Comparison of normalized operators and string equations. No algebraic-solution exclusion or full irreducibility theorem is imported."
        },
        "written_proof_coverage": [
            "Every algebraic-coefficient order-(3,4) pair reaches the normal form, including constant lower flows.",
            "The scalar commutator creates an algebraic time primitive with derivative one.",
            "All finite singular balances are exhausted, including unequal poles and bounded components.",
            "The first fractional Puiseux index is excluded by integer resonance roots.",
            "Bounded components give a regular first-order ODE and no fractional branching.",
            "Riemann-Hurwitz makes the algebraic coefficient functions rational in the created time.",
            "Infinity degree comparison excludes rational solutions.",
            "The exact-R h-regular prescribed Long pair with generic S order three requires generic T order at least five, using QSO2."
        ],
        "hashes": {name: sha(HERE/name) for name in [
            "CUBIC_QUARTIC_STRING_OBSTRUCTION.md", "check_cubic_quartic_string.py",
            "QUADRATIC_OPERATOR_MATE_CLASSIFICATION.md"]},
        "coverage_boundary": "Exact symbolic identities, finite exhaustive principal-coefficient balances and resonance factors. Universal normalization, all algebraic places, Puiseux induction, unramified curve argument, infinity exclusion and Long valuation consequence use written proofs. No generic order-(3,5) or arbitrary partner obstruction; no bounded full Long lift or DC1/DC2 solution.",
        "independent_replay": False, "formal_verification": False, "external_review": False
    }
    output.write_text(json.dumps(record, indent=2)+"\n", encoding="utf-8")
    print("QSO34 exact controls PASS")
    print("Verification output:", output)


if __name__ == "__main__":
    main()
