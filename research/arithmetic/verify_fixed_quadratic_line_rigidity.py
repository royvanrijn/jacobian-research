#!/usr/bin/env python3
"""Exact symbolic controls for FIXED_QUADRATIC_LINE_RIGIDITY.md.

Requires SymPy, as in the existing fixed-map ambient checker. General
elliptic descent and local-global inputs remain external written theorems.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, isqrt
from pathlib import Path
import platform

import sympy as sp


def run() -> dict:
    checks = []

    def equal(name: str, left, right=0) -> None:
        assert sp.cancel(sp.expand(left - right)) == 0, name
        checks.append(name)

    S, W, pi, b, c, u, v, t = sp.symbols("S W pi b c u v t")
    a, rpar = sp.symbols("a rpar")
    E = (pi**5 * S**5 - sp.Rational(3, 2) * pi**4 * S**4
         + sp.Rational(3, 2) * pi * S**3
         + (-sp.Rational(5, 4) * pi - sp.Rational(9, 32) * b) * S**2
         + sp.Rational(9, 16) * S - sp.Rational(9, 32) * c)
    A, C = sp.Rational(3, 2) / pi**2, sp.Rational(9, 16) / pi
    beta = -sp.Rational(5, 4) / pi - sp.Rational(9, 32) * b / pi**2
    delta = -sp.Rational(9, 32) * c
    P = W**5 - sp.Rational(3, 2) * W**4 + A * W**3 + beta * W**2 + C * W + delta
    equal("normalized_inverse_polynomial", E.subs(S, W / pi), P)

    lam, p0, p1, b0, b1, c0, c1 = sp.symbols("lam p0 p1 b0 b1 c0 c1")
    on_line = E.subs({pi: p0 + p1 * lam, b: b0 + b1 * lam, c: c0 + c1 * lam})
    equal("varying_pi_degree_five_obstruction", sp.Poly(on_line, lam).coeff_monomial(lam**5), p1**5 * S**5)
    db, dd = sp.symbols("db dd")
    equal("line_difference_centering", db * W**2 + dd - db * (W**2 + u * W + v),
          -db * u * W + dd - db * v)

    h = W**3 - sp.Rational(3, 2) * W**2 + (A - v) * W + t
    coefficient_identity = ((W**2 + v) * h
                            - (W**5 - sp.Rational(3, 2) * W**4 + A * W**3
                               + (t - sp.Rational(3, 2) * v) * W**2 + C * W + v * t))
    equal("centered_factor_condition", coefficient_identity, (v * (A - v) - C) * W)
    line_b = -sp.Rational(32, 9) * pi**2 * (t - sp.Rational(3, 2) * v + sp.Rational(5, 4) / pi)
    line_c = -sp.Rational(32, 9) * v * t
    constraint = 16 * pi**2 * v**2 + 9 * pi - 24 * v
    equal("line_factorization_modulo_constraint",
          P.subs({b: line_b, c: line_c}) - (W**2 + v) * h,
          constraint * W / (16 * pi**2))

    X, Y = 8 * v / 3, (32 * v**2 * pi + 9) / 9
    equal("elliptic_curve_reduction", Y**2 - X**3 - 1, 64 * v**2 * constraint / 81)
    xx, yy = sp.symbols("xx yy")
    equal("elliptic_inverse_constraint",
          constraint.subs({v: 3 * xx / 8, pi: 2 * (yy - 1) / xx**2}),
          9 * (yy**2 - xx**3 - 1) / xx**2)
    equal("rank_zero_curve_translation", yy**2 - (xx - 1)**3 - 1,
          yy**2 - xx**3 + 3 * xx**2 - 3 * xx)
    assert (-2 * (-3), (-3)**2 - 4 * 3) == (6, -3)
    assert 16 * 3**2 * ((-3)**2 - 4 * 3) != 0
    assert 16 * (-3)**2 * (6**2 - 4 * (-3)) != 0
    checks.append("two_isogenous_integral_coefficients_and_nonsingularity")

    squares_mod_nine = {n * n % 9 for n in range(9)}
    primitive_residue_controls = 0
    for m in range(9):
        for e in range(9):
            if gcd(gcd(m, e), 3) != 1:
                continue
            assert (-m**4 + 6 * m*m * e*e + 3 * e**4) % 9 not in squares_mod_nine
            primitive_residue_controls += 1
    assert primitive_residue_controls == 72
    checks.append("complete_primitive_mod_nine_descent_obstruction")

    ys = [0] + [sign * d for d in range(1, 6) if 27 % (d*d) == 0 for sign in (-1, 1)]
    affine_points = set()
    for y in ys:
        roots = [x for x in range(-1, 3) if x**3 == y*y - 1]
        assert len(roots) == 1
        affine_points.add((roots[0], y))
    expected = {(-1, 0), (0, 1), (0, -1), (2, 3), (2, -3)}
    assert affine_points == expected
    checks.append("complete_lutz_nagell_candidate_enumeration_after_rank_zero")
    cases = []
    for x, y in sorted(affine_points):
        if x == 0:
            continue
        vv = sp.Rational(3, 8) * x
        pp = sp.Rational(2 * (y - 1), x*x)
        equal(f"elliptic_point_{x}_{y}_line_constraint", constraint.subs({pi: pp, v: vv}))
        cases.append((pp, vv))
    assert set(cases) == {(sp.Integer(1), sp.Rational(3, 4)),
                          (sp.Integer(-2), sp.Rational(3, 4)),
                          (sp.Integer(-2), -sp.Rational(3, 8))}

    targets = [
        (sp.Integer(1), sp.Rational(3, 4), -a - sp.Rational(1, 8),
         sp.Rational(32, 9) * a, (8 * a + 1) / 3),
        (sp.Integer(-2), sp.Rational(3, 4), sp.Rational(7, 16) + 8 * rpar,
         sp.Rational(56, 3) - sp.Rational(1024, 9) * rpar,
         -sp.Rational(7, 6) - sp.Rational(64, 3) * rpar),
        (sp.Integer(-2), -sp.Rational(3, 8), -a - sp.Rational(1, 8),
         sp.Rational(8, 3) + sp.Rational(128, 9) * a,
         -sp.Rational(1, 6) - sp.Rational(4, 3) * a),
    ]
    for index, (pp, vv, tt, bb, cc) in enumerate(targets, 1):
        equal(f"line_{index}_b_formula", line_b.subs({pi: pp, v: vv, t: tt}), bb)
        equal(f"line_{index}_c_formula", line_c.subs({pi: pp, v: vv, t: tt}), cc)
        equal(f"line_{index}_factorization", P.subs({pi: pp, b: bb, c: cc}),
              ((W**2 + v) * h).subs({pi: pp, v: vv, t: tt}))

    for index, (c_weight, b_weight, expected) in enumerate(((12, -9, 4),
                                                          (48, -9, -224),
                                                          (96, 9, 8))):
        bb, cc = targets[index][3:5]
        equal(f"line_{index + 1}_integral_target_obstruction", c_weight * cc + b_weight * bb,
              expected)
        assert c_weight % 3 == b_weight % 3 == 0 and expected % 3 != 0

    quadratic = W**2 - sp.Rational(3, 8)
    cubic = (W - sp.Rational(1, 2))**3 - a
    equal("third_line_quadratic_discriminant", sp.discriminant(quadratic, W), sp.Rational(3, 2))
    equal("third_line_cubic_discriminant", sp.discriminant(cubic, W), -27 * a*a)
    equal("third_line_coprime_resultant", sp.resultant(quadratic, cubic, W),
          (a + sp.Rational(11, 16))**2 - sp.Rational(243, 512))
    assert isqrt(6)**2 != 6
    checks.append("third_line_wrong_quadratic_resolvent_squareclasses")

    rr = -sp.Rational(3, 2) - u
    ss = A - u * rr - v
    tt = (C - v * ss) / u
    chart_b = -sp.Rational(32, 9) * pi**2 * (tt + u * ss + v * rr + sp.Rational(5, 4) / pi)
    chart_c = -sp.Rational(32, 9) * v * tt
    equal("noncentered_quadratic_incidence_chart", P.subs({b: chart_b, c: chart_c}),
          (W**2 + u * W + v) * (W**3 + rr * W**2 + ss * W + tt))
    d = sp.symbols("d")
    equal("quadratic_K_minus_three_discriminant", (u*u - 4*v).subs(v, (u*u + 3*d*d) / 4), -3*d*d)

    # Reject candidates from nonpermitted elliptic points, a moving-pi line,
    # and the mistaken assertion that every centered factor has field -3.
    assert constraint.subs({pi: 2, v: sp.Rational(3, 4)}) != 0
    assert sp.Poly(on_line, lam).coeff_monomial(lam**5).subs({p1: 1, S: 1}) != 0
    assert sp.discriminant(quadratic, W) / (-3) < 0
    rejection_controls = 3

    here = Path(__file__).resolve().parent
    proof = here / "FIXED_QUADRATIC_LINE_RIGIDITY.md"
    context_path = here / "sources/dujella-number-theory-2021-context.json"
    context = json.loads(context_path.read_text())
    for entry in context["sources"]:
        assert sha256((context_path.parent / entry["path"]).read_bytes()).hexdigest() == entry["sha256"]
    return {
        "status": "pass", "python": platform.python_version(), "sympy": sp.__version__,
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "proof_sha256": sha256(proof.read_bytes()).hexdigest(),
        "source_context_sha256": sha256(context_path.read_bytes()).hexdigest(),
        "command": ".venv/bin/python research/arithmetic/verify_fixed_quadratic_line_rigidity.py",
        "symbolic_controls": checks,
        "primitive_mod_nine_controls": primitive_residue_controls,
        "rejection_controls": rejection_controls,
        "affine_auxiliary_points": [[x, y] for x, y in sorted(affine_points)],
        "classified_pi_v": [[str(pp), str(vv)] for pp, vv in cases],
        "source_pdfs": context["sources"],
        "scope": ("Exact polynomial, finite residue and finite Lutz-Nagell candidate controls. "
                  "The imported general two-isogeny descent, Mordell-Weil, Lutz-Nagell "
                  "and arithmetic local-global theorems are not formalized by this checker."),
        "formal_verification": False, "independent_whole_proof_replay": False,
        "external_review": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", action="store_true")
    args = parser.parse_args()
    result = run()
    output = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.write_receipt:
        Path(__file__).resolve().with_name("fixed_quadratic_line_rigidity_verification.json").write_text(
            output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()
