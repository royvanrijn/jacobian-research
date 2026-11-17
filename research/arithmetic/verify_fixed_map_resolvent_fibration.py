#!/usr/bin/env python3
"""Exact controls for FIXED_MAP_RESOLVENT_FIBRATION.md (SymPy)."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform

import sympy as sp


def run() -> dict:
    checks = []

    def equal(name, left, right=0):
        assert sp.cancel(sp.expand(left - right)) == 0, name
        checks.append(name)

    pi, ell, v, W, Y = sp.symbols("pi ell v W Y")
    A = sp.Rational(3, 2) / pi**2
    C = sp.Rational(9, 16) / pi
    r = -sp.Rational(3, 2) - ell
    s = A - ell*r - v
    t = (C-v*s)/ell
    h = W**3 + r*W**2 + s*W + t
    q = W**2 + ell*W + v
    beta, delta = t+ell*s+v*r, v*t
    P = W**5 - sp.Rational(3, 2)*W**4 + A*W**3 + beta*W**2 + C*W + delta
    equal("general_factorization_chart", P, q*h)
    D = ell**2-4*v
    Delta = r*r*s*s-4*s**3-4*r**3*t-27*t*t+18*r*s*t
    equal("quadratic_discriminant", sp.discriminant(q, W), D)
    equal("cubic_discriminant", sp.discriminant(h, W), Delta)
    F = sp.expand(D*Delta)
    assert sp.Poly(Delta, v).degree() == 4
    assert sp.Poly(F, v).degree() == 5
    equal("cubic_discriminant_leading_coefficient", sp.Poly(Delta, v).LC(), -27/ell**2)
    equal("resolvent_quintic_leading_coefficient", sp.Poly(F, v).LC(), 108/ell**2)
    z = sp.symbols("z")
    equal("matching_squareclass_forward", (Y*Y-D*Delta).subs(Y, D*z), D*D*(z*z-Delta/D))
    equal("rational_branch_point", F.subs(v, ell**2/4))
    equal("rational_branch_derivative", sp.diff(F, v).subs(v, ell**2/4),
          -4*Delta.subs(v, ell**2/4))

    b = -32*pi*pi/9*(beta+sp.Rational(5, 4)/pi)
    c = -32*delta/9
    assert sp.Poly(b, v).degree() == 2 and sp.Poly(c, v).degree() == 3
    equal("height_b_quadratic_leading_coefficient", sp.Poly(b, v).LC(), -32*pi*pi/(9*ell))
    equal("height_c_cubic_leading_coefficient", sp.Poly(c, v).LC(), -32/(9*ell))
    M, NN = sp.symbols("M NN")
    b2 = sp.cancel(NN**2*b.subs(v, M/NN))
    c3 = sp.cancel(NN**3*c.subs(v, M/NN))
    forms = [NN**3, -pi*NN**3, NN*b2, c3]
    assert all(sp.Poly(form, M, NN).total_degree() == 3 for form in forms)
    checks.append("height_map_homogeneous_degree_three")
    equal("height_map_no_basepoint_at_infinity", c3.subs(NN, 0), -32*M**3/(9*ell))
    qh_resultant = sp.resultant(q, h, W)
    assert sp.Poly(qh_resultant, v).degree() == 4
    equal("quadratic_cubic_resultant_leading_coefficient", sp.Poly(qh_resultant, v).LC(), 1/ell**2)

    aa = (2*ell+3)**2*((20*ell**2+12*ell-27)*pi*pi+54)-243*pi
    bb = (ell**3*(2*ell+3)**4*(5*ell+8)*pi**8
          -3*ell*(7*ell+12)*(23*ell**2+48*ell+18)*pi**7
          +(376*ell**6+2088*ell**5+4302*ell**4+3888*ell**3+1296*ell**2+729)*pi**6
          -162*ell*(7*ell+12)*pi**5
          +6*ell*(2*ell+3)*(76*ell**2+141*ell+36)*pi**4-1458*pi**3
          +18*ell*(71*ell+108)*pi**2+729)
    rr = (ell**5*(ell+2)*(25*ell**2+30*ell-36)*pi**6
          -6*ell*(5*ell**3-72*ell**2-108*ell+108)*pi**5
          +(260*ell**6+696*ell**5+144*ell**4-432*ell**3+729)*pi**4
          +324*ell*(ell+6)*pi**3+12*ell**2*(83*ell**2+132*ell-36)*pi**2
          +1152*ell**2)
    quartic_numerator = sp.cancel(256*ell**2*pi**6*Delta)
    assert sp.denom(quartic_numerator) == 1
    equal("exceptional_quartic_numerator_discriminant", sp.discriminant(quartic_numerator, v),
          -2**40*3*ell**4*pi**22*aa**3*bb)
    equal("exceptional_quartic_rational_discriminant", sp.discriminant(Delta, v),
          -3*aa**3*bb/(256*ell**8*pi**14))
    equal("exceptional_branch_quartic_resultant", sp.resultant(D, Delta, v),
          -3*rr/(ell**2*pi**6))
    equal("exceptional_quintic_discriminant", sp.discriminant(F, v),
          -27*aa**3*bb*rr**2/(256*ell**12*pi**26))
    equal("exceptional_collision_root_evaluation", 256*Delta.subs(v, ell**2/4),
          -3*rr/(ell**2*pi**6))
    assert tuple(expr.subs({pi: 1, ell: 1}) for expr in (aa, bb, rr)) == (1232, 22736, 7424)
    checks.append("exceptional_factors_at_smooth_specialization")

    x, yy = sp.symbols("x yy")
    a_x = x**3*(5*x-24)*pi**2-243*pi+54*x*x
    equal("first_exceptional_component_coordinate_change", aa.subs(ell, (x-3)/2), a_x)
    sextic = 59049-216*x**5*(5*x-24)
    y_map = 2*x**3*(5*x-24)*pi-243
    equal("first_exceptional_component_square_completion", y_map**2-sextic,
          4*x**3*(5*x-24)*a_x)
    pi_inverse = (yy+243)/(2*x**3*(5*x-24))
    equal("first_exceptional_component_inverse_map", y_map.subs(pi, pi_inverse), yy)
    equal("first_exceptional_component_sextic_derivative", sp.diff(sextic, x),
          -6480*x**4*(x-4))
    assert (sextic.subs(x, 0), sextic.subs(x, 4)) == (59049, 943785)
    assert sp.Poly(sextic, x, domain=sp.QQ).gcd(
        sp.Poly(sp.diff(sextic, x), x, domain=sp.QQ)).degree() == 0
    checks.append("first_exceptional_component_squarefree_sextic")
    equal("first_exceptional_component_excluded_case_x_zero", a_x.subs(x, 0), -243*pi)
    exceptional_base = {pi: sp.Rational(128, 25), ell: sp.Rational(9, 10)}
    equal("first_exceptional_component_linear_case", a_x.subs(x, sp.Rational(24, 5)),
          -243*(pi-sp.Rational(128, 25)))
    equal("first_exceptional_component_rational_base_pair", aa.subs(exceptional_base))
    equal("first_exceptional_component_singular_cover", sp.discriminant(F.subs(exceptional_base), v))
    degeneration_linear = 819200*v-243483
    degeneration_cubic = (400*v-81)*(20480000*v*v-110592000*v+143467857)
    degeneration_denominator = sp.Integer(41231686041600000000)
    equal("elliptic_degeneration_exact_factorization", F.subs(exceptional_base),
          degeneration_linear**2*degeneration_cubic/degeneration_denominator)
    degeneration_discriminant = 2**33*3**6*5**14*7**3*1889*2709173**2
    assert sp.discriminant(degeneration_cubic, v) == degeneration_discriminant != 0
    checks.append("elliptic_degeneration_nonzero_cubic_discriminant")
    equal("elliptic_degeneration_rational_branch", degeneration_cubic.subs(v, sp.Rational(81, 400)))
    assert degeneration_cubic.subs(v, sp.Rational(243483, 819200)) != 0
    checks.append("elliptic_degeneration_square_factor_is_separate")
    equal("elliptic_degeneration_deleted_value_is_inadmissible",
          Delta.subs(exceptional_base).subs(v, sp.Rational(243483, 819200)))
    ZZ, XX, TT = sp.symbols("ZZ XX TT")
    curve_a, curve_b = sp.Integer(-10758382), sp.Integer(6600759744)
    X_map = sp.Rational(6400, 3)*v-3984
    T_scale = sp.Rational(20971520000, 3)
    T_map = T_scale*ZZ
    equal("elliptic_integral_short_model_transformation",
          T_map**2-X_map**3-curve_a*X_map-curve_b,
          T_scale**2*(ZZ**2-degeneration_cubic/degeneration_denominator))
    equal("elliptic_integral_model_inverse_v", X_map.subs(v, 3*(XX+3984)/6400), XX)
    equal("elliptic_integral_model_inverse_Z", T_map.subs(ZZ, 3*TT/20971520000), TT)
    model_discriminant = -16*(4*curve_a**3+27*curve_b**2)
    assert model_discriminant == 2**9*5**2*7**3*1889*2709173**2 != 0
    checks.append("elliptic_integral_model_nonsingular_discriminant")
    point_x, point_y = sp.Rational(199200, 49), sp.Rational(59459616, 343)
    equal("elliptic_nontorsion_point_on_curve", point_y**2,
          point_x**3+curve_a*point_x+curve_b)
    assert sp.denom(point_x) == 49 != 1
    checks.append("elliptic_point_nonintegral_x_for_lutz_nagell")
    assert 3*(point_x+3984)/6400 == sp.Rational(73953, 19600)
    checks.append("elliptic_nontorsion_point_parameter")

    f = ((4*v-1)*(6912*v**4-67840*v**3+205216*v*v-194256*v+59043)/256)
    equal("genus_two_specialization", F.subs({pi: 1, ell: 1}), f)
    expected_discriminant = -2**24 * 3**3 * 7**5 * 11**3 * 29**3
    exact_discriminant = sp.discriminant(f, v)
    assert exact_discriminant == expected_discriminant != 0
    checks.append("exact_nonzero_specialization_discriminant")

    f5 = 3*v**5+3*v**4+4*v**3+3*v+2
    numerator = sp.expand(256*f)
    assert sp.Poly(numerator-f5, v, modulus=5).is_zero
    assert sp.Poly((v*v+2)*f5+(v**4-2*v-1)*sp.diff(f5, v)-1, v, modulus=5).is_zero
    assert sp.gcd(sp.Poly(f5, v, modulus=5), sp.Poly(sp.diff(f5, v), v, modulus=5)).degree() == 0
    checks.extend(["specialization_reduction_mod_five", "explicit_mod_five_bezout_identity",
                   "independent_finite_field_squarefreeness"])
    assert 2*(-2)+6 == 2 and (2+2)//2 == 2
    checks.append("six_branch_riemann_hurwitz_integer_arithmetic")

    # The centered branch is independent of the noncentered denominator ell.
    tt = sp.symbols("tt")
    centered_h = W**3-sp.Rational(3, 2)*W**2+(A-v)*W+tt
    centered_P = (W**5-sp.Rational(3, 2)*W**4+A*W**3
                  +(tt-sp.Rational(3, 2)*v)*W**2+C*W+v*tt)
    equal("centered_single_target_coefficient_constraint", (W*W+v)*centered_h-centered_P,
          (v*(A-v)-C)*W)
    constraint = 16*pi*pi*v*v+9*pi-24*v
    equal("centered_elliptic_equation", v*(A-v)-C, -constraint/(16*pi*pi))
    X = 8*v/3
    YY = (32*v*v*pi+9)/9
    equal("centered_elliptic_reduction", YY*YY-X**3-1, 64*v*v*constraint/81)
    centered_cases = [(1, sp.Rational(3, 4)), (-2, sp.Rational(3, 4)),
                      (-2, -sp.Rational(3, 8))]
    for index, (pp, vv) in enumerate(centered_cases, 1):
        equal(f"centered_case_{index}", constraint.subs({pi: pp, v: vv}))

    # A real Hasse example from the moving surface lifts to the cover, so
    # irreducibility/local checks are not confused with the algebraic open set.
    example = {pi: sp.Rational(17, 29), ell: -sp.Rational(18, 17),
               v: sp.Rational(4431, 1156)}
    value = sp.factor(F.subs(example))
    assert value > 0
    yn = sp.integer_nthroot(int(sp.numer(value)), 2)
    yd = sp.integer_nthroot(int(sp.denom(value)), 2)
    assert yn[1] and yd[1]
    y_value = sp.Rational(yn[0], yd[0])
    equal("moving_hasse_example_lift", y_value*y_value, value)
    assert sp.discriminant(P.subs(example), W) != 0
    checks.append("moving_hasse_example_has_squarefree_inverse")

    other_example = {pi: sp.Rational(128, 25), ell: sp.Rational(9, 10),
                     v: sp.Rational(73953, 19600)}
    other_q = q.subs(other_example)
    other_h = h.subs(other_example)
    equal("nonpure_example_quadratic_field", sp.discriminant(other_q, W),
          -6*sp.Rational(54, 35)**2)
    shifted_other_h = sp.Poly(sp.expand(other_h.subs(W, W+1)), W)
    equal("nonpure_example_shifted_cubic", shifted_other_h.as_expr(),
          W**3+sp.Rational(3, 5)*W**2-sp.Rational(134708181, 40140800)*W
          +sp.Rational(290236959237, 78675968000))
    for index in range(3):
        numerator, denominator = sp.fraction(shifted_other_h.nth(index))
        assert numerator % 3 == 0 and numerator % 9 != 0 and denominator % 3 != 0
    checks.append("nonpure_example_eisenstein_at_three")
    ratio = sp.factor(sp.discriminant(other_h, W)/sp.discriminant(other_q, W))
    numerator, denominator = sp.fraction(ratio)
    assert numerator > 0 and denominator > 0
    assert sp.integer_nthroot(int(numerator), 2)[1] and sp.integer_nthroot(int(denominator), 2)[1]
    assert sp.degree(sp.gcd(other_q, other_h), W) == 0
    checks.append("nonpure_example_matching_resolvent_full_rank_five")

    here = Path(__file__).resolve().parent
    faltings_path = here/"sources"/"faltings-1983-endlichkeitssaetze.pdf"
    assert faltings_path.is_file() and faltings_path.read_bytes().startswith(b"%PDF-")
    lutz_path = here/"sources"/"dujella-number-theory-2021-section15-3.pdf"
    assert lutz_path.is_file() and lutz_path.read_bytes().startswith(b"%PDF-")
    return {
        "status": "pass", "python": platform.python_version(), "sympy": sp.__version__,
        "command": ".venv/bin/python research/arithmetic/verify_fixed_map_resolvent_fibration.py",
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "proof_sha256": sha256((here/"FIXED_MAP_RESOLVENT_FIBRATION.md").read_bytes()).hexdigest(),
        "exact_controls": checks,
        "specialization": {"pi": 1, "quadratic_linear_coefficient": 1,
                           "discriminant": str(exact_discriminant),
                           "discriminant_factorization": "-2^24*3^3*7^5*11^3*29^3",
                           "squarefreeness_certificate_prime": 5},
        "moving_example_lift_Y": str(y_value),
        "external_primary_source": {
            "path": "sources/faltings-1983-endlichkeitssaetze.pdf",
            "sha256": sha256(faltings_path.read_bytes()).hexdigest(),
            "statement": "Satz 7, p. 365 (PDF p. 17): genus at least two implies finitely many rational points",
            "proof_imported_not_replayed": True,
        },
        "external_primary_lutz_nagell_source": {
            "path": "sources/dujella-number-theory-2021-section15-3.pdf",
            "sha256": sha256(lutz_path.read_bytes()).hexdigest(),
            "statement": "Theorem 15.5, p. 488: rational torsion on a nonsingular integral short Weierstrass model has integral coordinates",
            "proof_imported_not_replayed": True,
        },
        "scope": ("Exact coefficient/discriminant identities and one explicit squarefree "
                  "specialization prove the stated nonvanishing and exact exceptional base "
                  "with the written arguments. Height-map identities and branch controls "
                  "support the elementary height and two-torsion proofs. Faltings is an "
                  "explicit primary-source import; Lutz--Nagell is imported to infer "
                  "nontorsion from the exact on-curve nonintegral point. "
                  "The checker does not classify rational points of the genus-two fibres, "
                  "their singular specializations or the finite arithmetic local conditions."),
        "formal_verification": False, "independent_whole_proof_replay": False,
        "external_review": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", action="store_true")
    args = parser.parse_args()
    result = run()
    output = json.dumps(result, indent=2, sort_keys=True)+"\n"
    if args.write_receipt:
        Path(__file__).resolve().with_name("fixed_map_resolvent_fibration_checks.json").write_text(
            output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()
