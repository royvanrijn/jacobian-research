#!/usr/bin/env python3
"""Exact controls for FIXED_MAP_MOVING_QUADRATIC_FAMILY.md (SymPy).

The checker does not certify an asymptotic, a rational-point density,
the general Keller reconstruction, or an external review.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import platform

import sympy as sp


def run() -> dict:
    controls = []

    def equal(name, left, right=0):
        assert sp.cancel(sp.expand(left - right)) == 0, name
        controls.append(name)

    pi, alpha, d, W, T, Z = sp.symbols("pi alpha d W T Z")
    m, n = sp.symbols("m n")
    ell = 3 * (alpha - sp.Rational(1, 2))
    v = sp.Rational(3, 2) / pi**2 + 6 * alpha**2 - sp.Rational(9, 2) * alpha
    a = ((10 * alpha**4 - 8 * alpha**3 + 3 * alpha**2 / pi**2
          - sp.Rational(3, 8) / pi) / (2 * alpha - 1))
    h0 = -alpha**3 - a
    beta = h0 + 3 * ell * alpha**2 - 3 * alpha * v
    b = -sp.Rational(32, 9) * pi**2 * (beta + sp.Rational(5, 4) / pi)
    c = -sp.Rational(32, 9) * v * h0
    q = W**2 + ell * W + v
    h = (W - alpha)**3 - a
    P = (W**5 - sp.Rational(3, 2) * W**4 + sp.Rational(3, 2) * W**3 / pi**2
         + beta * W**2 + sp.Rational(9, 16) * W / pi + v * h0)
    equal("whole_factorization", P, q * h)
    equal("quadratic_discriminant_relation", ell**2 - 4*v,
          -3 * (5 * alpha**2 - 3 * alpha - sp.Rational(3, 4) + 2 / pi**2))
    quadric = Z**2 - 5*T**2 + 3*pi*T + sp.Rational(3, 4)*pi**2 - 2
    equal("quadric_from_discriminant", quadric.subs({T: pi*alpha, Z: pi*d}),
          pi**2 * (d*d - 5*alpha**2 + 3*alpha + sp.Rational(3, 4) - 2/pi**2))

    denominator = n*n - 5*m*m + 3*m + sp.Rational(3, 4)
    lam = -(2*n - 2*m + 3) / denominator
    par_pi, par_T, par_Z = 1 + lam, sp.Rational(1, 2) + m*lam, 1 + n*lam
    equal("quadric_rational_parametrization", quadric.subs({pi: par_pi, T: par_T, Z: par_Z}))
    matrix = sp.Matrix([[sp.Rational(3, 4), sp.Rational(3, 2), 0, 0],
                        [sp.Rational(3, 2), -5, 0, 0], [0, 0, 1, 0], [0, 0, 0, -2]])
    assert matrix.det() == 12
    controls.append("projective_quadric_nonsingularity")

    example = {m: 1, n: -sp.Rational(2, 3)}
    pp = sp.cancel(par_pi.subs(example))
    aa = sp.cancel((par_T / par_pi).subs(example))
    dd = sp.cancel((par_Z / par_pi).subs(example))
    assert (pp, aa, dd) == (sp.Rational(17, 29), sp.Rational(5, 34), sp.Rational(37, 17))
    example_subs = {pi: pp, alpha: aa}
    expected_a = sp.Rational(26263, 34**3)
    equal("example_pure_cubic_parameter", a.subs(example_subs), expected_a)
    equal("example_quadratic_linear_coefficient", ell.subs(example_subs), -sp.Rational(18, 17))
    equal("example_quadratic_constant", v.subs(example_subs), sp.Rational(4431, 1156))
    equal("example_b", b.subs(example_subs), sp.Rational(46972, 128673))
    equal("example_c", c.subs(example_subs), sp.Rational(12991692, 1419857))
    equal("example_discriminant", (ell*ell - 4*v).subs(example_subs), -3*dd*dd)
    assert -pp not in (-1, 2)
    controls.append("example_outside_both_old_target_lines")
    primality_bound = isqrt(26263)
    assert all(26263 % p for p in range(2, primality_bound + 1))
    assert 26263 % 9 == 1 and 34 % 3 != 0 and pow(34, 3, 9) == 1
    controls.append("example_prime_and_three_adic_cube_unit")
    assert 34 == 2 * 17 and 2 % 3 == 17 % 3 == 2
    controls.append("example_inert_denominator_exponents_are_multiples_of_three")
    assert sp.discriminant(h.subs(example_subs), W) == -27 * expected_a**2
    assert sp.resultant(q.subs(example_subs), h.subs(example_subs), W) != 0
    controls.append("example_squarefree_coprime_rank_five")

    # The claimed exceptional algebras are controlled by independent root
    # factors: a nonzero cube has a simple rational root, while a=0 is triple.
    rr, aa0 = sp.symbols("rr aa0")
    general_h = (W - aa0)**3 - rr**3
    equal("cube_factorization", general_h,
          (W - aa0 - rr) * ((W - aa0)**2 + rr*(W - aa0) + rr*rr))
    equal("cube_rational_root_derivative", sp.diff(general_h, W).subs(W, aa0 + rr), 3*rr*rr)
    equal("zero_cubic_is_triple", general_h.subs(rr, 0), (W - aa0)**3)

    # Exact identities and exhaustive residue sets used in the 3-adic proof.
    equal("quadratic_constant_divided_by_three", v/3,
          1/(2*pi*pi) + 2*alpha*alpha - sp.Rational(3, 2)*alpha)
    equal("coefficient_beta_from_target", -sp.Rational(5, 4)/pi - sp.Rational(9, 32)*b/pi**2,
          beta)
    equal("scaled_pure_cubic_parameter", (a*pi**3).subs(alpha, T/pi),
          (10*T**4 - 8*pi*T**3 + 3*T*T - sp.Rational(3, 8)*pi**3)/(2*T-pi))
    numerator = 10*T**4 - 8*pi*T**3 + 3*T*T - sp.Rational(3, 8)*pi**3
    equal("mod_nine_cube_obstruction_identity", numerator - 2*T**3*(2*T-pi),
          3*T*T*(2*T*T - 2*pi*T + 1) - sp.Rational(3, 8)*pi**3)
    zeta = sp.symbols("zeta")
    assert sp.rem((zeta-1)**2 + 3*zeta, zeta*zeta + zeta + 1, zeta) == 0
    controls.append("cyclotomic_root_difference_valuation_identity")
    unit_residue_controls = 0
    for pi_residue in (1, 2):
        for alpha_residue in range(3):
            assert (pow(2*pi_residue*pi_residue, -1, 3) + 2*alpha_residue**2) % 3 != 0
            unit_residue_controls += 1
    cube_residue_controls = 0
    for pi_residue in (0, 3, 6):
        for t_residue in (1, 2, 4, 5, 7, 8):
            val = ((10*t_residue**4 - 8*pi_residue*t_residue**3 + 3*t_residue*t_residue)
                   * pow((2*t_residue-pi_residue) % 9, -1, 9)) % 9
            assert val == 2*pow(t_residue, 3, 9) % 9 and val in (2, 7)
            assert val not in (1, 8)
            cube_residue_controls += 1

    centered_pi_constraint = sp.factor((3*v*alpha**2 - ell*alpha**3
                                        - sp.Rational(9, 16)/pi).subs(alpha, sp.Rational(1, 2)))
    equal("centered_exception_only_old_pi_values", centered_pi_constraint,
          -sp.Rational(9, 16)*(pi-1)*(pi+2)/pi**2)
    equal("centered_pi_one_quadratic_constant", v.subs({pi: 1, alpha: sp.Rational(1, 2)}),
          sp.Rational(3, 4))
    equal("centered_pi_minus_two_has_wrong_quadratic_field", v.subs({pi: -2, alpha: sp.Rational(1, 2)}),
          -sp.Rational(3, 8))

    here = Path(__file__).resolve().parent
    proof = here / "FIXED_MAP_MOVING_QUADRATIC_FAMILY.md"
    return {
        "status": "pass", "python": platform.python_version(), "sympy": sp.__version__,
        "command": ".venv/bin/python research/arithmetic/verify_fixed_map_moving_quadratic.py",
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "proof_sha256": sha256(proof.read_bytes()).hexdigest(),
        "symbolic_controls": controls, "unit_mod_three_controls": unit_residue_controls,
        "cube_mod_nine_controls": cube_residue_controls,
        "example": {"pi": str(pp), "alpha": str(aa), "d": str(dd),
                    "pure_cubic_parameter": str(expected_a),
                    "target": [str(-pp), "46972/128673", "12991692/1419857"],
                    "prime": 26263, "complete_trial_division_through": primality_bound},
        "scope": ("Exact coefficient, parametrization, primality, coprimality and residue controls. "
                  "The general derivative-unit Keller reconstruction and rational local criterion "
                  "are adopted from their canonical proofs; no density or counting extrapolation."),
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
        Path(__file__).resolve().with_name("fixed_map_moving_quadratic_checks.json").write_text(
            output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()
