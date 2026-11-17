"""Exact controls for algebraic Gamma prime-power parameter domains.

Written proof: ALGEBRAIC_GAMMA_PRIME_POWER_DOMAINS.tex.
The all-index conclusions use its coefficientwise congruence and closed-disk
proof. These finite controls verify the explicit obstruction and illustrative
residue calculations; they do not establish the unbounded quantifiers.
"""
from pathlib import Path
import hashlib
import json
import platform
import sympy as s


def replay():
    x, t = s.symbols("x t")
    checked = []

    def check(condition, label):
        assert condition, label
        checked.append(label)

    quadratic = s.Poly(x ** 2 + 7 * x + 15, x)
    determinant = s.expand((s.rf(x + 4, 4) - s.rf(x, 4)) / 4)
    check(s.expand(determinant - 2 * (2 * x + 7) * quadratic.as_expr()) == 0,
          "literal gap-four determinant factorization")
    check(s.Poly(determinant / 2, x, modulus=2)
          == s.Poly(x ** 2 + x + 1, x, modulus=2),
          "normalized determinant has nonconstant residue")
    check(s.rem(s.rf(x, 4), quadratic.as_expr(), x) == -45,
          "fourth rising factorial in the quadratic field")
    check(s.rem(s.rf(x, 8), quadratic.as_expr(), x) == 2025,
          "eighth rising factorial in the quadratic field")
    moment = 1 + 2 * s.rf(x, 4) * t + s.rf(x, 8) * t ** 2
    check(s.rem(moment - (1 - 45 * t) ** 2, quadratic.as_expr(), x) == 0,
          "actual moment is a nonconstant perfect square")
    check(s.rem(s.discriminant(moment, t), quadratic.as_expr(), x) == 0,
          "literal moment discriminant vanishes")
    check(s.discriminant(quadratic) == -11,
          "obstruction shape has no real embedding")
    for residue in range(2):
        check(int(quadratic.eval(residue)) % 2 != 0,
              "obstruction residue has degree two over F2: " + str(residue))
    for prime in (2, 3, 5, 7):
        degree_two_det = s.Poly(
            s.expand((s.rf(x + prime, prime) - s.rf(x, prime)) / prime), x)
        check(degree_two_det.trunc(prime) == s.Poly(-1, x).trunc(prime),
              "degree-two prime-gap coefficient congruence: p=" + str(prime))
    golden = s.Poly(x ** 2 - x - 1, x)
    for residue in range(3):
        check(int(golden.eval(residue)) % 3 != 0,
              "positive golden-ratio shape has non-prime-field residue: " + str(residue))
    assert len(checked) == 16
    return {
        "status": "pass",
        "arithmetic": "exact integer-polynomial and finite-field",
        "exact_controls": len(checked),
        "controls": checked,
        "counterexample": {
            "gap": 4, "degree": 2,
            "parameter_minimal_polynomial_ascending": [15, 7, 1],
            "parameter_is_two_adic_unit": True,
            "parameter_residue_degree_over_F2": 2,
            "moment_polynomial": "(1-45*t)^2",
            "normalized_determinant_residue": "x^2+x+1",
            "has_positive_real_shape": False,
        },
        "prime_domain": "all integral Cp gamma, by the written coefficientwise congruence",
        "proper_prime_power_domain": "Zp+p^s O_Cp, by the written closed-disk proof",
        "all_integral_proper_prime_power_extension": False,
        "all_positive_algebraic_shapes_claimed": False,
        "finite_controls_prove_unbounded_scope": False,
        "proof": "ALGEBRAIC_GAMMA_PRIME_POWER_DOMAINS.tex",
        "python_version": platform.python_version(),
        "sympy_version": s.__version__,
        "independent_replay": False,
        "formal_verification": False,
        "external_review": False,
        "priority_claim": False,
    }


if __name__ == "__main__":
    report = replay()
    here = Path(__file__).resolve().parent
    report["input_sha256"] = {
        name: hashlib.sha256((here / name).read_bytes()).hexdigest()
        for name in (
            "ALGEBRAIC_GAMMA_PRIME_POWER_DOMAINS.tex",
            "ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex",
            "PRIME_POWER_MARKED_RADIAL_SFC.tex",
            "algebraic_gamma_prime_power_domains.py",
        )
    }
    destination = here / "algebraic_gamma_prime_power_domains_verification.json"
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
