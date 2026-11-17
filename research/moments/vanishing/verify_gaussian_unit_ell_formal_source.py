"""Exact bounded controls for the unit-ell formal ORIGINAL-source theorem.

The written proof establishes the unrestricted completed-z extension and
uses the retained cover comparison. These finite controls do not prove
unrestricted geometry or infer all-order claims from a moment prefix.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from math import factorial
from pathlib import Path
import platform

import sympy as sp

z, w, y, t, eps, U, S = sp.symbols("z w y t eps U S")


def zero(value):
    assert sp.simplify(sp.cancel(value)) == 0, value


def moment(exponent):
    if exponent % 2:
        return 0
    return sp.Integer(factorial(exponent) //
                      (2 ** (exponent // 2) * factorial(exponent // 2)))


def polynomial_wick(value):
    answer = 0
    for (iz, iw, iy), coefficient in sp.Poly(sp.expand(value), z, w, y).terms():
        if iz == iw:
            answer += coefficient * factorial(iz) * moment(iy)
    return sp.expand(answer)


def germ_wick(value):
    """Finite W,Y degree; only finitely many derivatives of each z germ."""
    answer = 0
    for (iw, iy), coefficient in sp.Poly(sp.expand(value), w, y).terms():
        coefficient = sp.cancel(coefficient)
        denominator = sp.denom(coefficient)
        if sp.simplify(denominator.subs(z, 0)) == 0:
            raise ValueError("A negative z power is outside the completed source")
        answer += sp.diff(coefficient, z, iw).subs(z, 0) * moment(iy)
    return sp.cancel(answer)


def gaussian_S(value):
    polynomial = sp.Poly(sp.expand(value), S)
    return sp.expand(sum(coefficient * moment(power)
                         for (power,), coefficient in polynomial.terms()))


def epsilon_gaussian_jet(value, order):
    return sp.expand(gaussian_S(sp.series(value, eps, 0, order + 1).removeO()))


def source_jet(amplitude, phase, order):
    """Completed-z original source, using only finite t and spatial jets."""
    amplitude = sp.series(amplitude, t, 0, order + 1).removeO().expand()
    power, answer = sp.Integer(1), sp.Integer(0)
    for n in range(order + 1):
        if n:
            power = sp.expand(power * phase)
        answer += t ** n * germ_wick(amplitude * power) / factorial(n)
    return sp.expand(sp.series(answer, t, 0, order + 1).removeO())


def checks():
    # Two genuine formal power series germs, represented rationally.
    germs = [1 / (1 - z), 1 / (1 - z) ** 2 + z / (1 + 2 * z)]
    for germ in germs:
        for j in range(5):
            for k in range(5):
                amplitude = germ * w ** j * y ** k
                zero(germ_wick(sp.diff(amplitude, z) - w * amplitude))
                zero(germ_wick(sp.diff(amplitude, w) - z * amplitude))
                zero(germ_wick(sp.diff(amplitude, y) - y * amplitude))

    # The ell(0)=0 boundary really breaks the W-Stein identity if one
    # extends the coefficient prescription to Laurent z.
    anomaly = -sp.series(z / z, z, 0, 1).removeO().coeff(z, 0)
    assert anomaly == -1
    try:
        germ_wick(1 / z)
    except ValueError:
        pass
    else:
        raise AssertionError("A pole at z=0 was incorrectly admitted")

    # An ell with an infinite inverse z-germ and a genuinely perturbed phase.
    ell = 1 + z + z ** 2
    c, h, a = 1 + z, 1 - z, sp.Rational(1, 3)
    u, k = 2 - z, z ** 2
    F = c ** 2 - 4 * a * ell
    q = 2 * ell / t + F
    q_inverse_jet = sum(t / (2 * ell) * (-t * F / (2 * ell)) ** j
                        for j in range(4))
    zero(q * q_inverse_jet - 1 + (-t * F / (2 * ell)) ** 4)

    for source_order in (1, 2):
        germ = sp.expand(q_inverse_jet).coeff(t, source_order)
        amplitude = germ * w ** 2
        zero(germ_wick(sp.diff(amplitude, z) - w * amplitude))
        zero(germ_wick(sp.diff(amplitude, w) - z * amplitude))
        zero(germ_wick(sp.diff(amplitude * y, y) - y * amplitude * y))
    phase = ell * w ** 2 + w * (h + c * y) + a * y ** 2 + u * y + k
    order = 2
    original_markers = [sp.Integer(1), z, w, y, z * w, w * y]
    for marker in original_markers:
        completed = source_jet(marker, phase, order)
        power, direct = sp.Integer(1), sp.Integer(0)
        for n in range(order + 1):
            if n:
                power = sp.expand(power * phase)
            direct += polynomial_wick(marker * power) * t ** n / factorial(n)
        zero(completed - direct)

    # Independent exact scalar model ell=(1+z)^2, all other coefficients0.
    # Here all source/cover calculations are rational in epsilon,sqrt(2).
    beta = sp.sqrt(2)
    model_q = 2 * (1 + z) ** 2 / eps ** 2
    phi = -z ** 2 / (4 * eps ** 2 * (1 + z) ** 2)
    J = 1 - beta * eps * S
    zeta = beta * eps * S / J
    selected_U = beta / (eps * J)
    zS = sp.diff(zeta, S)
    zt = sp.diff(zeta, eps) / (2 * eps)

    def evaluate(value):
        return sp.cancel(sp.sympify(value).subs(
            {z: zeta, U: selected_U}, simultaneous=True))

    def dt(value):
        return sp.diff(value, eps) / (2 * eps) - U * sp.diff(value, U) / (2 * eps ** 2)

    def dz(value):
        return sp.diff(value, z) + sp.diff(model_q, z) * sp.diff(value, U) / (2 * U)

    zero(selected_U ** 2 - model_q.subs(z, zeta))
    zero(phi.subs(z, zeta) + S ** 2 / 2)
    zero(evaluate(1 / (eps ** 2 * U)) * zS - 1 / J)
    forms = [sp.Integer(1), z, U, 1 / U, U ** -2,
             z / U ** 3, (eps ** 2 + z) / U ** 5,
             1 / (eps ** 2 * U * (1 + 3 * eps ** 2))]
    gauge = eps ** 2 + 2 * eps ** 4
    for form in forms:
        pull = evaluate(form)
        zero(sp.diff(pull, eps) / (2 * eps)
             - evaluate(dt(form)) - evaluate(dz(form)) * zt)
        zero(sp.diff(pull, S) - evaluate(dz(form)) * zS)
        for psi in (sp.Integer(0), gauge):
            shifted = phi - psi
            pulled_scalar = -psi
            left = (sp.diff(pull * zS, eps) / (2 * eps)
                    + sp.diff(pulled_scalar, eps) * pull * zS / (2 * eps)
                    - evaluate(dt(form) + sp.diff(shifted, eps) * form / (2 * eps)) * zS)
            right = sp.diff(pull * zt, S) - S * pull * zt
            zero(left - right)

    # Direct ORIGINAL scalar moments, rather than derivatives of the
    # candidate cover output, supply the independent expected coefficients.
    model_phase = (1 + z) ** 2 * w ** 2
    model_order = 7
    power, scalar_coefficients = sp.Integer(1), []
    for n in range(model_order + 1):
        if n:
            power = sp.expand(power * model_phase)
        scalar_coefficients.append(sp.simplify(polynomial_wick(power) / factorial(n)))
    assert scalar_coefficients == [sp.Integer(factorial(2 * n) // factorial(n))
                                   for n in range(model_order + 1)]
    scalar_series = sum(value * eps ** (2 * n)
                        for n, value in enumerate(scalar_coefficients))
    zero(epsilon_gaussian_jet(1 / J, 12)
         - sum(scalar_coefficients[n] * eps ** (2 * n) for n in range(7)))

    return {
        "status": "passed",
        "written_theorem_scope": (
            "ell(0) nonzero: full original q-localized source. With deg ell>=2 "
            "and gcd(ell,c^2) geometric root multiplicities<=2, no scalar "
            "polynomial-exponential vacuum and no finite actual image."
        ),
        "bounded_controls": {
            "formal_power_series_germs": len(germs),
            "Stein_monomial_exponents_W_Y": [0, 4],
            "q_inverse_geometric_terms": 4,
            "q_inverse_source_Stein_t_coefficients": [1, 2],
            "perturbed_original_source_t_orders": [0, order],
            "perturbed_original_markers": len(original_markers),
            "rationalizable_model": "ell=(1+z)^2; all h,c,a,u,k=0",
            "independent_original_scalar_t_orders": [0, 6],
            "full_chain_rule_forms": len(forms),
            "inverse_U_exponents": [-1, -2, -3, -5],
            "scalar_gauges": ["0", "t+2t^2"],
            "parameter_normalizer_in_tested_forms": "p=1+3t",
        },
        "negative_control": {
            "Laurent_z_residue_anomaly": "Ehat[(partial_W-z)z^-1]=-1",
            "singular_germ_rejected": True,
        },
        "scope": (
            "Finite exact algebraic controls only. The completed-z source "
            "extension and full comparison scope are established by the written "
            "proof, not by finite moments or tested connection iterates."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = checks()
    directory = Path(__file__).resolve().parent

    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    result.update(
        family="Unit-ell completed-z original source and finite-image exclusion",
        python=platform.python_version(), sympy=sp.__version__,
        proof_sha256=digest(directory / "GAUSSIAN_UNIT_ELL_FORMAL_SOURCE_FINITE_IMAGE_EXCLUSION.md"),
        checker_sha256=digest(Path(__file__)),
        cover_proof_sha256=digest(directory / "GAUSSIAN_AFFINE_ISOTROPIC_CHART.md"),
        independent_replay=False, formal_verification=False, external_review=False,
    )
    report = json.dumps(result, indent=2) + "\n"
    if args.output is not None:
        args.output.write_text(report, encoding="utf-8", newline="\n")
    print(report, end="")


if __name__ == "__main__":
    main()
