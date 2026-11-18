#!/usr/bin/env python3
"""Exact bounded controls for the uniform scalar-secant obstruction.

The all-h theorem is the proof in SCALAR_SECANT_PROPORTIONAL_SCALE_OBSTRUCTION.md.
This checker verifies both obstructions, symbolic identities and finite Sturm controls; it
does not turn a bounded control into a uniform formal verification.
"""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    t, x, a, c, Z, u, nu, M, s0, alpha, z2 = sp.symbols(
        "t x a c Z u nu M s0 alpha z2", nonzero=True
    )

    # Divide nu=a(u-cM) by u, and substitute M/u=(nu/u)^2.
    scalar_eq = nu - a * (u - c * M)
    substituted = scalar_eq.subs({nu: Z * u, M: Z**2 * u})
    require(sp.expand(substituted / u) == a * c * Z**2 + Z - a,
            "constant-slope quadratic")

    # Vieta gives t+t^-1=(z1+z2)^2/(z1*z2)-2.
    zsum, zprod = -1 / (a * c), -1 / c
    trace = sp.factor(zsum**2 / zprod - 2)
    require(sp.cancel(trace + 2 + 1 / (a**2 * c)) == 0, "real scalar trace")

    # Value/derivative jets for H=N(t^2*u)-t*N(u), at a slope-z1 root.
    # z1=t*z2; N(alpha)=z1*alpha, N(t^2 alpha)=z2*t^2 alpha.
    value = z2 * t**2 * alpha - t * (t * z2 * alpha)
    derivative = t**2 * z2 / 2 - t * (t * z2 / 2)
    require(sp.expand(value) == 0 and sp.expand(derivative) == 0,
            "paired critical double-zero jets")

    # Proposition2: both scalar elements, over complex coefficients.
    # Euler coefficient comparison is an identity for symbolic degree m.
    m = sp.symbols("m", integer=True, positive=True)
    rho = sp.symbols("rho", nonzero=True)
    Nm = u**m + (2*m-1)*a*u/(2*m-2) - (2*m-1)*a*rho/(2*m)
    Sm = u**m + a*u/(2*m-2) + a*rho/(2*m)
    require(sp.simplify(2*u*sp.diff(Nm, u)-Nm-(2*m-1)*Sm) == 0,
            "uniform both-scalar Euler identity")
    require(sp.simplify(Nm-a*u+a*rho-Sm) == 0,
            "uniform both-scalar monic identity")

    controls = []
    both_scalar_controls = []
    witnesses = {}
    for h in range(1, 13):
        L = 4 * h - 1
        center = 2 * h - 1
        mm = 2*h
        scalar_a = 2*mm*rho**(mm-1)
        Sscalar = sp.Poly(Sm.subs({m: mm, a: scalar_a}), u)
        reversed_S = sp.Poly(sp.cancel(
            u**mm / rho**mm * Sscalar.as_expr().subs(u, rho**2/u)), u)
        require(Sscalar.nth(0) == rho**mm, f"paired norm constant h={h}")
        if mm >= 4:
            require(Sscalar.nth(mm-1) == 0,
                    f"sparse direct coefficient h={h}")
            require(sp.cancel(reversed_S.nth(mm-1)-mm*rho/(mm-1)) == 0,
                    f"nonzero reciprocal coefficient h={h}")
            both_scalar_controls.append({
                "h": h, "numerator_degree": mm,
                "direct_coefficient_u_to_m_minus_1": "0",
                "reciprocal_coefficient_u_to_m_minus_1":
                    str(sp.factor(reversed_S.nth(mm-1))),
                "contradiction": "rho is a unit, so reciprocal coefficient is nonzero"
            })
        else:
            require(sp.expand(Sscalar.as_expr()-(u+rho)**2) == 0,
                    "degree2 repeated critical root")
            require(sp.discriminant(Sscalar.as_expr(), u) == 0,
                    "degree2 zero critical discriminant")
            both_scalar_controls.append({
                "h": h, "numerator_degree": mm,
                "forced_critical_polynomial": "(u+rho)^2",
                "critical_discriminant": "0"
            })
        Q = sp.Poly(sum(t**j for j in range(L)) - L * t**center, t)

        # Constants of H=(t^(4h)-t)S1^2 and N0=-L*t^(2h)*S10^2.
        constant_residual = ((1 - t) * (-L * t**(2*h) * s0**2)
                             - (t**(4*h) - t) * s0**2)
        require(sp.expand(constant_residual
                          + s0**2 * t * (t - 1) * Q.as_expr()) == 0,
                f"constant identity h={h}")

        T = [sp.Integer(2), x]
        while len(T) <= center:
            T.append(sp.expand(x * T[-1] - T[-2]))
        numerator = sp.Poly(1 + sum(T[1:center+1]) - L, x, domain=sp.ZZ)
        B, remainder = sp.div(numerator, sp.Poly(x - 2, x))
        require(remainder.is_zero, f"trace divisibility h={h}")
        require(B.LC() == 1 and B.degree() == 2*h - 2, f"trace degree h={h}")
        require(B.eval(2) == sum(j*j for j in range(1, center+1)),
                f"value at repeated-root trace h={h}")
        laurent_residual = sp.cancel(
            Q.as_expr() - t**center * (t + 1/t - 2) * B.as_expr().subs(x, t + 1/t)
        )
        require(laurent_residual == 0, f"trace Laurent identity h={h}")

        # This is only a bounded check; real positivity in all h is proved
        # in the note by positive/negative real ratio and unit-circle cases.
        roots = int(B.count_roots(-sp.oo, sp.oo)) if B.degree() else 0
        require(roots == 0 and B.eval(0) > 0, f"Sturm positivity h={h}")

        # The degree-drop t^L=1 cannot satisfy Q=0 except at t=1.
        common = sp.gcd(Q, sp.Poly(t**L - 1, t)).monic()
        require(common.as_expr() == t - 1, f"degree-drop control h={h}")
        witnesses[str(h)] = [int(v) for v in B.all_coeffs()]
        controls.append({
            "h": h, "pencil_degree": 4*h, "witness_degree": B.degree(),
            "B_h_at_2": int(B.eval(2)),
            "exact_real_root_count": roots,
            "B_h_at_0": int(B.eval(0)),
            "degree_drop_gcd": "t-1",
            "constant_identity": True,
            "trace_laurent_identity": True
        })

    require(witnesses["1"] == [1], "h1 witness")
    require(witnesses["2"] == [1, 3, 4], "h2 witness")
    require(witnesses["3"] == [1, 3, 2, 1, 5], "degree12 witness")

    inputs = [
        "SCALAR_SECANT_PROPORTIONAL_SCALE_OBSTRUCTION.md",
        "PAIRED_CRITICAL_TRACE_CHART.md",
        "verify_scalar_secant_proportional_scale.py"
    ]
    report = {
        "result": "PASS",
        "claims": [
            "No good real paired even-numerator pencil has scalar a and rho=c*M.",
            "No good complex paired even-numerator pencil has both a and rho scalar."
        ],
        "uniform_degree_scope": "Every h>=1 (pencil degree4h); supplied by the written proof.",
        "bounded_control_scope": "h=1,...,12, exact identities and Sturm counts.",
        "symbolic_controls": {
            "constant_slope_quadratic": True,
            "double_zero_value_and_derivative_jets": True,
            "real_scalar_trace": str(trace),
            "both_scalar_uniform_Euler_identity": True,
            "both_scalar_uniform_monic_identity": True
        },
        "witness_coefficients_descending": witnesses,
        "controls": controls,
        "both_scalar_complex_controls": both_scalar_controls,
        "input_sha256": {
            str((HERE / name).relative_to(ROOT)): digest(HERE / name)
            for name in inputs
        },
        "assurance": {
            "formal_verification": False,
            "independent_external_replay": False,
            "external_human_review": False
        }
    }
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
