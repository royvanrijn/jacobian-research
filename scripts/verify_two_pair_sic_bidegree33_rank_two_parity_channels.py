#!/usr/bin/env python3
"""Exact parity-channel exclusion on a direct rank-two cubic chart.

The family uses the reversal-even and reversal-odd channel lines

    U = [[1, 1], [b, d], [b, -d], [1, -1]],
    W = [[p, q, q, p], [r, s, -s, -r]],
    C = U W.

The one-parameter SL2 centralizer of reversal moves every pair of nonzero
even/odd U channels into the nonzero-endpoint chart.  There the two channel
scalings are quotiented by normalizing the top entries of the two columns
of U.  The projective W channel row is covered by r=1, r=0,s=1, and the
rank-one boundary r=s=0.  Exact coefficient rank two is imposed by
complete two-open covers rather than by a selected minor.

All polynomial construction is over QQ.  Two independent variable orders
are passed to msolve over characteristic zero.  The semistable exact-rank-
two opens are unit ideals through moment six.  On the remaining theta=0
slice, Singular gives exactly two minimal primes through moment six; both
are certified over QQ(q) as fixed-flag one-sided families.  Their relative
periods give the recurrence nu_(m+1)=0 and an exact mixed-tail bound.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from math import factorial
from pathlib import Path
import shutil
import subprocess
import tempfile

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_parity_channels.json"
)


def primitive_integer_polynomial(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> tuple[sp.Poly, int, int]:
    """Return primitive P and integers numerator,denominator with f=n/d P."""

    polynomial = sp.Poly(sp.together(expression), *variables, domain=sp.QQ)
    denominator, integral = polynomial.clear_denoms(convert=True)
    content, primitive = integral.primitive()
    return primitive, int(content), int(denominator)


def msolve_expression(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


def msolve_source(
    variable_order: tuple[str, ...],
    moments: list[sp.Poly],
    rabinowitsch: sp.Expr,
) -> str:
    equations = [msolve_expression(moment.as_expr()) for moment in moments]
    equations.append(msolve_expression(rabinowitsch))
    return (
        ",".join(variable_order)
        + "\n0\n"
        + ",\n".join(equations)
        + "\n"
    )


def run_msolve(source: str, threads: int) -> dict[str, object]:
    executable = shutil.which("msolve")
    if executable is None:
        raise RuntimeError("msolve is required for the exact QQ replay")
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".ms", prefix="sic33-parity-", delete=True
    ) as input_handle, tempfile.NamedTemporaryFile(
        mode="r", suffix=".out", prefix="sic33-parity-", delete=True
    ) as output_handle:
        input_handle.write(source)
        input_handle.flush()
        process = subprocess.run(
            [
                executable,
                "-f",
                input_handle.name,
                "-o",
                output_handle.name,
                "-t",
                str(threads),
                "-v",
                "1",
                "--random-seed",
                "0",
            ],
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        solver_output = output_handle.read().strip()
    if process.returncode != 0:
        raise RuntimeError(
            "msolve failed with return code "
            f"{process.returncode}: {process.stderr[-1000:]}"
        )
    diagnostic = process.stdout + process.stderr
    if solver_output == "[-1]:":
        status = "unit_ideal"
    elif "The ideal has positive dimension" in diagnostic:
        status = "positive_dimensional"
    else:
        status = "zero_dimensional"
    return {
        "status": status,
        "output": solver_output,
        "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
        "solver_summary": (
            "positive-dimensional Groebner basis"
            if status == "positive_dimensional"
            else "Groebner basis has a single element; no solution"
            if status == "unit_ideal"
            else "zero-dimensional output"
        ),
    }


def run_singular_min_ass(source: str) -> dict[str, str]:
    """Run an exact characteristic-zero minimal-prime decomposition."""

    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required for the exact QQ replay")
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".sing", prefix="sic33-parity-", delete=True
    ) as input_handle:
        input_handle.write(source)
        input_handle.flush()
        process = subprocess.run(
            [executable, "-q", input_handle.name],
            check=False,
            capture_output=True,
            text=True,
            timeout=300,
        )
    if process.returncode != 0:
        raise RuntimeError(
            "Singular failed with return code "
            f"{process.returncode}: {process.stderr[-1000:]}"
        )
    return {
        "output": process.stdout.strip(),
        "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    if arguments.threads < 1:
        raise ValueError("--threads must be positive")

    h, b, d, p, q, r, s, x, y = sp.symbols(
        "h b d p q r s x y"
    )
    chart_variables = (b, d, q, s)

    # The centralizer of reversal supplies a global orbit cover of the
    # normalized-U chart.  In eigen-coordinates X=W+V,Y=W-V, put
    # X -> lambda*X and Y -> lambda^(-1)*Y.  This is induced by an SL2
    # matrix commuting with reversal.  The endpoint coefficients of a
    # nonzero even and odd binary cubic become nonzero affine-linear
    # polynomials in z=lambda^4, so one z avoids both exceptional roots.
    a0, b0, c0, d0, lam = sp.symbols("a0 b0 c0 d0 lam")
    capital_w0, capital_v0 = sp.symbols("W0 V0")
    alpha = (lam + lam**-1) / 2
    beta = (lam - lam**-1) / 2
    assert sp.simplify(alpha**2 - beta**2) == 1
    centralizer_matrix = sp.Matrix([[alpha, beta], [beta, alpha]])
    reversal_matrix = sp.Matrix([[0, 1], [1, 0]])
    assert centralizer_matrix.det().simplify() == 1
    assert (
        centralizer_matrix * reversal_matrix
        == reversal_matrix * centralizer_matrix
    )
    even_cubic = (
        a0 * (capital_w0**3 + capital_v0**3)
        + b0
        * (
            capital_w0**2 * capital_v0
            + capital_w0 * capital_v0**2
        )
    )
    odd_cubic = (
        c0 * (capital_w0**3 - capital_v0**3)
        + d0
        * (
            capital_w0**2 * capital_v0
            - capital_w0 * capital_v0**2
        )
    )
    centralizer_substitution = {
        capital_w0: alpha * capital_w0 + beta * capital_v0,
        capital_v0: beta * capital_w0 + alpha * capital_v0,
    }
    transformed_even = sp.Poly(
        sp.expand(even_cubic.subs(centralizer_substitution, simultaneous=True)),
        capital_w0,
        capital_v0,
    )
    transformed_odd = sp.Poly(
        sp.expand(odd_cubic.subs(centralizer_substitution, simultaneous=True)),
        capital_w0,
        capital_v0,
    )
    even_endpoint = sp.factor(
        transformed_even.coeff_monomial(capital_w0**3)
    )
    odd_endpoint = sp.factor(
        transformed_odd.coeff_monomial(capital_w0**3)
    )
    expected_even_endpoint = (
        (a0 + b0) * lam**4 + 3 * a0 - b0
    ) / (4 * lam)
    expected_odd_endpoint = (
        (3 * c0 + d0) * lam**4 + c0 - d0
    ) / (4 * lam**3)
    assert sp.cancel(even_endpoint - expected_even_endpoint) == 0
    assert sp.cancel(odd_endpoint - expected_odd_endpoint) == 0
    even_coefficient_matrix = sp.Matrix([[1, 1], [3, -1]])
    odd_coefficient_matrix = sp.Matrix([[3, 1], [1, -1]])
    assert even_coefficient_matrix.det() == -4
    assert odd_coefficient_matrix.det() == -4

    even_dual = (1, b, b, 1)
    odd_dual = (1, d, -d, -1)
    even_target = (p, q, q, p)
    odd_target = (r, s, -s, -r)
    u_matrix = sp.Matrix.hstack(
        sp.Matrix(even_dual), sp.Matrix(odd_dual)
    )
    w_matrix = sp.Matrix([even_target, odd_target])
    coefficient_matrix = u_matrix * w_matrix

    coefficient_polynomial = sum(
        coefficient_matrix[i, j] * x**i * y**j
        for i in range(4)
        for j in range(4)
    )
    power = sp.Poly(1, x, y)
    raw_moments: list[sp.Expr] = []
    for order in range(1, 7):
        power *= sp.Poly(coefficient_polynomial, x, y)
        raw_moments.append(
            sp.expand(
                sum(
                    factorial(index)
                    * factorial(3 * order - index)
                    * power.coeff_monomial(x**index * y**index)
                    for index in range(3 * order + 1)
                )
            )
        )

    expected_first = 4 * (b * q + d * s + 3 * p + 3 * r)
    assert sp.expand(raw_moments[0] - expected_first) == 0

    projective_substitution = {
        r: sp.Integer(1),
        p: -(b * q + d * s + 3) / 3,
    }
    assert sp.expand(raw_moments[0].subs(projective_substitution)) == 0

    chart_matrix = coefficient_matrix.subs(projective_substitution)
    apolar_diagonal = sp.diag(6, 2, 2, 6)
    trace_square = sp.factor(
        sp.trace((chart_matrix * apolar_diagonal) ** 2)
    )
    theta = d * s + 3
    assert sp.expand(trace_square - 32 * theta**2) == 0

    u_minor = sp.det(u_matrix.extract((0, 1), (0, 1)))
    w_minor = sp.det(w_matrix.extract((0, 1), (0, 1)))
    c_minor = sp.det(coefficient_matrix.extract((0, 1), (0, 1)))
    assert sp.expand(c_minor - u_minor * w_minor) == 0
    localized_minor = sp.factor(
        3 * (u_minor * w_minor).subs(projective_substitution)
    )
    delta = (b - d) * (b * q * s + d * s**2 + 3 * q + 3 * s)
    assert sp.expand(localized_minor - delta) == 0

    primitive_moments: list[sp.Poly] = []
    scalings: list[dict[str, int]] = []
    for order, raw in enumerate(raw_moments[1:], start=2):
        primitive, numerator, denominator = primitive_integer_polynomial(
            raw.subs(projective_substitution), chart_variables
        )
        primitive_moments.append(primitive)
        scalings.append(
            {
                "order": order,
                "raw_equals_primitive_times_numerator": numerator,
                "raw_equals_primitive_divided_by_denominator": denominator,
            }
        )

    expected_profiles = {
        2: (29, 6),
        3: (93, 9),
        4: (233, 12),
        5: (487, 15),
        6: (914, 18),
    }
    profiles: dict[str, dict[str, int]] = {}
    for order, polynomial in enumerate(primitive_moments, start=2):
        profile = (len(polynomial.terms()), polynomial.total_degree())
        assert profile == expected_profiles[order]
        profiles[str(order)] = {
            "terms": profile[0],
            "total_degree": profile[1],
        }

    # U is always rank two.  On r=1 the odd W row is nonzero, and the
    # even W row is nonzero exactly on the two-open cover q!=0 or p!=0.
    fixed_u_minor = sp.det(u_matrix.extract((0, 3), (0, 1)))
    assert fixed_u_minor == -2
    p_numerator = b * q + d * s + 3
    rank_cover = {
        "q_nonzero": q,
        "p_nonzero": p_numerator,
    }
    variable_orders = (
        ("h", "b", "d", "q", "s"),
        ("h", "s", "q", "d", "b"),
    )
    semistable_results: dict[str, dict[str, list[dict[str, object]]]] = {}
    for open_name, rank_factor in rank_cover.items():
        cutoff_results: dict[str, list[dict[str, object]]] = {}
        for cutoff in (5, 6):
            runs = []
            rabinowitsch = sp.expand(h * rank_factor * theta - 1)
            for variable_order in variable_orders:
                source = msolve_source(
                    variable_order,
                    primitive_moments[: cutoff - 1],
                    rabinowitsch,
                )
                result = run_msolve(source, arguments.threads)
                result["variable_order"] = list(variable_order)
                runs.append(result)
            expected_status = (
                "positive_dimensional" if cutoff == 5 else "unit_ideal"
            )
            if any(run["status"] != expected_status for run in runs):
                raise RuntimeError(
                    "unexpected exact r=1 semistable status on "
                    f"{open_name} through mu_{cutoff}: {runs}"
                )
            cutoff_results[str(cutoff)] = runs
        semistable_results[open_name] = cutoff_results

    # The complementary projective chart r=0,s!=0 is normalized by s=1.
    # Its exact-rank-two locus is covered by q!=0 and p!=0.  Both opens are
    # zero-dimensional through mu_5 and empty after adjoining mu_6.
    r0_variables = (b, d, q)
    r0_substitution = {
        r: sp.Integer(0),
        s: sp.Integer(1),
        p: -(b * q + d) / 3,
    }
    assert sp.expand(raw_moments[0].subs(r0_substitution)) == 0
    r0_moments = [
        primitive_integer_polynomial(
            raw.subs(r0_substitution), r0_variables
        )[0]
        for raw in raw_moments[1:]
    ]
    r0_matrix = coefficient_matrix.subs(r0_substitution)
    r0_trace = sp.factor(
        sp.trace((r0_matrix * apolar_diagonal) ** 2)
    )
    assert r0_trace == 32 * d**2
    r0_rank_cover = {
        "q_nonzero": q,
        "p_nonzero": b * q + d,
    }
    r0_variable_orders = (
        ("h", "b", "d", "q"),
        ("h", "q", "d", "b"),
    )
    r0_results: dict[str, dict[str, list[dict[str, object]]]] = {}
    for open_name, rank_factor in r0_rank_cover.items():
        cutoff_results = {}
        for cutoff, expected_status in (
            (5, "zero_dimensional"),
            (6, "unit_ideal"),
        ):
            runs = []
            for variable_order in r0_variable_orders:
                source = msolve_source(
                    variable_order,
                    r0_moments[: cutoff - 1],
                    sp.expand(h * rank_factor - 1),
                )
                result = run_msolve(source, arguments.threads)
                result["variable_order"] = list(variable_order)
                runs.append(result)
            if any(run["status"] != expected_status for run in runs):
                raise RuntimeError(
                    "unexpected exact r=0 status on "
                    f"{open_name} through mu_{cutoff}: {runs}"
                )
            cutoff_results[str(cutoff)] = runs
        r0_results[open_name] = cutoff_results

    # On theta=0 one has d=-3/s and s!=0.  Exact rank two is precisely
    # q!=0, because p=-b*q/3 and U remains rank two.  Compute the complete
    # minimal-prime decomposition through mu_6 over QQ on that full open.
    theta_variables = (b, q, s)
    theta_substitution = {
        r: sp.Integer(1),
        d: -3 / s,
        p: -b * q / 3,
    }
    theta_moments: list[sp.Poly] = []
    for raw in raw_moments[1:]:
        numerator = sp.together(
            sp.cancel(raw.subs(theta_substitution))
        ).as_numer_denom()[0]
        theta_moments.append(
            primitive_integer_polynomial(numerator, theta_variables)[0]
        )
    singular_equations = [
        msolve_expression(moment.as_expr()) for moment in theta_moments
    ]
    singular_equations.append(msolve_expression(h * s * q - 1))
    singular_source = "\n".join(
        [
            "ring R=0,(h,b,q,s),dp;",
            "ideal I=" + ",\n".join(singular_equations) + ";",
            'LIB "primdec.lib";',
            "ideal G=std(I);",
            'print("DIM"); print(dim(G));',
            "list L=minAssGTZ(I);",
            'print("COUNT"); print(size(L));',
            "for (int i=1; i<=size(L); i++)",
            "{",
            '  print("COMPONENT"); print(i); print(std(L[i]));',
            "}",
            "quit;",
        ]
    )
    theta_decomposition = run_singular_min_ass(singular_source)
    normalized_decomposition = "".join(
        theta_decomposition["output"].split()
    )
    for signature in (
        "DIM1COUNT2",
        "s-1,b+1,hq-1",
        "s+3,b-3,3hq+1",
    ):
        if signature not in normalized_decomposition:
            raise RuntimeError(
                "unexpected theta=0 minimal-prime decomposition: "
                + theta_decomposition["output"]
            )

    plus_line = {
        r: 1,
        b: -1,
        d: -3,
        s: 1,
        p: q / 3,
    }
    minus_line = {
        r: 1,
        b: 3,
        d: 1,
        s: -3,
        p: -q,
    }
    for line in (plus_line, minus_line):
        for raw in raw_moments:
            assert sp.expand(raw.subs(line)) == 0

    # Exhibit both components in one-sided fixed-flag coordinates.  This
    # also certifies their relative periods and the function-field
    # recurrence without specializing q.
    capital_w, capital_v, capital_z, capital_y = sp.symbols("W V Z Y")
    wp, vp, zp, yp = sp.symbols("Wp Vp Zp Yp")
    homogeneous = sum(
        coefficient_matrix[i, j]
        * capital_w ** (3 - i)
        * capital_v**i
        * capital_z ** (3 - j)
        * capital_y**j
        for i in range(4)
        for j in range(4)
    )
    plus_transformed = sp.expand(
        homogeneous.subs(plus_line).subs(
            {
                capital_w: wp,
                capital_v: vp + wp,
                capital_z: zp - yp,
                capital_y: yp,
            }
        )
    )
    minus_transformed = sp.expand(
        homogeneous.subs(minus_line).subs(
            {
                capital_w: wp,
                capital_v: vp - wp,
                capital_z: zp + yp,
                capital_y: yp,
            }
        )
    )
    expected_plus = vp**2 * zp**2 * (
        q * (2 * wp + vp) * zp / 3 - vp * (zp - 2 * yp)
    )
    expected_minus = vp**2 * zp**2 * (
        -q * vp * (zp + 2 * yp) + (2 * wp - vp) * zp
    )
    assert sp.expand(plus_transformed - expected_plus) == 0
    assert sp.expand(minus_transformed - expected_minus) == 0

    plus_entries = {
        (2, 0): 2 * q / 3,
        (3, 0): q / 3 - 1,
        (3, 1): 2,
    }
    minus_entries = {
        (2, 0): 2,
        (3, 0): -q - 1,
        (3, 1): -2 * q,
    }
    u, t = sp.symbols("u t")

    def period_laurent(entries: dict[tuple[int, int], sp.Expr]) -> sp.Expr:
        return sp.expand(
            sum(
                value
                * u ** (target_index - dual_index)
                * t**target_index
                * (1 - t) ** (3 - target_index)
                for (dual_index, target_index), value in entries.items()
            )
        )

    mixed_values: dict[str, sp.Expr] = {}
    periods: dict[str, sp.Expr] = {}
    for name, entries, expected_mixed in (
        ("b=-1,s=1", plus_entries, 8 * (q + 3)),
        ("b=3,s=-3", minus_entries, 24 * (1 - q)),
    ):
        assert min(i - j for i, j in entries) == 2
        period = period_laurent(entries)
        periods[name] = period
        normalized_mixed = sp.integrate(
            sp.expand(u**2 * t**2 * period).coeff(u, 0),
            (t, 0, 1),
        )
        mixed_value = sp.factor(factorial(6) * normalized_mixed)
        assert sp.expand(mixed_value - expected_mixed) == 0
        assert mixed_value != 0
        mixed_values[name] = mixed_value

    payload = {
        "format": "two-pair-sic-bidegree33-rank-two-parity-channels-v2",
        "field": "QQ",
        "status": (
            "the complete exact-rank-two reversal-parity factor family is "
            "SIC-safe: centralizer orbit coverage reduces it to the "
            "normalized-U chart, its semistable opens are empty through "
            "mu_6, and its two all-order theta=0 components are one-sided"
        ),
        "centralizer_orbit_cover": {
            "matrix": "[[alpha,beta],[beta,alpha]]",
            "alpha": str(alpha),
            "beta": str(beta),
            "determinant": "alpha^2-beta^2=1",
            "eigen_coordinate_action": "X->lambda*X, Y->lambda^(-1)*Y",
            "even_endpoint": str(expected_even_endpoint),
            "odd_endpoint": str(expected_odd_endpoint),
            "even_bad_coefficient_matrix_determinant": -4,
            "odd_bad_coefficient_matrix_determinant": -4,
            "conclusion": (
                "for every nonzero even/odd U-channel pair, the two bad "
                "conditions are proper affine-linear equations in "
                "z=lambda^4; over an algebraic closure one nonzero z "
                "avoids both, after which diagonal internal gauge gives "
                "U[0,*]=(1,1)"
            ),
        },
        "family": {
            "coefficient_factorization": "C=U*W",
            "U": [[str(value) for value in row] for row in u_matrix.tolist()],
            "W": [[str(value) for value in row] for row in w_matrix.tolist()],
            "internal_gauge": (
                "reversal parity fixes the two channel lines; their two "
                "scalings are quotiented by U[0,*]=(1,1)"
            ),
            "factor_parameters": ["b", "d", "p", "q", "r", "s"],
            "projective_cover": ["r=1", "r=0,s=1", "r=s=0"],
            "rank_fact": "det(U[rows 0,3])=-2, so U is always rank two",
            "r_equals_s_equals_zero": "W has rank at most one",
        },
        "r_equals_one_chart": {
            "mu_1_elimination": "p=-(b*q+d*s+3)/3",
            "quotient_coordinates": ["b", "d", "q", "s"],
            "rank_two_open_cover": {
                "q_nonzero": "q!=0",
                "p_nonzero": "b*q+d*s+3!=0",
            },
            "selected_minor_times_3_for_comparison": str(delta),
            "semistability_factor": str(theta),
            "trace_L_squared": str(trace_square),
            "semistable_exact_schemes": semistable_results,
            "semistable_status": (
                "on both exact-rank-two opens the scheme is positive-"
                "dimensional through mu_5 and the unit ideal through mu_6"
            ),
        },
        "primitive_moment_profiles": profiles,
        "primitive_moment_scalings": scalings,
        "primitive_moments": {
            str(order): str(polynomial.as_expr())
            for order, polynomial in enumerate(primitive_moments, start=2)
        },
        "r_equals_zero_chart": {
            "normalization": "s=1",
            "mu_1_elimination": "p=-(b*q+d)/3",
            "rank_two_open_cover": {
                "q_nonzero": "q!=0",
                "p_nonzero": "b*q+d!=0",
            },
            "trace_L_squared": str(r0_trace),
            "exact_schemes": r0_results,
            "status": (
                "on both exact-rank-two opens the scheme is zero-dimensional "
                "through mu_5 and the unit ideal through mu_6"
            ),
        },
        "theta_equals_zero_components": {
            "substitution": "d=-3/s, p=-b*q/3, r=1",
            "full_exact_rank_two_localization": "s*q!=0",
            "ideal": "(mu_2,...,mu_6,h*s*q-1)",
            "dimension": 1,
            "minimal_primes": [
                "(s-1,b+1,h*q-1)",
                "(s+3,b-3,3*h*q+1)",
            ],
            "singular_certificate": theta_decomposition,
            "initial_vanishing": (
                "mu_1,...,mu_6 vanish identically in QQ(q) on both lines"
            ),
        },
        "all_order_component_certificate": {
            "normalized_moment": (
                "mu_m/(3*m+1)! = CT_u integral_0^1 P(u,t)^m dt"
            ),
            "fixed_flag_factorizations": {
                "b=-1,s=1": str(expected_plus),
                "b=3,s=-3": str(expected_minus),
            },
            "period_laurent_polynomials": {
                name: str(period) for name, period in periods.items()
            },
            "support_certificate": (
                "every transformed monomial has i-j>=2, so every monomial "
                "of P has u-degree at most -2"
            ),
            "creative_telescoping_recurrence": (
                "nu_(m+1)=0 for m>=0 over QQ(q), where "
                "nu_m=mu_m/(3*m+1)!; the valuation telescoper has forward "
                "coefficient 1 and no singular step"
            ),
            "initial_condition": "mu_1=0",
            "mixed_multiplier": (
                "in fixed-flag coordinates take e=2,a=0,b=2, whose period "
                "numerator is u^2*t^2"
            ),
            "mixed_values_at_m_equals_one": {
                name: str(value) for name, value in mixed_values.items()
            },
            "function_field_nonvanishing": (
                "8*(q+3) and 24*(1-q) are nonzero in QQ(q)"
            ),
            "mixed_tail": (
                "for every degree-e multiplier the mixed sequence vanishes "
                "when 2*m>e"
            ),
        },
        "scope": (
            "This classifies the complete exact-rank-two reversal-parity "
            "factor family up to pair-preserving SL2 orbit and internal "
            "channel gauge.  The full direct rank-two factor space remains "
            "unclassified."
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=2) + "\n")

    print("PASS centralizer orbit cover of every parity U chart")
    print("PASS complete exact QQ parity-channel rank-two family")
    print("semistable r=1 and r=0 charts: unit ideals through mu_6")
    print("theta=0: exactly two one-sided QQ(q) components through mu_6")
    print("all-order recurrence: nu_(m+1)=0; mixed tail: 2*m>e")
    print(arguments.output)


if __name__ == "__main__":
    main()
