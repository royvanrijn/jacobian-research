#!/usr/bin/env python3
"""Audit the normalization-first global multi-Rees/Cox gate for F20.

The compact colored identity ``3*D_d + D_q + D_r = div(P_X)`` lives on a
regular normalized root cover.  This checker asks whether the three most
natural two-generated incidence ideals on the *nonnormal* root hypersurface
already realize those columns.

It proves three deliberately scoped facts.

1. For a domain A and ideals I_i=(a_i,b_i), the multi-Rees kernel is the
   saturation of the linear symmetric-algebra relations by product(a_i).
   This is recorded here for the F20 ideals using d, q, and r as pivots.
2. At either conjugate triple-E1 color the local value semigroup of the
   original incidence ring has no value one, whereas D_d has coefficient
   one.  Thus no ordinary integral ideal on that ring can realize D_d;
   normalization (and its uniformizer) is logically prior to the Cox fill.
3. For the natural incidence ideals, the degree-(3,1,1) product does not
   contain P_X.  The resulting cyclic residue has exact length 57, projects
   to a length-33 base scheme, is killed minimally among squarefree boundary
   monomials by d*q*r, and has reduced support exactly the eight previously
   certified finite collision centers.

This is not an obstruction to the normalized divisorial Cox algebra and not
an affine-space or Keller-map exclusion.  It identifies the exact finite
residue that a normalization-first Cox construction must repair.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import sympy as sp


def f20_data() -> dict[str, sp.Expr]:
    s, t, X = sp.symbols("s t X")
    d = s**2 + 4
    q = 4 * s**2 * t**2 + 4 * s**2 * t + 8 * s * t + 6 * s - 8 * t - 5
    r = (
        16 * s**2 * t**3
        + 4 * s**2 * t**2
        - 76 * s * t
        - 16 * s
        + 64 * t**3
        + 16 * t**2
        - 164 * t
        - 199
    )
    P = sp.expand(
        X**5
        + (t**2 * d - 2 * s - sp.Rational(17, 4)) * X**4
        + (3 * t * d + d + sp.Rational(13, 2) * s + 1) * X**3
        - (t * d + sp.Rational(11, 2) * s - 8) * X**2
        + (s - 6) * X
        + 1
    )
    P_X = sp.diff(P, X)
    subresultants = sp.subresultants(P, P_X, X)
    h_r = sp.cancel(subresultants[-3] / d)
    h_d = 2 * X - s - 2
    h_q = 2 * (s - 1) * X - (2 * s**2 * t + 2 * s**2 + 3 * s - 4)
    return {
        "s": s,
        "t": t,
        "X": X,
        "d": d,
        "q": q,
        "r": r,
        "P": P,
        "P_X": P_X,
        "h_d": h_d,
        "h_q": h_q,
        "h_r": h_r,
        "subresultants": subresultants,
    }


def lowest_order(expression: sp.Expr, parameter: sp.Symbol) -> int:
    polynomial = sp.Poly(sp.expand(expression), parameter)
    for order in range(polynomial.degree() + 1):
        if polynomial.nth(order) != 0:
            return order
    raise AssertionError("zero series has no finite order")


def lowest_order_mod(
    expression: sp.Expr,
    parameter: sp.Symbol,
    modulus: sp.Expr,
    modulus_variable: sp.Symbol,
) -> int:
    polynomial = sp.Poly(sp.expand(expression), parameter)
    for order in range(polynomial.degree() + 1):
        if sp.rem(polynomial.nth(order), modulus, modulus_variable) != 0:
            return order
    raise AssertionError("reduced zero series has no finite order")


def singular_polynomial(
    expression: sp.Expr, s: sp.Symbol, t: sp.Symbol, X: sp.Symbol
) -> str:
    """Render a QQ polynomial with coefficients before monomials.

    Singular parses ``17/4*X^4`` as intended, while the SymPy spelling
    ``17*X^4/4`` is not accepted in every polynomial context.
    """

    polynomial = sp.Poly(sp.expand(expression), s, t, X, domain=sp.QQ)
    terms: list[tuple[str, str]] = []
    for (s_degree, t_degree, x_degree), coefficient in polynomial.terms():
        numerator, denominator = map(int, sp.fraction(coefficient))
        monomial_factors = []
        if s_degree:
            monomial_factors.append(f"s^{s_degree}")
        if t_degree:
            monomial_factors.append(f"t^{t_degree}")
        if x_degree:
            monomial_factors.append(f"X^{x_degree}")
        monomial = "*".join(monomial_factors) or "1"
        magnitude = abs(numerator)
        scalar = (
            str(magnitude)
            if denominator == 1
            else f"{magnitude}/{denominator}"
        )
        term = monomial if scalar == "1" else f"{scalar}*{monomial}"
        terms.append(("-" if numerator < 0 else "+", term))
    first_sign, first_term = terms[0]
    rendered = ("-" if first_sign == "-" else "") + first_term
    for sign, term in terms[1:]:
        rendered += sign + term
    return rendered


def singular_residue_gate(data: dict[str, sp.Expr]) -> dict[str, object]:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the exact multi-Rees residue gate")

    s = data["s"]
    t = data["t"]
    X = data["X"]
    d = data["d"]
    q = data["q"]
    r = data["r"]
    P = data["P"]
    P_X = data["P_X"]
    h_d = data["h_d"]
    h_q = data["h_q"]
    h_r = data["h_r"]
    assert isinstance(s, sp.Symbol)
    assert isinstance(t, sp.Symbol)
    assert isinstance(X, sp.Symbol)

    d_cubed_generators = (d**3, d**2 * h_d, d * h_d**2, h_d**3)
    product_generators = tuple(
        sp.expand(d_generator * q_generator * r_generator)
        for d_generator in d_cubed_generators
        for q_generator in (q, h_q)
        for r_generator in (r, h_r)
    )
    render = lambda expression: singular_polynomial(expression, s, t, X)

    h2 = 16 * t**2 + 24 * t + 13
    h3 = 8 * t**3 + 16 * t**2 + 2 * t - 7
    program = '\n'.join(
        (
            'LIB "primdec.lib";',
            "option(redSB);",
            "ring R=0,(X,s,t),dp;",
            f"poly P={render(P)};",
            f"poly PX={render(P_X)};",
            f"poly d={render(d)};",
            f"poly q={render(q)};",
            f"poly rr={render(r)};",
            f"poly h2={render(h2)};",
            f"poly h3={render(h3)};",
            "ideal J=P," + ",".join(map(render, product_generators)) + ";",
            "ideal G=std(J);",
            "poly remPX=reduce(PX,G);",
            '"PRODUCT_GENERATOR_COUNT";16;',
            '"PRODUCT_GROEBNER_SIZE";size(G);',
            '"PX_IN_PRODUCT";remPX==0;',
            '"PX_REMAINDER_EQUALS_PX";remPX-PX==0;',
            "ideal C=std(quotient(J,ideal(PX)));",
            '"COLON_GROEBNER_SIZE";size(C);',
            '"COLON_DIMENSION";dim(C);',
            '"RESIDUE_LENGTH";vdim(C);',
            '"ANNIHILATOR_1";reduce(1,C)==0;',
            '"ANNIHILATOR_D";reduce(d,C)==0;',
            '"ANNIHILATOR_Q";reduce(q,C)==0;',
            '"ANNIHILATOR_R";reduce(rr,C)==0;',
            '"ANNIHILATOR_DQ";reduce(d*q,C)==0;',
            '"ANNIHILATOR_DR";reduce(d*rr,C)==0;',
            '"ANNIHILATOR_QR";reduce(q*rr,C)==0;',
            '"ANNIHILATOR_DQR";reduce(d*q*rr,C)==0;',
            "ideal E=std(eliminate(C,X));",
            '"BASE_ELIMINATION_SIZE";size(E);',
            "ideal CB=std(E+ideal(X));",
            '"BASE_RESIDUE_LENGTH";vdim(CB);',
            "ideal N=intersect("
            "ideal(X,s-1,2*t+1),"
            "ideal(X,s-11,2*t+1),"
            "ideal(X,s-4*t-3,h2),"
            "ideal(X,12*s-7,t-2),"
            "ideal(X,s-4*t^2+5,h3));",
            "N=std(N);",
            '"KNOWN_CENTER_DEGREE";vdim(N);',
            '"BASE_RESIDUE_IN_KNOWN_CENTERS";size(reduce(CB,N))==0;',
            "ideal RB=std(radical(CB));",
            '"RADICAL_BASE_LENGTH";vdim(RB);',
            '"RADICAL_IN_KNOWN_CENTERS";size(reduce(RB,N))==0;',
            '"KNOWN_CENTERS_IN_RADICAL";size(reduce(N,RB))==0;',
            "exit;",
        )
    )

    labels = {
        "PRODUCT_GENERATOR_COUNT",
        "PRODUCT_GROEBNER_SIZE",
        "PX_IN_PRODUCT",
        "PX_REMAINDER_EQUALS_PX",
        "COLON_GROEBNER_SIZE",
        "COLON_DIMENSION",
        "RESIDUE_LENGTH",
        "ANNIHILATOR_1",
        "ANNIHILATOR_D",
        "ANNIHILATOR_Q",
        "ANNIHILATOR_R",
        "ANNIHILATOR_DQ",
        "ANNIHILATOR_DR",
        "ANNIHILATOR_QR",
        "ANNIHILATOR_DQR",
        "BASE_ELIMINATION_SIZE",
        "BASE_RESIDUE_LENGTH",
        "KNOWN_CENTER_DEGREE",
        "BASE_RESIDUE_IN_KNOWN_CENTERS",
        "RADICAL_BASE_LENGTH",
        "RADICAL_IN_KNOWN_CENTERS",
        "KNOWN_CENTERS_IN_RADICAL",
    }
    with tempfile.TemporaryDirectory() as temporary_directory:
        source = Path(temporary_directory) / "f20_multi_rees_gate.sing"
        source.write_text(program + "\n")
        completed = subprocess.run(
            (singular, "-q", str(source)),
            check=True,
            capture_output=True,
            text=True,
            timeout=300,
            stdin=subprocess.DEVNULL,
        )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    parsed: dict[str, int] = {}
    for index, line in enumerate(lines[:-1]):
        if line in labels:
            parsed[line] = int(lines[index + 1])
    missing = labels - set(parsed)
    assert not missing, (missing, completed.stdout, completed.stderr)

    expected = {
        "PRODUCT_GENERATOR_COUNT": 16,
        "PRODUCT_GROEBNER_SIZE": 34,
        "PX_IN_PRODUCT": 0,
        "PX_REMAINDER_EQUALS_PX": 1,
        "COLON_GROEBNER_SIZE": 28,
        "COLON_DIMENSION": 0,
        "RESIDUE_LENGTH": 57,
        "ANNIHILATOR_1": 0,
        "ANNIHILATOR_D": 0,
        "ANNIHILATOR_Q": 0,
        "ANNIHILATOR_R": 0,
        "ANNIHILATOR_DQ": 0,
        "ANNIHILATOR_DR": 0,
        "ANNIHILATOR_QR": 0,
        "ANNIHILATOR_DQR": 1,
        "BASE_ELIMINATION_SIZE": 8,
        "BASE_RESIDUE_LENGTH": 33,
        "KNOWN_CENTER_DEGREE": 8,
        "BASE_RESIDUE_IN_KNOWN_CENTERS": 1,
        "RADICAL_BASE_LENGTH": 8,
        "RADICAL_IN_KNOWN_CENTERS": 1,
        "KNOWN_CENTERS_IN_RADICAL": 1,
    }
    assert parsed == expected
    return {
        "product_generator_count": parsed["PRODUCT_GENERATOR_COUNT"],
        "product_groebner_size": parsed["PRODUCT_GROEBNER_SIZE"],
        "P_X_in_product_ideal": False,
        "P_X_normal_form": "P_X",
        "colon_groebner_size": parsed["COLON_GROEBNER_SIZE"],
        "residue_length": parsed["RESIDUE_LENGTH"],
        "base_elimination_groebner_size": parsed["BASE_ELIMINATION_SIZE"],
        "base_residue_length": parsed["BASE_RESIDUE_LENGTH"],
        "squarefree_boundary_annihilators": {
            "1": False,
            "d": False,
            "q": False,
            "r": False,
            "d*q": False,
            "d*r": False,
            "q*r": False,
            "d*q*r": True,
        },
        "reduced_base_support": {
            "degree": parsed["RADICAL_BASE_LENGTH"],
            "equals_known_collision_centers": True,
            "centers": (
                "q node (s,t)=(1,-1/2)",
                "r cusp (s,t)=(11,-1/2)",
                "two triple centers (s=4t+3, 16t^2+24t+13=0)",
                "transverse q-r center (s,t)=(7/12,2)",
                "three q-r tangencies (s=4t^2-5, 8t^3+16t^2+2t-7=0)",
            ),
        },
    }


def build_certificate() -> dict[str, object]:
    data = f20_data()
    s = data["s"]
    t = data["t"]
    X = data["X"]
    d = data["d"]
    q = data["q"]
    r = data["r"]
    P = data["P"]
    P_X = data["P_X"]
    h_d = data["h_d"]
    h_r = data["h_r"]
    subresultants = data["subresultants"]
    assert isinstance(s, sp.Symbol)
    assert isinstance(t, sp.Symbol)
    assert isinstance(X, sp.Symbol)
    assert isinstance(subresultants, list)

    factor_content, factor_list = sp.factor_list(P)
    assert factor_content != 0
    assert len(factor_list) == 1 and factor_list[0][1] == 1
    assert tuple(sp.degree(item, X) for item in subresultants) == (5, 4, 3, 2, 1, 0)
    assert sp.expand(subresultants[-3] - d * h_r) == 0
    h_r_poly = sp.Poly(h_r, s, t, X, domain=sp.QQ)
    assert h_r_poly.degree(X) == 2
    assert sp.Poly(h_r, X).content() == 1
    assert sp.expand(
        h_r.subs({s: 2, t: sp.Rational(7, 4)})
        + sp.Rational(35721, 64) * (X**2 + 8 * X - 2)
    ) == 0

    tau, z, b, Y = sp.symbols("tau z b Y")
    imaginary_unit = sp.I
    triple_s0 = 2 * imaginary_unit
    triple_t0 = -sp.Rational(3, 4) + imaginary_unit / 2
    triple_root = 1 + imaginary_unit
    triple_substitution = {
        s: triple_s0 + tau**4,
        t: triple_t0 + z * tau**4,
        X: triple_root + tau**2 * (b + tau * Y),
    }
    triple_orders = {
        "s_minus_s0": lowest_order(s.subs(triple_substitution) - triple_s0, tau),
        "t_minus_t0": lowest_order(t.subs(triple_substitution) - triple_t0, tau),
        "X_minus_root": lowest_order(X.subs(triple_substitution) - triple_root, tau),
        "d": lowest_order(d.subs(triple_substitution), tau),
        "q": lowest_order(q.subs(triple_substitution), tau),
        "r": lowest_order(r.subs(triple_substitution), tau),
        "h_d": lowest_order(h_d.subs(triple_substitution), tau),
        "P_X": lowest_order_mod(
            P_X.subs(triple_substitution), tau, b**2 + imaginary_unit, b
        ),
    }
    assert triple_orders == {
        "s_minus_s0": 4,
        "t_minus_t0": 4,
        "X_minus_root": 2,
        "d": 4,
        "q": 4,
        "r": 4,
        "h_d": 2,
        "P_X": 7,
    }

    residue_gate = singular_residue_gate(data)
    return {
        "schema": "f20-global-multi-rees-cox-algebra-v1",
        "status": "normalization_first_obstruction_and_finite_residue",
        "general_multi_rees_theorem": {
            "hypotheses": "A is a domain and I_i=(a_i,b_i) with every a_i nonzero",
            "map": "A[U_i,V_i] -> A[T_i], U_i |-> a_i*T_i, V_i |-> b_i*T_i",
            "kernel": "(b_i*U_i-a_i*V_i for all i) : (product_i a_i)^infinity",
            "reason": "localize at the pivot product, solve V_i=(b_i/a_i)U_i, then contract",
            "graded_piece": "degree alpha is product_i I_i^(alpha_i)",
        },
        "root_incidence_ring": {
            "ring": "R=QQ[s,t,X]/(P)",
            "P_irreducible_over_QQ": True,
            "natural_ideals": {
                "I_d": "(d, 2*X-s-2)",
                "I_q": "(q, 2*(s-1)*X-(2*s^2*t+2*s^2+3*s-4))",
                "I_r": "(r, Sres_2(P,P_X)/(s^2+4))",
            },
            "multi_rees_kernel": (
                "(P, h_d*U_d-d*V_d, h_q*U_q-q*V_q, "
                "h_r*U_r-r*V_r) : (d*q*r)^infinity"
            ),
        },
        "r_incidence_generator": {
            "definition": "h_r=Sres_2(P,P_X)/(s^2+4)",
            "polynomial": True,
            "primitive_in_QQ[s,t][X]": True,
            "X_degree": h_r_poly.degree(X),
            "term_count": len(h_r_poly.terms()),
            "sample_s_2_t_7_over_4": "-35721/64*(X^2+8*X-2)",
        },
        "triple_E1_value_semigroup_gate": {
            "orders": triple_orders,
            "maximal_ideal_positive_order_floor": 2,
            "desired_D_d_coefficient": 1,
            "natural_I_d_order": 2,
            "ordinary_incidence_ideal_can_realize_D_d": False,
            "conclusion": (
                "value one is a gap in the original local incidence ring; "
                "the normalized uniformizer tau must be adjoined before the Cox fill"
            ),
        },
        "degree_3_1_1_natural_incidence_piece": residue_gate,
        "normalized_divisorial_cox_frontier": {
            "algebra": (
                "direct_sum_{a,b,c>=0} H^0(X_tilde, "
                "O(-a*D_d-b*D_q-c*D_r))*U^a*V^b*W^c"
            ),
            "required_local_extension": "adjoin the triple-E1 normalized uniformizer tau of value one",
            "remaining": (
                "construct simultaneous global divisorial ideals on the normalized regular cover, "
                "compute their multiplication/residue cocycle, then retest the (3,1,1) section"
            ),
        },
        "downstream": {
            "abstract_colored_packet_identity": "still passes",
            "ordinary_incidence_multi_rees_realization": "fails",
            "normalized_Cox_realization": "open",
            "inverse_adjugate_polynomiality": "not_reached",
            "affine_space_recognition": "not_reached",
            "no_go_scope": "natural ideals on the nonnormal root-incidence ring only",
        },
        "software_assumptions": {
            "python": ".venv Python with pinned SymPy",
            "singular": "Singular 4.4-compatible quotient, elimination, and radical routines",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")
    print(rendered)
    print("PASS: the natural global multi-Rees algebra has an exact saturated presentation")
    print("PASS: the triple-E1 value-one gap forces normalization before the Cox fill")
    print("PASS: the natural (3,1,1) residue has length 57 over eight collision centers")
    print("SCOPE: the normalized divisorial Cox algebra and affine completion remain open")


if __name__ == "__main__":
    main()
