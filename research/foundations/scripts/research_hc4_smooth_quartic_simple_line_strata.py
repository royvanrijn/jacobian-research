#!/usr/bin/env python3
"""Explore explicit strata of the first HC4NHM16 exceptional divisor."""

from __future__ import annotations

import argparse
import subprocess

import sympy as sp

from research_hc4_smooth_quartic_simple_line import (
    build_equations,
    singular,
    unknown_degree,
)


def stratum_substitution(name: str) -> dict[sp.Symbol, sp.Expr]:
    tau, p, q, r = sp.symbols("tau p q r")
    c, k, m = sp.symbols("c k m")
    substitutions = {
        "h2-zero": {p: 0, q: 0, r: 0},
        "tau0-delta-p-nonzero": {tau: 0, p: 3 * c, q: c * m**2, r: c * m},
        "tau0-delta-p-zero": {tau: 0, p: 0, r: 0},
        "tau0-delta-m3": {tau: 0, p: 3 * c, q: c * m**2, r: c * m},
        "h2-zero-tau-cubic": {p: 0, q: 0, r: 0},
        "taum1-p-equals-r": {tau: -1, r: p},
        "taum1-p-equals-r-quadratic": {tau: -1, r: p, q: k * p},
        "taum1-linear-factor": {tau: -1, q: -3 * p - 3 * r},
        "taum1-linear-pivot-zero": {
            tau: -1,
            r: sp.Rational(7, 33) * p,
            q: -sp.Rational(40, 11) * p,
        },
    }
    return substitutions[name]


def stratum_relation(name: str) -> sp.Expr | None:
    tau, k, m = sp.symbols("tau k m")
    relations = {
        "tau0-delta-m3": m**3 - 48,
        "h2-zero-tau-cubic": 6 * tau**3 + 1,
        "taum1-p-equals-r-quadratic": k**2 + 3 * k + 8,
    }
    return relations.get(name)


def deformation_substitution(name: str) -> dict[sp.Symbol, sp.Expr]:
    if name != "tau0-delta-m3":
        return {}
    c, m = sp.symbols("c m")
    bs = sp.symbols("b0:18")
    u, v = sp.symbols("u v")
    b10 = bs[10]
    b13 = (18 * c * m * v - 36 * c * u) / 5
    b12 = 2 * c * m**2 * v - 2 * c * m * u
    b9 = -m * b10
    b7 = c * (-3 * b10 * m**2 + 2 * m**2 * u - 6 * v)
    return {
        bs[0]: 0,
        bs[1]: m * b10,
        bs[3]: -3 * b10 + 2 * u,
        bs[4]: -2 * v,
        bs[6]: -24 * c * b10 - sp.Rational(2, 3) * b13 + 6 * c * u,
        bs[7]: b7,
        bs[9]: b9,
        bs[12]: b12,
        bs[13]: b13,
    }


def specialize_tau0_delta_m3(
    equations: list[sp.Expr],
    substitution: dict[sp.Symbol, sp.Expr],
    relation: sp.Expr,
    active_substitution: dict[sp.Symbol, sp.Expr],
) -> list[sp.Expr]:
    relation_symbol = next(iter(relation.free_symbols))
    specialized: list[sp.Expr] = []
    seen: set[str] = set()
    for equation in equations:
        value = sp.expand(equation.subs(substitution).subs(active_substitution))
        if value == 0:
            continue
        value = sp.rem(
            sp.Poly(value, relation_symbol), sp.Poly(relation, relation_symbol)
        ).as_expr()
        if value == 0:
            continue
        value_symbols = sorted(value.free_symbols, key=str)
        _, cleared = sp.Poly(value, *value_symbols).clear_denoms()
        value = sp.expand(cleared.as_expr())
        key = sp.srepr(value)
        if key not in seen:
            seen.add(key)
            specialized.append(value)
    return specialized


def build_tau0_delta_m3_program(
    equations: list[sp.Expr],
    unknowns: tuple[sp.Symbol, ...],
    stage: str,
) -> tuple[str, int, dict[int, int]]:
    substitution = stratum_substitution("tau0-delta-m3")
    relation = stratum_relation("tau0-delta-m3")
    assert relation is not None
    active_substitution = deformation_substitution("tau0-delta-m3")
    raw = specialize_tau0_delta_m3(equations, substitution, relation, {})
    specialized = specialize_tau0_delta_m3(
        equations, substitution, relation, active_substitution
    )
    active_unknowns = tuple(
        symbol for symbol in unknowns if symbol not in active_substitution
    )
    degree_counts: dict[int, int] = {}
    for equation in specialized:
        degree = unknown_degree(equation, active_unknowns)
        degree_counts[degree] = degree_counts.get(degree, 0) + 1

    relation_symbol = next(iter(relation.free_symbols))
    parameter_set: set[sp.Symbol] = set()
    for equation in (*raw, *specialized):
        parameter_set.update(equation.free_symbols - set(unknowns) - {relation_symbol})
    parameters = sorted(parameter_set, key=str)
    field = ",".join(map(str, parameters))
    variables = ",".join(map(str, (relation_symbol, *unknowns)))
    raw_generators = ",\n".join(singular(equation) for equation in raw)
    linear = [equation for equation in raw if unknown_degree(equation, unknowns) <= 1]
    linear_generators = ",\n".join(singular(equation) for equation in linear)
    substitution_targets: list[sp.Expr] = []
    for symbol, expression in active_substitution.items():
        value = sp.expand(symbol - expression)
        value_symbols = sorted(value.free_symbols, key=str)
        _, cleared = sp.Poly(value, *value_symbols).clear_denoms()
        substitution_targets.append(sp.expand(cleared.as_expr()))
    target_generators = ",\n".join(
        singular(equation) for equation in substitution_targets
    )

    lines = [
        "option(redSB);",
        f"ring ambient=(0,{field}),({variables}),dp;",
        f"ideal quotient_relation={singular(relation)};",
        "qring rr=std(quotient_relation);",
        f"ideal Iraw={raw_generators};",
        f"ideal L={linear_generators};",
        "timer=1;",
        "ideal GL=std(L);",
        'print("LINEAR_BEGIN"); print(GL); print("LINEAR_END");',
        f"ideal T={target_generators};",
        "ideal RT=reduce(T,GL);",
        'print("SUBSTITUTION_REMAINDER_BEGIN"); print(RT); '
        'print("SUBSTITUTION_REMAINDER_END");',
    ]
    if stage == "full":
        active_set = set(active_unknowns)
        core_set = set(sp.symbols("b10 u v"))
        core = [
            equation
            for equation in specialized
            if equation.free_symbols & active_set <= core_set
        ]
        specialized_generators = ",\n".join(
            singular(equation) for equation in specialized
        )
        core_generators = ",\n".join(singular(equation) for equation in core)
        support_generators = ",".join(map(str, unknowns))
        lines.extend(
            [
                f"ideal I={specialized_generators};",
                f"ideal Icore={core_generators};",
                "ideal Gcore=std(Icore);",
                'print("CORE_BEGIN"); print(Gcore); print("CORE_END");',
                "ideal J=reduce(I,Gcore);",
                "ideal H=std(J);",
                'print("REDUCED_BEGIN"); print(H); print("REDUCED_END");',
                "ideal G=std(GL+Gcore+H);",
                f"ideal S={support_generators};",
                "ideal RS=reduce(S,G);",
                'print("SUPPORT_REMAINDER_BEGIN"); print(RS); '
                'print("SUPPORT_REMAINDER_END");',
                'print("core_size="+string(size(Gcore)));',
                'print("reduced_size="+string(size(H)));',
            ]
        )
    lines.extend(['print("elapsed_ticks="+string(timer));', 'print("DONE");', "quit;"])
    return "\n".join(lines) + "\n", len(specialized), degree_counts


def build_program(stratum: str, stage: str) -> tuple[str, int, dict[int, int]]:
    equations, unknowns, _ = build_equations("squarefree-line", False)
    assert len(equations) == 81
    if stratum == "tau0-delta-m3":
        return build_tau0_delta_m3_program(equations, unknowns, stage)
    substitution = stratum_substitution(stratum)
    relation = stratum_relation(stratum)
    active_substitution = deformation_substitution(stratum)
    active_unknowns = tuple(symbol for symbol in unknowns if symbol not in active_substitution)
    specialized: list[sp.Expr] = []
    seen: set[str] = set()
    for equation in equations:
        value = sp.expand(equation.subs(substitution).subs(active_substitution))
        if value == 0:
            continue
        if relation is not None and stratum == "tau0-delta-m3":
            relation_symbol = next(iter(relation.free_symbols))
            value = sp.rem(
                sp.Poly(sp.expand(value), relation_symbol),
                sp.Poly(relation, relation_symbol),
            ).as_expr()
            if value == 0:
                continue
        # Clear rational scalar denominators.  In a Singular parameter field,
        # an unparenthesized term such as sigma^2/33 can otherwise be parsed
        # as a nonintegral power of the coefficient-field element sigma.
        value_symbols = sorted(value.free_symbols, key=str)
        if value_symbols:
            _, cleared = sp.Poly(value, *value_symbols).clear_denoms()
            value = sp.expand(cleared.as_expr())
        key = sp.srepr(value)
        if key not in seen:
            seen.add(key)
            specialized.append(value)

    unknowns = active_unknowns
    unknown_set = set(unknowns)
    parameter_set: set[sp.Symbol] = set()
    for equation in specialized:
        parameter_set.update(equation.free_symbols - unknown_set)
    relation_symbol: sp.Symbol | None = None
    if relation is not None:
        relation_symbols = sorted(relation.free_symbols, key=str)
        assert len(relation_symbols) == 1
        relation_symbol = relation_symbols[0]
        parameter_set.remove(relation_symbol)
    parameters = sorted(parameter_set, key=str)
    degree_counts: dict[int, int] = {}
    for equation in specialized:
        degree = unknown_degree(equation, unknowns)
        degree_counts[degree] = degree_counts.get(degree, 0) + 1

    field = ",".join(map(str, parameters))
    ring_variables = unknowns if relation_symbol is None else (relation_symbol, *unknowns)
    variables = ",".join(map(str, ring_variables))
    generators = ",\n".join(singular(equation) for equation in specialized)
    linear = [equation for equation in specialized if unknown_degree(equation, unknowns) <= 1]
    linear_generators = ",\n".join(singular(equation) for equation in linear) or "0"
    lines = ["option(redSB);"]
    if relation is None:
        lines.append(f"ring rr=(0,{field}),({variables}),dp;")
    else:
        lines.extend(
            [
                f"ring ambient=(0,{field}),({variables}),dp;",
                f"ideal quotient_relation={singular(relation)};",
                "qring rr=std(quotient_relation);",
            ]
        )
    lines.extend(
        [
        f"ideal I={generators};",
        f"ideal L={linear_generators};",
        "timer=1;",
        "ideal GL=std(L);",
        'print("LINEAR_BEGIN"); print(GL); print("LINEAR_END");',
        ]
    )
    if stage == "full":
        basis_algorithm = "std" if relation is not None else "slimgb"
        support_targets = [
            symbol**2
            if stratum == "taum1-p-equals-r-quadratic" and str(symbol) == "b12"
            else symbol
            for symbol in unknowns
        ]
        support_generators = ",".join(singular(target) for target in support_targets)
        lines.extend(
            [
                "ideal J=reduce(I,GL);",
                f"ideal H={basis_algorithm}(J);",
                'print("REDUCED_BEGIN"); print(H); print("REDUCED_END");',
                "ideal G=std(GL+H);",
                f"ideal S={support_generators};",
                "ideal RS=reduce(S,G);",
                'print("SUPPORT_REMAINDER_BEGIN"); print(RS); '
                'print("SUPPORT_REMAINDER_END");',
                'print("linear_size="+string(size(GL)));',
                'print("reduced_size="+string(size(H)));',
                'print("reduced_dimension="+string(dim(H)));',
            ]
        )
    lines.extend(['print("elapsed_ticks="+string(timer));', 'print("DONE");', "quit;"])
    return "\n".join(lines) + "\n", len(specialized), degree_counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stratum",
        choices=(
            "h2-zero",
            "tau0-delta-p-nonzero",
            "tau0-delta-p-zero",
            "tau0-delta-m3",
            "h2-zero-tau-cubic",
            "taum1-p-equals-r",
            "taum1-p-equals-r-quadratic",
            "taum1-linear-factor",
            "taum1-linear-pivot-zero",
        ),
        required=True,
    )
    parser.add_argument("--stage", choices=("linear", "full"), default="linear")
    parser.add_argument("--timeout", type=int, default=300)
    args = parser.parse_args()

    program, equation_count, degree_counts = build_program(args.stratum, args.stage)
    result = subprocess.run(
        ["Singular", "--no-tty", "--quiet"],
        input=program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=args.timeout,
        check=False,
    )
    print(f"stratum={args.stratum}")
    print(f"equation_count={equation_count}")
    print(f"unknown_degree_counts={degree_counts}")
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
