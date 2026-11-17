#!/usr/bin/env python3
"""Verify the global F20 normalization and its Cox-conductor frontier.

This is the normalization-first continuation of
``verify_f20_global_multi_rees_cox_algebra.py``.  It does not duplicate the
F20 polynomial or polynomial renderer; both are imported from that checker.

The exact calculation:

* computes the conductor and the finite R-module structure of the global
  normalization of R=QQ[s,t,X]/(P);
* decomposes the previously found length-57 incidence residue and its
  length-33 base projection by collision packet;
* identifies the conductor square with the rational q-crossing double cover
  w -> y and factors d, r, and the derivative slope residue on that cover;
* proves that the two integral normalization generators still have order two
  at triple E1, so normalization alone supplies no value-one Cox generator;
* retests the natural degree-(3,1,1) product after normalization and proves
  that it still does not contain P_X.

Thus the next required algebra is the exceptional/controlled-transform Cox
algebra of the resolved normalization.  No inverse-adjugate or affine-space
claim is made here.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import sympy as sp

from verify_f20_global_multi_rees_cox_algebra import (
    f20_data,
    lowest_order_mod,
    singular_polynomial,
)


CENTER_NAMES = {
    1: "q_node",
    2: "r_cusp",
    3: "triple_orbit",
    4: "q_r_transverse",
    5: "q_r_tangency_orbit",
}


def conductor_polynomials(data: dict[str, sp.Expr]) -> tuple[sp.Expr, ...]:
    s = data["s"]
    t = data["t"]
    X = data["X"]
    return (
        4 * X * t - 2 * s * t + 2 * X - 2 * s + 1,
        X**2 - X * s - 2 * X + s,
        2 * s**2 * t - 2 * X * s + 2 * s**2 + 2 * X + 3 * s - 4,
    )


def parse_singular_output(stdout: str) -> dict[str, object]:
    lines = [
        line.strip()
        for line in stdout.splitlines()
        if line.strip() and not line.lstrip().startswith("//")
    ]
    scalar_labels = {
        "CONDUCTOR_EQUAL",
        "NORMAL_MODULE_EQUAL",
        "NORMAL_VARIABLE_COUNT",
        "NORMAL_RELATION_COUNT",
        "NORMAL_PRODUCT_GROEBNER_SIZE",
        "PX_IN_NORMALIZED_PRODUCT",
    }
    result: dict[str, object] = {}
    root_components: list[tuple[int, int, int]] = []
    base_components: list[tuple[int, int, int]] = []
    index = 0
    while index < len(lines):
        line = lines[index]
        if line in scalar_labels:
            result[line] = int(lines[index + 1])
            index += 2
            continue
        if line == "ROOT_COMPONENT":
            root_components.append(
                (int(lines[index + 1]), int(lines[index + 2]), int(lines[index + 3]))
            )
            index += 4
            continue
        if line == "BASE_COMPONENT":
            base_components.append(
                (int(lines[index + 1]), int(lines[index + 2]), int(lines[index + 3]))
            )
            index += 4
            continue
        index += 1
    result["root_components"] = root_components
    result["base_components"] = base_components
    return result


def singular_normalization_gate(
    data: dict[str, sp.Expr], conductor: tuple[sp.Expr, ...]
) -> dict[str, object]:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the F20 normalization gate")

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
    render = lambda expression: singular_polynomial(expression, s, t, X)

    product_generators = tuple(
        sp.expand(d_generator * q_generator * r_generator)
        for d_generator in (d**3, d**2 * h_d, d * h_d**2, h_d**3)
        for q_generator in (q, h_q)
        for r_generator in (r, h_r)
    )
    c1, c2, c3 = conductor
    module_generators = (sp.expand(d * c1), sp.expand(d * c3), c2)

    common_lines = [
        "option(redSB);",
        "ring R=0,(X,s,t),dp;",
        "proc eq(ideal A,ideal B)"
        "{A=std(A);B=std(B);return(size(reduce(A,B))==0 && size(reduce(B,A))==0);}",
        f"poly P={render(P)};",
        f"poly PX={render(P_X)};",
        f"poly d={render(d)};",
        f"poly q={render(q)};",
        f"poly rr={render(r)};",
        f"poly c1={render(c1)};",
        f"poly c2={render(c2)};",
        f"poly c3={render(c3)};",
        "ideal I=P;",
    ]
    decomposition_lines = [
        'LIB "primdec.lib";',
        *common_lines,
        "ideal J=P," + ",".join(map(render, product_generators)) + ";",
        "ideal C=std(quotient(J,ideal(PX)));",
        "ideal N1=X,s-1,2*t+1;",
        "ideal N2=X,s-11,2*t+1;",
        "ideal N3=X,s-4*t-3,16*t^2+24*t+13;",
        "ideal N4=X,12*s-7,t-2;",
        "ideal N5=X,s-4*t^2+5,8*t^3+16*t^2+2*t-7;",
        "list LR=primdecGTZ(C);",
        "int i;",
        "for(i=1;i<=size(LR);i++)"
        "{ideal QP=std(LR[i][1]);ideal PP=std(LR[i][2]);"
        "ideal EP=std(radical(eliminate(PP,X)+ideal(X)));"
        '"ROOT_COMPONENT";'
        "eq(EP,N1)+2*eq(EP,N2)+3*eq(EP,N3)+4*eq(EP,N4)+5*eq(EP,N5);"
        "vdim(QP);vdim(PP);kill QP,PP,EP;}",
        "ideal EB=std(eliminate(C,X)+ideal(X));",
        "list LB=primdecGTZ(EB);",
        "for(i=1;i<=size(LB);i++)"
        "{ideal QB=std(LB[i][1]);ideal PB=std(LB[i][2]);"
        '"BASE_COMPONENT";'
        "eq(PB,N1)+2*eq(PB,N2)+3*eq(PB,N3)+4*eq(PB,N4)+5*eq(PB,N5);"
        "vdim(QB);vdim(PB);kill QB,PB;}",
        "exit;",
    ]
    normalization_lines = [
        'LIB "normal.lib";',
        *common_lines,
        "ideal CE=std(ideal(c1,c2,c3));",
        "ideal CA=std(normalConductor(I));",
        '"CONDUCTOR_EQUAL";eq(CE,CA);',
        'list NN=normal(I,"useRing","prim","wd");',
        "ideal ME=" + ",".join(map(render, module_generators)) + ";",
        '"NORMAL_MODULE_EQUAL";eq(ME,NN[2][1]);',
        "def RN=NN[1][1];setring RN;",
        '"NORMAL_VARIABLE_COUNT";nvars(basering);',
        '"NORMAL_RELATION_COUNT";size(norid);',
        f"poly PX={render(P_X)};",
        "ideal JN=norid," + ",".join(map(render, product_generators)) + ";",
        "ideal GN=std(JN);",
        '"NORMAL_PRODUCT_GROEBNER_SIZE";size(GN);',
        '"PX_IN_NORMALIZED_PRODUCT";reduce(PX,GN)==0;',
        "exit;",
    ]
    with tempfile.TemporaryDirectory() as temporary_directory:
        decomposition_source = Path(temporary_directory) / "f20_residue_packets.sing"
        normalization_source = Path(temporary_directory) / "f20_normalization.sing"
        decomposition_source.write_text("\n".join(decomposition_lines) + "\n")
        normalization_source.write_text("\n".join(normalization_lines) + "\n")
        decomposition_completed = subprocess.run(
            (singular, "-q", str(decomposition_source)),
            check=True,
            capture_output=True,
            text=True,
            timeout=180,
            stdin=subprocess.DEVNULL,
        )
        normalization_completed = subprocess.run(
            (singular, "-q", str(normalization_source)),
            check=True,
            capture_output=True,
            text=True,
            timeout=120,
            stdin=subprocess.DEVNULL,
        )
    parsed = parse_singular_output(decomposition_completed.stdout)
    normalization_parsed = parse_singular_output(normalization_completed.stdout)
    parsed.update(
        {
            key: value
            for key, value in normalization_parsed.items()
            if key not in {"root_components", "base_components"}
        }
    )
    expected_scalars = {
        "CONDUCTOR_EQUAL": 1,
        "NORMAL_MODULE_EQUAL": 1,
        "NORMAL_VARIABLE_COUNT": 5,
        "NORMAL_RELATION_COUNT": 12,
        "NORMAL_PRODUCT_GROEBNER_SIZE": 43,
        "PX_IN_NORMALIZED_PRODUCT": 0,
    }
    for label, expected in expected_scalars.items():
        assert parsed.get(label) == expected, (
            label,
            parsed,
            decomposition_completed.stdout,
            normalization_completed.stdout,
        )

    root_by_center: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for center, primary_length, prime_degree in parsed["root_components"]:
        root_by_center[center].append((primary_length, prime_degree))
    base_by_center: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for center, primary_length, prime_degree in parsed["base_components"]:
        base_by_center[center].append((primary_length, prime_degree))
    for values in root_by_center.values():
        values.sort()
    for values in base_by_center.values():
        values.sort()
    assert dict(root_by_center) == {
        1: [(2, 1), (4, 2)],
        2: [(12, 1)],
        3: [(20, 2)],
        4: [(1, 1)],
        5: [(6, 3), (12, 3)],
    }
    assert dict(base_by_center) == {
        1: [(3, 1)],
        2: [(6, 1)],
        3: [(14, 2)],
        4: [(1, 1)],
        5: [(9, 3)],
    }
    assert sum(length for values in root_by_center.values() for length, _ in values) == 57
    assert sum(length for values in base_by_center.values() for length, _ in values) == 33

    return {
        "normalization_verified": True,
        "normalization_variable_count": parsed["NORMAL_VARIABLE_COUNT"],
        "normalization_relation_count": parsed["NORMAL_RELATION_COUNT"],
        "normalization_module_generators": (
            "d*c1/c2",
            "d*c3/c2",
            "1",
        ),
        "conductor_verified": True,
        "root_residue_packets": {
            CENTER_NAMES[center]: [
                {"primary_length": length, "prime_degree": degree}
                for length, degree in values
            ]
            for center, values in sorted(root_by_center.items())
        },
        "base_residue_packets": {
            CENTER_NAMES[center]: [
                {"primary_length": length, "prime_degree": degree}
                for length, degree in values
            ]
            for center, values in sorted(base_by_center.items())
        },
        "natural_normalized_product_groebner_size": parsed[
            "NORMAL_PRODUCT_GROEBNER_SIZE"
        ],
        "P_X_in_natural_normalized_degree_3_1_1_piece": False,
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
    h_q = data["h_q"]
    assert isinstance(s, sp.Symbol)
    assert isinstance(t, sp.Symbol)
    assert isinstance(X, sp.Symbol)

    c1, c2, c3 = conductor_polynomials(data)
    assert sp.expand(c3 + h_q) == 0
    singular_gate = singular_normalization_gate(data, (c1, c2, c3))

    # Exact conductor normalization and connected slope cover.
    y, w, epsilon, Y = sp.symbols("y w epsilon Y")
    q_normal_t = (y**2 - 9) / 8
    q_normal_s = 4 * (y + 2) / ((y + 1) * (y + 3))
    q_repeated_root = 2 / (y + 3)
    y_of_w = (5 + 3 * w**2) / (1 - w**2)
    assert sp.cancel(q.subs({s: q_normal_s, t: q_normal_t})) == 0
    assert all(
        sp.cancel(
            conductor_generator.subs(
                {s: q_normal_s, t: q_normal_t, X: q_repeated_root}
            )
        )
        == 0
        for conductor_generator in (c1, c2, c3)
    )
    assert sp.cancel(((y - 5) / (y + 3) - w**2).subs(y, y_of_w)) == 0

    A_minus = w**2 - 2 * w + 5
    A_plus = w**2 + 2 * w + 5
    B_minus = w**3 - 3 * w**2 - w - 5
    B_plus = w**3 + 3 * w**2 - w + 5
    node_packet = w**4 + 10 * w**2 + 5
    d_on_y = sp.factor(d.subs({s: q_normal_s, t: q_normal_t}))
    r_on_y = sp.factor(r.subs({s: q_normal_s, t: q_normal_t}))
    d_on_w = sp.factor(sp.cancel(d_on_y.subs(y, y_of_w)))
    r_on_w = sp.factor(sp.cancel(r_on_y.subs(y, y_of_w)))
    expected_d_on_w = A_minus**2 * A_plus**2 / (16 * (w**2 + 3) ** 2)
    expected_r_on_w = (
        4
        * w**2
        * A_minus
        * A_plus
        * B_minus**2
        * B_plus**2
        / ((w - 1) ** 6 * (w + 1) ** 6 * (w**2 + 3) ** 2)
    )
    assert sp.cancel(d_on_w - expected_d_on_w) == 0
    assert sp.cancel(r_on_w - expected_r_on_w) == 0
    assert sp.factor(
        (y**2 - 5).subs(y, y_of_w)
        - 4
        * node_packet
        / ((w - 1) ** 2 * (w + 1) ** 2)
    ) == 0

    slope_deformation = sp.cancel(
        P.subs(
            {
                s: q_normal_s + epsilon,
                t: q_normal_t,
                X: q_repeated_root + epsilon * Y,
            }
        )
    )
    slope_residual = sp.factor(sp.Poly(slope_deformation, epsilon).nth(2))
    slope_discriminant = sp.factor(sp.discriminant(slope_residual, Y))
    derivative_residue = (
        w
        * (3 * w**2 + 5)
        * A_minus
        * A_plus
        / (32 * (w - 1) * (w + 1))
    )
    assert sp.cancel(
        slope_discriminant.subs(y, y_of_w) - derivative_residue**2
    ) == 0
    assert sp.cancel(derivative_residue.subs(w, -w) + derivative_residue) == 0

    # The normalization module generators still have value two at triple E1.
    tau, z, b = sp.symbols("tau z b")
    imaginary_unit = sp.I
    triple_substitution = {
        s: 2 * imaginary_unit + tau**4,
        t: -sp.Rational(3, 4) + imaginary_unit / 2 + z * tau**4,
        X: 1 + imaginary_unit + tau**2 * (b + tau * Y),
    }
    modulus = b**2 + imaginary_unit
    conductor_orders = tuple(
        lowest_order_mod(
            generator.subs(triple_substitution), tau, modulus, b
        )
        for generator in (c1, c2, c3)
    )
    normalization_generator_orders = (
        4 + conductor_orders[0] - conductor_orders[1],
        4 + conductor_orders[2] - conductor_orders[1],
    )
    assert conductor_orders == (2, 4, 2)
    assert normalization_generator_orders == (2, 2)

    unit_pullback = sp.Matrix(
        (
            (-1, -1, -2),
            (-1, -1, -2),
            (1, 0, 0),
            (0, 0, 1),
        )
    )
    selector_completion = unit_pullback.row_join(sp.Matrix((1, 0, 0, 0)))
    assert selector_completion.det() == -1

    return {
        "schema": "f20-normalized-cox-conductor-v1",
        "status": "global_normalization_certified_exceptional_cox_fill_required",
        "global_normalization": {
            "ring": "R=QQ[s,t,X]/(P)",
            "conductor": (str(c1), str(c2), str(c3)),
            "c3_equals_minus_h_q": True,
            "module": "R + R*(d*c1/c2) + R*(d*c3/c2)",
            "module_generator_orders_at_triple_E1": normalization_generator_orders,
            "triple_E1_value_one_present": False,
            "singular_gate": singular_gate,
        },
        "conductor_equalizer": {
            "square": "R = Rbar fiber_product_(Rbar/c) R/c",
            "base_conductor_parameter": "y",
            "normal_conductor_parameter": "w",
            "cover": "y=(5+3*w^2)/(1-w^2)",
            "involution": "w -> -w",
            "selector": "w-1",
            "unit_completion_determinant": -1,
            "collision_packets_on_w_line": {
                "q_node": str(node_packet),
                "triple_orbit": (str(A_minus), str(A_plus)),
                "q_r_transverse": "w",
                "q_r_tangency_orbit": (str(B_minus), str(B_plus)),
                "r_cusp": "off the q conductor",
            },
            "boundary_restrictions": {
                "d": str(expected_d_on_w),
                "r": str(expected_r_on_w),
            },
            "derivative_slope_residue": str(derivative_residue),
            "derivative_residue_parity": "anti_invariant",
            "individual_D_d_D_q_D_r_residue_cocycle": "not_yet_distributed",
        },
        "exceptional_cox_gate": {
            "normalization_alone_suffices": False,
            "reason": (
                "both integral normalization generators have triple-E1 order two; "
                "the desired D_d coefficient is one"
            ),
            "required_algebra": (
                "normalized multi-Rees/total-coordinate algebra of the resolved "
                "exceptional divisors, with a primitive exceptional Cox variable"
            ),
        },
        "downstream": {
            "natural_normalized_degree_3_1_1_contains_P_X": False,
            "abstract_divisorial_degree_3_1_1_contains_P_X": (
                "formal_on_a_certified_regular_resolution"
            ),
            "inverse_adjugate_polynomiality": "not_reached",
            "affine_space_recognition": "not_reached",
        },
        "software_assumptions": {
            "python": ".venv Python with pinned SymPy",
            "singular": (
                "Singular 4.4-compatible normal.lib and primdec.lib exact routines"
            ),
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
    print("PASS: the global F20 normalization module and conductor are exact")
    print("PASS: the length-57 residue splits into five certified center packets")
    print("PASS: the conductor derivative residue is anti-invariant under w -> -w")
    print("OBSTRUCTION: normalization alone does not repair the (3,1,1) product")
    print("SCOPE: exceptional Cox variables on the resolved normalization remain open")


if __name__ == "__main__":
    main()
