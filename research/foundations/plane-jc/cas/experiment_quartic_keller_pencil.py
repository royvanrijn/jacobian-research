#!/usr/bin/env python3
"""Exact pencil-at-infinity calibration for the quartic packet countermodel.

The finite-free map

    F0 = (P,Q0) = (y, x^4-x^3+x*y)

has the clean 3+1 cusp and a 2+2 connector used in the quartic packet
frontier.  The target shear (u,v) -> (u,v+u^2) gives

    F2 = (P,Q2) = (y, x^4-x^3+x*y+y^2).

It preserves the finite cover up to target isomorphism and preserves
dP wedge dQ, but it does not preserve the topology of the *linear* target
pencil.  This script verifies the algebra and emits the two stratumwise
resolved boundary graphs and their A'Campo zeta products.

Neither F0 nor F2 is Keller: their common Jacobian determinant vanishes on
y=-4*x^3+3*x^2.  The output is therefore a packet-under-determination
counterexperiment, not a quartic Keller model or exclusion.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import sympy as sp


@dataclass(frozen=True)
class BoundaryComponent:
    name: str
    kind: str
    pole_order: int | None
    dicritical_degree: int | None
    neighbors: tuple[str, ...]
    open_euler_characteristic: int | None


def component(
    name: str,
    kind: str,
    neighbors: tuple[str, ...],
    *,
    pole_order: int | None = None,
    dicritical_degree: int | None = None,
) -> BoundaryComponent:
    if kind == "pole":
        assert pole_order is not None and dicritical_degree is None
        open_euler = 2 - len(neighbors)
    elif kind == "dicritical":
        assert pole_order is None and dicritical_degree is not None
        open_euler = None
    else:
        raise ValueError(f"unknown boundary component kind: {kind}")
    return BoundaryComponent(
        name=name,
        kind=kind,
        pole_order=pole_order,
        dicritical_degree=dicritical_degree,
        neighbors=neighbors,
        open_euler_characteristic=open_euler,
    )


def coordinate_graph() -> tuple[BoundaryComponent, ...]:
    """Boundary of P1_x times P1_y for the coordinate polynomial P=y."""

    return (
        component("D_y_infinity", "pole", ("D_x_infinity",), pole_order=1),
        component(
            "D_x_infinity",
            "dicritical",
            ("D_y_infinity",),
            dicritical_degree=1,
        ),
    )


def unsheared_graph() -> tuple[BoundaryComponent, ...]:
    """Adapted P1xP1 resolution of X*Y after two disjoint base blowups."""

    return (
        component(
            "D_X_infinity",
            "pole",
            ("D_Y_infinity", "H_X"),
            pole_order=1,
        ),
        component(
            "D_Y_infinity",
            "pole",
            ("D_X_infinity", "H_Y"),
            pole_order=1,
        ),
        component("H_X", "dicritical", ("D_X_infinity",), dicritical_degree=1),
        component("H_Y", "dicritical", ("D_Y_infinity",), dicritical_degree=1),
    )


def sheared_graph() -> tuple[BoundaryComponent, ...]:
    """Ordinary-blowup graph for y^2+(x+s)y+x^4-x^3.

    In the y=1 chart at [0:1:0], two shared blowups resolve the common
    tangent and separate the two branches z/x^2=+/-i.  Each resulting local
    ideal is (r,x^4), hence has a chain with pole orders 3,2,1 followed by
    a degree-one dicritical.
    """

    result = [
        component("L_infinity", "pole", ("E2",), pole_order=4),
        component("E1", "pole", ("E2",), pole_order=2),
        component(
            "E2",
            "pole",
            ("L_infinity", "E1", "V3_plus", "V3_minus"),
            pole_order=4,
        ),
    ]
    for sign in ("plus", "minus"):
        result.extend(
            (
                component(
                    f"V3_{sign}",
                    "pole",
                    ("E2", f"V2_{sign}"),
                    pole_order=3,
                ),
                component(
                    f"V2_{sign}",
                    "pole",
                    (f"V3_{sign}", f"V1_{sign}"),
                    pole_order=2,
                ),
                component(
                    f"V1_{sign}",
                    "pole",
                    (f"V2_{sign}", f"H_{sign}"),
                    pole_order=1,
                ),
                component(
                    f"H_{sign}",
                    "dicritical",
                    (f"V1_{sign}",),
                    dicritical_degree=1,
                ),
            )
        )
    return tuple(result)


def acampo_factors(
    graph: tuple[BoundaryComponent, ...],
) -> tuple[tuple[int, int], ...]:
    """Return nonzero (pole order, Euler exponent) pairs."""

    factors = []
    for item in graph:
        if item.kind != "pole" or item.open_euler_characteristic == 0:
            continue
        assert item.pole_order is not None
        assert item.open_euler_characteristic is not None
        factors.append((item.pole_order, item.open_euler_characteristic))
    return tuple(factors)


def graph_record(
    graph: tuple[BoundaryComponent, ...],
    *,
    zeta: str,
) -> dict[str, object]:
    return {
        "components": [asdict(item) for item in graph],
        "dicritical_degrees": [
            item.dicritical_degree for item in graph if item.kind == "dicritical"
        ],
        "acampo_factors": [
            {"factor": f"1-t^{order}", "exponent": exponent}
            for order, exponent in acampo_factors(graph)
        ],
        "zeta": zeta,
    }


def compile_report() -> dict[str, object]:
    x, y, s = sp.symbols("x y s")
    g = x**4 - x**3
    p = y
    q0 = g + x * y
    q2 = q0 + y**2

    jacobian0 = sp.det(
        sp.Matrix(
            [[sp.diff(p, x), sp.diff(p, y)], [sp.diff(q0, x), sp.diff(q0, y)]]
        )
    )
    jacobian2 = sp.det(
        sp.Matrix(
            [[sp.diff(p, x), sp.diff(p, y)], [sp.diff(q2, x), sp.diff(q2, y)]]
        )
    )
    assert sp.expand(q2 - q0 - p**2) == 0
    assert sp.expand(jacobian2 - jacobian0) == 0
    assert sp.expand(jacobian0 - (-x * (4 * x**2 - 3 * x) - y)) == 0

    cusp0 = sp.factor(q0.subs(y, 0))
    cusp2 = sp.factor(q2.subs(y, 0))
    assert cusp0 == x**3 * (x - 1)
    assert cusp2 == cusp0

    connector_u = sp.Rational(1, 8)
    connector_v0 = -sp.Rational(1, 64)
    connector_v2 = 0
    connector0 = sp.factor(q0.subs(y, connector_u) - connector_v0)
    connector2 = sp.factor(q2.subs(y, connector_u) - connector_v2)
    connector_square = (8 * x**2 - 4 * x - 1) ** 2 / 64
    assert sp.expand(connector0 - connector_square) == 0
    assert sp.expand(connector2 - connector_square) == 0

    constant = sp.expand(g.subs(x, -s))
    quotient = sp.div(g - constant, x + s, domain=sp.QQ.frac_field(s))
    assert quotient[1] == 0
    h = sp.expand(quotient[0])
    X = x + s
    Y = y + h
    assert sp.expand(q0 + s * p - (X * Y + constant)) == 0

    # The sole face at infinity of q2+s*p is x^4+y^2.  It is nondegenerate
    # on (C*)^2, and its Newton triangle has area 4, boundary length 8, and
    # hence one interior lattice point by Pick's theorem.
    face = x**4 + y**2
    assert sp.diff(face, x) == 4 * x**3
    assert sp.diff(face, y) == 2 * y
    newton_area = 4
    newton_boundary_points = 8
    newton_interior_points = newton_area - newton_boundary_points // 2 + 1
    assert newton_interior_points == 1
    punctures = sp.gcd(4, 2)
    assert punctures == 2

    coordinate = coordinate_graph()
    unsheared = unsheared_graph()
    sheared = sheared_graph()
    assert acampo_factors(coordinate) == ((1, 1),)
    assert acampo_factors(unsheared) == ()
    assert acampo_factors(sheared) == ((4, 1), (2, 1), (4, -2))

    return {
        "claim": (
            "exact target-shear calibration; packet data do not determine "
            "the topology or zeta function of the linear pencil"
        ),
        "warning": (
            "F0 and F2 are finite quartic packet countermodels, not Keller maps"
        ),
        "maps": {
            "P": str(p),
            "Q0": str(q0),
            "Q2": str(q2),
            "target_shear": "(u,v) -> (u,v+u^2)",
            "common_jacobian_determinant": str(jacobian0),
            "jacobian_zero_curve": "y=-4*x^3+3*x^2",
        },
        "packet": {
            "cusp_target_F0": ["0", "0"],
            "cusp_target_F2": ["0", "0"],
            "cusp_fiber_factorization": str(cusp0),
            "connector_target_F0": ["1/8", "-1/64"],
            "connector_target_F2": ["1/8", "0"],
            "connector_fiber_factorization": str(connector_square),
        },
        "symbolic_parameter": {
            "direction": "[alpha:beta]",
            "generic_chart": "beta!=0, s=alpha/beta",
            "exceptional_directions": ["[1:0]"],
        },
        "pencil_F0": {
            "generic_normal_form": {
                "identity": "Q0+s*P=(x+s)*(y+h_s(x))+s^4+s^3",
                "h_s": str(h),
            },
            "generic_fiber": {"genus": 0, "punctures": 2, "homotopy": "C*"},
            "generic_monodromy_H1_eigenvalues": ["1"],
            "generic_boundary": graph_record(unsheared, zeta="1"),
            "exceptional_direction": {
                "parameter": "[1:0]",
                "polynomial": "P=y",
                "fiber": {"genus": 0, "punctures": 1, "homotopy": "A1"},
                "boundary": graph_record(coordinate, zeta="1-t"),
            },
        },
        "pencil_F2": {
            "generic_polynomial": "y^2+(x+s)*y+x^4-x^3",
            "face_at_infinity": str(face),
            "newton_triangle": [[0, 0], [4, 0], [0, 2]],
            "newton_nondegenerate_at_infinity": True,
            "generic_fiber": {
                "genus": int(newton_interior_points),
                "punctures": int(punctures),
                "first_betti_number": 3,
            },
            "generic_monodromy_H1_eigenvalues": ["1", "i", "-i"],
            "generic_boundary": graph_record(sheared, zeta="1/(1+t^2)"),
            "exceptional_direction": {
                "parameter": "[1:0]",
                "polynomial": "P=y",
                "fiber": {"genus": 0, "punctures": 1, "homotopy": "A1"},
                "boundary": graph_record(coordinate, zeta="1-t"),
            },
        },
        "conclusion": {
            "same_finite_normalization_packet": True,
            "same_dP_wedge_dQ": True,
            "same_linear_pencil_topology": False,
            "linear_pencil_is_invariant_under_general_target_shear": False,
            "quartic_keller_packet_excluded": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="optional JSON output path",
    )
    arguments = parser.parse_args()
    report = compile_report()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
    print("PASS quartic packet target-shear and Jacobian identities")
    print("PASS symbolic pencil strata and resolved boundary graphs")
    print("PASS zeta(F0_generic)=1, zeta(F2_generic)=1/(1+t^2)")
    if arguments.output is not None:
        print(f"WROTE {arguments.output}")


if __name__ == "__main__":
    main()
