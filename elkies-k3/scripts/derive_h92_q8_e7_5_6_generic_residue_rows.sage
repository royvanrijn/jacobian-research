#!/usr/bin/env sage -python
"""Exact q=8 generic residue rows on the actual conic components E7_5,E7_6.

Unlike a rational parametrization, this evaluator works directly in the
component coordinate fields

    QQ(U,Y)/(F(0,U,Y)).

It clears one exact common denominator per normal-order group and reduces in
the actual conic coordinate ring. This keeps the q=8 residue calculation
exact while avoiding artificial coefficient expansion from a chosen conic
parameter. At E7_5 it retains the certified cancellation in x-x(P1) and
y-y(P1); at E7_6 the strict inequalities give m=y/x.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, sage_eval


ROOT = Path(__file__).resolve().parents[2]
FRAMES = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-generic-component-chart-frames.json"
CONDITIONS = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-conditions.json"
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-5-6-generic-residue-rows.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def path_label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def poly_from_coefficients(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def invert_base(rational_u):
    numerator, denominator = rational_u.numerator(), rational_u.denominator()
    t_ring = PolynomialRing(QQ, "t")
    t = t_ring.gen()
    field = t_ring.fraction_field()
    return field(
        t**(denominator.degree()-numerator.degree())
        * t_ring(list(reversed(numerator.list())))
        / t_ring(list(reversed(denominator.list())))
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--frames", type=Path, default=FRAMES)
parser.add_argument("--conditions", type=Path, default=CONDITIONS)
parser.add_argument("--p1", type=Path, default=P1)
parser.add_argument("--component", action="append", choices=("E7_5", "E7_6"))
parser.add_argument("--trace", action="store_true")
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

frames = json.loads(args.frames.read_text())
conditions = json.loads(args.conditions.read_text())
p1 = json.loads(args.p1.read_text())
assert frames["status"] == "PASS_EXACT_Q8_GENERIC_COMPONENT_CHART_FRAMES"
assert conditions["status"] == "PASS_EXACT_Q8_ALL_COMPONENT_GENERIC_CONDITION_TEMPLATE"
assert p1["status"] == "PASS_EXACT_H92_P1"

frame_by_component = {entry["component"]: entry for entry in frames["component_frames"]}
condition_by_component = {entry["component"]: entry for entry in conditions["component_conditions"]}
ambient_ref = conditions["inputs"]["endpoint_ambient"]
ambient_data = json.loads((ROOT / ambient_ref["path"]).read_text())
assert digest(ROOT / ambient_ref["path"]) == ambient_ref["sha256"]
basis = ambient_data["ambient_basis"]

u_base_ring = PolynomialRing(QQ, "u")
u_base = u_base_ring.fraction_field()
x_p_t = invert_base(
    u_base(poly_from_coefficients(u_base_ring, p1["x_entrance_base"]["numerator_coefficients"]))
    / u_base(poly_from_coefficients(u_base_ring, p1["x_entrance_base"]["denominator_coefficients"]))
)
y_p_t = invert_base(
    u_base(poly_from_coefficients(u_base_ring, p1["y_entrance_base"]["numerator_coefficients"]))
    / u_base(poly_from_coefficients(u_base_ring, p1["y_entrance_base"]["denominator_coefficients"]))
)
t_base = x_p_t.parent().gen()
x_p2 = QQ((x_p_t/t_base**2)(0))
x_p3 = QQ((x_p_t/t_base**2).derivative()(0))
y_p3 = QQ((y_p_t/t_base**3)(0))
h_leading = QQ(p1["structured_denominator"]["Z4_coefficients"][-1])
assert h_leading

chart_ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = chart_ring.gens()


def initial_z(poly, ring, u, y):
    """Return the first Z coefficient as a polynomial in U,Y."""
    terms = poly.dict()
    order = min(exponents[0] for exponents in terms)
    value = sum(
        coefficient*u**exponents[1]*y**exponents[2]
        for exponents, coefficient in terms.items() if exponents[0] == order
    )
    return order, ring(value)


def restriction_z_zero(poly, ring, u, y):
    return ring(sum(
        coefficient*u**exponents[1]*y**exponents[2]
        for exponents, coefficient in poly.dict().items() if exponents[0] == 0
    ))


component_payloads = []
selected_components = tuple(args.component or ("E7_5", "E7_6"))
for component in selected_components:
    g_y_power = {"E7_5": 5, "E7_6": 6}[component]
    if args.trace:
        print("H92Q8E756|component={}|stage=chart".format(component), flush=True)
    frame = frame_by_component[component]
    assert frame["component_equation"] == "Z"
    surface = chart_ring(sage_eval(
        frame["surface_equation"], locals={"Z": Z, "U": U, "Y": Y}
    ))
    old = {
        name: chart_ring(sage_eval(expression, locals={"Z": Z, "U": U, "Y": Y}))
        for name, expression in frame["old_coordinate_pullback"].items()
    }
    component_ring = PolynomialRing(QQ, names=("U", "Y"))
    uu, yy = component_ring.gens()
    relation = restriction_z_zero(surface, component_ring, uu, yy)
    assert relation.degree(yy) == 2
    groebner = component_ring.ideal((relation,)).groebner_basis()
    component_field = component_ring.fraction_field()

    t_order, t0 = initial_z(old["t"], component_ring, uu, yy)
    x_order, x0 = initial_z(old["x"], component_ring, uu, yy)
    y_order, y0 = initial_z(old["y"], component_ring, uu, yy)
    expected = condition_by_component[component]["orders"]
    chart_orders = frame["old_coordinate_weight_orders"]
    assert (t_order, x_order, y_order) == (
        expected["t"], expected["x"], int(chart_orders["y"])
    )
    if component == "E7_5":
        _, lx = initial_z(old["x"] - x_p2*old["t"]**2 - x_p3*old["t"]**3, component_ring, uu, yy)
        _, ly = initial_z(old["y"] - y_p3*old["t"]**3, component_ring, uu, yy)
        m0 = component_field(ly)/component_field(lx)
        # The audit proves these are the exact Z^3 leading numerator and
        # denominator, so their quotient has the certified order zero.
        assert expected["m"] == 0
    else:
        m0 = component_field(y0)/component_field(x0)
        assert y_order-x_order == expected["m"]
    g0 = component_field(yy**g_y_power)

    rows = []
    for group in condition_by_component[component]["negative_order_groups"]:
        if len(group["basis_indices"]) == 1:
            continue
        residues = []
        for index in group["basis_indices"]:
            entry = basis[index]
            degree = 4*int(entry["h_power"])-int(entry["u_power"])
            residue = (
                g0 * component_field(t0)**(degree-9)
                * component_field(x0)**int(entry["x_power"])
                * m0**int(entry["m_power"])
                / h_leading**int(entry["h_power"])
            )
            residues.append((int(index), residue))
        if args.trace:
            print("H92Q8E756|component={}|order={}|stage=denominator".format(
                component, group["residual_order"]
            ), flush=True)
        if component == "E7_5":
            # The audited chord is ly/lx and every source monomial has
            # m-exponent at most nine.  This explicit common denominator
            # avoids costly factorization/lcm work on the large H92 conic
            # coefficients while clearing every t and chord denominator.
            max_t_pole = max(
                9 - (4*int(basis[index]["h_power"])-int(basis[index]["u_power"]))
                for index in group["basis_indices"]
            )
            denominator = component_ring(t0)**max_t_pole * component_ring(lx)**9
        else:
            denominator = component_ring.one()
            for _, value in residues:
                denominator = denominator.lcm(component_ring(value.denominator()))
        coefficient_rows = {}
        for index, value in residues:
            numerator = component_ring(value.numerator())
            value_denominator = component_ring(value.denominator())
            quotient, remainder = denominator.quo_rem(value_denominator)
            assert not remainder
            reduced = component_ring(numerator*quotient).reduce(groebner)
            for exponents, coefficient in reduced.dict().items():
                coefficient_rows.setdefault((int(exponents[0]), int(exponents[1])), []).append({
                    "basis_index": index, "coefficient": str(QQ(coefficient)),
                })
        for (u_power, y_power), entries in sorted(coefficient_rows.items()):
            rows.append({
                "residual_order": int(group["residual_order"]),
                "component_monomial": {"U_power": u_power, "Y_power": y_power},
                "entries": entries,
                "cleared_denominator": str(denominator),
                "interpretation": "the displayed coordinate-ring coefficient must vanish",
            })
    component_payloads.append({
        "component": component,
        "chart": frame["actual_edge_chart"],
        "component_relation": str(relation),
        "orders": {"t": int(t_order), "x": int(x_order), "y": int(y_order), "m": int(expected["m"])},
        "non_singleton_residue_rows": rows,
    })

payload = {
    "schema": "elkies-k3.h92-q8-e7-5-6-generic-residue-rows.v1",
    "status": "PASS_EXACT_Q8_E7_5_6_GENERIC_RESIDUE_ROWS",
    "inputs": {
        "component_frames": {"path": path_label(args.frames), "sha256": digest(args.frames)},
        "generic_conditions": {"path": path_label(args.conditions), "sha256": digest(args.conditions)},
        "p1": {"path": path_label(args.p1), "sha256": digest(args.p1)},
    },
    "ambient_basis_sha256": conditions["ambient_basis_sha256"],
    "components": component_payloads,
    "boundary": (
        "This is an exact generic component-ring calculation on E7_5 and E7_6 "
        "only. It does not cover E7_1--E7_4/E7_7, nodes, marked branch, "
        "overlaps, a global kernel, or a pencil."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q8E756GENERICRESIDUES|components={}|rows={}|status="
    "PASS_EXACT_Q8_E7_5_6_GENERIC_RESIDUE_ROWS".format(
        len(component_payloads),
        sum(len(item["non_singleton_residue_rows"]) for item in component_payloads)
    ),
    flush=True,
)
