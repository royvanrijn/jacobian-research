#!/usr/bin/env sage -python
"""Exact q=8 generic residue rows on actual Y-branch components E7_1..E7_3.

For each selected H92 chart, substitute the certified normal weights
``Y=e`` and the adjacent quadratic coordinate ``e^2*A``.  The leading actual
surface equation is linear in A, so it solves an exact rational component
normalization A=A(s).  This uses the blow-up chart itself, including the
nonconstant E7_2/E7_3 normal coefficient, rather than a Kodaira proxy.
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
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-1-3-generic-residue-rows.json"


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
    return field(t**(denominator.degree()-numerator.degree())
                 * t_ring(list(reversed(numerator.list())))
                 / t_ring(list(reversed(denominator.list()))))


def e_leading(value):
    def leading_polynomial(poly):
        terms = poly.dict()
        order = min(exponents[0] for exponents in terms)
        coefficient = sum(
            scalar*s**exponents[1]
            for exponents, scalar in terms.items() if exponents[0] == order
        )
        return order, coefficient
    numerator_order, numerator = leading_polynomial(value.numerator())
    denominator_order, denominator = leading_polynomial(value.denominator())
    return numerator_order-denominator_order, numerator/denominator


def coefficient_rows(residues):
    s_ring = PolynomialRing(QQ, "s")
    s_only = s_ring.gen()
    projection = series_ring.hom([s_ring(0), s_only], s_ring)
    field = s_ring.fraction_field()
    converted = [
        (index, field(projection(value.numerator()))/field(projection(value.denominator())))
        for index, value in residues
    ]
    denominator = s_ring.one()
    for _, value in converted:
        denominator = denominator.lcm(s_ring(value.denominator()))
    rows = {}
    for index, value in converted:
        quotient, remainder = denominator.quo_rem(s_ring(value.denominator()))
        assert not remainder
        polynomial = s_ring(value.numerator())*quotient
        for exponent, coefficient in polynomial.dict().items():
            rows.setdefault(int(exponent), []).append({
                "basis_index": int(index), "coefficient": str(QQ(coefficient)),
            })
    return denominator, rows


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--frames", type=Path, default=FRAMES)
parser.add_argument("--conditions", type=Path, default=CONDITIONS)
parser.add_argument("--p1", type=Path, default=P1)
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
x_p_t = invert_base(u_base(poly_from_coefficients(u_base_ring, p1["x_entrance_base"]["numerator_coefficients"])) / u_base(poly_from_coefficients(u_base_ring, p1["x_entrance_base"]["denominator_coefficients"])))
t_base = x_p_t.parent().gen()
x_p2 = QQ((x_p_t/t_base**2)(0))
h_leading = QQ(p1["structured_denominator"]["Z4_coefficients"][-1])

chart_ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = chart_ring.gens()
series_ring = PolynomialRing(QQ, names=("e", "s"))
e, s = series_ring.gens()
series_field = series_ring.fraction_field()
aux_ring = PolynomialRing(QQ, names=("e", "s", "A"))
ee, ss, A = aux_ring.gens()


def normalization(surface, weights):
    if tuple(weights) == (2, 0, 1):
        substituted = aux_ring(surface.subs({Z: ee**2*A, U: ss, Y: ee}))
    else:
        assert tuple(weights) == (0, 2, 1)
        substituted = aux_ring(surface.subs({Z: ss, U: ee**2*A, Y: ee}))
    terms = substituted.dict()
    order = min(exponents[0] for exponents in terms)
    leading = sum(
        coefficient*ss**exponents[1]*A**exponents[2]
        for exponents, coefficient in terms.items() if exponents[0] == order
    )
    polynomial = leading.polynomial(A)
    assert polynomial.degree() == 1
    a_value = -polynomial[0]/polynomial[1]
    def into_series(value):
        numerator = aux_ring(value.numerator())
        denominator = aux_ring(value.denominator())
        def map_polynomial(poly):
            return sum(
                coefficient*s**exponents[1]
                for exponents, coefficient in poly.dict().items()
                if not exponents[0] and not exponents[2]
            )
        mapped_numerator = map_polynomial(numerator)
        mapped_denominator = map_polynomial(denominator)
        assert mapped_numerator and mapped_denominator
        return series_field(mapped_numerator/mapped_denominator)
    if tuple(weights) == (2, 0, 1):
        return {Z: e**2*into_series(a_value), U: s, Y: e}
    return {Z: s, U: e**2*into_series(a_value), Y: e}


component_payloads = []
for component, g_power in (("E7_1", (2, 4)), ("E7_2", (5, 6)), ("E7_3", (6, 4))):
    frame = frame_by_component[component]
    assert frame["component_equation"] == "Y"
    surface = chart_ring(sage_eval(frame["surface_equation"], locals={"Z": Z, "U": U, "Y": Y}))
    old = {name: chart_ring(sage_eval(expression, locals={"Z": Z, "U": U, "Y": Y}))
           for name, expression in frame["old_coordinate_pullback"].items()}
    values = normalization(surface, frame["normal_weights_Z_U_Y"])
    pullback = {name: series_field(value.subs(values)) for name, value in old.items()}
    t_value, x_value, y_value = (pullback[name] for name in ("t", "x", "y"))
    t_order, _ = e_leading(t_value)
    x_order, _ = e_leading(x_value)
    y_order, _ = e_leading(y_value)
    expected = condition_by_component[component]["orders"]
    assert (t_order, x_order, y_order) == (
        expected["t"], expected["x"], frame["old_coordinate_weight_orders"]["y"]
    )
    m_value = y_value/(x_value-x_p2*t_value**2) if component == "E7_2" else y_value/x_value
    m_order, _ = e_leading(m_value)
    assert m_order == expected["m"]
    g = series_field(e**g_power[0]*s**g_power[1])

    rows = []
    for group in condition_by_component[component]["negative_order_groups"]:
        if len(group["basis_indices"]) == 1:
            continue
        residues = []
        for index in group["basis_indices"]:
            entry = basis[index]
            degree = 4*int(entry["h_power"])-int(entry["u_power"])
            value = (g*t_value**(degree-9)*x_value**int(entry["x_power"])
                     *m_value**int(entry["m_power"])/h_leading**int(entry["h_power"]))
            order, residue = e_leading(value)
            assert order == int(group["residual_order"])
            residues.append((int(index), residue))
        denominator, grouped = coefficient_rows(residues)
        for s_power, entries in sorted(grouped.items()):
            rows.append({
                "residual_order": int(group["residual_order"]),
                "component_parameter_power": int(s_power),
                "entries": entries,
                "cleared_denominator": str(denominator),
                "interpretation": "the displayed component-parameter coefficient must vanish",
            })
    component_payloads.append({
        "component": component,
        "chart": frame["actual_edge_chart"],
        "normalization": {str(key): str(value) for key, value in values.items()},
        "orders": {"t": int(t_order), "x": int(x_order), "y": int(y_order), "m": int(m_order)},
        "non_singleton_residue_rows": rows,
    })

payload = {
    "schema": "elkies-k3.h92-q8-e7-1-3-generic-residue-rows.v1",
    "status": "PASS_EXACT_Q8_E7_1_3_GENERIC_RESIDUE_ROWS",
    "inputs": {"component_frames": {"path": path_label(args.frames), "sha256": digest(args.frames)}, "generic_conditions": {"path": path_label(args.conditions), "sha256": digest(args.conditions)}, "p1": {"path": path_label(args.p1), "sha256": digest(args.p1)}},
    "ambient_basis_sha256": conditions["ambient_basis_sha256"],
    "components": component_payloads,
    "boundary": "This is the generic residue layer on E7_1,E7_2,E7_3 only; nodes, marked branch, overlaps, a global kernel, and a pencil are not evaluated.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H92Q8E713GENERICRESIDUES|components=3|rows={}|status=PASS_EXACT_Q8_E7_1_3_GENERIC_RESIDUE_ROWS".format(sum(len(item["non_singleton_residue_rows"]) for item in component_payloads)), flush=True)
