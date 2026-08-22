#!/usr/bin/env sage -python
"""Derive exact non-singleton generic residue rows on actual E7 components.

This is the first chart-function-field layer beyond component valuations.  On
the actual H92 charts, E7_4 and E7_7 have rational normalizations

    E7_4: Z=e, U=s^2, Y=s;
    E7_7: U=e, Z=gamma*s^2, Y=s,

where gamma is solved from the *actual* restricted surface equation.  These
are the default audited components. Optional E7_5/E7_6 conic charts are kept
behind an explicit selection while their larger exact residue reductions are
being optimized. The script records every exact leading residue equation from
a non-singleton negative-order group. It does not replace E7_1--E7_3, nodes,
marked branch, or overlaps.
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
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-4-7-generic-residue-rows.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def path_label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def e_leading(value, e, s):
    """Return (e-order, leading QQ(s) coefficient) of a nonzero rational map."""
    value = value.parent()(value)

    def leading_polynomial(poly):
        terms = poly.dict()
        order = min(exponents[0] for exponents in terms)
        coefficient = sum(
            scalar * s**exponents[1]
            for exponents, scalar in terms.items() if exponents[0] == order
        )
        return order, coefficient

    numerator_order, numerator = leading_polynomial(value.numerator())
    denominator_order, denominator = leading_polynomial(value.denominator())
    return numerator_order-denominator_order, numerator/denominator


def residue_coefficient_rows(residues, s):
    """Clear a common QQ[s] denominator and return coefficient equations."""
    s_ring = PolynomialRing(QQ, "s")
    s_only = s_ring.gen()
    projection = series_ring.hom([s_ring(0), s_only], s_ring)
    s_field = s_ring.fraction_field()
    converted = [
        (index, s_field(projection(value.numerator())) / s_field(projection(value.denominator())))
        for index, value in residues
    ]
    common_denominator = s_ring.one()
    for _, value in converted:
        common_denominator = common_denominator.lcm(s_ring(value.denominator()))
    rows = {}
    for index, value in converted:
        numerator = s_ring(value.numerator())
        denominator = s_ring(value.denominator())
        quotient, remainder = common_denominator.quo_rem(denominator)
        assert not remainder
        scaled = numerator*quotient
        for exponent, coefficient in scaled.dict().items():
            rows.setdefault(int(exponent), []).append({
                "basis_index": int(index), "coefficient": str(QQ(coefficient)),
            })
    return rows


def solve_linear_coordinate(surface, coordinate, substitution, ring):
    """Solve a linear restricted chart equation after the chosen substitution."""
    value = ring(surface.subs(substitution))
    polynomial = value.polynomial(coordinate)
    assert polynomial.degree() == 1
    return -polynomial[0]/polynomial[1]


def conic_z_normalization(surface):
    """Parametrize F(0,U,Y)=0 from (U,Y)=(0,0) by Y=s*U."""
    restricted = chart_ring(surface.subs({Z: 0}))
    polynomial = restricted.polynomial(U)
    assert polynomial.degree() == 2
    constant, linear, quadratic = polynomial[0], polynomial[1], polynomial[2]
    quotient, remainder = constant.quo_rem(Y**2)
    assert not remainder and linear and quadratic and quotient
    linear = series_field(chart_ring(linear).subs({Y: s}))
    quadratic = series_field(chart_ring(quadratic).subs({Y: s}))
    quotient = series_field(chart_ring(quotient).subs({Y: s}))
    return -linear/(quadratic + quotient*s**2)


def polynomial_from_coefficients(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def invert_base(rational_u):
    """Rewrite QQ(u) as QQ(t), with t=1/u, preserving the exact germ at t=0."""
    numerator = rational_u.numerator()
    denominator = rational_u.denominator()
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
parser.add_argument(
    "--component", action="append", choices=("E7_4", "E7_5", "E7_6", "E7_7"),
    help="component to evaluate; repeat to select several (default: E7_4 and E7_7)",
)
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
ambient = conditions["inputs"]["endpoint_ambient"]
ambient_data = json.loads((ROOT / ambient["path"]).read_text())
assert digest(ROOT / ambient["path"]) == ambient["sha256"]
basis = ambient_data["ambient_basis"]

chart_ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = chart_ring.gens()
series_ring = PolynomialRing(QQ, names=("e", "s"))
e, s = series_ring.gens()
series_field = series_ring.fraction_field()
h_leading = QQ(p1["structured_denominator"]["Z4_coefficients"][-1])
assert h_leading
u_ring = PolynomialRing(QQ, "u")
u_field = u_ring.fraction_field()
x_p_t = invert_base(
    u_field(polynomial_from_coefficients(
        u_ring, p1["x_entrance_base"]["numerator_coefficients"]
    )) / u_field(polynomial_from_coefficients(
        u_ring, p1["x_entrance_base"]["denominator_coefficients"]
    ))
)
y_p_t = invert_base(
    u_field(polynomial_from_coefficients(
        u_ring, p1["y_entrance_base"]["numerator_coefficients"]
    )) / u_field(polynomial_from_coefficients(
        u_ring, p1["y_entrance_base"]["denominator_coefficients"]
    ))
)
t_base = x_p_t.parent().gen()
x_p2 = QQ((x_p_t/t_base**2)(0))
x_p3 = QQ((x_p_t/t_base**2).derivative()(0))
y_p3 = QQ((y_p_t/t_base**3)(0))


def chart_pullbacks(component):
    frame = frame_by_component[component]
    surface = chart_ring(sage_eval(
        frame["surface_equation"], locals={"Z": Z, "U": U, "Y": Y}
    ))
    old = {
        name: chart_ring(sage_eval(expression, locals={"Z": Z, "U": U, "Y": Y}))
        for name, expression in frame["old_coordinate_pullback"].items()
    }
    if component == "E7_4":
        assert frame["component_equation"] == "Z"
        # At Z=0 the actual equation is Y^2-U, not a guessed E7 normal form.
        u_value = solve_linear_coordinate(surface, U, {Z: 0}, chart_ring)(Y=s)
        values = {Z: e, U: series_field(u_value), Y: s}
        g = series_field(e**4 * s**6)
    elif component == "E7_7":
        assert frame["component_equation"] == "U"
        z_value = solve_linear_coordinate(surface, Z, {U: 0}, chart_ring)(Y=s)
        values = {Z: series_field(z_value), U: e, Y: s}
        g = series_field(e**5 * s**6)
    elif component in ("E7_5", "E7_6"):
        assert frame["component_equation"] == "Z"
        u_value = conic_z_normalization(surface)
        values = {Z: e, U: u_value, Y: s*u_value}
        g = series_field(
            e**6 * (s*u_value)**5 if component == "E7_5"
            else e**3 * (s*u_value)**6
        )
    else:
        raise ValueError("unsupported component")
    pullback = {
        name: series_field(value.subs(values)) for name, value in old.items()
    }
    return pullback, g, values


component_payloads = []
selected_components = tuple(args.component or ("E7_4", "E7_7"))
for component in selected_components:
    pullback, g, normalization = chart_pullbacks(component)
    t_value, x_value, y_value = (pullback[name] for name in ("t", "x", "y"))
    t_order, _ = e_leading(t_value, e, s)
    x_order, _ = e_leading(x_value, e, s)
    y_order, _ = e_leading(y_value, e, s)
    expected = condition_by_component[component]["orders"]
    chart_orders = frame_by_component[component]["old_coordinate_weight_orders"]
    assert (t_order, x_order, y_order) == (
        expected["t"], expected["x"], int(chart_orders["y"])
    )
    # At E7_5 both P1 entrance leading orders meet the old coordinate orders;
    # retain the first two/one exact P1 coefficients. The other selected
    # components have strict inequalities and therefore use m=y/x.
    m_value = (
        (y_value-y_p3*t_value**3)
        / (x_value-x_p2*t_value**2-x_p3*t_value**3)
        if component == "E7_5" else y_value/x_value
    )
    m_order, _ = e_leading(m_value, e, s)
    assert m_order == expected["m"]

    rows = []
    for group in condition_by_component[component]["negative_order_groups"]:
        if len(group["basis_indices"]) == 1:
            continue
        residues = []
        for index in group["basis_indices"]:
            entry = basis[index]
            degree = 4*int(entry["h_power"])-int(entry["u_power"])
            value = (g * t_value**(degree-9) * x_value**int(entry["x_power"])
                     * m_value**int(entry["m_power"]) / h_leading**int(entry["h_power"]))
            order, residue = e_leading(value, e, s)
            assert order == int(group["residual_order"])
            residues.append((index, residue))
        for s_power, entries in sorted(residue_coefficient_rows(residues, s).items()):
            if len(entries) > 1:
                rows.append({
                    "residual_order": int(group["residual_order"]),
                    "component_parameter_power": int(s_power),
                    "entries": entries,
                    "interpretation": "the displayed coefficient-weighted sum must vanish",
                })
            else:
                # The all-component valuation compiler did not call this a
                # singleton because other terms had the same normal order,
                # but their distinct component residues cannot cancel it.
                rows.append({
                    "residual_order": int(group["residual_order"]),
                    "component_parameter_power": int(s_power),
                    "entries": entries,
                    "interpretation": "a distinct component-parameter power; this coefficient must vanish",
                })
    component_payloads.append({
        "component": component,
        "chart": frame_by_component[component]["actual_edge_chart"],
        "normalization": {str(name): str(value) for name, value in normalization.items()},
        "cartier_factor": str(g),
        "orders": {"t": int(t_order), "x": int(x_order), "y": int(y_order), "m": int(m_order)},
        "non_singleton_residue_rows": rows,
    })

payload = {
    "schema": "elkies-k3.h92-q8-e7-4-7-generic-residue-rows.v1",
    "status": "PASS_EXACT_Q8_E7_4_7_GENERIC_RESIDUE_ROWS",
    "inputs": {
        "component_frames": {"path": path_label(args.frames), "sha256": digest(args.frames)},
        "generic_conditions": {"path": path_label(args.conditions), "sha256": digest(args.conditions)},
        "p1": {"path": path_label(args.p1), "sha256": digest(args.p1)},
    },
    "ambient_basis_sha256": conditions["ambient_basis_sha256"],
    "components": component_payloads,
    "compiler_instruction": (
        "Stack these residue rows with the singleton coordinate rows. Continue "
        "with E7_1,E7_2,E7_3,E7_5,E7_6, edge nodes, the marked branch, and "
        "overlap conditions before claiming a complete E7 cover."
    ),
    "boundary": (
        "This is an exact generic residue calculation on two actual rational "
        "E7 components only. It is not a complete q8 chart cover, global kernel, "
        "pencil, child equation, or rank statement."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q8E747GENERICRESIDUES|components=2|rows={}|status="
    "PASS_EXACT_Q8_E7_4_7_GENERIC_RESIDUE_ROWS".format(
        sum(len(item["non_singleton_residue_rows"]) for item in component_payloads)
    ),
    flush=True,
)
