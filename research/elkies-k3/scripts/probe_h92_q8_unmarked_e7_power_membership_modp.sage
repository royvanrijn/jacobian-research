#!/usr/bin/env sage -python
"""Test cleared q=8 E7 power-module membership modulo a good prime.

For a source q8 term f=u^i*x^a*m^b/h^k, with
Lx=x-x(P1), Ly=y-y(P1), the actual E7 comparison is

    Lx^9*g*f in t^9*(Lx,Ly)^9.

On an actual edge chart all P1 and reversed-h denominators are units at the
edge origin.  This program clears those units and all negative t powers, then
uses a Groebner membership calculation in the corresponding localized chart
ring.  It is a finite-field chart calculation, not a characteristic-zero
q8-pencil certificate.
"""

import argparse
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
GLUING = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-actual-e7-gluing.json"
KERNEL = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra4.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-power-membership-mod-43-extra4.json"


def coefficient(field, value):
    value = QQ(value)
    denominator = field(ZZ(value.denominator()))
    if not denominator:
        raise ValueError("the chosen prime divides an input denominator")
    return field(ZZ(value.numerator())) / denominator


def polynomial(ring, field, values):
    return ring([coefficient(field, value) for value in values])


def strict(ring, value, substitutions, exceptional):
    transformed = ring(value(*substitutions))
    quotient, remainder = transformed.quo_rem(exceptional**2)
    assert not remainder
    return quotient


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--kernel", type=Path, default=KERNEL)
parser.add_argument("--gluing", type=Path, default=GLUING)
parser.add_argument(
    "--chart", action="append", default=[],
    help="actual chart name; repeat to select a subset (default: all five unmarked charts)",
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if not ZZ(args.prime).is_prime() or args.prime in (2, 3):
    raise ValueError("prime must be odd and different from 3")

finite = GF(args.prime)
kernel = json.loads(args.kernel.read_text())
gluing = json.loads(args.gluing.read_text())
p1 = json.loads(P1.read_text())
assert kernel["status"] == "EXPERIMENTAL_MODULAR_SMOOTH_BLOCK_RANK"
assert kernel["prime"] == args.prime
assert "kernel_basis_rows" in kernel
assert gluing["status"] == "PASS_EXACT_Q8_ACTUAL_E7_GLUING"
assert p1["status"] == "PASS_EXACT_H92_P1"

# The input kernel may be a different declared endpoint enlargement, but it
# must still be a degree-one/degree-zero generic q8 ambient with h powers.
basis = kernel["ambient_basis"]
for entry in basis:
    assert int(entry["x_power"]) in (0, 1)
    assert 0 <= int(entry["m_power"]) <= 9
    assert int(entry["h_power"]) >= 0

# Reconstruct the six exact H92 edge charts over the selected finite field.
ring = PolynomialRing(finite, names=("Z", "U", "Y", "S"), order="degrevlex")
Z, U, Y, S = ring.gens()
anchor = SourceFileLoader("h92_q8_power_module_anchor", str(ANCHOR)).load_module()
r92, s92 = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (coefficient(finite, value(r92, s92)) for value in formulas)
f0 = Y**2-U**3-(A1*Z**3+A*Z**4)*U-(B1*Z**5+B*Z**6+B2*Z**7)
f1_z = strict(ring, f0, (Z, Z*U, Z*Y, S), Z)
f2_u = strict(ring, f1_z, (U*Z, U, U*Y, S), U)
f3_u = strict(ring, f2_u, (U*Z, U, U*Y, S), U)
f3_z = strict(ring, f2_u, (Z, Z*U, Z*Y, S), Z)
second = -A1/B1
third = -finite.one()/A1
charts = {
    "E7_2--E7_5": (
        strict(ring, ring(f2_u(Z+second, U, Y, S)), (Z, Z*U, Z*Y, S), Z),
        (Z*U*(Z+second), Z**2*U**2*(Z+second), Z**3*U**2*(Z+second)*Y),
    ),
    "E7_1--E7_4": (
        strict(ring, f3_u, (U*Z, U, U*Y, S), U),
        (U**3*Z, U**4*Z, U**6*Z*Y),
    ),
    "E7_4--E7_3": (
        strict(ring, f3_u, (Z, Z*U, Z*Y, S), Z),
        (Z**3*U**2, Z**4*U**3, Z**6*U**4*Y),
    ),
    "E7_3--E7_7": (
        strict(ring, f3_z, (U*Z, U, U*Y, S), U),
        (U**3*Z**2, U**5*Z**3, U**7*Z**4*Y),
    ),
    "E7_7--E7_2": (
        strict(ring, f3_z, (Z, Z*U, Z*Y, S), Z),
        (Z**3*U, Z**5*U**2, Z**7*U**2*Y),
    ),
    "E7_3--E7_6": (
        strict(ring, ring(f3_u(Z+third, U, Y, S)), (Z, Z*U, Z*Y, S), Z),
        (Z**2*U**2*(Z+third), Z**3*U**3*(Z+third), Z**5*U**4*(Z+third)*Y),
    ),
}
assert set(charts) == {entry["name"] for entry in gluing["actual_edge_chart_gluing"]}

# Turn the P1 input into exact rational functions of t=1/u.  At t=0 their
# denominators are units; the chart-localization below explicitly inverts all
# such denominators rather than assuming a polynomial representative.
t_ring = PolynomialRing(finite, "t")
t = t_ring.gen()


def reversed_fraction(numerator_values, denominator_values):
    numerator = polynomial(t_ring, finite, numerator_values)
    denominator = polynomial(t_ring, finite, denominator_values)
    rev_numerator = sum(numerator[index]*t**(numerator.degree()-index) for index in range(numerator.degree()+1))
    rev_denominator = sum(denominator[index]*t**(denominator.degree()-index) for index in range(denominator.degree()+1))
    shift = denominator.degree()-numerator.degree()
    assert shift >= 0
    return t**shift*rev_numerator, rev_denominator


xp_num, xp_den = reversed_fraction(
    p1["x_entrance_base"]["numerator_coefficients"],
    p1["x_entrance_base"]["denominator_coefficients"],
)
yp_num, yp_den = reversed_fraction(
    p1["y_entrance_base"]["numerator_coefficients"],
    p1["y_entrance_base"]["denominator_coefficients"],
)
h = polynomial(t_ring, finite, p1["structured_denominator"]["Z4_coefficients"])
h_reverse = sum(h[index]*t**(h.degree()-index) for index in range(h.degree()+1))
assert xp_den(0) and yp_den(0) and h_reverse(0)

gluing_by_name = {entry["name"]: entry for entry in gluing["actual_edge_chart_gluing"]}
selected = tuple(args.chart) if args.chart else tuple(name for name in charts if name != "E7_2--E7_5")
if not set(selected) <= set(charts):
    raise ValueError("unknown chart")


def g_factor(edge):
    value = ring.one()
    for component in edge["components"]:
        variable = {"Z": Z, "U": U, "Y": Y}[component["equation"]]
        value *= variable**int(component["w_coefficient"])
    return value


results = []
for name in selected:
    surface, (t_map, x_map, y_map) = charts[name]
    xp_n = ring(xp_num(t_map))
    xp_d = ring(xp_den(t_map))
    yp_n = ring(yp_num(t_map))
    yp_d = ring(yp_den(t_map))
    h_r = ring(h_reverse(t_map))
    assert xp_d(0, 0, 0, 0) and yp_d(0, 0, 0, 0) and h_r(0, 0, 0, 0)

    # Nx/xp_d=x-x(P1) and Ny/yp_d=y-y(P1); denominators are chart units.
    nx = x_map*xp_d-xp_n
    ny = y_map*yp_d-yp_n
    assert nx and ny
    g = g_factor(gluing_by_name[name])
    nonzero_labels = [
        index for index in range(len(basis))
        if any(int(row[index]) % args.prime for row in kernel["kernel_basis_rows"])
    ]
    common_h_power = max(int(basis[index]["h_power"]) for index in nonzero_labels)
    common_t_clear = max(max(int(basis[index]["u_power"])-4*int(basis[index]["h_power"]), 0) for index in nonzero_labels)
    unit_denominator = h_r**common_h_power * xp_d**9 * yp_d**9
    assert unit_denominator(0, 0, 0, 0)

    # Localize exactly at the known unit denominator.  The surface relation
    # and the ten power-ideal generators are all retained in one ideal.
    ideal_generators = [surface, S*unit_denominator-ring.one()]
    ideal_generators.extend(
        t_map**(common_t_clear+9) * nx**(9-j) * ny**j
        for j in range(10)
    )
    local_ideal = ring.ideal(ideal_generators)
    groebner = local_ideal.groebner_basis()
    vector_results = []
    for row in kernel["kernel_basis_rows"]:
        cleared = ring.zero()
        for coefficient_value, entry in zip(row, basis):
            c = finite(int(coefficient_value))
            if not c:
                continue
            a, b = int(entry["x_power"]), int(entry["m_power"])
            i, k = int(entry["u_power"]), int(entry["h_power"])
            exponent = common_t_clear + 4*k-i
            assert exponent >= 0
            cleared += c * g * x_map**a * nx**(9-b) * ny**b * t_map**exponent * h_r**(common_h_power-k) * xp_d**b * yp_d**(9-b)
        remainder = cleared.reduce(groebner)
        vector_results.append({
            "cleared_membership": remainder == 0,
            "remainder_terms": len(remainder.dict()),
            "remainder_total_degree": int(remainder.total_degree()) if remainder else -1,
        })
    results.append({
        "chart": name,
        "common_h_power": common_h_power,
        "common_t_clear": common_t_clear,
        "localized_ideal_generators": 12,
        "kernel_vectors": vector_results,
    })
    print(
        "H92Q8E7POWERPROBE|chart={}|vectors={}|pass={}|status=EXPERIMENTAL_MODULAR_E7_POWER_MEMBERSHIP".format(
            name, len(vector_results), sum(result["cleared_membership"] for result in vector_results)
        ),
        flush=True,
    )

payload = {
    "schema": "elkies-k3.h92-q8-unmarked-e7-power-membership-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_E7_POWER_MEMBERSHIP",
    "prime": args.prime,
    "inputs": {
        "smooth_kernel": str(args.kernel.relative_to(ROOT)),
        "q8_gluing": str(args.gluing.relative_to(ROOT)),
        "p1": str(P1.relative_to(ROOT)),
    },
    "clearing_identity": "(x-xP1)^9*g*f belongs to t^9*(x-xP1,y-yP1)^9",
    "localization": "invert the product of reversed h and P1 denominators, each a unit at the selected edge origin",
    "results": results,
    "boundary": "This checks selected finite-field chart germs only. It does not establish compatible overlaps, a complete E7 condition matrix, a characteristic-zero q8 kernel, a pencil, or any rank claim.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
