#!/usr/bin/env sage-python
"""Compile section-first Tate charts and replay the Golay/NS0031 controls.

The reusable polynomial chart is

    y^2 + a1*x*y + a3*y = x^3 + a2*x^2,

with P=(0,0).  For a prescribed affine intersection divisor h, put

    Q=(h*r, h^2*s),
    R=h*(r^3-h*s^2-a1*r*s).

If alpha*s+beta*r^2=1, then

    a3=alpha*R+k*r^2,  a2=-beta*R+k*s

makes Q a section identically.  The coprimality conditions
gcd(r,h)=gcd(r,s)=1 make the affine intersection divisor of P and Q exactly
h.  Fibre searches can therefore work only with discriminant jets, while
resolved-component, minimality, and saturation gates remain separate.

The script also translates the exact Golay QQ pair and the marked NS0031
GF(7) pair into this chart.  These are algebraic positive controls; it does
not promote either source to a saturated characteristic-zero foundry model.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_GOLAY = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-3a5-source-qq-v1.json"
)
DEFAULT_NS0031 = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-marked-gf7-hensel-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-section-first-normal-form-controls-v1.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialized_polynomial(poly):
    return [str(value) for value in poly.list()]


def rational_signature(value):
    value = value.parent()(value)
    numerator = value.numerator()
    denominator = value.denominator()
    return {
        "numerator_degree": int(numerator.degree()),
        "denominator_degree": int(denominator.degree()),
        "sha256": hashlib.sha256(str(value).encode()).hexdigest(),
    }


def tate_invariants(a1, a2, a3):
    """Return the standard invariants for a4=a6=0."""
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3
    b6 = a3**2
    b8 = a2 * a3**2
    c4 = b2**2 - 24 * b4
    c6 = -(b2**3) + 36 * b2 * b4 - 216 * b6
    discriminant = -(b2**2) * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    return {
        "b2": b2,
        "b4": b4,
        "b6": b6,
        "b8": b8,
        "c4": c4,
        "c6": c6,
        "discriminant": discriminant,
    }


def tate_residual(a1, a2, a3, x_coordinate, y_coordinate):
    return (
        y_coordinate**2
        + a1 * x_coordinate * y_coordinate
        + a3 * y_coordinate
        - x_coordinate**3
        - a2 * x_coordinate**2
    )


def compile_two_section_chart(a1, h, r, s, k):
    """Compile a polynomial two-section chart from coprime input data."""
    ring = a1.parent()
    inputs = [ring(value) for value in (h, r, s, k)]
    h, r, s, k = inputs
    if r.gcd(h).degree() != 0 or r.gcd(s).degree() != 0:
        raise ValueError("two-section chart requires gcd(r,h)=gcd(r,s)=1")
    gcd_value, alpha, beta = s.xgcd(r**2)
    if gcd_value.degree() != 0:
        raise ValueError("two-section chart requires gcd(s,r^2)=1")
    inverse_unit = gcd_value[0] ** (-1)
    alpha *= inverse_unit
    beta *= inverse_unit
    if alpha * s + beta * r**2 != 1:
        raise ArithmeticError("normalized Bezout identity failed")
    relation = h * (r**3 - h * s**2 - a1 * r * s)
    a3 = alpha * relation + k * r**2
    a2 = -beta * relation + k * s
    point_q = (h * r, h**2 * s)
    if tate_residual(a1, a2, a3, 0, 0) != 0:
        raise ArithmeticError("the first marked point was not built in")
    if tate_residual(a1, a2, a3, *point_q) != 0:
        raise ArithmeticError("the second marked point was not built in")
    if point_q[0].gcd(point_q[1]).monic() != h.monic():
        raise ArithmeticError("the compiled affine intersection divisor changed")
    return {
        "a1": a1,
        "a2": a2,
        "a3": a3,
        "P": (ring.zero(), ring.zero()),
        "Q": point_q,
        "h": h.monic(),
        "r": r,
        "s": s,
        "k": k,
        "alpha": alpha,
        "beta": beta,
        "relation": relation,
        "invariants": tate_invariants(a1, a2, a3),
    }


def as_polynomial(value, ring):
    value = ring.fraction_field()(value)
    denominator = value.denominator()
    if denominator.degree() > 0:
        raise ArithmeticError("expected invariant cancellation to a polynomial")
    return ring(value.numerator() / denominator[0])


def exact_order_at(poly, support):
    ring = poly.parent()
    t = ring.gen()
    shifted = ring(poly(t + support))
    for index, coefficient in enumerate(shifted.list()):
        if coefficient:
            return index
    return None


def semistable_fibre_conditions(
    delta, c4, support, order, total_delta_degree=24, total_c4_degree=8
):
    """Compile the closed jets and open gates for a prescribed I_n fibre."""
    if order < 1:
        raise ValueError("a semistable fibre order must be positive")
    if support == "infinity":
        closed = [
            delta[index]
            for index in range(total_delta_degree - order + 1, total_delta_degree + 1)
        ]
        exact_order_gate = delta[total_delta_degree - order]
        c4_unit_gate = c4[total_c4_degree]
    else:
        ring = delta.parent()
        t = ring.gen()
        shifted_delta = ring(delta(t + delta.base_ring()(support)))
        closed = [shifted_delta[index] for index in range(order)]
        exact_order_gate = shifted_delta[order]
        c4_unit_gate = c4(delta.base_ring()(support))
    return {
        "closed_discriminant_jets": closed,
        "exact_order_open_gate": exact_order_gate,
        "c4_unit_open_gate": c4_unit_gate,
    }


def semistable_fibre_gate(delta, c4, expected, total_delta_degree=24, total_c4_degree=8):
    """Verify exact I_n discriminant orders and c4-unit open gates."""
    rows = []
    for support, order in expected:
        conditions = semistable_fibre_conditions(
            delta,
            c4,
            support,
            order,
            total_delta_degree=total_delta_degree,
            total_c4_degree=total_c4_degree,
        )
        closed_zero = not any(conditions["closed_discriminant_jets"])
        exact_order = bool(conditions["exact_order_open_gate"])
        c4_unit = bool(conditions["c4_unit_open_gate"])
        actual = (
            total_delta_degree - delta.degree()
            if support == "infinity"
            else exact_order_at(delta, delta.base_ring()(support))
        )
        if not closed_zero or not exact_order or actual != order or not c4_unit:
            raise ArithmeticError(
                f"semistable fibre gate failed at {support}: order={actual}, c4_unit={c4_unit}"
            )
        rows.append(
            {
                "support": str(support),
                "kodaira": f"I{order}",
                "discriminant_order": int(actual),
                "closed_discriminant_jet_count": int(order),
                "exact_order_open_gate": True,
                "c4_unit": True,
            }
        )
    return rows


def translate_short_to_tate(A, B, first_point, second_point):
    """Put a non-2-torsion first point at (0,0) with tangent y=0."""
    ring = A.parent()
    field = ring.fraction_field()
    xp, yp = [field(value) for value in first_point]
    xq, yq = [field(value) for value in second_point]
    A = field(A)
    B = field(B)
    if yp == 0:
        raise ValueError("the selected first section is 2-torsion in the function field")
    if yp**2 != xp**3 + A * xp + B:
        raise ArithmeticError("first point is not on the short Weierstrass model")
    if yq**2 != xq**3 + A * xq + B:
        raise ArithmeticError("second point is not on the short Weierstrass model")
    slope = (3 * xp**2 + A) / (2 * yp)
    a1 = 2 * slope
    a3 = 2 * yp
    a2 = 3 * xp - slope**2
    u = xq - xp
    v = yq - yp - slope * u
    if tate_residual(a1, a2, a3, 0, 0) != 0:
        raise ArithmeticError("translated first point failed")
    if tate_residual(a1, a2, a3, u, v) != 0:
        raise ArithmeticError("translated second point failed")
    invariants = tate_invariants(a1, a2, a3)
    if invariants["c4"] != -48 * A:
        raise ArithmeticError("c4 changed under a unit Weierstrass translation")
    if invariants["c6"] != -864 * B:
        raise ArithmeticError("c6 changed under a unit Weierstrass translation")
    if invariants["discriminant"] != -16 * (4 * A**3 + 27 * B**2):
        raise ArithmeticError("discriminant changed under a unit Weierstrass translation")
    return {
        "a1": a1,
        "a2": a2,
        "a3": a3,
        "P": (field.zero(), field.zero()),
        "Q": (u, v),
        "slope": slope,
        "invariants": invariants,
    }


def affine_intersection_polynomial(first_point, second_point, ring):
    field = ring.fraction_field()
    xp, yp = [field(value) for value in first_point]
    xq, yq = [field(value) for value in second_point]
    dx = xq - xp
    dy = yq - yp
    h = ring(dx.numerator()).gcd(ring(dy.numerator())).monic()
    if h.degree() < 0:
        raise ArithmeticError("the marked sections coincide generically")
    if h.gcd(ring(dx.denominator())).degree() > 0:
        raise ArithmeticError("x-coordinate pole contaminates the intersection divisor")
    if h.gcd(ring(dy.denominator())).degree() > 0:
        raise ArithmeticError("y-coordinate pole contaminates the intersection divisor")
    return h


def control_record(name, base_field, ring, A, B, first_point, second_point, expected_fibres):
    field = ring.fraction_field()
    D = ring(4 * A**3 + 27 * B**2)
    if field(first_point[1]) ** 2 != field(first_point[0]) ** 3 + A * first_point[0] + B:
        raise ArithmeticError(f"{name}: first short-Weierstrass section failed")
    if field(second_point[1]) ** 2 != field(second_point[0]) ** 3 + A * second_point[0] + B:
        raise ArithmeticError(f"{name}: second short-Weierstrass section failed")
    h = affine_intersection_polynomial(first_point, second_point, ring)
    if h.gcd(D).degree() > 0:
        raise ArithmeticError(f"{name}: marked intersection meets a singular fibre")
    translated = translate_short_to_tate(A, B, first_point, second_point)
    delta = as_polynomial(translated["invariants"]["discriminant"], ring)
    c4 = as_polynomial(translated["invariants"]["c4"], ring)
    fibre_rows = semistable_fibre_gate(delta, c4, expected_fibres)
    u, v = translated["Q"]
    r = u / field(h)
    s = v / field(h**2)
    if ring(r.numerator()).gcd(h).degree() > 0:
        raise ArithmeticError(f"{name}: x-direction is not a unit along the intersection divisor")
    for value, label in ((r, "r"), (s, "s"), (translated["slope"], "tangent slope")):
        if ring(value.denominator()).gcd(h).degree() > 0:
            raise ArithmeticError(f"{name}: {label} has a pole along the intersection divisor")
    relation = field(h) * (
        r**3 - field(h) * s**2 - translated["a1"] * r * s
    )
    # Over the function field choose alpha=1/s, beta=0.  The recovered k
    # specializes the polynomial Bezout compiler and must reproduce a2,a3.
    k = (translated["a3"] - relation / s) / r**2
    recovered_a2 = k * s
    recovered_a3 = relation / s + k * r**2
    if recovered_a2 != translated["a2"] or recovered_a3 != translated["a3"]:
        raise ArithmeticError(f"{name}: two-section compiler failed to recover Tate coefficients")
    return {
        "name": name,
        "base_field": base_field,
        "short_weierstrass_section_identities": "PASS",
        "tate_equation": "y^2+a1*x*y+a3*y=x^3+a2*x^2",
        "first_marked_section": "P=(0,0)",
        "second_marked_section": "Q=(h*r,h^2*s)",
        "intersection_polynomial_monic": serialized_polynomial(h),
        "smooth_affine_intersection_degree": int(h.degree()),
        "intersection_disjoint_from_discriminant": True,
        "tangent_chart_regular_along_intersection": True,
        "tate_invariant_replay": {
            "c4_equals_short_c4": True,
            "c6_equals_short_c6": True,
            "delta_equals_minus_16_discriminant_core": True,
        },
        "two_section_parameterization_recovered": True,
        "translated_data_signatures": {
            "a1": rational_signature(translated["a1"]),
            "a2": rational_signature(translated["a2"]),
            "a3": rational_signature(translated["a3"]),
            "r": rational_signature(r),
            "s": rational_signature(s),
            "k": rational_signature(k),
        },
        "semistable_fibre_gates": fibre_rows,
    }


def load_golay(path: Path):
    payload = json.loads(path.read_text())
    if payload.get("schema") != "elkies-k3.golay-det720-3a5-source-qq.v1":
        raise ValueError("unexpected Golay control schema")
    if payload.get("status") != "PASS_EXACT_QQ_3I6_MW2_RANK19_SUBLATTICE_DET720":
        raise ValueError("Golay control did not pass its source certificate")
    ring = PolynomialRing(QQ, "t")
    model = payload["weierstrass_model"]
    A = ring(model["A_coefficients_low_to_high"])
    B = ring(model["B_coefficients_low_to_high"])
    points = []
    for section in payload["marked_sections"]:
        points.append(
            (
                ring(section["X_coefficients_low_to_high"]),
                ring(section["Y_coefficients_low_to_high"]),
            )
        )
    record = control_record(
        "Golay G720-S0128 rational specialization",
        "QQ",
        ring,
        A,
        B,
        points[0],
        points[1],
        [(0, 6), (1, 6), ("infinity", 6)],
    )
    if record["intersection_polynomial_monic"] != payload["section_pair"][
        "smooth_intersection_polynomial_monic"
    ]:
        raise ArithmeticError("Golay intersection polynomial changed")
    return record


def load_ns0031(path: Path):
    payload = json.loads(path.read_text())
    if payload.get("schema") != "elkies-k3.lattice-foundry-ns0031-marked-gf7-hensel.v1":
        raise ValueError("unexpected NS0031 control schema")
    if payload.get("prime") != 7 or payload.get("status") != (
        "PASS_ONE_DIMENSIONAL_MARKED_TANGENT_AND_EXPLICIT_Z7_LIFT_TO_REQUESTED_PRECISION"
    ):
        raise ValueError("NS0031 control did not pass its marked GF(7) gate")
    field = GF(7)
    ring = PolynomialRing(field, "t")
    t = ring.gen()
    function_field = ring.fraction_field()
    values = payload["seed"]["coordinates_mod_7"]
    if len(values) != 52:
        raise ArithmeticError("NS0031 marked seed changed size")
    if payload["seed"].get("fibre_example_index") != 157:
        raise ArithmeticError("NS0031 positive-control model index changed")
    A = ring(values[0:9])
    B = ring(values[9:22])
    X_P = ring(values[22:27])
    Y_P = ring(values[27:34])
    C = t + field(values[34])
    N_R = ring(values[35:42])
    M_R = ring(values[42:52])
    first_point = (function_field(X_P), function_field(Y_P))
    second_point = (
        function_field(N_R) / C**2,
        function_field(M_R) / C**3,
    )
    record = control_record(
        "NS0031-S001 marked model 157",
        "GF(7)",
        ring,
        A,
        B,
        first_point,
        second_point,
        [(0, 2), (1, 8), ("infinity", 8)],
    )
    if record["smooth_affine_intersection_degree"] != 2:
        raise ArithmeticError("NS0031 marked pair lost its required intersection degree")
    if record["intersection_polynomial_monic"] != ["2", "3", "1"]:
        raise ArithmeticError("NS0031 marked intersection polynomial changed")
    return record


def polynomial_compiler_smoke_test():
    ring = PolynomialRing(QQ, "t")
    t = ring.gen()
    chart = compile_two_section_chart(
        a1=t**2 + t + 1,
        h=t**2 + t + 1,
        r=1,
        s=1,
        k=t**4 + 2,
    )
    delta = chart["invariants"]["discriminant"]
    if delta == 0:
        raise ArithmeticError("polynomial compiler smoke test produced a singular generic fibre")
    degree_bounds = {"a1": 2, "a2": 4, "a3": 6, "discriminant": 24}
    actual_degrees = {
        "a1": int(chart["a1"].degree()),
        "a2": int(chart["a2"].degree()),
        "a3": int(chart["a3"].degree()),
        "discriminant": int(delta.degree()),
    }
    if any(actual_degrees[name] > bound for name, bound in degree_bounds.items()):
        raise ArithmeticError("polynomial compiler smoke test exceeded elliptic-K3 degree bounds")
    fibre_rows = semistable_fibre_gate(
        delta, chart["invariants"]["c4"], [("infinity", 4)]
    )
    return {
        "status": "PASS_EXACT_POLYNOMIAL_IDENTITY",
        "input_degrees": {name: int(chart[name].degree()) for name in ("h", "r", "s", "k")},
        "elliptic_K3_degree_bounds": degree_bounds,
        "output_degrees": actual_degrees,
        "intersection_polynomial_monic": serialized_polynomial(chart["h"]),
        "first_section_residual": "0",
        "second_section_residual": "0",
        "compiled_fibre_condition_example": fibre_rows[0],
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--golay", type=Path, default=DEFAULT_GOLAY)
parser.add_argument("--ns0031", type=Path, default=DEFAULT_NS0031)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

golay_path = arguments.golay.resolve()
ns0031_path = arguments.ns0031.resolve()
output_path = arguments.output.resolve()
payload = {
    "schema": "elkies-k3.section-first-normal-form-controls.v1",
    "status": "PASS_SECTION_FIRST_MW1_MW2_IDENTITIES_AND_POSITIVE_CONTROLS",
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/compile_section_first_normal_forms.sage"
    ),
    "inputs": {
        relative(golay_path): digest(golay_path),
        relative(ns0031_path): digest(ns0031_path),
    },
    "one_marked_section_template": {
        "equation": "y^2+a1*x*y+a3*y=x^3+a2*x^2",
        "marked_section": "P=(0,0)",
        "section_equations_remaining": 0,
        "search_variables": "degree-bounded a1,a2,a3 coefficients",
        "fibre_tuning": "impose discriminant jets; retain exact-order, c4-unit, minimality, splitness and component gates",
    },
    "two_marked_section_template": {
        "equation": "y^2+a1*x*y+a3*y=x^3+a2*x^2",
        "marked_sections": ["P=(0,0)", "Q=(h*r,h^2*s)"],
        "assumptions": ["gcd(r,h)=1", "gcd(r,s)=1", "alpha*s+beta*r^2=1"],
        "relation": "R=h*(r^3-h*s^2-a1*r*s)",
        "compiled_coefficients": ["a3=alpha*R+k*r^2", "a2=-beta*R+k*s"],
        "section_equations_remaining": 0,
        "affine_intersection_divisor": "h",
        "fibre_tuning": "impose discriminant jets after coefficient compilation",
    },
    "polynomial_compiler_smoke_test": polynomial_compiler_smoke_test(),
    "positive_controls": [load_golay(golay_path), load_ns0031(ns0031_path)],
    "proof_boundary": {
        "proved": (
            "The displayed one- and two-marked-section formulas are exact algebraic identities. "
            "The Golay QQ and NS0031 GF(7) marked pairs translate into the two-section chart with "
            "their smooth intersection divisor and semistable discriminant orders preserved."
        ),
        "not_proved": (
            "Built-in points do not prove exact Mordell-Weil rank, independence, resolved component "
            "profiles, global minimality, Neron-Severi saturation, rational algebraization of NS0031, "
            "or membership of the rational Golay specialization in the determinant-720 lattice."
        ),
    },
}

rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != rendered:
        raise SystemExit("section-first normal-form control artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(rendered)

print(
    "SECTIONFIRST|mw1=P00|mw2=P00,Q=hr_h2s|"
    "controls=Golay_QQ,NS0031_GF7|status=PASS",
    flush=True,
)
