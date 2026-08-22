#!/usr/bin/env sage -python
"""Exhaust polynomial sections of the exact H92 q=6 child modulo a prime.

The q=8-to-D13 lattice class has old-fibre degree two *on the q=6 child*.
Its two smallest Mordell--Weil directions have height 8/3, which is
compatible with a section disjoint from the zero section.  On a minimal K3
Weierstrass model such a section can be represented by polynomial coordinates
of degrees at most four and six.  This script exhausts that small ansatz over
``GF(p)`` as reconnaissance for an equation-level q=8 construction.

It is deliberately an experiment.  A modular polynomial section neither
proves that it lifts to characteristic zero nor identifies its resolved E8/E6
components, its Mordell--Weil class, or the selected q=8 divisor.
"""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-polynomial-sections-mod-7.json"


def coefficient(field, value):
    value = QQ(value)
    numerator = field(ZZ(value.numerator()))
    denominator = field(ZZ(value.denominator()))
    if not denominator:
        raise ValueError("chosen prime divides a child-model coefficient denominator")
    return numerator / denominator


def polynomial(ring, field, coefficients):
    return ring([coefficient(field, value) for value in coefficients])


def coefficient_list(value):
    return [int(entry) for entry in value.list()]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def vanishing_order(value, factor):
    if not value:
        return None
    result = 0
    while value % factor == 0:
        value //= factor
        result += 1
    return result


def has_order_at_least(value, factor, order):
    value_order = vanishing_order(value, factor)
    return value_order is None or value_order >= order


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--prime", type=int, default=7)
parser.add_argument("--max-x-degree", type=int, default=4)
parser.add_argument(
    "--require-iv-star-singular", action="store_true",
    help=(
        "restrict to the singular IV* branch x divisible by f^2 and "
        "y divisible by f^2; this is a bounded experimental filter"
    ),
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

prime = ZZ(args.prime)
if not prime.is_prime() or prime in (2, 3):
    raise ValueError("prime must be an odd prime different from 3")
if args.max_x_degree < 0:
    raise ValueError("max-x-degree must be nonnegative")

child = json.loads(args.child.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
finite = GF(prime)
ring = PolynomialRing(finite, "T")
T = ring.gen()
source_ring = PolynomialRing(QQ, "T")
coefficient_a = polynomial(
    ring, finite, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"]
)
coefficient_b = polynomial(
    ring, finite, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"]
)

# Reject a bad reduction before searching.  In particular, a reduction at
# which the additive factors or their Kodaira valuations collapse produces
# many irrelevant polynomial solutions of a different surface.
additive = {}
for item in child["finite_fibres"]:
    if item["kodaira"] not in ("II*", "IV*"):
        continue
    source_factor = source_ring(item["factor"])
    factor = polynomial(ring, finite, source_factor.list()).monic()
    if factor.degree() != 1:
        raise ValueError("additive fibre factor did not retain degree one modulo the chosen prime")
    additive[item["kodaira"]] = factor
if set(additive) != {"II*", "IV*"}:
    raise ArithmeticError("child input does not provide exactly the expected II* and IV* factors")
discriminant = 4 * coefficient_a**3 + 27 * coefficient_b**2
expected_orders = {"II*": (4, 5, 10), "IV*": (3, 4, 8)}
for kind, factor in additive.items():
    orders = (
        vanishing_order(coefficient_a, factor),
        vanishing_order(coefficient_b, factor),
        vanishing_order(discriminant, factor),
    )
    if orders != expected_orders[kind]:
        raise ValueError(
            "chosen prime does not retain the {} fibre valuations: {}".format(
                kind, orders
            )
        )

# The degree-six y bound is the standard K3 polynomial-section bound.  The
# x degree may be tightened for a small reconnaissance run, but enlarging it
# is an explicit change of the bounded experiment rather than a proof step.
max_y_degree = 6
sections = []
iv_factor = additive["IV*"]
if args.require_iv_star_singular:
    x_degree_after_iv = args.max_x_degree - 2
    if x_degree_after_iv < 0:
        raise ValueError("IV* singular-branch filter requires max-x-degree at least two")
    x_values = (
        iv_factor**2 * ring(entries)
        for entries in itertools.product(finite, repeat=x_degree_after_iv + 1)
    )
else:
    x_values = (
        ring(entries)
        for entries in itertools.product(finite, repeat=args.max_x_degree + 1)
    )
for x_value in x_values:
    right_side = x_value**3 + coefficient_a * x_value + coefficient_b
    if not right_side.is_square():
        continue
    y_value = right_side.sqrt()
    if y_value.degree() > max_y_degree:
        continue
    signs = (y_value,) if not y_value else (y_value, -y_value)
    for signed_y in signs:
        if signed_y and signed_y.degree() > max_y_degree:
            continue
        if args.require_iv_star_singular and not (
            has_order_at_least(x_value, iv_factor, 2)
            and has_order_at_least(signed_y, iv_factor, 2)
        ):
            continue
        assert signed_y**2 == right_side
        sections.append({
            "x_coefficients_low_to_high": coefficient_list(x_value),
            "y_coefficients_low_to_high": coefficient_list(signed_y),
        })

sections.sort(key=lambda entry: (
    entry["x_coefficients_low_to_high"], entry["y_coefficients_low_to_high"]
))

for entry in sections:
    x_value = ring(entry["x_coefficients_low_to_high"])
    y_value = ring(entry["y_coefficients_low_to_high"])
    entry["additive_fibre_vanishing_orders"] = {
        kind: {
            "x": vanishing_order(x_value, factor),
            "y": vanishing_order(y_value, factor),
        }
        for kind, factor in sorted(additive.items())
    }

payload = {
    "schema": "elkies-k3.h92-q6-child-polynomial-sections-modp.v1",
    "status": "EXPERIMENTAL_EXHAUSTIVE_MODULAR_ANSATZ",
    "inputs": {
        "child_jacobian": str(args.child.relative_to(ROOT)),
        "child_jacobian_sha256": digest(args.child),
        "child_status": child["status"],
    },
    "prime": int(prime),
    "ansatz": {
        "x_degree_at_most": int(args.max_x_degree),
        "y_degree_at_most": max_y_degree,
        "iv_star_singular_filter": bool(args.require_iv_star_singular),
        "x_space_size": int(
            prime ** (args.max_x_degree - 1)
            if args.require_iv_star_singular else prime ** (args.max_x_degree + 1)
        ),
        "complete_over_this_field_and_ansatz": True,
    },
    "sections": sections,
    "section_count": len(sections),
    "boundary": (
        "This exhausts only polynomial sections in the stated finite-field "
        "ansatz. It does not establish a characteristic-zero lift, resolved "
        "component labels, a q=8 pencil, a rootless bisection, an extension "
        "collision, or generic rank 18 or 19."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDPOLYMOD|prime={}|x_degree={}|iv_singular={}|x_space={}|sections={}|"
    "status=EXPERIMENTAL_EXHAUSTIVE_MODULAR_ANSATZ".format(
        prime, args.max_x_degree, int(args.require_iv_star_singular),
        prime ** (args.max_x_degree - 1)
        if args.require_iv_star_singular else prime ** (args.max_x_degree + 1),
        len(sections)
    ),
    flush=True,
)
