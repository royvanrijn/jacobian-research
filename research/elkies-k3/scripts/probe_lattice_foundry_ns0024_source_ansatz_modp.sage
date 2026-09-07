#!/usr/bin/env sage -python
"""Construct the A6+A4+A3 source fibre ansatz over a small finite field.

status: ACTIVE_SEARCH
claim: exact modular feasibility of the semistable I7+I5+I4+8I1 stratum
inputs: an exact NS0024 root-rank-13 source certificate
output: caller-selected bounded-search JSON artifact

This is the first equation gate, not an identification of the NS0024 source.
Four independent Mordell--Weil sections and the marked lattice embedding still
have to be imposed before a modular surface belongs to that Picard-19 locus.
"""

import argparse
import hashlib
import json
import random
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ, binomial, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-source-hunt-r13.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-source-ansatz-mod11.json"
)


def display_path(path):
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def local_square_root(unit_coefficients, root0):
    """Square root of a truncated unit series with prescribed constant sign."""
    field = root0.parent()
    answer = [field.zero()] * len(unit_coefficients)
    answer[0] = root0
    for degree in range(1, len(answer)):
        known = sum(
            answer[left] * answer[degree - left]
            for left in range(1, degree)
        )
        answer[degree] = (unit_coefficients[degree] - known) / (2 * root0)
    return answer


def truncated_product(left, right, precision):
    field = left[0].parent()
    answer = [field.zero()] * precision
    for i, left_value in enumerate(left[:precision]):
        for j, right_value in enumerate(right[: precision - i]):
            answer[i + j] += left_value * right_value
    return answer


def multiplicative_branch(a_series, sign=1):
    """Return B=2*(-A/3)^(3/2) to the supplied precision."""
    field = a_series[0].parent()
    u = [-value / field(3) for value in a_series]
    if not u[0] or not u[0].is_square():
        return None
    root0 = u[0].sqrt()
    if sign == -1:
        root0 = -root0
    h = local_square_root(u, root0)
    h2 = truncated_product(h, h, len(h))
    h3 = truncated_product(h2, h, len(h))
    return [2 * value for value in h3]


def order_at(poly, point):
    if not poly:
        return None
    shifted = poly(poly.parent().gen() + point)
    return min(index for index, value in enumerate(shifted.list()) if value)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--prime", type=int, default=11)
parser.add_argument("--seed", type=int, default=20260901)
parser.add_argument("--max-samples", type=int, default=100000)
parser.add_argument("--examples", type=int, default=3)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

SOURCE = args.source.resolve()
OUTPUT = args.output.resolve()
source = json.loads(SOURCE.read_text())
assert source["status"] == "PASS_EXACT_NEW_K3_ROOTFUL_MW4_SOURCE_AND_NIEMEIER_CERTIFICATE"
assert source["source"]["root_type"] == "A3+A4+A6"
assert source["source"]["mw_rank_for_rho_19"] == 4

field = GF(args.prime)
assert field.characteristic() not in (2, 3, 5, 7)
ring = PolynomialRing(field, "t")
t = ring.gen()
rng = random.Random(args.seed)

# Rows map the thirteen B coefficients to the seven jets at zero, five at
# one, and four at infinity.  The rank is thirteen, so the right-hand branch
# data obey exactly three compatibility equations.
rows = []
for jet in range(7):
    rows.append([field(index == jet) for index in range(13)])
for jet in range(5):
    rows.append([
        field(binomial(index, jet)) if index >= jet else field.zero()
        for index in range(13)
    ])
for jet in range(4):
    rows.append([field(index == 12 - jet) for index in range(13)])
hermite = matrix(field, rows)
assert hermite.nrows() == 16 and hermite.ncols() == hermite.rank() == 13

examples = []
branch_eligible = 0
compatible = 0
exact_orders = 0
for sample in range(1, args.max_samples + 1):
    a_coefficients = [field(-3)] + [field(rng.randrange(args.prime)) for _ in range(8)]
    if not a_coefficients[8]:
        continue
    A = ring(a_coefficients)

    at_zero = a_coefficients[:7]
    at_one = [
        sum(
            a_coefficients[index] * field(binomial(index, jet))
            for index in range(jet, 9)
        )
        for jet in range(5)
    ]
    at_infinity = [a_coefficients[8 - jet] for jet in range(4)]
    branches = (
        multiplicative_branch(at_zero),
        multiplicative_branch(at_one),
        multiplicative_branch(at_infinity),
    )
    if any(branch is None for branch in branches):
        continue
    branch_eligible += 1
    target = vector(field, branches[0] + branches[1] + branches[2])
    try:
        b_coefficients = list(hermite.solve_right(target))
    except ValueError:
        continue
    compatible += 1
    B = ring(b_coefficients)
    discriminant_core = 4 * A**3 + 27 * B**2
    order_zero = order_at(discriminant_core, field.zero())
    order_one = order_at(discriminant_core, field.one())
    order_infinity = 24 - discriminant_core.degree()
    if (order_zero, order_one, order_infinity) != (7, 5, 4):
        continue
    exact_orders += 1
    divisor = t**7 * (t - 1)**5
    residual, remainder = discriminant_core.quo_rem(divisor)
    assert not remainder and residual.degree() == 8
    if residual(0) == 0 or residual(1) == 0:
        continue
    if residual.gcd(residual.derivative()).degree() != 0:
        continue
    factorization = [
        {"degree": int(factor.degree()), "multiplicity": int(multiplicity)}
        for factor, multiplicity in residual.factor()
    ]
    examples.append(
        {
            "sample_index": sample,
            "A_coefficients_low_to_high": [int(value) for value in a_coefficients],
            "B_coefficients_low_to_high": [int(value) for value in b_coefficients],
            "discriminant_orders": {"0": 7, "1": 5, "infinity": 4},
            "residual_discriminant_coefficients_low_to_high": [
                int(value) for value in residual
            ],
            "residual_factorization": factorization,
            "residual_squarefree_over_algebraic_closure": True,
            "geometric_fibre_profile": "I7+I5+I4+8I1",
        }
    )
    if len(examples) >= args.examples:
        break

payload = {
    "schema": "elkies-k3.lattice-foundry-ns0024-source-ansatz-modp.v1",
    "status": (
        "PASS_EXACT_MODULAR_SOURCE_FIBRE_ANSATZ"
        if examples else "PASS_BOUNDED_NO_MODULAR_SOURCE_FIBRE_ANSATZ"
    ),
    "prime": args.prime,
    "seed": args.seed,
    "search_bound": args.max_samples,
    "samples_consumed": sample,
    "accounting": {
        "branch_eligible": branch_eligible,
        "hermite_compatible": compatible,
        "exact_prescribed_orders": exact_orders,
        "accepted_squarefree_examples": len(examples),
    },
    "ansatz": {
        "short_weierstrass": "y^2=x^3+A(t)x+B(t)",
        "degree_bounds": {"A": 8, "B": 12},
        "normalized_reducible_supports": ["0:I7", "1:I5", "infinity:I4"],
        "hermite_conditions": 16,
        "B_coefficient_rank": 13,
        "compatibility_equations_on_A": 3,
        "fibre_stratum_dimension_after_scaling_and_base_normalization": 5,
        "expected_NS0024_MW4_locus_dimension": 1,
        "expected_MW_conditions_still_missing": 4,
    },
    "examples": examples,
    "source_certificate": display_path(SOURCE),
    "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    "proof_boundary": {
        "proved": (
            "Every emitted example is an exact short Weierstrass K3 model over the "
            "displayed finite field with geometric fibre profile I7+I5+I4+8I1."
        ),
        "not_proved": (
            "The four MW sections, the NS0024 lattice embedding, characteristic-zero "
            "lifting, and equation-side marking of the neighbour route are not proved."
        ),
    },
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if args.check:
    assert OUTPUT.read_text() == serialized
else:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(serialized)

print(
    "FOUNDRYANSATZ|p={}|samples={}|compatible={}|examples={}|profile=I7+I5+I4+8I1|status={}".format(
        args.prime,
        sample,
        compatible,
        len(examples),
        "PASS" if examples else "BOUNDED_NEGATIVE",
    ),
    flush=True,
)
