#!/usr/bin/env sage-python
"""Construct the NS0011 A2+A6+A8 source fibre ansatz over a finite field.

The selected exact source has semistable-compatible fibre profile
``I3+I7+I9+5I1``, MW rank one, and minimum section pole order two.  This first
equation gate imposes only the three reducible fibres.  The default exhaustive
``GF(5)`` run scans every normalized degree-eight ``A`` polynomial, solves the
overdetermined Hermite system for ``B``, and retains exact squarefree residual
discriminants.

An emitted model proves modular feasibility of the two-dimensional fibre
stratum.  It does not yet impose the pole-two MW section, the NS0011 marking,
or a characteristic-zero lift.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, binomial, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-a-v1.json"
)
DEFAULT_POLES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-rank1-section-poles-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0011-source-ansatz-mod5.json"
)


def display_path(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def local_square_root(unit_coefficients, root0):
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
    """Return ``B=2*(-A/3)^(3/2)`` to the supplied precision."""
    field = a_series[0].parent()
    u = [-value / field(3) for value in a_series]
    if not u[0] or not u[0].is_square():
        return None
    root0 = u[0].sqrt()
    if sign == -1:
        root0 = -root0
    h = local_square_root(u, root0)
    return [
        2 * value
        for value in truncated_product(
            truncated_product(h, h, len(h)), h, len(h)
        )
    ]


def order_at(poly, point):
    if not poly:
        return None
    shifted = poly(poly.parent().gen() + point)
    return min(index for index, value in enumerate(shifted.list()) if value)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
parser.add_argument("--source-id", default="NS0011-S005")
parser.add_argument("--section-poles", type=Path, default=DEFAULT_POLES)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--prime", type=int, default=5)
parser.add_argument("--examples", type=int, default=3)
parser.add_argument(
    "--max-samples",
    type=int,
    default=0,
    help="truncate the deterministic coefficient scan; zero exhausts it",
)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

source_path = arguments.source.resolve()
pole_path = arguments.section_poles.resolve()
output_path = arguments.output.resolve()
source_payload = json.loads(source_path.read_text())
source_entry = next(
    row
    for row in source_payload["sources"]
    if row["ns_id"] == "NS0011" and row["source_id"] == arguments.source_id
)
source = source_entry["source"]
assert source["root_type"] == "A2+A6+A8"
assert source["root_rank"] == 16
assert source["mw_rank_for_rho_19"] == 1
assert source["root_lattice_primitive"] and source["torsion"] == 1
pole_payload = json.loads(pole_path.read_text())
pole_row = next(
    row
    for row in pole_payload["sources"]
    if row["source_artifact"] == display_path(source_path)
    and row["source_id"] == arguments.source_id
)
assert pole_row["minimum_section_pole_order"] == 2

field = GF(arguments.prime)
assert field.characteristic() not in (2, 3)
ring = PolynomialRing(field, "t")
t = ring.gen()

# The nineteen branch jets are nine at zero, seven at one, and three at
# infinity.  A degree-twelve B polynomial has thirteen coefficients, leaving
# six compatibility equations on the normalized degree-eight A polynomial.
rows = []
for jet in range(9):
    rows.append([field(index == jet) for index in range(13)])
for jet in range(7):
    rows.append(
        [
            field(binomial(index, jet)) if index >= jet else field.zero()
            for index in range(13)
        ]
    )
for jet in range(3):
    rows.append([field(index == 12 - jet) for index in range(13)])
hermite = matrix(field, rows)
assert hermite.nrows() == 19 and hermite.ncols() == hermite.rank() == 13
compatibility_matrix = hermite.left_kernel().basis_matrix()
assert compatibility_matrix.nrows() == 6

examples = []
branch_eligible = 0
compatible = 0
exact_orders = 0
squarefree = 0
sample = 0
exhaustive_total = arguments.prime**8
for digits in itertools.product(range(arguments.prime), repeat=8):
    sample += 1
    if arguments.max_samples and sample > arguments.max_samples:
        sample -= 1
        break
    a_coefficients = [field(-3)] + [field(value) for value in digits]
    if not a_coefficients[8]:
        continue
    A = ring(a_coefficients)
    at_zero = a_coefficients[:9]
    at_one = [
        sum(
            a_coefficients[index] * field(binomial(index, jet))
            for index in range(jet, 9)
        )
        for jet in range(7)
    ]
    at_infinity = [a_coefficients[8 - jet] for jet in range(3)]
    branches_positive = (
        multiplicative_branch(at_zero, 1),
        multiplicative_branch(at_one, 1),
        multiplicative_branch(at_infinity, 1),
    )
    if any(branch is None for branch in branches_positive):
        continue
    for sign_one in (1, -1):
        for sign_infinity in (1, -1):
            branches = (
                branches_positive[0],
                [sign_one * value for value in branches_positive[1]],
                [sign_infinity * value for value in branches_positive[2]],
            )
            branch_eligible += 1
            target = vector(field, branches[0] + branches[1] + branches[2])
            if compatibility_matrix * target:
                continue
            b_coefficients = list(hermite.solve_right(target))
            compatible += 1
            B = ring(b_coefficients)
            discriminant_core = 4 * A**3 + 27 * B**2
            orders = (
                order_at(discriminant_core, field.zero()),
                order_at(discriminant_core, field.one()),
                24 - discriminant_core.degree(),
            )
            if orders != (9, 7, 3):
                continue
            exact_orders += 1
            divisor = t**9 * (t - 1) ** 7
            residual, remainder = discriminant_core.quo_rem(divisor)
            assert not remainder and residual.degree() == 5
            if residual(0) == 0 or residual(1) == 0:
                continue
            if residual.gcd(residual.derivative()).degree() != 0:
                continue
            squarefree += 1
            if len(examples) < arguments.examples:
                examples.append(
                    {
                        "sample_index": sample,
                        "branch_signs_at_one_and_infinity": [
                            sign_one,
                            sign_infinity,
                        ],
                        "A_coefficients_low_to_high": [
                            int(value) for value in a_coefficients
                        ],
                        "B_coefficients_low_to_high": [
                            int(value) for value in b_coefficients
                        ],
                        "discriminant_orders": {
                            "0": 9,
                            "1": 7,
                            "infinity": 3,
                        },
                        "residual_discriminant_coefficients_low_to_high": [
                            int(value) for value in residual
                        ],
                        "residual_factorization": [
                            {
                                "degree": int(factor.degree()),
                                "multiplicity": int(multiplicity),
                            }
                            for factor, multiplicity in residual.factor()
                        ],
                        "geometric_fibre_profile": "I9+I7+I3+5I1",
                    }
                )

exhausted = sample == exhaustive_total and not arguments.max_samples
payload = {
    "schema": "elkies-k3.lattice-foundry-ns0011-source-ansatz-modp.v1",
    "status": (
        "PASS_EXACT_EXHAUSTIVE_MODULAR_SOURCE_FIBRE_ANSATZ"
        if exhausted and examples
        else (
            "PASS_BOUNDED_MODULAR_SOURCE_FIBRE_ANSATZ"
            if examples
            else "PASS_BOUNDED_NO_MODULAR_SOURCE_FIBRE_ANSATZ"
        )
    ),
    "prime": arguments.prime,
    "scan": {
        "normalized_A_polynomials": exhaustive_total,
        "samples_consumed": sample,
        "exhausted": exhausted,
    },
    "accounting": {
        "branch_eligible_with_signs": branch_eligible,
        "hermite_compatible_with_signs": compatible,
        "exact_prescribed_orders": exact_orders,
        "squarefree_examples_with_signs": squarefree,
        "stored_examples": len(examples),
    },
    "ansatz": {
        "short_weierstrass": "y^2=x^3+A(t)x+B(t)",
        "degree_bounds": {"A": 8, "B": 12},
        "normalization": "A(0)=-3; supports at 0,1,infinity",
        "normalized_reducible_supports": ["0:I9", "1:I7", "infinity:I3"],
        "hermite_conditions": 19,
        "B_coefficient_rank": 13,
        "compatibility_equations_on_A": 6,
        "expected_fibre_stratum_dimension": 2,
        "expected_NS0011_MW1_locus_dimension": 1,
        "expected_MW_conditions_still_missing": 1,
        "minimum_section_pole_order_to_impose": 2,
    },
    "examples": examples,
    "source": {
        "artifact": display_path(source_path),
        "artifact_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "source_id": arguments.source_id,
        "source_gram_sha256": source["gram_sha256"],
        "section_pole_artifact": display_path(pole_path),
        "section_pole_artifact_sha256": hashlib.sha256(
            pole_path.read_bytes()
        ).hexdigest(),
    },
    "proof_boundary": {
        "proved": (
            "Every stored example is an exact short Weierstrass K3 model over "
            "the displayed finite field with fibre profile I9+I7+I3+5I1."
        ),
        "not_proved": (
            "The MW1 pole-two section, NS0011 lattice marking, rational "
            "parameterization, characteristic-zero lifting, and neighbour route "
            "are not proved."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/probe_lattice_foundry_ns0011_source_ansatz_modp.sage"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("NS0011 modular source-ansatz artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "FOUNDRYNS0011ANSATZ|"
    f"p={arguments.prime}|samples={sample}|compatible={compatible}|"
    f"squarefree={squarefree}|exhausted={int(exhausted)}|"
    f"status={'PASS' if examples else 'BOUNDED_NEGATIVE'}",
    flush=True,
)
