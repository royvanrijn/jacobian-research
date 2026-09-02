#!/usr/bin/env sage-python
"""Exhaust the normalized NS0028 ``I3+I7+I8`` fibre ansatz modulo p.

The selected exact source ``NS0028-S001`` has root type ``A2+A6+A7``,
Mordell--Weil rank two, and two independent pole-zero generators.  This first
gate imposes only the three semistable reducible fibres, normalized at zero,
one, and infinity.  It scans every degree-eight short-Weierstrass A polynomial
when ``--max-samples`` is zero and solves the overdetermined Hermite problem
for B exactly.

Stored models prove only finite-field feasibility of the fibre stratum.  The
two polynomial sections, the full NS0028 marking, characteristic-zero descent,
and an elliptic-neighbour corridor remain separate gates.
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
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-ns0028-source-ansatz-mod5.json"
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
parser.add_argument("--source-id", default="NS0028-S001")
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--prime", type=int, default=5)
parser.add_argument(
    "--examples",
    type=int,
    default=20,
    help="number of exact models to retain; zero retains every model",
)
parser.add_argument(
    "--max-samples",
    type=int,
    default=0,
    help="truncate the deterministic coefficient scan; zero exhausts it",
)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
if arguments.examples < 0 or arguments.max_samples < 0:
    parser.error("examples and max-samples must be nonnegative")

source_path = arguments.source.resolve()
output_path = arguments.output.resolve()
source_payload = json.loads(source_path.read_text())
source_entry = next(
    row
    for row in source_payload["sources"]
    if row["ns_id"] == "NS0028" and row["source_id"] == arguments.source_id
)
source = source_entry["source"]
assert source["root_type"] == "A2+A6+A7"
assert source["root_rank"] == 15
assert source["mw_rank_for_rho_19"] == 2
assert source["root_lattice_primitive"] and source["torsion"] == 1
assert source["mw_height_gram"] in (
    [["52/21", "-1"], ["-1", "25/8"]],
    [["52/21", "1"], ["1", "25/8"]],
)

field = GF(arguments.prime)
assert field.characteristic() not in (2, 3)
ring = PolynomialRing(field, "t")
t = ring.gen()

# Three jets at zero, seven at one, and eight at infinity impose the
# I3+I7+I8 profile.  B has thirteen coefficients, so five compatibility
# equations remain on the normalized degree-eight A polynomial.
rows = []
for jet in range(3):
    rows.append([field(index == jet) for index in range(13)])
for jet in range(7):
    rows.append(
        [
            field(binomial(index, jet)) if index >= jet else field.zero()
            for index in range(13)
        ]
    )
for jet in range(8):
    rows.append([field(index == 12 - jet) for index in range(13)])
hermite = matrix(field, rows)
assert hermite.nrows() == 18 and hermite.ncols() == hermite.rank() == 13
compatibility_matrix = hermite.left_kernel().basis_matrix()
assert compatibility_matrix.nrows() == 5

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
    at_zero = a_coefficients[:3]
    at_one = [
        sum(
            a_coefficients[index] * field(binomial(index, jet))
            for index in range(jet, 9)
        )
        for jet in range(7)
    ]
    at_infinity = [a_coefficients[8 - jet] for jet in range(8)]
    positive = (
        multiplicative_branch(at_zero, 1),
        multiplicative_branch(at_one, 1),
        multiplicative_branch(at_infinity, 1),
    )
    if any(branch is None for branch in positive):
        continue
    for sign_one in (1, -1):
        for sign_infinity in (1, -1):
            branches = (
                positive[0],
                [sign_one * value for value in positive[1]],
                [sign_infinity * value for value in positive[2]],
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
            if orders != (3, 7, 8):
                continue
            exact_orders += 1
            divisor = t**3 * (t - 1) ** 7
            residual, remainder = discriminant_core.quo_rem(divisor)
            assert not remainder and residual.degree() == 6
            if residual(0) == 0 or residual(1) == 0:
                continue
            if residual.gcd(residual.derivative()).degree() != 0:
                continue
            squarefree += 1
            if not arguments.examples or len(examples) < arguments.examples:
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
                            "0": 3,
                            "1": 7,
                            "infinity": 8,
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
                        "geometric_fibre_profile": "I3+I7+I8+6I1",
                    }
                )

exhausted = sample == exhaustive_total and not arguments.max_samples
payload = {
    "schema": "elkies-k3.lattice-foundry-ns0028-source-ansatz-modp.v1",
    "status": (
        "PASS_EXACT_EXHAUSTIVE_MODULAR_SOURCE_FIBRE_ANSATZ"
        if exhausted and squarefree
        else (
            "PASS_BOUNDED_MODULAR_SOURCE_FIBRE_ANSATZ"
            if squarefree
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
        "normalized_reducible_supports": ["0:I3", "1:I7", "infinity:I8"],
        "hermite_conditions": 18,
        "B_coefficient_rank": 13,
        "compatibility_equations_on_A": 5,
        "expected_fibre_stratum_dimension": 3,
        "expected_NS0028_MW2_locus_dimension": 1,
        "expected_MW_conditions_still_missing": 2,
        "section_marking": {
            "P": "pole zero; depth one at I3 and I7; identity at I8",
            "Q": "pole zero; identity at I3 and I7; depth one at I8",
            "height_gram": source["mw_height_gram"],
        },
    },
    "examples": examples,
    "source": {
        "artifact": display_path(source_path),
        "artifact_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "source_id": arguments.source_id,
        "source_gram_sha256": source["gram_sha256"],
    },
    "proof_boundary": {
        "proved": (
            "Every stored example is an exact short Weierstrass K3 model over "
            "the displayed finite field with fibre profile I3+I7+I8+6I1."
        ),
        "not_proved": (
            "The two pole-zero MW sections, full NS0028 lattice marking, rational "
            "parameterization, characteristic-zero lifting, and neighbour route "
            "are not proved."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/probe_lattice_foundry_ns0028_source_ansatz_modp.sage "
        f"--prime {arguments.prime} --examples {arguments.examples}"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("NS0028 modular source-ansatz artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "FOUNDRYNS0028ANSATZ|"
    f"p={arguments.prime}|samples={sample}|compatible={compatible}|"
    f"squarefree={squarefree}|stored={len(examples)}|exhausted={int(exhausted)}|"
    f"status={'PASS' if squarefree else 'BOUNDED_NEGATIVE'}",
    flush=True,
)
