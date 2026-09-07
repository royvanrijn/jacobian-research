#!/usr/bin/env sage-python
"""Exhaust the normalized G720-S0260 ``A11+A4/MW2`` fibre ansatz mod p.

The canonical determinant-720 source has semistable fibres ``I12+I5``, two
reducible supports, trivial torsion, and a complete pole-zero MW basis.  This
first equation gate fixes the supports at zero and infinity and normalizes
``A(0)=-3`` in a short Weierstrass model.  It exhausts every degree-eight A
polynomial, solves the 17-jet Hermite problem for B, and retains exactly the
``I12+I5+7I1`` models.  The two polynomial sections are imposed separately.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-octad-det720-prescribed-root-sources-v1.json"
)
DEFAULT_POLES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-octad-det720-source-poles-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-a11-a4-source-ansatz-mod5-v1.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    unit = [-value / field(3) for value in a_series]
    if not unit[0] or not unit[0].is_square():
        return None
    root0 = unit[0].sqrt()
    if sign == -1:
        root0 = -root0
    square_root = local_square_root(unit, root0)
    return [
        2 * value
        for value in truncated_product(
            truncated_product(square_root, square_root, len(square_root)),
            square_root,
            len(square_root),
        )
    ]


def order_at_zero(poly):
    return next(
        (index for index, value in enumerate(poly.list()) if value), None
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
parser.add_argument("--poles", type=Path, default=DEFAULT_POLES)
parser.add_argument("--source-id", default="G720-S0260")
parser.add_argument("--prime", type=int, default=5)
parser.add_argument("--examples", type=int, default=0)
parser.add_argument("--max-samples", type=int, default=0)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
if arguments.examples < 0 or arguments.max_samples < 0:
    parser.error("examples and max-samples must be nonnegative")
if arguments.prime in (2, 3):
    parser.error("short Weierstrass scan requires characteristic other than 2 or 3")

sources_path = arguments.sources.resolve()
poles_path = arguments.poles.resolve()
output_path = arguments.output.resolve()
sources = json.loads(sources_path.read_text())
poles = json.loads(poles_path.read_text())
source_entry = next(
    row for row in sources["sources"] if row["source_id"] == arguments.source_id
)
source = source_entry["source"]
if (
    source["root_type"] != "A11+A4"
    or source["root_rank"] != 15
    or source["mw_rank_for_rho_19"] != 2
    or source["support_count"] != 2
    or not source["root_lattice_primitive"]
    or source["torsion"] != 1
):
    raise ValueError("selected source is not the primitive A11+A4/MW2 target")
pole_entry = next(
    row
    for row in poles["audits"]
    if row["source_id"] == arguments.source_id
    and row["source_gram_sha256"] == source["gram_sha256"]
)
audit = pole_entry["audit"]
if (
    audit["basis_sorted_pole_profile"] != [0, 0]
    or audit["torsion_order"] != 1
    or len(audit["basis"]) != 2
    or any(any(section["simple_root_pairings"]) for section in audit["basis"])
):
    raise ValueError("selected source lost its identity-component pole-zero basis")

field = GF(arguments.prime)
ring = PolynomialRing(field, "t")
t = ring.gen()

# B has thirteen coefficients.  Twelve jets at zero and five at infinity give
# a rank-thirteen Hermite system and four compatibility equations on A.
rows = []
for jet in range(12):
    rows.append([field(index == jet) for index in range(13)])
for jet in range(5):
    rows.append([field(index == 12 - jet) for index in range(13)])
hermite = matrix(field, rows)
if hermite.nrows() != 17 or hermite.ncols() != hermite.rank() != 13:
    raise ArithmeticError("unexpected I12+I5 Hermite rank")
compatibility = hermite.left_kernel().basis_matrix()
if compatibility.nrows() != 4:
    raise ArithmeticError("unexpected I12+I5 compatibility codimension")

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
    at_zero = a_coefficients[:12] + [field.zero()] * max(0, 12 - len(a_coefficients))
    at_zero = at_zero[:12]
    at_infinity = [a_coefficients[8 - jet] for jet in range(5)]
    zero_branch = multiplicative_branch(at_zero, 1)
    infinity_branch = multiplicative_branch(at_infinity, 1)
    if zero_branch is None or infinity_branch is None:
        continue
    for sign_infinity in (1, -1):
        branch_eligible += 1
        target = vector(
            field,
            zero_branch + [sign_infinity * value for value in infinity_branch],
        )
        if compatibility * target:
            continue
        b_coefficients = list(hermite.solve_right(target))
        compatible += 1
        B = ring(b_coefficients)
        discriminant_core = 4 * A**3 + 27 * B**2
        order_zero = order_at_zero(discriminant_core)
        order_infinity = (
            None if not discriminant_core else 24 - discriminant_core.degree()
        )
        if (order_zero, order_infinity) != (12, 5):
            continue
        exact_orders += 1
        residual, remainder = discriminant_core.quo_rem(t**12)
        if remainder or residual.degree() != 7 or residual[0] == 0:
            raise ArithmeticError("exact fibre orders gave the wrong residual")
        if residual.gcd(residual.derivative()).degree() != 0:
            continue
        squarefree += 1
        if not arguments.examples or len(examples) < arguments.examples:
            examples.append(
                {
                    "sample_index": sample,
                    "branch_sign_at_infinity": sign_infinity,
                    "A_coefficients_low_to_high": [
                        int(value) for value in a_coefficients
                    ],
                    "B_coefficients_low_to_high": [
                        int(value) for value in b_coefficients
                    ],
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
                }
            )

exhausted = sample == exhaustive_total and not arguments.max_samples
output = {
    "schema": "elkies-k3.golay-det720-a11-a4-source-ansatz-modp.v1",
    "status": (
        "PASS_EXACT_EXHAUSTIVE_MODULAR_FIBRE_ANSATZ"
        if exhausted and squarefree
        else (
            "PASS_BOUNDED_MODULAR_FIBRE_ANSATZ"
            if squarefree
            else "PASS_BOUNDED_NO_MODULAR_FIBRE_ANSATZ"
        )
    ),
    "prime": arguments.prime,
    "inputs": {
        relative(sources_path): digest(sources_path),
        relative(poles_path): digest(poles_path),
    },
    "source": {
        "source_id": arguments.source_id,
        "source_gram_sha256": source["gram_sha256"],
        "root_type": source["root_type"],
        "mw_rank": source["mw_rank_for_rho_19"],
        "mw_height_gram": source["mw_height_gram"],
        "basis_sorted_pole_profile": audit["basis_sorted_pole_profile"],
        "basis_on_identity_components": True,
    },
    "ansatz": {
        "short_weierstrass": "y^2=x^3+A(t)x+B(t)",
        "degree_bounds": {"A": 8, "B": 12},
        "normalization": "A(0)=-3; I12 at zero; I5 at infinity",
        "geometric_fibre_profile": "I12+I5+7I1",
        "hermite_conditions": 17,
        "B_coefficient_rank": 13,
        "compatibility_equations_on_A": 4,
        "expected_fibre_stratum_dimension_before_residual_coordinate_quotients": 4,
        "expected_MW_conditions_still_missing": 2,
    },
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
    "examples": examples,
    "proof_boundary": {
        "proved": (
            "Every stored example is an exact short-Weierstrass model over the "
            "displayed finite field with fibre profile I12+I5+7I1; with an "
            "exhaustive scan every normalized A polynomial and both relative "
            "local B branches are covered."
        ),
        "not_proved": (
            "The two pole-zero sections, full determinant-720 source marking, "
            "characteristic-zero lift, rational parameterization, target "
            "multisection spectrum, and neighbour corridor are not proved."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/probe_golay_det720_a11_a4_source_ansatz_modp.sage "
        f"--prime {arguments.prime} --examples {arguments.examples}"
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("Golay-720 A11+A4 fibre artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "GOLAY720A11A4FIBRE|"
    f"p={arguments.prime}|samples={sample}|compatible={compatible}|"
    f"squarefree={squarefree}|stored={len(examples)}|exhausted={int(exhausted)}|"
    f"status={'PASS' if squarefree else 'EMPTY'}",
    flush=True,
)
