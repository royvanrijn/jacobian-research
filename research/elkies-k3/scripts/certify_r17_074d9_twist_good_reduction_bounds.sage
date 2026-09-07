#!/usr/bin/env sage-python
"""Certify Frobenius rank bounds for the four 074d9 singleton twists.

The regular model ``q(u)y^2=x^3+A(u)x+B(u)`` is a nondegenerate toric
hypersurface at each declared prime.  Toric controlled reduction returns its
degree-28 primitive ``det(1-T*Frob | H^2)`` polynomial and Hodge vector
``(2,24,2)``.  The six split ambient toric divisor classes restore full
degree 34.

Cyclotomic roots after scaling by p, together with the six-dimensional toric
complement, bound the geometric Mordell--Weil rank after subtracting the
geometric degree-10 trivial lattice ``U+D4+D4``.  No arithmetic Frobenius
action is assigned to the toric complement or the nonsplit fibre components.
"""

from __future__ import annotations

import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, divisors, factor
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
MODEL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
)
LOCAL_ROOT = (
    ROOT / "artifacts/local/elkies-k3/r17-074d9-singleton-toric-frobenius"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-074d9-twist-good-reduction-bounds-v1.json"
)
LABELS = (
    "074d9-orbit-04b07",
    "074d9-orbit-11a44",
    "074d9-orbit-11279",
    "074d9-orbit-080fa",
)
PRIMES = (131, 137, 151, 157, 167, 173, 181, 193)
TORIC_CONTROLLED_REDUCTION_COMMIT = "74cda9e8148cd8e9a3928fc15a558c9a70b67cc1"


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def coefficients_text(polynomial) -> list[str]:
    return [str(value) for value in polynomial.list()]


def parse_driver_output(path: Path, expected_label: str, prime: int):
    fields = path.read_text().strip().split(":")
    if len(fields) != 8:
        raise ValueError(f"{path}: unexpected controlled-reduction output")
    label = fields[0]
    output_prime = int(fields[5])
    hodge = ast.literal_eval(fields[6])
    zeta = ast.literal_eval(fields[7])
    if label != expected_label or output_prime != prime:
        raise ValueError(f"{path}: output label or prime mismatch")
    if hodge != [2, 24, 2] or len(zeta) != 29:
        raise ArithmeticError(f"{path}: primitive cohomology dimensions changed")
    ring = PolynomialRing(ZZ, "T")
    return hodge, ring(list(reversed(zeta)))


def cyclotomic_scaled_degree(polynomial, prime: int) -> int:
    """Multiplicity of eigenvalues p*zeta with zeta a root of unity."""

    ring = polynomial.parent()
    total = 0
    for irreducible, multiplicity in factor(polynomial):
        scaled_coefficients = []
        for degree, coefficient in enumerate(irreducible.list()):
            divisor = ZZ(prime) ** degree
            if coefficient % divisor:
                break
            scaled_coefficients.append(coefficient // divisor)
        else:
            scaled = ring(scaled_coefficients)
            if scaled.is_cyclotomic():
                total += int(multiplicity * irreducible.degree())
    return total


def factor_rows(polynomial):
    return [
        {
            "coefficients_low_to_high": coefficients_text(irreducible),
            "degree": int(irreducible.degree()),
            "multiplicity": int(multiplicity),
        }
        for irreducible, multiplicity in factor(polynomial)
    ]


def component_cycles(model, input_record, prime: int):
    representative = model["representative"]
    q_mod_p = list(map(int, input_record["quadratic_coefficients_low_to_high_mod_p"]))
    fixed_counts = {}
    for extension_degree in range(1, 7):
        constants = GF(prime**extension_degree, name=f"a{extension_degree}")
        base_ring = PolynomialRing(constants, "u")
        u = base_ring.gen()
        q = base_ring(q_mod_p)
        A = base_ring(
            [constants(QQ(value)) for value in representative["A_coefficients_low_to_high"]]
        )
        B = base_ring(
            [constants(QQ(value)) for value in representative["B_coefficients_low_to_high"]]
        )
        cubic_ring = PolynomialRing(constants, "x")
        x = cubic_ring.gen()
        fixed = 0
        for alpha, multiplicity in q.roots():
            if multiplicity != 1:
                raise ArithmeticError("the twist branch is not squarefree")
            # Central D4 component plus the outer components indexed by E[2]-{O}.
            fixed += 1 + len((x**3 + A(alpha) * x + B(alpha)).roots())
        fixed_counts[extension_degree] = fixed

    cycles = {}
    for length in range(1, 7):
        residual = fixed_counts[length] - sum(
            divisor * cycles[divisor]
            for divisor in divisors(length)
            if divisor < length
        )
        if residual % length:
            raise ArithmeticError("component fixed counts do not define a permutation")
        cycles[length] = residual // length
    cycles = {length: count for length, count in cycles.items() if count}
    if sum(length * count for length, count in cycles.items()) != 8:
        raise ArithmeticError("the two D4 component permutations do not have degree eight")
    return fixed_counts, cycles


def multiplicity_of(polynomial, divisor):
    multiplicity = 0
    quotient = polynomial
    while True:
        next_quotient, remainder = quotient.quo_rem(divisor)
        if remainder:
            return multiplicity
        multiplicity += 1
        quotient = next_quotient


def build_payload():
    model = json.loads(MODEL.read_text())
    if model.get("status") != "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
        raise ValueError("unexpected lineage-model status")

    twist_records = []
    all_raw_inputs = {relative(MODEL): digest(MODEL)}
    for label in LABELS:
        tag = label.removeprefix("074d9-orbit-")
        reduction_rows = []
        for prime in PRIMES:
            directory = LOCAL_ROOT / tag / f"p{prime}"
            input_certificate_path = directory / "input-certificate.json"
            raw_input_path = directory / "toric-controlled-reduction.input"
            raw_output_path = directory / "toric-controlled-reduction.output"
            input_record = json.loads(input_certificate_path.read_text())
            if (
                input_record.get("status") != "PASS_EXACT_TORIC_FROBENIUS_INPUT_EXPORT"
                or input_record.get("label") != label
                or int(input_record.get("prime")) != prime
            ):
                raise ValueError("stale toric input certificate")
            if input_record["toric_input"]["sha256"] != digest(raw_input_path):
                raise ArithmeticError("toric input hash mismatch")

            expected_driver_label = f"r17-074d9-singleton-{tag}-p{prime}"
            hodge, primitive = parse_driver_output(
                raw_output_path, expected_driver_label, prime
            )
            if primitive.degree() != 28 or primitive[0] != 1:
                raise ArithmeticError("unexpected primitive Frobenius polynomial")
            ring = primitive.parent()
            T = ring.gen()
            fixed_counts, cycles = component_cycles(model, input_record, prime)
            trivial = (1 - prime * T) ** 2
            for length, count in cycles.items():
                trivial *= (1 - (prime * T) ** length) ** count
            if trivial.degree() != 10:
                raise ArithmeticError("the elliptic trivial lattice changed rank")
            primitive_geometric = cyclotomic_scaled_degree(primitive, prime)
            primitive_arithmetic = multiplicity_of(primitive, 1 - prime * T)
            trivial_geometric = 10
            trivial_arithmetic = multiplicity_of(trivial, 1 - prime * T)
            geometric_bound = primitive_geometric + 6 - trivial_geometric
            if geometric_bound < 0:
                raise ArithmeticError("cohomological upper bound became negative")
            reduction_rows.append(
                {
                    "prime": prime,
                    "good_reduction": True,
                    "toric_nondegeneracy_certified_by_controlled_reduction": True,
                    "primitive_hodge_numbers": hodge,
                    "primitive_frobenius_coefficients_low_to_high": coefficients_text(
                        primitive
                    ),
                    "primitive_factorization": factor_rows(primitive),
                    "toric_complement_dimension": 6,
                    "toric_complement_geometric_algebraicity": 6,
                    "toric_complement_plus_p_multiplicity_upper_bound": 6,
                    "D4_component_fixed_counts_frobenius_powers_1_through_6": [
                        int(fixed_counts[index]) for index in range(1, 7)
                    ],
                    "D4_component_cycle_counts": {
                        str(length): int(count) for length, count in cycles.items()
                    },
                    "trivial_lattice_frobenius_coefficients_low_to_high": coefficients_text(
                        trivial
                    ),
                    "primitive_geometric_algebraic_multiplicity": primitive_geometric,
                    "primitive_plus_p_multiplicity": primitive_arithmetic,
                    "trivial_lattice_geometric_rank": trivial_geometric,
                    "trivial_lattice_plus_p_multiplicity": trivial_arithmetic,
                    "geometric_MW_rank_upper_bound": geometric_bound,
                    "arithmetic_component_cycle_diagnostic_only": True,
                    "raw_input": {
                        "path": relative(raw_input_path),
                        "sha256": digest(raw_input_path),
                    },
                    "raw_output": {
                        "path": relative(raw_output_path),
                        "sha256": digest(raw_output_path),
                    },
                    "input_certificate": {
                        "path": relative(input_certificate_path),
                        "sha256": digest(input_certificate_path),
                    },
                }
            )
            all_raw_inputs.update(
                {
                    relative(raw_input_path): digest(raw_input_path),
                    relative(raw_output_path): digest(raw_output_path),
                    relative(input_certificate_path): digest(input_certificate_path),
                }
            )
        twist_records.append(
            {
                "label": label,
                "reductions": reduction_rows,
                "best_geometric_MW_rank_upper_bound": min(
                    row["geometric_MW_rank_upper_bound"] for row in reduction_rows
                ),
            }
        )

    return {
        "schema": "elkies-k3.r17-074d9-twist-good-reduction-bounds.v1",
        "status": "PASS_EXACT_FOUR_TWIST_GOOD_REDUCTION_RANK_BOUNDS",
        "claim": (
            "Six exact good-reduction Frobenius polynomials are converted into "
            "geometric function-field Mordell--Weil rank upper bounds for each "
            "of the four record-specific singleton twists."
        ),
        "method": {
            "surface_model": "q(u)*y^2=x^3+A(u)*x+B(u)",
            "full_H2": (
                "primitive degree 28 plus a six-dimensional algebraic toric "
                "complement; its exact arithmetic Frobenius action is not assumed"
            ),
            "trivial_lattice": "geometrically U plus two D4 root lattices",
            "geometric_bound": (
                "primitive p-times-root-of-unity multiplicity plus six, minus "
                "the degree-10 geometric trivial lattice"
            ),
            "arithmetic_diagnostic_boundary": (
                "the recorded +p multiplicities and component cycles are not used "
                "as rank bounds because the arithmetic action on the six toric "
                "complement classes has not been identified"
            ),
        },
        "twists": twist_records,
        "proof_boundary": (
            "Cycle classes give unconditional geometric upper bounds; equality "
            "with analytic orders is not assumed. Exact characteristic-zero rank "
            "still requires matching known-section lower bounds and saturation."
        ),
        "generation": {
            "command": (
                "sage -python elkies-k3/scripts/"
                "certify_r17_074d9_twist_good_reduction_bounds.sage"
            ),
            "toric_controlled_reduction_repository": (
                "https://github.com/edgarcosta/ToricControlledReduction"
            ),
            "toric_controlled_reduction_commit": TORIC_CONTROLLED_REDUCTION_COMMIT,
            "sage_version": SAGE_VERSION,
        },
        "inputs": all_raw_inputs,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit("stale 074d9 good-reduction bound certificate")
        terminal = "PASS"
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
        terminal = "WROTE"
    print(
        "R17074D9FROBENIUS|bounds="
        + ",".join(
            f"{row['label']}:{row['best_geometric_MW_rank_upper_bound']}"
            for row in payload["twists"]
        )
        + f"|status={terminal}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
