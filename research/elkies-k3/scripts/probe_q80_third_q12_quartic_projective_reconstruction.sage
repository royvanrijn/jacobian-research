#!/usr/bin/env sage -python
"""Probe a common-projective-scale reconstruction of the third-q12 quartic.

status: ACTIVE_COMPILER
claim: fail-closed reconstruction probe, accepted only after three held-out replays
inputs: exact pencil, exact descent data, p-adic factor lift, and H candidate
outputs: diagnostic JSON only when --output is supplied

The monic quartic has four nonleading coefficients N_i(V)/H(V).  This probe
places the sixteen rational coordinates of the four degree-one N_i in one
projective lattice, in either the omega or denominator-integral delta basis.
It is deliberately not a proof unless an LLL row replays direct exact-pencil
factorizations at all held-out primes.  Even a three-prime pass would still
require characteristic-zero Q^2 division before promotion.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys

from fpylll import IntegerMatrix, LLL
from sage.all import CRT_list, GF, PolynomialRing, QQ, ZZ, inverse_mod


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
DEFAULT_PENCIL = RESULTS / "q80-third-q12-um2-biquadratic-resolved-pencil-qq.json"
DEFAULT_OPERANDS = (
    RESULTS / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json"
)
DEFAULT_LIFT = (
    ROOT
    / "artifacts/local/elkies-k3/"
    / "q80-third-q12-discriminant-factors-p19-adic-precision12288.json"
)
DEFAULT_H = RESULTS / "elkies-k3-q80-third-q12-quartic-denominator-candidate-v1.json"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/local/elkies-k3/"
    / "q80-third-q12-quartic-projective-reconstruction-probe.json"
)
HELD_OUT_PRIMES = (163, 191, 199)
COEFFICIENT = re.compile(r"^(-?\d+)/(\d+)\*theta\^2 ([+-]) (\d+)/(\d+)$")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational(record):
    return QQ(ZZ(record["numerator"])) / ZZ(record["denominator"])


def rational_record(value):
    value = QQ(value)
    return {"numerator": str(value.numerator()), "denominator": str(value.denominator())}


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--pencil", type=Path, default=DEFAULT_PENCIL)
parser.add_argument("--operands", type=Path, default=DEFAULT_OPERANDS)
parser.add_argument("--factor-lift", type=Path, default=DEFAULT_LIFT)
parser.add_argument("--H-candidate", type=Path, default=DEFAULT_H)
parser.add_argument("--basis", choices=("omega", "delta", "both"), default="both")
parser.add_argument("--output", type=Path)
args = parser.parse_args()
for key in ("pencil", "operands", "factor_lift", "H_candidate"):
    setattr(args, key, getattr(args, key).resolve())
if args.output:
    args.output = args.output.resolve()

pencil = json.loads(args.pencil.read_text())
operands = json.loads(args.operands.read_text())
lift = json.loads(args.factor_lift.read_text())
H_artifact = json.loads(args.H_candidate.read_text())
q1 = rational(operands["biquadratic_field"]["q1"])
q2 = rational(operands["biquadratic_field"]["q2"])
product = q1 * q2
product_root = ZZ(product.numerator()).isqrt()
if product_root**2 != product.numerator():
    raise ArithmeticError("q1*q2 numerator is not square")
omega_to_delta = QQ(4 * product_root) / product.denominator()
omega_square = QQ(16) * product
delta_square = ZZ(product.denominator())
modulus = ZZ(lift["specialization"]["modulus"])
omega_to_delta_modulus = (
    ZZ(omega_to_delta.numerator())
    * inverse_mod(ZZ(omega_to_delta.denominator()), modulus)
) % modulus
h0_rational = rational(H_artifact["candidate"]["h0_rational"])
h0_delta = rational(H_artifact["candidate"]["h0_delta"])

q_records = lift["factorization"]["Q"]["coefficients_low_to_high_W"]
if len(q_records) != 5 or q_records[-1]["degrees_numerator_denominator"] != [0, 0]:
    raise ArithmeticError("unexpected monic quartic record")


def numerator_residues(basis):
    values = []
    for record in q_records[:4]:
        numerator = record["numerator_coefficients_low_to_high_U_1_omega"]
        if len(numerator) != 2:
            raise ArithmeticError("quartic numerator lost degree-one representation")
        for pair in numerator:
            values.append(ZZ(pair[0]) % modulus)
            second = ZZ(pair[1]) % modulus
            if basis == "delta":
                second = second * omega_to_delta_modulus % modulus
            values.append(second)
    return values


def projective_rows(residues):
    dimension = len(residues) + 1
    basis = IntegerMatrix(dimension, dimension)
    for index in range(dimension - 1):
        basis[index, index] = int(modulus)
    for index, value in enumerate(residues):
        basis[dimension - 1, index] = int(value)
    basis[dimension - 1, dimension - 1] = 1
    LLL.reduction(basis, delta=0.99)
    rows = []
    for row_index in range(dimension):
        row = [ZZ(basis[row_index, column]) for column in range(dimension)]
        if row[-1] < 0:
            row = [-value for value in row]
        rows.append(row)
    return sorted(rows, key=lambda row: sum(value * value for value in row))


# Parse the exact pencil once; each held-out replay below is a direct reduction
# of these characteristic-zero coefficients, not a replay of the p-adic lift.
terms = []
for v_degree, w_degree, x_degree, encoded in pencil["moving_equation"][
    "terms_T_W_x_coefficient_1_r"
]:
    match = COEFFICIENT.fullmatch(encoded[0])
    if match is None:
        raise ArithmeticError("unexpected exact pencil coefficient encoding")
    theta2 = QQ(ZZ(match[1])) / ZZ(match[2])
    sign = 1 if match[3] == "+" else -1
    constant = sign * QQ(ZZ(match[4])) / ZZ(match[5])
    terms.append(
        (v_degree, w_degree, x_degree, constant + theta2 * (q1 + q2), theta2 / 2)
    )


def held_out_quartics():
    answer = {}
    for prime in HELD_OUT_PRIMES:
        constants = GF(prime)
        local_omega_square = constants(omega_square)
        if local_omega_square.is_square():
            raise ArithmeticError(f"quadratic field splits at held-out prime {prime}")
        z_ring = PolynomialRing(constants, "z")
        z = z_ring.gen()
        finite = GF(prime**2, "omega", modulus=z**2 - local_omega_square)
        omega = finite.gen()
        v_ring = PolynomialRing(finite, "V")
        V = v_ring.gen()
        v_field = v_ring.fraction_field()
        w_ring = PolynomialRing(v_field, "W")
        W = w_ring.gen()
        x_ring = PolynomialRing(w_ring, "old_x")
        coefficients = [w_ring.zero() for _ in range(4)]
        for v_degree, w_degree, x_degree, rational_part, omega_part in terms:
            coefficient = finite(rational_part) + finite(omega_part) * omega
            coefficients[x_degree] += coefficient * v_field(V) ** v_degree * W**w_degree
        cubic = x_ring([value / coefficients[3] for value in coefficients])
        b, c, d = cubic[2], cubic[1], cubic[0]
        discriminant = b**2 * c**2 - 4 * c**3 - 4 * b**3 * d - 27 * d**2 + 18 * b * c * d
        factors = discriminant.factor()
        shape = sorted((int(factor.degree()), int(exponent)) for factor, exponent in factors)
        if shape != [(1, 3), (4, 1), (4, 2)]:
            raise ArithmeticError(f"held-out factor shape changed at p={prime}")
        answer[prime] = next(factor.monic() for factor, exponent in factors if int(exponent) == 2)
    return answer


finite_quartics = held_out_quartics()


def rational_mod(value, prime):
    value = QQ(value)
    if value.denominator() % prime == 0:
        raise ZeroDivisionError
    return int(value.numerator() * inverse_mod(value.denominator(), prime) % prime)


def finite_numerator_coordinates(prime, basis):
    expected = finite_quartics[prime]
    finite = expected.base_ring().base_ring()
    omega = finite.gen()
    V = expected.base_ring().gen()
    local_H = V + finite(h0_rational) + finite(h0_delta / omega_to_delta) * omega
    coordinates = []
    for w_degree in range(4):
        numerator = expected[w_degree] * local_H
        if numerator.denominator().degree() != 0:
            raise ArithmeticError(
                f"held-out Q*H is not polynomial at p={prime}: "
                f"denominator={numerator.denominator()} H={local_H} "
                f"Qden={expected[w_degree].denominator()}"
            )
        coefficients = [
            value / numerator.denominator()[0]
            for value in numerator.numerator().list()
        ]
        coefficients += [finite.zero()] * (2 - len(coefficients))
        if len(coefficients) != 2:
            raise ArithmeticError(f"held-out Q*H changed degree at p={prime}")
        for coefficient in coefficients:
            pair = list(coefficient.polynomial())
            pair += [finite.zero()] * (2 - len(pair))
            if basis == "delta":
                pair[1] *= finite(omega_to_delta)
            coordinates.extend([int(pair[0]), int(pair[1])])
    return coordinates


def asymmetric_scalar_candidates(residue, expected_by_prime):
    held_out_modulus = ZZ(math.prod(expected_by_prime))
    expected_crt = ZZ(
        CRT_list(
            [ZZ(expected_by_prime[prime]) for prime in expected_by_prime],
            [ZZ(prime) for prime in expected_by_prime],
        )
    )
    old_remainder, remainder = modulus, ZZ(residue) % modulus
    old_cofactor, cofactor = ZZ(0), ZZ(1)
    old_remainder_held, remainder_held = int(modulus % held_out_modulus), int(
        residue % held_out_modulus
    )
    old_cofactor_held, cofactor_held = 0, 1
    tested = 0
    accepted = []
    while remainder:
        if cofactor and int(cofactor % 19):
            tested += 1
            if math.gcd(cofactor_held, int(held_out_modulus)) == 1 and (
                remainder_held
                * pow(cofactor_held, -1, int(held_out_modulus))
                - int(expected_crt)
            ) % int(held_out_modulus) == 0:
                accepted.append(QQ(remainder) / QQ(cofactor))
        quotient = old_remainder // remainder
        old_remainder, remainder = remainder, old_remainder - quotient * remainder
        old_cofactor, cofactor = cofactor, old_cofactor - quotient * cofactor
        quotient_held = int(quotient % held_out_modulus)
        old_remainder_held, remainder_held = (
            remainder_held,
            (old_remainder_held - quotient_held * remainder_held)
            % int(held_out_modulus),
        )
        old_cofactor_held, cofactor_held = (
            cofactor_held,
            (old_cofactor_held - quotient_held * cofactor_held)
            % int(held_out_modulus),
        )
    unique = {str(value): value for value in accepted}
    ordered = sorted(
        unique.values(),
        key=lambda value: (
            max(abs(value.numerator()).nbits(), value.denominator().nbits()),
            abs(value.numerator()).nbits() + value.denominator().nbits(),
        ),
    )
    return ordered, tested


def decode_row(row, basis):
    scale = row[-1]
    if not scale or math.gcd(int(scale), int(modulus)) != 1:
        return None
    if any((row[index] - scale * numerator_residues(basis)[index]) % modulus for index in range(16)):
        raise ArithmeticError("LLL row lost its modular congruence")
    values = [QQ(value) / scale for value in row[:-1]]
    numerators = []
    position = 0
    for unused_w in range(4):
        coefficient = []
        for unused_v in range(2):
            coefficient.append((values[position], values[position + 1]))
            position += 2
        numerators.append(coefficient)
    return numerators


def replay_candidate(numerators, basis):
    passed = []
    for prime, expected in finite_quartics.items():
        finite = expected.base_ring().base_ring()
        omega = finite.gen()
        V = expected.base_ring().gen()
        if basis == "omega":
            second_scale = QQ(1)
        else:
            second_scale = QQ(1) / omega_to_delta
        local_H = V + finite(h0_rational) + finite(h0_delta / omega_to_delta) * omega
        observed = expected.parent().gen() ** 4
        W = expected.parent().gen()
        for w_degree, coefficient in enumerate(numerators):
            local_numerator = expected.base_ring().zero()
            for v_degree, (first, second) in enumerate(coefficient):
                local_value = finite(first) + finite(second * second_scale) * omega
                local_numerator += local_value * V**v_degree
            observed += (local_numerator / local_H) * W**w_degree
        if observed != expected:
            return passed
        passed.append(prime)
    return passed


selected_bases = ("omega", "delta") if args.basis == "both" else (args.basis,)
diagnostics = {}
accepted = None
for selected_basis in selected_bases:
    residues = numerator_residues(selected_basis)
    finite_coordinates = {
        prime: finite_numerator_coordinates(prime, selected_basis)
        for prime in HELD_OUT_PRIMES
    }
    scalar_values = []
    scalar_diagnostics = []
    for coordinate_index, residue in enumerate(residues):
        candidates, tested = asymmetric_scalar_candidates(
            residue,
            {
                prime: coordinates[coordinate_index]
                for prime, coordinates in finite_coordinates.items()
            },
        )
        scalar_diagnostics.append(
            {
                "coordinate_index": coordinate_index,
                "euclidean_convergents_tested": tested,
                "held_out_valid_candidate_count": len(candidates),
                "candidate_height_bits": [
                    [abs(value.numerator()).nbits(), value.denominator().nbits()]
                    for value in candidates
                ],
            }
        )
        scalar_values.append(candidates[0] if len(candidates) == 1 else None)
    scalar_numerators = None
    if all(value is not None for value in scalar_values):
        scalar_numerators = []
        position = 0
        for unused_w in range(4):
            coefficient = []
            for unused_v in range(2):
                coefficient.append((scalar_values[position], scalar_values[position + 1]))
                position += 2
            scalar_numerators.append(coefficient)
        if replay_candidate(scalar_numerators, selected_basis) != list(HELD_OUT_PRIMES):
            raise ArithmeticError("coordinatewise held-out candidates fail assembled replay")
        accepted = {
            "basis": selected_basis,
            "row": None,
            "numerators": scalar_numerators,
            "diagnostic": {
                "method": "asymmetric_scalar_convergents",
                "held_out_primes_replayed": list(HELD_OUT_PRIMES),
            },
        }
    coefficient_values = []
    coefficient_diagnostics = []
    for w_degree in range(4):
        block = residues[4 * w_degree : 4 * w_degree + 4]
        block_rows = projective_rows(block)
        passing = []
        block_records = []
        for row_index, row in enumerate(block_rows):
            scale = row[-1]
            replayed = []
            values = None
            if scale and math.gcd(int(scale), int(modulus)) == 1:
                if any(
                    (row[index] - scale * block[index]) % modulus
                    for index in range(4)
                ):
                    raise ArithmeticError("coefficient LLL row lost modular congruence")
                values = [QQ(value) / scale for value in row[:-1]]
                try:
                    for prime in HELD_OUT_PRIMES:
                        expected_block = finite_coordinates[prime][
                            4 * w_degree : 4 * w_degree + 4
                        ]
                        if [rational_mod(value, prime) for value in values] != expected_block:
                            break
                        replayed.append(prime)
                except ZeroDivisionError:
                    replayed = []
            block_records.append(
                {
                    "row_index_by_norm": row_index,
                    "maximum_primitive_coordinate_bits": max(
                        abs(value).nbits() for value in row
                    ),
                    "scale_bits": abs(scale).nbits(),
                    "held_out_primes_replayed": replayed,
                }
            )
            if replayed == list(HELD_OUT_PRIMES):
                passing.append(values)
        coefficient_diagnostics.append(
            {
                "w_degree": w_degree,
                "lattice_dimension": 5,
                "random_lattice_boundary_bits": int(
                    math.ceil(modulus.nbits() * 4 / 5)
                ),
                "passing_row_count": len(passing),
                "rows": block_records,
            }
        )
        coefficient_values.append(passing[0] if len(passing) == 1 else None)
    if accepted is None and all(value is not None for value in coefficient_values):
        coefficient_numerators = [
            [(value[0], value[1]), (value[2], value[3])]
            for value in coefficient_values
        ]
        if replay_candidate(coefficient_numerators, selected_basis) != list(HELD_OUT_PRIMES):
            raise ArithmeticError("coefficient-projective candidates fail assembled replay")
        accepted = {
            "basis": selected_basis,
            "row": None,
            "numerators": coefficient_numerators,
            "diagnostic": {
                "method": "per_coefficient_projective_lll",
                "held_out_primes_replayed": list(HELD_OUT_PRIMES),
            },
        }
    pair_values = []
    pair_diagnostics = []
    for pair_index in range(8):
        block = residues[2 * pair_index : 2 * pair_index + 2]
        block_rows = projective_rows(block)
        passing = []
        block_records = []
        for row_index, row in enumerate(block_rows):
            scale = row[-1]
            replayed = []
            values = None
            if scale and math.gcd(int(scale), int(modulus)) == 1:
                if any(
                    (row[index] - scale * block[index]) % modulus
                    for index in range(2)
                ):
                    raise ArithmeticError("pair LLL row lost modular congruence")
                values = [QQ(value) / scale for value in row[:-1]]
                try:
                    for prime in HELD_OUT_PRIMES:
                        expected_block = finite_coordinates[prime][
                            2 * pair_index : 2 * pair_index + 2
                        ]
                        if [rational_mod(value, prime) for value in values] != expected_block:
                            break
                        replayed.append(prime)
                except ZeroDivisionError:
                    replayed = []
            block_records.append(
                {
                    "row_index_by_norm": row_index,
                    "maximum_primitive_coordinate_bits": max(
                        abs(value).nbits() for value in row
                    ),
                    "scale_bits": abs(scale).nbits(),
                    "held_out_primes_replayed": replayed,
                }
            )
            if replayed == list(HELD_OUT_PRIMES):
                passing.append(values)
        pair_diagnostics.append(
            {
                "pair_index": pair_index,
                "lattice_dimension": 3,
                "random_lattice_boundary_bits": int(
                    math.ceil(modulus.nbits() * 2 / 3)
                ),
                "passing_row_count": len(passing),
                "rows": block_records,
            }
        )
        pair_values.append(passing[0] if len(passing) == 1 else None)
    if accepted is None and all(value is not None for value in pair_values):
        pair_numerators = [
            [pair_values[2 * w_degree], pair_values[2 * w_degree + 1]]
            for w_degree in range(4)
        ]
        if replay_candidate(pair_numerators, selected_basis) != list(HELD_OUT_PRIMES):
            raise ArithmeticError("pair-projective candidates fail assembled replay")
        accepted = {
            "basis": selected_basis,
            "row": None,
            "numerators": pair_numerators,
            "diagnostic": {
                "method": "per_field_coefficient_projective_lll",
                "held_out_primes_replayed": list(HELD_OUT_PRIMES),
            },
        }
    rows = projective_rows(residues)
    row_diagnostics = []
    for row_index, row in enumerate(rows):
        numerators = decode_row(row, selected_basis)
        replayed = [] if numerators is None else replay_candidate(numerators, selected_basis)
        record = {
            "row_index_by_norm": row_index,
            "maximum_primitive_coordinate_bits": max(abs(value).nbits() for value in row),
            "scale_bits": abs(row[-1]).nbits(),
            "held_out_primes_replayed": replayed,
        }
        row_diagnostics.append(record)
        if accepted is None and replayed == list(HELD_OUT_PRIMES):
            accepted = {
                "basis": selected_basis,
                "row": row,
                "numerators": numerators,
                "diagnostic": record,
            }
            break
    diagnostics[selected_basis] = {
        "lattice_dimension": 17,
        "random_lattice_boundary_bits": int(math.ceil(modulus.nbits() * 16 / 17)),
        "asymmetric_scalar_convergents": scalar_diagnostics,
        "per_coefficient_projective_lll": coefficient_diagnostics,
        "per_field_coefficient_projective_lll": pair_diagnostics,
        "rows": row_diagnostics,
    }
    if accepted is not None:
        break

payload = {
    "schema": "elkies-k3-q80-third-q12-quartic-projective-reconstruction-probe-v1",
    "status": (
        "PASS_THREE_HELDOUT_PROJECTIVE_QUARTIC_CANDIDATE_EXACT_DIVISION_OPEN"
        if accepted is not None
        else "FAIL_NO_PROJECTIVE_LLL_ROW_REPLAYS_THREE_HELDOUTS"
    ),
    "modulus_bits": int(modulus.nbits()),
    "bases_tested": list(selected_bases),
    "diagnostics": diagnostics,
    "accepted": (
        {
            "basis": accepted["basis"],
            "diagnostic": accepted["diagnostic"],
            "numerators_low_to_high_W_then_V_1_basis": [
                [[rational_record(first), rational_record(second)] for first, second in coefficient]
                for coefficient in accepted["numerators"]
            ],
        }
        if accepted is not None
        else None
    ),
    "claim_boundary": (
        "A passing row has three untouched-prime replays only; exact characteristic-zero "
        "Q^2 division remains mandatory. A failure rules out only this common rational "
        "projective scale and the two displayed field bases at the current modulus."
    ),
    "inputs": {
        "pencil": {"path": str(args.pencil.relative_to(ROOT)), "sha256": sha256(args.pencil)},
        "operands": {"path": str(args.operands.relative_to(ROOT)), "sha256": sha256(args.operands)},
        "factor_lift": {
            "path": str(args.factor_lift.relative_to(ROOT)),
            "sha256": sha256(args.factor_lift),
        },
        "H_candidate": {
            "path": str(args.H_candidate.relative_to(ROOT)),
            "sha256": sha256(args.H_candidate),
        },
        "worker": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256(Path(__file__).resolve()),
        },
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/probe_q80_third_q12_quartic_projective_reconstruction.sage"
    ),
}
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if args.output:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded)
else:
    print(encoded, end="")
print(f"Q80Q12QPROJECTIVE|status={payload['status']}", file=sys.stderr)
