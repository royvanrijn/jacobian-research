#!/usr/bin/env sage
"""Reconstruct the aligned alternate-Q80 rootless Weierstrass model over Q.

Inputs are prime-by-prime endpoint artifacts with a common canonical marking
identifier and fixed short-Weierstrass gauge.  The final ``--heldout-count``
primes are never used for CRT or LLL.  The script reconstructs the coefficient
vectors, replays them at every training and held-out prime, classifies the
exact rational fibres, and emits the reduced j-map.

This is deliberately an endpoint consumer.  It does not attempt to align
prime-specific signs, base PGL2 gauges, or Weierstrass scalings: producers
must do that using retained parent/child maps and record the same nonempty
``canonical_marking_id`` and ``route_map_chain_id``.  Every modular record
must also retain a digest for the same ordered set of parent/child map slots.
Refusing unaligned or unattested residues prevents a plausible but meaningless
coefficientwise CRT.
"""

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

from sage.all import CRT_list, PolynomialRing, QQ, ZZ, gcd, lcm, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
DEFAULT_OUTPUT = GENERATED / "q80-alternate-rootless-mw17-qq.json"
INPUT_SCHEMA = "elkies-k3.q80-alternate-rootless-aligned-modp.v1"
INPUT_STATUS = "PASS_ALIGNED_Q80_ALTERNATE_ROOTLESS_MODP"
ROUTE_ID = "Q80-q4-q4-q12-q12-alternate-q4-q6"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("inputs", nargs="*", type=Path)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--heldout-count", type=int, default=2)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
inputs = args.inputs or sorted(GENERATED.glob("q80-alternate-rootless-aligned-modp-*.json"))
inputs = [path.resolve() for path in inputs]
output_path = args.output.resolve()
if args.heldout_count < 1:
    raise ValueError("at least one prime must be held out")
if len(inputs) < args.heldout_count + 2:
    raise ValueError("need at least two training primes plus the held-out ensemble")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def path_label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def pad(values, size):
    if len(values) > size:
        raise ValueError("modular coefficient vector exceeds declared K3 degree")
    return [int(value) for value in values] + [0] * (size - len(values))


def reduce_rational(value, prime):
    value = QQ(value)
    prime = ZZ(prime)
    denominator = ZZ(value.denominator()) % prime
    if denominator == 0:
        return None
    return int((ZZ(value.numerator()) % prime) * denominator.inverse_mod(prime) % prime)


def primitive_integer_pair(numerator, denominator):
    coefficients = list(numerator) + list(denominator)
    common_denominator = ZZ.one()
    for coefficient in coefficients:
        common_denominator = lcm(common_denominator, QQ(coefficient).denominator())
    integers = [
        ZZ(QQ(coefficient).numerator())
        * (common_denominator // QQ(coefficient).denominator())
        for coefficient in coefficients
    ]
    content = ZZ.zero()
    for coefficient in integers:
        content = gcd(content, abs(coefficient))
    if not content:
        raise ArithmeticError("zero rational-function pair")
    integers = [coefficient // content for coefficient in integers]
    split = len(numerator)
    numerator_integers = integers[:split]
    denominator_integers = integers[split:]
    if denominator_integers[-1] < 0:
        numerator_integers = [-value for value in numerator_integers]
        denominator_integers = [-value for value in denominator_integers]
    return numerator_integers, denominator_integers


def kodaira_data(ord_a, ord_b, ord_delta):
    """Return root rank/count/determinant, Euler number, and Kodaira name."""
    if ord_delta <= 0:
        return 0, 0, 1, 0, "smooth"
    if ord_a == 0 or ord_b == 0:
        if ord_delta == 1:
            return 0, 0, 1, 1, "I1"
        n = int(ord_delta)
        return n - 1, n * (n - 1), n, n, f"I{n}"
    if ord_delta == 2:
        return 0, 0, 1, 2, "II"
    if ord_delta == 3:
        return 1, 2, 2, 3, "III"
    if ord_delta == 4:
        return 2, 6, 3, 4, "IV"
    if ord_delta == 6 and ord_a >= 2 and ord_b >= 3:
        return 4, 24, 4, 6, "I0*"
    if ord_delta >= 7 and ord_a == 2 and ord_b == 3:
        n = int(ord_delta - 6)
        rank = n + 4
        return rank, 2 * rank * (rank - 1), 4, n + 6, f"I{n}*"
    if ord_delta == 8:
        return 6, 72, 3, 8, "IV*"
    if ord_delta == 9:
        return 7, 126, 2, 9, "III*"
    if ord_delta == 10:
        return 8, 240, 1, 10, "II*"
    raise ArithmeticError(f"unclassified minimal orders {(ord_a, ord_b, ord_delta)}")


def classify_exact_model(A, B, delta):
    root_rank = root_count = euler = 0
    root_determinant = 1
    finite = []
    for factor, exponent in delta.factor():
        ord_a = int(A.valuation(factor))
        ord_b = int(B.valuation(factor))
        ord_delta = int(exponent)
        scaling = min(ord_a // 4, ord_b // 6, ord_delta // 12)
        orders = (
            ord_a - 4 * scaling,
            ord_b - 6 * scaling,
            ord_delta - 12 * scaling,
        )
        rank, count, determinant, local_euler, name = kodaira_data(*orders)
        degree = int(factor.degree())
        root_rank += degree * rank
        root_count += degree * count
        root_determinant *= determinant**degree
        euler += degree * local_euler
        finite.append(
            {
                "factor": str(factor.monic()),
                "degree": degree,
                "multiplicity": int(exponent),
                "minimal_orders_A_B_Delta": list(orders),
                "kodaira": name,
                "geometric_places": degree,
            }
        )

    infinity_raw = (8 - A.degree(), 12 - B.degree(), 24 - delta.degree())
    infinity_scaling = min(
        infinity_raw[0] // 4, infinity_raw[1] // 6, infinity_raw[2] // 12
    )
    infinity_orders = (
        infinity_raw[0] - 4 * infinity_scaling,
        infinity_raw[1] - 6 * infinity_scaling,
        infinity_raw[2] - 12 * infinity_scaling,
    )
    rank, count, determinant, local_euler, infinity_name = kodaira_data(
        *infinity_orders
    )
    root_rank += rank
    root_count += count
    root_determinant *= determinant
    euler += local_euler
    return {
        "finite": finite,
        "infinity": {
            "minimal_orders_A_B_Delta": list(map(int, infinity_orders)),
            "kodaira": infinity_name,
        },
        "root_rank": int(root_rank),
        "root_count": int(root_count),
        "root_determinant": int(root_determinant),
        "euler_number": int(euler),
        "rootless": root_rank == 0,
        "MW_rank_at_Picard_19": int(17 - root_rank),
    }


started = time.monotonic()
raw_inputs = [path.read_bytes() for path in inputs]
records = [json.loads(raw) for raw in raw_inputs]
for path, record in zip(inputs, records):
    if record.get("schema") != INPUT_SCHEMA or record.get("status") != INPUT_STATUS:
        raise ValueError(f"unaccepted modular endpoint artifact: {path}")
    if record.get("route_id") != ROUTE_ID:
        raise ValueError(f"wrong Q80 route in {path}")
    marking_id = record.get("canonical_marking_id")
    if not isinstance(marking_id, str) or not marking_id:
        raise ValueError(f"missing canonical marking identifier in {path}")
    transport = record.get("marking_transport")
    if not isinstance(transport, dict):
        raise ValueError(f"missing marking transport certificate in {path}")
    chain_id = transport.get("route_map_chain_id")
    if not isinstance(chain_id, str) or not chain_id:
        raise ValueError(f"missing route-map chain identifier in {path}")
    map_digests = transport.get("parent_child_map_digests_modp")
    if (
        not isinstance(map_digests, dict)
        or not map_digests
        or any(not isinstance(key, str) or not key for key in map_digests)
        or any(not isinstance(value, str) or not value for value in map_digests.values())
    ):
        raise ValueError(f"missing parent/child map digests in {path}")

primes = [ZZ(record["prime"]) for record in records]
if len(set(primes)) != len(primes) or not all(prime.is_prime() for prime in primes):
    raise ValueError("input primes must be distinct primes")
marking_ids = {record["canonical_marking_id"] for record in records}
if len(marking_ids) != 1:
    raise ValueError("cross-prime canonical marking identifiers disagree")
route_map_chain_ids = {
    record["marking_transport"]["route_map_chain_id"] for record in records
}
if len(route_map_chain_ids) != 1:
    raise ValueError("cross-prime route-map chain identifiers disagree")
map_slot_sets = {
    tuple(sorted(record["marking_transport"]["parent_child_map_digests_modp"]))
    for record in records
}
if len(map_slot_sets) != 1:
    raise ValueError("cross-prime parent/child map slots disagree")
normalizations = {
    json.dumps(record.get("normalization"), sort_keys=True, separators=(",", ":"))
    for record in records
}
if len(normalizations) != 1:
    raise ValueError("cross-prime base/Weierstrass normalizations disagree")
if any(record["model"].get("A_degree") != 8 for record in records):
    raise ValueError("aligned endpoint must retain degree-eight A")
if any(record["model"].get("B_degree") != 12 for record in records):
    raise ValueError("aligned endpoint must retain degree-twelve B")
if any(record["model"].get("root_rank") != 0 for record in records):
    raise ValueError("a modular endpoint is not rootless")

residue_vectors = []
for prime, record in zip(primes, records):
    A_values = pad(record["model"]["A_coefficients_low_to_high"], 9)
    B_values = pad(record["model"]["B_coefficients_low_to_high"], 13)
    residue_vectors.append([value % int(prime) for value in A_values + B_values])

training_count = len(records) - args.heldout_count
training_records = records[:training_count]
training_primes = primes[:training_count]
training_vectors = residue_vectors[:training_count]
heldout_records = records[training_count:]
heldout_primes = primes[training_count:]
heldout_vectors = residue_vectors[training_count:]
modulus = math.prod(training_primes)


def reconstruct_independently():
    values = []
    failures = []
    for index in range(22):
        residue = ZZ(
            CRT_list([row[index] for row in training_vectors], training_primes)
        )
        try:
            value = QQ(residue.rational_reconstruction(modulus))
        except (ArithmeticError, ValueError, ZeroDivisionError) as error:
            failures.append({"index": index, "reason": str(error)})
            value = None
        values.append(value)
    return values, failures


def simultaneous_projective_candidates():
    residues = [
        ZZ(CRT_list([row[index] for row in training_vectors], training_primes))
        for index in range(22)
    ]
    dimension = 23
    basis = matrix(ZZ, dimension, dimension)
    for index in range(22):
        basis[index, index] = modulus
    basis[-1] = vector(ZZ, residues + [1])
    reduced = basis.LLL(delta=0.99)
    candidates = []
    for row in sorted(reduced.rows(), key=lambda value: value * value):
        scale = ZZ(row[-1])
        if not scale or gcd(scale, modulus) != 1:
            continue
        if not all((row[index] - scale * residues[index]) % modulus == 0 for index in range(22)):
            continue
        candidates.append(
            {
                "values": tuple(QQ(row[index]) / scale for index in range(22)),
                "scale": scale,
                "max_vector_bits": max(abs(ZZ(value)).nbits() for value in row),
                "norm_squared": str(row * row),
            }
        )
    return candidates


def replay(values, selected_primes, selected_vectors):
    failures = []
    for prime, expected in zip(selected_primes, selected_vectors):
        actual = [reduce_rational(value, prime) for value in values]
        for index, (left, right) in enumerate(zip(actual, expected)):
            if left is None or left != right:
                failures.append({"prime": int(prime), "index": index})
    return failures


independent_values, independent_failures = reconstruct_independently()
candidates = []
if not independent_failures:
    candidates.append(
        {
            "method": "coefficientwise_Sage_rational_reconstruction",
            "values": tuple(independent_values),
            "scale": ZZ.one(),
            "max_vector_bits": max(
                max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())
                for value in independent_values
            ),
            "norm_squared": None,
        }
    )
for candidate in simultaneous_projective_candidates():
    candidate["method"] = "simultaneous_projective_LLL"
    candidates.append(candidate)

deduplicated = {}
for candidate in candidates:
    deduplicated.setdefault(tuple(candidate["values"]), candidate)
candidates = list(deduplicated.values())

accepted = []
for candidate in candidates:
    values = candidate["values"]
    training_failures = replay(values, training_primes, training_vectors)
    heldout_failures = replay(values, heldout_primes, heldout_vectors)
    if not training_failures and not heldout_failures:
        candidate["training_failures"] = training_failures
        candidate["heldout_failures"] = heldout_failures
        accepted.append(candidate)
if len(accepted) != 1:
    raise ArithmeticError(
        f"expected one cross-prime reconstruction, found {len(accepted)}; "
        f"independent failures={len(independent_failures)}, LLL candidates={len(candidates)}"
    )
selected = accepted[0]
values = selected["values"]

ring = PolynomialRing(QQ, "u")
u = ring.gen()
A = ring(values[:9])
B = ring(values[9:])
if A.degree() != 8 or B.degree() != 12:
    raise ArithmeticError("reconstructed model lost its canonical leading degrees")
delta = ring(-16 * (4 * A**3 + 27 * B**2))
if not delta:
    raise ArithmeticError("reconstructed model is identically singular")
classification = classify_exact_model(A, B, delta)
if not classification["rootless"] or classification["euler_number"] != 24:
    raise ArithmeticError(
        "reconstructed endpoint is not a rootless elliptic K3: "
        + json.dumps(classification, sort_keys=True)
    )

c4 = ring(-48 * A)
j_numerator_raw = ring(c4**3)
j_denominator_raw = delta
common = j_numerator_raw.gcd(j_denominator_raw)
j_numerator = j_numerator_raw // common
j_denominator = j_denominator_raw // common
j_num_ZZ, j_den_ZZ = primitive_integer_pair(j_numerator.list(), j_denominator.list())

payload = {
    "schema": "q80-alternate-rootless-mw17-qq-v1",
    "status": "PASS_Q80_ALTERNATE_ROOTLESS_MW17_CRT_HELDOUT",
    "route_id": ROUTE_ID,
    "canonical_marking_id": next(iter(marking_ids)),
    "marking_transport": {
        "route_map_chain_id": next(iter(route_map_chain_ids)),
        "parent_child_map_slots": list(next(iter(map_slot_sets))),
        "prime_specific_map_digests_replayed": {
            str(prime): record["marking_transport"]["parent_child_map_digests_modp"]
            for prime, record in zip(primes, records)
        },
    },
    "normalization": records[0]["normalization"],
    "reconstruction": {
        "method": selected["method"],
        "training_primes": list(map(int, training_primes)),
        "heldout_primes": list(map(int, heldout_primes)),
        "training_modulus": str(modulus),
        "training_modulus_bits": int(ZZ(modulus).nbits()),
        "coefficient_count": 22,
        "independent_reconstruction_failures": independent_failures,
        "candidate_count_before_replay": len(candidates),
        "unique_candidate_after_training_and_heldout_replay": True,
        "selected_projective_scale": str(selected["scale"]),
        "selected_vector_max_bits": int(selected["max_vector_bits"]),
        "all_training_coefficients_replayed": True,
        "all_heldout_coefficients_replayed": True,
    },
    "rootless_family": {
        "equation": "y^2=x^3+A(u)x+B(u)",
        "A_coefficients_low_to_high": list(map(str, A.list())),
        "B_coefficients_low_to_high": list(map(str, B.list())),
        "A_degree": int(A.degree()),
        "B_degree": int(B.degree()),
        "Delta_coefficients_low_to_high": list(map(str, delta.list())),
        "Delta_degree": int(delta.degree()),
        "classification": classification,
        "generic_MW_rank_at_Picard_19": 17,
        "reduced_j": {
            "formula": "(-48*A(u))^3/(-16*(4*A(u)^3+27*B(u)^2))",
            "cancelled_gcd_degree": int(common.degree()),
            "primitive_integer_numerator_coefficients_low_to_high": list(map(str, j_num_ZZ)),
            "primitive_integer_denominator_coefficients_low_to_high": list(map(str, j_den_ZZ)),
            "numerator_degree": len(j_num_ZZ) - 1,
            "denominator_degree": len(j_den_ZZ) - 1,
        },
    },
    "inputs": {
        path_label(path): {
            "sha256": hashlib.sha256(raw).hexdigest(),
            "prime": int(prime),
            "role": "training" if index < training_count else "heldout",
        }
        for index, (path, raw, prime) in enumerate(zip(inputs, raw_inputs, primes))
    },
    "proof_boundary": (
        "The displayed QQ Weierstrass equation, fibre classification, rootlessness, "
        "Euler number, and j-map are exact. Its identification with the selected "
        "generic Q80 neighbour chain uses the common canonical marking and exact "
        "replay at every supplied training and held-out prime; the modular producers "
        "must separately retain and certify the parent/child rational maps."
    ),
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/reconstruct_q80_alternate_rootless_crt_qq.sage "
        + " ".join(path_label(path) for path in inputs)
        + f" --heldout-count {args.heldout_count}"
    ),
    "runtime_seconds": time.monotonic() - started,
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if args.check:
    if not output_path.exists():
        raise SystemExit(f"missing reconstructed endpoint artifact: {output_path}")
    existing = json.loads(output_path.read_text())
    existing.pop("runtime_seconds", None)
    payload.pop("runtime_seconds", None)
    if existing != payload:
        raise SystemExit("stale Q80 rootless CRT artifact")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "Q80ROOTLESSCRT|"
    f"training={len(training_primes)}|heldout={len(heldout_primes)}|"
    f"modulus_bits={ZZ(modulus).nbits()}|method={selected['method']}|"
    f"root_rank={classification['root_rank']}|MW=17|"
    f"j_degrees={len(j_num_ZZ)-1},{len(j_den_ZZ)-1}|"
    "status=PASS_Q80_ALTERNATE_ROOTLESS_MW17_CRT_HELDOUT",
    flush=True,
)
