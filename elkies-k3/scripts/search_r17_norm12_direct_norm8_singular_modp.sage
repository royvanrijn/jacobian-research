#!/usr/bin/env sage-python
"""Certify norm-eight singular-pencil exclusions with modular trace arithmetic.

For a minimum norm-eight parity class, let ``m`` be the number of minimum
representatives up to sign.  Those representatives give ``m`` distinct split
members of the associated genus-one pencil; the regular chord gauge puts one
at infinity, where ``q_infinity=h^2``.  Hence the finite pencil discriminant
has an a-priori square divisor of degree ``2*(m-1)``.

At a good prime this script requires the full discriminant to retain degree
22 and its polynomial squareclass to have the complementary degree
``22-2*(m-1)``.  Equality proves that the known split members exhaust the
even-multiplicity part in characteristic zero.  If the modular squareclass
has no projective root, its characteristic-zero counterpart has no rational
root.  Together these two checks exclude both nodal and even-multiplicity
nonsplit rational normalizations without forming the large exact trace.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
from pathlib import Path
import runpy

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
DEFAULT_MODEL = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
DEFAULT_PRIORITY = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-singular-bisection-search-modp-v1.json"
)
INVERSION_SCRIPT = SCRIPTS / "search_r17_norm12_11952_product_bisection_inversion.sage"
CHORD_SCRIPT = SCRIPTS / "construct_elkies_2026_bisections.sage"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def relative(path):
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def parity_mask(word):
    return sum((int(entry) & 1) << index for index, entry in enumerate(word))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--priority-table", type=Path, default=DEFAULT_PRIORITY)
parser.add_argument("--primes", default="131,137,151,157,167,173,181,191,193")
parser.add_argument("--max-l1", type=int, default=44)
parser.add_argument("--start", type=int, default=0)
parser.add_argument("--trace-limit", type=int)
parser.add_argument("--pari-stack-gb", type=int, default=2)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
primes = tuple(int(value) for value in args.primes.split(",") if value)
if not primes or any(not ZZ(prime).is_prime() or prime <= 22 for prime in primes):
    parser.error("--primes must contain primes greater than 22")
if args.max_l1 <= 0 or args.start < 0 or args.pari_stack_gb <= 0:
    parser.error("invalid bound")
if args.trace_limit is not None and args.trace_limit <= 0:
    parser.error("--trace-limit must be positive")

model_path = args.model.resolve()
priority_path = args.priority_table.resolve()
model = json.loads(model_path.read_text())
if model.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    raise ValueError("expected a certified direct norm-12 model")

gram = matrix(ZZ, model["sections"]["height_gram"])
pari.allocatemem(args.pari_stack_gb * 1024**3)
short_vectors = matrix(ZZ, pari(gram).qfminim(8)[2])
minimum_norm_by_mask = {}
best_norm_eight_by_mask = {}
for column in range(short_vectors.ncols()):
    value = short_vectors.column(column)
    norm = int(value * gram * value)
    mask = parity_mask(value)
    minimum_norm_by_mask[mask] = min(norm, minimum_norm_by_mask.get(mask, norm))
    if norm != 8:
        continue
    oriented = min(tuple(value), tuple(-value))
    candidate = vector(ZZ, oriented)
    score = (
        int(sum(abs(entry) for entry in candidate)),
        sum(bool(entry) for entry in candidate),
        int(max(abs(entry) for entry in candidate)),
        oriented,
    )
    if mask not in best_norm_eight_by_mask or score < best_norm_eight_by_mask[mask][0]:
        best_norm_eight_by_mask[mask] = (score, candidate)
selected = [
    item[1]
    for mask, item in best_norm_eight_by_mask.items()
    if minimum_norm_by_mask[mask] == 8 and item[0][0] <= args.max_l1
]
selected.sort(key=lambda value: (sum(abs(entry) for entry in value), tuple(value)))
available_in_prefix = len(selected)
selected = selected[args.start :]
if args.trace_limit is not None:
    selected = selected[: args.trace_limit]
if not selected:
    raise ValueError("selected interval is empty")

minimum_count_by_mask = {}
with priority_path.open(newline="") as stream:
    for row in csv.DictReader(stream, delimiter="\t"):
        word = tuple(map(int, row["section_basis_w"].split()))
        mask = parity_mask(word)
        if mask in minimum_count_by_mask:
            raise ValueError("duplicate parity mask in priority table")
        minimum_count_by_mask[mask] = int(row["minimal_unoriented_count"])
expected_class_count = len(best_norm_eight_by_mask)
if len(minimum_count_by_mask) != expected_class_count:
    raise ValueError(
        "priority-table class count does not match the exact norm-eight "
        f"enumeration: {len(minimum_count_by_mask)} != {expected_class_count}"
    )
if set(minimum_count_by_mask) != set(best_norm_eight_by_mask):
    raise ValueError("priority-table parity classes do not match the exact enumeration")

inverse = runpy.run_path(str(INVERSION_SCRIPT))
helper = runpy.run_path(str(CHORD_SCRIPT))
maximum_coefficient = max(abs(int(entry)) for word in selected for entry in word)
contexts = {}
context_errors = {}
for prime in primes:
    try:
        contexts[prime] = inverse["modular_context"](
            prime, model, helper, maximum_coefficient
        )
    except (ArithmeticError, ValueError, ZeroDivisionError) as error:
        context_errors[str(prime)] = type(error).__name__
if not contexts:
    raise ValueError("none of the displayed primes gives a good modular model")


def modular_discriminant_record(context, word, expected_odd_degree):
    trace = inverse["trace_from_word"](context, word, context["multiples"])
    if trace.is_zero():
        raise ArithmeticError("trace reduced to zero")
    pencil = inverse["trace_pencil_family"](context, trace)
    coefficient_field = context["coefficient_field"]
    parameter_ring = PolynomialRing(coefficient_field, "lambda")
    parameter = parameter_ring.gen()
    bivariate_ring = PolynomialRing(parameter_ring, "u")
    q_symbolic = sum(
        (
            bivariate_ring(pencil["q_lambda_family"][power]) * parameter**power
            for power in range(5)
        ),
        bivariate_ring.zero(),
    )
    discriminant = parameter_ring(q_symbolic.discriminant())
    odd_part = parameter_ring(discriminant.squarefree_part())
    root_count = sum(odd_part(value) == 0 for value in coefficient_field)
    good = (
        discriminant.degree() == 22
        and odd_part.degree() == expected_odd_degree
        and root_count == 0
    )
    return {
        "prime": int(context["prime"]),
        "chart": pencil["chart"],
        "pencil_discriminant_degree": int(discriminant.degree()),
        "odd_multiplicity_pencil_discriminant_degree": int(odd_part.degree()),
        "odd_part_projective_root_count": int(root_count),
        "full_degree_and_split_exhaustion_verified": bool(
            discriminant.degree() == 22
            and odd_part.degree() == expected_odd_degree
        ),
        "excludes_rational_nonsplit_singular_parameter": bool(good),
    }


trace_records = []
excluding_prime_histogram = Counter()
bad_reduction_histogram = Counter()
unresolved = []
for offset, word in enumerate(selected):
    mask = parity_mask(word)
    minimum_count = minimum_count_by_mask[mask]
    expected_odd_degree = 22 - 2 * (minimum_count - 1)
    reductions = []
    excluding_prime = None
    for prime, context in contexts.items():
        try:
            reduction = modular_discriminant_record(context, word, expected_odd_degree)
        except (ArithmeticError, ValueError, ZeroDivisionError) as error:
            bad_reduction_histogram[f"p{prime}:{type(error).__name__}"] += 1
            reductions.append(
                {
                    "prime": int(prime),
                    "status": "BAD_REDUCTION",
                    "error_type": type(error).__name__,
                }
            )
            continue
        reductions.append(reduction)
        if reduction["excludes_rational_nonsplit_singular_parameter"]:
            excluding_prime = prime
            excluding_prime_histogram[str(prime)] += 1
            break
    record = {
        "trace_index": args.start + offset,
        "basis_coordinates": list(map(int, word)),
        "coefficient_l1": int(sum(abs(entry) for entry in word)),
        "finite_pole_degree": 2,
        "minimum_unoriented_split_member_count": minimum_count,
        "pencil_discriminant_degree": 22,
        "odd_multiplicity_pencil_discriminant_degree": expected_odd_degree,
        "modular_reductions": reductions,
        "rational_singular_parameters": [],
    }
    if excluding_prime is None:
        record["status"] = "UNRESOLVED_MODULAR_SINGULAR_PARAMETER"
        unresolved.append(record["trace_index"])
    else:
        record["status"] = "PASS_MODULAR_FULL_DISCRIMINANT_SPLIT_EXHAUSTION"
        record["excluding_prime"] = int(excluding_prime)
    trace_records.append(record)

status = (
    "PASS_BOUNDED_NO_NONSPLIT_RATIONAL_SINGULAR_MEMBER"
    if not unresolved
    else "UNRESOLVED_MODULAR_SINGULAR_PARAMETERS"
)
output = args.output if args.output.is_absolute() else ROOT / args.output
payload = {
    "schema": "elkies-k3.r17-norm12-direct-singular-bisection-search.v1",
    "status": status,
    "search": {
        "norm": 8,
        "coefficient_l1_bound": args.max_l1,
        "trace_limit": args.trace_limit,
        "start": args.start,
        "minimum_norm_eight_translation_classes": expected_class_count,
        "available_in_l1_prefix": available_in_prefix,
        "processed_trace_count": len(selected),
        "processed_half_open_range": [args.start, args.start + len(selected)],
        "finite_pole_trace_count": len(selected),
        "modular_trace_arithmetic_primes": list(primes),
        "excluding_prime_histogram": dict(sorted(excluding_prime_histogram.items())),
        "bad_reduction_histogram": dict(sorted(bad_reduction_histogram.items())),
        "discarded_model_primes": context_errors,
        "unresolved_trace_indices": unresolved,
    },
    "candidate_count": 0,
    "smooth_atlas_match_count": 0,
    "candidate_collision_count": 0,
    "candidate_collisions": {},
    "candidates": [],
    "trace_records": trace_records,
    "proof_boundary": (
        "For each resolved trace, one displayed good prime retains full finite "
        "discriminant degree 22 and squareclass degree 22-2*(m-1), where m is "
        "the exact minimum-representative count. The m known split members, one "
        "at infinity, therefore exhaust the even-multiplicity part over QQ. The "
        "remaining modular odd part has no projective root, excluding a rational "
        "nonsplit singular parameter."
    ),
    "inputs": {
        relative(path): digest(path)
        for path in (model_path, priority_path, INVERSION_SCRIPT, CHORD_SCRIPT)
    },
    "reproducing_command": (
        "sage -python elkies-k3/scripts/search_r17_norm12_direct_norm8_singular_modp.sage "
        f"--model {relative(model_path)} --priority-table {relative(priority_path)} "
        f"--primes {args.primes} --max-l1 {args.max_l1} --start {args.start}"
        + ("" if args.trace_limit is None else f" --trace-limit {args.trace_limit}")
        + f" --output {relative(output)}"
    ),
}
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "R17NORM8SINGULARMODP"
    f"|traces={len(selected)}|unresolved={len(unresolved)}"
    f"|output={relative(output)}|status={status}",
    flush=True,
)
