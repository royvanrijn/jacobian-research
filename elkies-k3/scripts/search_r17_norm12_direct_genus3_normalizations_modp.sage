#!/usr/bin/env sage-python
"""Exhaust a bounded mod-p search for rational genus-3 chord normalizations.

A norm-four trace section on the rootless direct R17 fibration has ``P.O=0``.
Its regular residual-chord slopes are the quadratic polynomials

    M = m0 + m1*u + m2*u^2,

and the chord branch has degree at most eight.  The normalization is rational
in the finite three-node chart when the branch squareclass has degree two,
equivalently

    q(u) = S(u)^2 * Q2(u),   deg(S)=3, deg(Q2)=2.

For every selected exact trace this script exhausts all ``p^3`` affine slope
parameters and tests the polynomial squareclass exactly over ``GF(p)``.  The
result is a complete finite-field discovery sieve for the displayed traces,
not a characteristic-zero existence or nonexistence theorem.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
DEFAULT_MODEL = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-genus3-normalization-modp-v1.json"
)
CHORD_SCRIPT = SCRIPTS / "construct_elkies_2026_bisections.sage"


def load_script(name, path):
    loader = SourceFileLoader(name, str(path))
    spec = spec_from_loader(name, loader)
    if spec is None:
        raise ImportError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    loader.exec_module(module)
    return module


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def relative(path):
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def rational_text(value):
    value = QQ(value)
    return (
        str(value.numerator())
        if value.denominator() == 1
        else f"{value.numerator()}/{value.denominator()}"
    )


def rational_function(record, ring, field):
    numerator = ring(
        [QQ(value) for value in record["numerator_coefficients_low_to_high"]]
    )
    denominator = ring(
        [QQ(value) for value in record["denominator_coefficients_low_to_high"]]
    )
    return field(numerator) / field(denominator)


def mod_p(value, field):
    value = QQ(value)
    if int(value.denominator()) % field.characteristic() == 0:
        raise ZeroDivisionError
    return field(value.numerator()) / field(value.denominator())


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--prime", type=int, required=True)
parser.add_argument("--max-l1", type=int, default=1)
parser.add_argument("--start", type=int, default=0)
parser.add_argument("--trace-limit", type=int)
parser.add_argument(
    "--trace-index",
    type=int,
    action="append",
    help="process this global ordered trace index (repeatable)",
)
parser.add_argument("--pari-stack-gb", type=int, default=2)
parser.add_argument(
    "--modular-trace-arithmetic",
    action="store_true",
    help="form trace sections directly on a minimal mod-p model",
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if not ZZ(args.prime).is_prime():
    parser.error("--prime must be prime")
if args.max_l1 <= 0 or args.start < 0 or args.pari_stack_gb <= 0:
    parser.error("invalid bound")
if args.trace_limit is not None and args.trace_limit <= 0:
    parser.error("--trace-limit must be positive")

chord = load_script("r17_direct_genus3_chord", CHORD_SCRIPT)
model_path = args.model.resolve()
model = json.loads(model_path.read_text())
if model.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    raise ValueError("expected a certified direct norm-12 model")

gram = matrix(ZZ, model["sections"]["height_gram"])
pari.allocatemem(args.pari_stack_gb * 1024**3)
short_vectors = matrix(ZZ, pari(gram).qfminim(4)[2])
best_norm_four = {}
exact_norm_four_count = 0
for column in range(short_vectors.ncols()):
    value = short_vectors.column(column)
    if int(value * gram * value) != 4:
        continue
    exact_norm_four_count += 1
    oriented = min(tuple(value), tuple(-value))
    candidate = vector(ZZ, oriented)
    mask = sum((int(entry) & 1) << index for index, entry in enumerate(candidate))
    score = (
        int(sum(abs(entry) for entry in candidate)),
        sum(bool(entry) for entry in candidate),
        int(max(abs(entry) for entry in candidate)),
        oriented,
    )
    if mask not in best_norm_four or score < best_norm_four[mask][0]:
        best_norm_four[mask] = (score, candidate)
l1_histogram = Counter(item[0][0] for item in best_norm_four.values())
selected = [
    (mask, item[1])
    for mask, item in best_norm_four.items()
    if item[0][0] <= args.max_l1
]
selected.sort(key=lambda item: (sum(abs(entry) for entry in item[1]), tuple(item[1])))
available_in_prefix = len(selected)
indexed_selected = list(enumerate(selected))
if args.trace_index:
    requested_indices = set(args.trace_index)
    indexed_selected = [
        item for item in indexed_selected if item[0] in requested_indices
    ]
    found_indices = {item[0] for item in indexed_selected}
    if found_indices != requested_indices:
        raise ValueError(
            f"trace indices outside selected prefix: {sorted(requested_indices - found_indices)}"
        )
else:
    indexed_selected = indexed_selected[args.start :]
    if args.trace_limit is not None:
        indexed_selected = indexed_selected[: args.trace_limit]

ring = PolynomialRing(QQ, "u")
u = ring.gen()
function_field = ring.fraction_field()
weierstrass = model["weierstrass_model"]
A = ring([QQ(value) for value in weierstrass["A_coefficients_low_to_high"]])
B = ring([QQ(value) for value in weierstrass["B_coefficients_low_to_high"]])
curve = EllipticCurve(function_field, [A, B])
basis = [
    curve(
        rational_function(record["X"], ring, function_field),
        rational_function(record["Y"], ring, function_field),
    )
    for record in model["sections"]["records"]
]

prime = int(args.prime)
field = GF(prime)
modp_ring = PolynomialRing(field, "u")
up = modp_ring.gen()
model_scale_valuation = 0
if args.modular_trace_arithmetic:
    minimum_a_valuation = min(value.valuation(prime) for value in A if value)
    minimum_b_valuation = min(value.valuation(prime) for value in B if value)
    model_scale_valuation = min(
        minimum_a_valuation // 4, minimum_b_valuation // 6
    )
    model_scale = QQ(prime) ** model_scale_valuation
    Ap = modp_ring([mod_p(value / model_scale**4, field) for value in A])
    Bp = modp_ring([mod_p(value / model_scale**6, field) for value in B])
    modp_function_field = modp_ring.fraction_field()
    curve_modp = EllipticCurve(modp_function_field, [Ap, Bp])

    def modular_rational_function(record, weight):
        numerator_values = [
            QQ(value) / model_scale**weight
            for value in record["numerator_coefficients_low_to_high"]
        ]
        denominator_values = [
            QQ(value)
            for value in record["denominator_coefficients_low_to_high"]
        ]
        common_valuation = min(
            value.valuation(prime)
            for value in numerator_values + denominator_values
            if value
        )
        common_scale = QQ(prime) ** common_valuation
        numerator = modp_ring(
            [mod_p(value / common_scale, field) for value in numerator_values]
        )
        denominator = modp_ring(
            [mod_p(value / common_scale, field) for value in denominator_values]
        )
        if not denominator:
            raise ZeroDivisionError("section denominator vanishes modulo p")
        return modp_function_field(numerator) / modp_function_field(denominator)

    basis_modp = [
        curve_modp(
            modular_rational_function(record["X"], 2),
            modular_rational_function(record["Y"], 3),
        )
        for record in model["sections"]["records"]
    ]
    trace_curve, trace_basis, trace_ring = curve_modp, basis_modp, modp_ring
else:
    trace_curve, trace_basis, trace_ring = curve, basis, ring
survivors = []
trace_records = []
profile_histogram = Counter()
total_specializations = 0
for trace_index, (mask, trace_vector) in indexed_selected:
    trace = sum(
        (
            coefficient * point
            for coefficient, point in zip(trace_vector, trace_basis)
            if coefficient
        ),
        trace_curve(0),
    )
    frame = chord.trace_chord_frame(trace[0], trace[1], trace_ring)
    h, Nx, Ny = (frame[key] for key in ("h", "Nx", "Ny"))
    record = {
        "trace_index": trace_index,
        "translation_orbit_mask": int(mask),
        "basis_coordinates": list(map(int, trace_vector)),
        "coefficient_l1": int(sum(abs(entry) for entry in trace_vector)),
        "finite_pole_degree": int(h.degree()),
        "trace_arithmetic": (
            "minimal_model_mod_p" if args.modular_trace_arithmetic else "exact_QQ_then_reduce"
        ),
    }
    if h.degree() != 0:
        record["status"] = "SKIPPED_NONPOLYNOMIAL_TRACE"
        trace_records.append(record)
        continue
    if args.modular_trace_arithmetic:
        Nxp, Nyp = Nx, Ny
    else:
        try:
            Nxp = modp_ring([mod_p(value, field) for value in Nx])
            Nyp = modp_ring([mod_p(value, field) for value in Ny])
            Ap = modp_ring([mod_p(value, field) for value in A])
        except ZeroDivisionError:
            record["status"] = "SKIPPED_BAD_REDUCTION_DENOMINATOR"
            trace_records.append(record)
            continue

    trace_survivor_count = 0
    for m0 in field:
        for m1 in field:
            for m2 in field:
                total_specializations += 1
                M = m0 + m1 * up + m2 * up**2
                q = M**4 - 6 * M**2 * Nxp - 8 * M * Nyp - 3 * Nxp**2 - 4 * Ap
                if not q:
                    profile_histogram["zero_polynomial"] += 1
                    continue
                odd_part = modp_ring(q.squarefree_part())
                profile = f"q{q.degree()}_odd{odd_part.degree()}"
                profile_histogram[profile] += 1
                if odd_part.degree() != 2:
                    continue
                factorization = q.factor()
                square_part = modp_ring.one()
                reduced = modp_ring(factorization.unit())
                for factor, exponent in factorization:
                    square_part *= factor ** (int(exponent) // 2)
                    if int(exponent) % 2:
                        reduced *= factor
                if square_part**2 * reduced != q or reduced.degree() != 2:
                    raise ArithmeticError("modular squareclass decomposition failed")
                trace_survivor_count += 1
                survivors.append(
                    {
                        "trace_index": trace_index,
                        "translation_orbit_mask": int(mask),
                        "basis_coordinates": list(map(int, trace_vector)),
                        "m0_m1_m2": [int(m0), int(m1), int(m2)],
                        "branch_coefficients_low_to_high": [int(value) for value in q],
                        "removed_square_factor_coefficients_low_to_high": [
                            int(value) for value in square_part
                        ],
                        "reduced_quadratic_coefficients_low_to_high": [
                            int(value) for value in reduced
                        ],
                    }
                )
    record["status"] = "PASS_COMPLETE_AFFINE_SLOPE_CENSUS_MODP"
    record["survivor_count"] = trace_survivor_count
    trace_records.append(record)

output_path = args.output if args.output.is_absolute() else ROOT / args.output
output_path.parent.mkdir(parents=True, exist_ok=True)
payload = {
    "schema": "elkies-k3.r17-norm12-direct-genus3-normalization-modp-search.v1",
    "status": (
        "PASS_COMPLETE_MODP_RATIONAL_NORMALIZATION_SURVIVORS"
        if survivors
        else "PASS_COMPLETE_MODP_NO_RATIONAL_NORMALIZATION"
    ),
    "prime": prime,
    "search": {
        "trace_norm": 4,
        "coefficient_l1_bound": args.max_l1,
        "start": args.start,
        "trace_limit": args.trace_limit,
        "exact_norm_four_representatives_up_to_sign": exact_norm_four_count,
        "minimum_norm_four_translation_classes": len(best_norm_four),
        "available_in_l1_prefix": available_in_prefix,
        "selected_trace_count": len(indexed_selected),
        "explicit_trace_indices": args.trace_index,
        "modular_trace_arithmetic": args.modular_trace_arithmetic,
        "minimal_model_scale_prime_valuation": int(model_scale_valuation),
        "processed_trace_count": sum(
            record["status"] == "PASS_COMPLETE_AFFINE_SLOPE_CENSUS_MODP"
            for record in trace_records
        ),
        "total_specialization_count": total_specializations,
        "l1_histogram_best_representative_by_translation_class": {
            str(key): value for key, value in sorted(l1_histogram.items())
        },
        "branch_profile_histogram": dict(sorted(profile_histogram.items())),
    },
    "survivor_count": len(survivors),
    "survivor_trace_count": len({item["trace_index"] for item in survivors}),
    "survivors": survivors,
    "trace_records": trace_records,
    "proof_boundary": (
        "For every processed trace all p^3 affine quadratic chord slopes were "
        "tested exactly. Survivors have branch squareclass degree two modulo p. "
        "This is a finite-field discovery sieve only; exact QQ reconstruction and "
        "cover-section verification are required. Skipped traces are out of scope."
    ),
    "inputs": {
        relative(path): digest(path) for path in (model_path, CHORD_SCRIPT)
    },
    "reproducing_command": (
        "sage -python "
        "elkies-k3/scripts/search_r17_norm12_direct_genus3_normalizations_modp.sage "
        f"--prime {prime} --max-l1 {args.max_l1} --start {args.start}"
        + ("" if args.trace_limit is None else f" --trace-limit {args.trace_limit}")
        + (
            ""
            if not args.trace_index
            else " " + " ".join(f"--trace-index {index}" for index in args.trace_index)
        )
        + (" --modular-trace-arithmetic" if args.modular_trace_arithmetic else "")
    ),
}
output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"R17GENUS3MODP|p={prime}|selected={len(indexed_selected)}"
    f"|processed={payload['search']['processed_trace_count']}"
    f"|specializations={total_specializations}|survivors={len(survivors)}"
    f"|survivor_traces={payload['survivor_trace_count']}"
    f"|output={relative(output_path)}|status={payload['status']}",
    flush=True,
)
