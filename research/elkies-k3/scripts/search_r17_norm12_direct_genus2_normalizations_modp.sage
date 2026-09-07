#!/usr/bin/env sage-python
"""Screen norm-six residual chords and their quadratic covers over GF(p).

For a norm-six trace section on the direct rootless R17 fibration, the
regular residual-chord slopes are

    M = M0 + (l0 + l1*u)*h^2.

The branch has degree at most six.  Its normalization is rational in this
finite two-node chart precisely when its odd squareclass has degree two.
This script exhausts all p^2 affine slopes for a deterministic chunk of the
26,645 minimum norm-six translation classes, using the reciprocal base chart
when the unique pole is at infinity.  The missing projective slope parameter
has branch ``(l0+l1*u)^4*h^2`` and hence is split.  The script also compares
the full finite-field binary-quadratic squareclass, including its scalar atom,
with the complete smooth rational-bisection atlas.

The result is an exact finite-field sieve for the displayed chunk.  It is not
a characteristic-zero existence or nonexistence theorem: a rational slope can
still have bad reduction at a displayed prime.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
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
DEFAULT_SMOOTH = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-alternate-bisections-full-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-genus2-normalization-modp-v1.json"
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


def smooth_branch_records(path):
    """Stream the 364-MiB pretty JSON without materializing lifted sections."""

    pending = None
    records = []
    with path.open() as stream:
        iterator = iter(stream)
        for line in iterator:
            if '"numerator_coefficients": [' in line:
                values = []
                for coefficient_line in iterator:
                    text = coefficient_line.strip().rstrip(",")
                    if text == "]":
                        break
                    values.append(QQ(json.loads(text)))
                if len(values) == 3:
                    pending = values
            elif pending is not None and '"label": ' in line:
                label = json.loads(line.split(":", 1)[1].strip().rstrip(","))
                records.append((label, pending))
                pending = None
    return records


def cover_key(polynomial, field):
    """Binary quadratic atom plus the squareclass of its first nonzero scalar."""

    values = [field(polynomial[index]) for index in range(3)]
    pivot = next((value for value in values if value), None)
    if pivot is None:
        return None
    scalar_sign = int(pivot ** ((field.cardinality() - 1) // 2))
    return tuple(int(value / pivot) for value in values), scalar_sign


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--smooth-covers", type=Path, default=DEFAULT_SMOOTH)
parser.add_argument("--prime", type=int, required=True)
parser.add_argument("--max-l1", type=int, default=42)
parser.add_argument("--start", type=int, default=0)
parser.add_argument("--trace-limit", type=int)
parser.add_argument(
    "--trace-index",
    type=int,
    action="append",
    help="process this global ordered trace index (repeatable)",
)
parser.add_argument("--pari-stack-gb", type=int, default=2)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if not ZZ(args.prime).is_prime():
    parser.error("--prime must be prime")
if args.max_l1 <= 0 or args.start < 0 or args.pari_stack_gb <= 0:
    parser.error("invalid bound")
if args.trace_limit is not None and args.trace_limit <= 0:
    parser.error("--trace-limit must be positive")

chord = load_script("r17_direct_genus2_screen_chord", CHORD_SCRIPT)
model_path = args.model.resolve()
smooth_path = args.smooth_covers.resolve()
model = json.loads(model_path.read_text())
if model.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    raise ValueError("expected a certified direct norm-12 model")

gram = matrix(ZZ, model["sections"]["height_gram"])
pari.allocatemem(args.pari_stack_gb * 1024**3)
short_vectors = matrix(ZZ, pari(gram).qfminim(6)[2])
minimum_norm_by_mask = {}
best_norm_six = {}
exact_norm_six_count = 0
for column in range(short_vectors.ncols()):
    value = short_vectors.column(column)
    norm = int(value * gram * value)
    mask = sum((int(entry) & 1) << index for index, entry in enumerate(value))
    minimum_norm_by_mask[mask] = min(norm, minimum_norm_by_mask.get(mask, norm))
    if norm != 6:
        continue
    exact_norm_six_count += 1
    oriented = min(tuple(value), tuple(-value))
    candidate = vector(ZZ, oriented)
    score = (
        int(sum(abs(entry) for entry in candidate)),
        sum(bool(entry) for entry in candidate),
        int(max(abs(entry) for entry in candidate)),
        oriented,
    )
    if mask not in best_norm_six or score < best_norm_six[mask][0]:
        best_norm_six[mask] = score, candidate
best_norm_six = {
    mask: item
    for mask, item in best_norm_six.items()
    if minimum_norm_by_mask[mask] == 6
}
l1_histogram = Counter(item[0][0] for item in best_norm_six.values())
selected = [
    (mask, item[1])
    for mask, item in best_norm_six.items()
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
qq_function_field = ring.fraction_field()
weierstrass = model["weierstrass_model"]
A = ring([QQ(value) for value in weierstrass["A_coefficients_low_to_high"]])
B = ring([QQ(value) for value in weierstrass["B_coefficients_low_to_high"]])
prime = int(args.prime)
field = GF(prime)
modp_ring = PolynomialRing(field, "u")
up = modp_ring.gen()
modp_function_field = modp_ring.fraction_field()
minimum_a_valuation = min(value.valuation(prime) for value in A if value)
minimum_b_valuation = min(value.valuation(prime) for value in B if value)
model_scale_valuation = min(minimum_a_valuation // 4, minimum_b_valuation // 6)
model_scale = QQ(prime) ** model_scale_valuation
Ap = modp_ring([mod_p(value / model_scale**4, field) for value in A])
Bp = modp_ring([mod_p(value / model_scale**6, field) for value in B])
curve = EllipticCurve(modp_function_field, [Ap, Bp])
Ap_inverted = chord.reciprocal_with_bound(Ap, 8, modp_ring)
Bp_inverted = chord.reciprocal_with_bound(Bp, 12, modp_ring)
curve_inverted = EllipticCurve(modp_function_field, [Ap_inverted, Bp_inverted])


def modularize_function(function, scale_weight):
    function = qq_function_field(function)
    numerator_values = [
        QQ(value) / model_scale**scale_weight for value in function.numerator()
    ]
    denominator_values = [QQ(value) for value in function.denominator()]
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


def section_coordinate(record):
    return rational_function(record, ring, qq_function_field)


basis = [
    curve(
        modularize_function(section_coordinate(record["X"]), 2),
        modularize_function(section_coordinate(record["Y"]), 3),
    )
    for record in model["sections"]["records"]
]
basis_inverted = [
    curve_inverted(
        modularize_function(
            chord.invert_rational(
                section_coordinate(record["X"]), 4, ring, qq_function_field
            ),
            2,
        ),
        modularize_function(
            chord.invert_rational(
                section_coordinate(record["Y"]), 6, ring, qq_function_field
            ),
            3,
        ),
    )
    for record in model["sections"]["records"]
]
smooth_by_key = defaultdict(list)
smooth_records = smooth_branch_records(smooth_path)
for label, coefficients in smooth_records:
    try:
        polynomial = modp_ring([mod_p(value, field) for value in coefficients])
    except ZeroDivisionError:
        continue
    key = cover_key(polynomial, field)
    if key is not None:
        smooth_by_key[key].append(int(label.rsplit("-", 1)[1], 16))

survivors = []
trace_records = []
profile_histogram = Counter()
total_specializations = 0
for trace_index, (mask, trace_vector) in indexed_selected:
    trace = sum(
        (
            coefficient * point
            for coefficient, point in zip(trace_vector, basis)
            if coefficient
        ),
        curve(0),
    )
    frame = chord.trace_chord_frame(trace[0], trace[1], modp_ring)
    h, Nx, Ny, M0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
    construction_chart = "finite_u"
    active_A = Ap
    if h.degree() != 1:
        trace = sum(
            (
                coefficient * point
                for coefficient, point in zip(trace_vector, basis_inverted)
                if coefficient
            ),
            curve_inverted(0),
        )
        frame = chord.trace_chord_frame(trace[0], trace[1], modp_ring)
        h, Nx, Ny, M0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
        construction_chart = "reciprocal_u"
        active_A = Ap_inverted
    record = {
        "trace_index": trace_index,
        "translation_orbit_mask": int(mask),
        "basis_coordinates": list(map(int, trace_vector)),
        "coefficient_l1": int(sum(abs(entry) for entry in trace_vector)),
        "finite_pole_degree_mod_p": int(h.degree()),
        "construction_chart": construction_chart,
    }
    if h.degree() != 1:
        record["status"] = "SKIPPED_DEGENERATE_TRACE_MODP"
        trace_records.append(record)
        continue
    hp, Nxp, Nyp, M0p = h, Nx, Ny, M0

    trace_survivors = 0
    trace_smooth_matches = set()
    for l0 in field:
        for l1 in field:
            total_specializations += 1
            M = M0p + (l0 + l1 * up) * hp**2
            numerator = (
                M**4
                - 6 * M**2 * Nxp
                - 8 * M * Nyp
                - 3 * Nxp**2
                - 4 * active_A * hp**4
            )
            q, remainder = numerator.quo_rem(hp**6)
            if remainder:
                raise ArithmeticError("modular genus-two branch division failed")
            if not q:
                profile_histogram["zero_polynomial"] += 1
                continue
            odd_part = modp_ring(q.squarefree_part())
            profile_histogram[f"q{q.degree()}_odd{odd_part.degree()}"] += 1
            if odd_part.degree() not in (1, 2):
                continue
            factorization = q.factor()
            square_part = modp_ring.one()
            reduced = modp_ring(factorization.unit())
            for factor, exponent in factorization:
                square_part *= factor ** (int(exponent) // 2)
                if int(exponent) % 2:
                    reduced *= factor
            if square_part**2 * reduced != q or reduced.degree() not in (1, 2):
                raise ArithmeticError("modular squareclass decomposition failed")
            reduced_original = (
                reduced
                if construction_chart == "finite_u"
                else chord.reciprocal_with_bound(reduced, 2, modp_ring)
            )
            key = cover_key(reduced_original, field)
            matches = [] if key is None else smooth_by_key.get(key, [])
            trace_survivors += 1
            trace_smooth_matches.update(matches)
            survivors.append(
                {
                    "trace_index": trace_index,
                    "translation_orbit_mask": int(mask),
                    "basis_coordinates": list(map(int, trace_vector)),
                    "l0_l1": [int(l0), int(l1)],
                    "construction_chart": construction_chart,
                    "branch_coefficients_low_to_high": [int(value) for value in q],
                    "removed_square_factor_coefficients_low_to_high": [
                        int(value) for value in square_part
                    ],
                    "reduced_quadratic_coefficients_low_to_high": [
                        int(reduced_original[index]) for index in range(3)
                    ],
                    "smooth_cover_match_masks": matches,
                }
            )
    record["status"] = "PASS_COMPLETE_AFFINE_SLOPE_CENSUS_MODP"
    record["survivor_count"] = trace_survivors
    record["smooth_cover_match_count"] = len(trace_smooth_matches)
    trace_records.append(record)

output_path = args.output if args.output.is_absolute() else ROOT / args.output
output_path.parent.mkdir(parents=True, exist_ok=True)
payload = {
    "schema": "elkies-k3.r17-norm12-direct-genus2-normalization-modp-search.v1",
    "status": (
        "PASS_COMPLETE_MODP_RATIONAL_NORMALIZATION_SURVIVORS"
        if survivors
        else "PASS_COMPLETE_MODP_NO_RATIONAL_NORMALIZATION"
    ),
    "prime": prime,
    "search": {
        "trace_norm": 6,
        "coefficient_l1_bound": args.max_l1,
        "start": args.start,
        "trace_limit": args.trace_limit,
        "explicit_trace_indices": args.trace_index,
        "exact_norm_six_representatives_up_to_sign": exact_norm_six_count,
        "minimum_norm_six_translation_classes": len(best_norm_six),
        "available_in_l1_prefix": available_in_prefix,
        "selected_trace_count": len(indexed_selected),
        "processed_trace_count": sum(
            record["status"] == "PASS_COMPLETE_AFFINE_SLOPE_CENSUS_MODP"
            for record in trace_records
        ),
        "total_specialization_count": total_specializations,
        "l1_histogram_best_representative_by_translation_class": {
            str(key): value for key, value in sorted(l1_histogram.items())
        },
        "branch_profile_histogram": dict(sorted(profile_histogram.items())),
        "smooth_branch_record_count": len(smooth_records),
        "smooth_good_reduction_cover_key_count": len(smooth_by_key),
        "minimal_model_scale_prime_valuation": int(model_scale_valuation),
    },
    "survivor_count": len(survivors),
    "survivor_trace_count": len({item["trace_index"] for item in survivors}),
    "smooth_cover_match_survivor_count": sum(
        bool(item["smooth_cover_match_masks"]) for item in survivors
    ),
    "smooth_cover_match_trace_count": len(
        {item["trace_index"] for item in survivors if item["smooth_cover_match_masks"]}
    ),
    "survivors": survivors,
    "trace_records": trace_records,
    "proof_boundary": (
        "For every processed norm-six trace all p^2 affine linear chord slopes "
        "were tested exactly, in the finite or reciprocal base chart as needed. "
        "The remaining projective slope parameter has split branch "
        "(l0+l1*u)^4*h^2. Survivors have binary-quadratic branch squareclass "
        "modulo p, and smooth-cover matches retain its scalar squareclass. This "
        "is a finite-field sieve only; rational slopes with bad reduction at p "
        "are not excluded."
    ),
    "inputs": {
        relative(path): digest(path) for path in (model_path, smooth_path, CHORD_SCRIPT)
    },
    "reproducing_command": (
        "sage -python "
        "elkies-k3/scripts/search_r17_norm12_direct_genus2_normalizations_modp.sage "
        f"--prime {prime} --max-l1 {args.max_l1} --start {args.start}"
        + ("" if args.trace_limit is None else f" --trace-limit {args.trace_limit}")
        + (
            ""
            if not args.trace_index
            else " " + " ".join(f"--trace-index {index}" for index in args.trace_index)
        )
    ),
}
output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"R17GENUS2MODP|p={prime}|selected={len(indexed_selected)}"
    f"|processed={payload['search']['processed_trace_count']}"
    f"|specializations={total_specializations}|survivors={len(survivors)}"
    f"|survivor_traces={payload['survivor_trace_count']}"
    f"|smooth_matches={payload['smooth_cover_match_survivor_count']}"
    f"|output={relative(output_path)}|status={payload['status']}",
    flush=True,
)
