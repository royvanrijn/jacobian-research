#!/usr/bin/env sage-python
"""CRT-reconstruct and exactly verify norm-four three-node normalizations."""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
from itertools import product
import json
from pathlib import Path
import sys

from sage.all import CRT_list, EllipticCurve, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
DEFAULT_MODEL = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
DEFAULT_SMOOTH = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisection-collisions-full-v1.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-genus3-normalization-reconstruction-v1.json"
CHORD_SCRIPT = SCRIPTS / "construct_elkies_2026_bisections.sage"
HASH_SCRIPT = SCRIPTS / "hash_bisection_extensions.py"


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
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def coefficients(polynomial):
    return [rational_text(polynomial[index]) for index in range(polynomial.degree() + 1)]


def rational_function(record, ring, field):
    numerator = ring([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = ring([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    return field(numerator) / field(denominator)


def squareclass_decomposition(polynomial, ring):
    factorization = polynomial.factor()
    square_part = ring.one()
    reduced = ring(factorization.unit())
    for factor, exponent in factorization:
        square_part *= factor ** (int(exponent) // 2)
        if int(exponent) % 2:
            reduced *= factor
    if square_part**2 * reduced != polynomial:
        raise ArithmeticError("squareclass decomposition failed")
    return square_part, reduced


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--screen", type=Path, action="append", required=True)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--smooth-collisions", type=Path, default=DEFAULT_SMOOTH)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if len(args.screen) < 2:
    parser.error("at least two --screen inputs are required")

chord = load_script("r17_direct_genus3_reconstruct_chord", CHORD_SCRIPT)
hasher = load_script("r17_direct_genus3_reconstruct_hasher", HASH_SCRIPT)
screen_paths = [path.resolve() for path in args.screen]
screens = [json.loads(path.read_text()) for path in screen_paths]
if any(screen.get("schema") != "elkies-k3.r17-norm12-direct-genus3-normalization-modp-search.v1" for screen in screens):
    raise ValueError("unexpected modular screen schema")
primes = [ZZ(screen["prime"]) for screen in screens]
if len(set(primes)) != len(primes):
    raise ValueError("screen primes must be distinct")

model_path = args.model.resolve()
model = json.loads(model_path.read_text())
if model.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    raise ValueError("expected a certified direct model")
smooth_path = args.smooth_collisions.resolve()
smooth = json.loads(smooth_path.read_text())
if smooth.get("status") != "PASS_EXTENSION_CANONICALIZATION" or not smooth.get("compact_output"):
    raise ValueError("expected compact smooth collision artifact")
smooth_by_digest = {item["extension_sha256"]: item for item in smooth["extension_manifest"]}

by_prime_and_trace = []
for screen in screens:
    current = {}
    for survivor in screen["survivors"]:
        current.setdefault(int(survivor["trace_index"]), []).append(survivor)
    by_prime_and_trace.append(current)
common_trace_indices = sorted(set.intersection(*(set(current) for current in by_prime_and_trace)))

ring = PolynomialRing(QQ, "u")
u = ring.gen()
field = ring.fraction_field()
weierstrass = model["weierstrass_model"]
A = ring([QQ(value) for value in weierstrass["A_coefficients_low_to_high"]])
B = ring([QQ(value) for value in weierstrass["B_coefficients_low_to_high"]])
surface_discriminant = ring(-16 * (4 * A**3 + 27 * B**2))
curve = EllipticCurve(field, [A, B])
basis = [curve(rational_function(record["X"], ring, field), rational_function(record["Y"], ring, field)) for record in model["sections"]["records"]]

modulus = ZZ.prod(primes)
crt_tuple_count = 0
reconstruction_count = 0
reconstructed_keys = set()
candidates = []
for trace_index in common_trace_indices:
    survivor_lists = [current[trace_index] for current in by_prime_and_trace]
    basis_coordinates = survivor_lists[0][0]["basis_coordinates"]
    if any(item["basis_coordinates"] != basis_coordinates for items in survivor_lists for item in items):
        raise ArithmeticError("trace mismatch across screens")
    trace = sum((coefficient * point for coefficient, point in zip(basis_coordinates, basis) if coefficient), curve(0))
    frame = chord.trace_chord_frame(trace[0], trace[1], ring)
    h, Nx, Ny = (frame[key] for key in ("h", "Nx", "Ny"))
    if h.degree() != 0:
        raise ArithmeticError("screened trace is not polynomial")

    for modular_tuple in product(*survivor_lists):
        crt_tuple_count += 1
        reconstructed = []
        failed = False
        for coordinate_index in range(3):
            residue = ZZ(CRT_list([ZZ(item["m0_m1_m2"][coordinate_index]) for item in modular_tuple], primes))
            try:
                reconstructed.append(QQ(residue.rational_reconstruction(modulus)))
            except (ArithmeticError, ValueError):
                failed = True
                break
        if failed:
            continue
        reconstruction_count += 1
        m0, m1, m2 = reconstructed
        key = (trace_index, m0, m1, m2)
        if key in reconstructed_keys:
            continue
        reconstructed_keys.add(key)

        M = ring(m0 + m1 * u + m2 * u**2)
        q = ring(M**4 - 6 * M**2 * Nx - 8 * M * Ny - 3 * Nx**2 - 4 * A)
        square_part, q_reduced = squareclass_decomposition(q, ring)
        if q_reduced.degree() != 2:
            continue
        sum_x = ring(M**2 - Nx)
        product_x = field(((M * Nx + Ny) ** 2 - B) / Nx)
        if product_x.denominator() != 1:
            raise ArithmeticError("residual x-product is not polynomial")
        product_x = ring(product_x)
        if sum_x**2 - 4 * product_x != q:
            raise ArithmeticError("residual discriminant identity failed")
        x0 = ring(sum_x / 2)
        x1 = ring(square_part / 2)
        intercept = ring(-(Ny + M * Nx))
        y0 = ring(M * x0 + intercept)
        y1 = ring(M * square_part / 2)
        if y0**2 + y1**2 * q_reduced != x0**3 + 3 * x0 * x1**2 * q_reduced + A * x0 + B:
            raise ArithmeticError("constant lift identity failed")
        if 2 * y0 * y1 != 3 * x0**2 * x1 + x1**3 * q_reduced + A * x1:
            raise ArithmeticError("linear lift identity failed")
        if q_reduced.gcd(q_reduced.derivative()).degree() != 0:
            raise ArithmeticError("reduced quadratic is singular")
        branch = {"numerator_coefficients": coefficients(q_reduced), "denominator_coefficients": ["1"]}
        sympy_u = hasher.sp.Symbol("u")
        extension = hasher.extension_key(branch, sympy_u)
        extension_digest = hasher.key_digest(extension)
        candidates.append({
            "label": f"genus3-normalization-trace-{trace_index:05d}",
            "trace_index": trace_index,
            "trace_basis_coordinates": basis_coordinates,
            "trace_translation_orbit_mask": int(survivor_lists[0][0]["translation_orbit_mask"]),
            "slope_coefficients_m0_m1_m2": list(map(rational_text, (m0, m1, m2))),
            "raw_branch_coefficients_low_to_high": coefficients(q),
            "removed_square_factor_coefficients_low_to_high": coefficients(square_part),
            "branch": branch,
            "extension_squareclass": extension,
            "extension_sha256": extension_digest,
            "smooth_atlas_match": smooth_by_digest.get(extension_digest),
            "branch_coprime_to_surface_discriminant": q_reduced.gcd(surface_discriminant).degree() == 0,
            "lifted_section": {
                "field": "QQ(u,r), r^2=q_reduced(u)",
                "x0_coefficients_low_to_high": coefficients(x0),
                "x1_coefficients_low_to_high": coefficients(x1),
                "y0_coefficients_low_to_high": coefficients(y0),
                "y1_coefficients_low_to_high": coefficients(y1),
                "constant_and_linear_identities_verified": True,
            },
        })

collisions = {}
for candidate in candidates:
    collisions.setdefault(candidate["extension_sha256"], []).append(candidate["label"])
collisions = {key: labels for key, labels in collisions.items() if len(labels) >= 2}
output_path = args.output if args.output.is_absolute() else ROOT / args.output
output_path.parent.mkdir(parents=True, exist_ok=True)
payload = {
    "schema": "elkies-k3.r17-norm12-direct-genus3-normalization-reconstruction.v1",
    "status": "PASS_EXACT_RATIONAL_NORMALIZATION_CANDIDATES" if candidates else "PASS_BOUNDED_CRT_NO_EXACT_RATIONAL_NORMALIZATION",
    "search": {
        "primes": list(map(int, primes)),
        "crt_modulus": int(modulus),
        "common_trace_indices": common_trace_indices,
        "crt_tuple_count": crt_tuple_count,
        "successful_coordinatewise_rational_reconstruction_count": reconstruction_count,
        "distinct_reconstructed_parameter_count": len(reconstructed_keys),
    },
    "candidate_count": len(candidates),
    "smooth_atlas_match_count": sum(candidate["smooth_atlas_match"] is not None for candidate in candidates),
    "candidate_collision_count": len(collisions),
    "candidate_collisions": collisions,
    "candidates": candidates,
    "proof_boundary": (
        "All Cartesian tuples of the displayed modular survivors were CRT-combined "
        "and subjected to standard coordinatewise rational reconstruction. Accepted "
        "candidates are factored and their cover sections verified exactly over QQ. "
        "A no-hit result is bounded by the reconstruction modulus and excludes neither "
        "parameters outside that height nor parameters nonintegral at an input prime."
    ),
    "inputs": {relative(path): digest(path) for path in screen_paths + [model_path, smooth_path, CHORD_SCRIPT, HASH_SCRIPT]},
    "reproducing_command": "sage -python elkies-k3/scripts/reconstruct_r17_norm12_direct_genus3_normalizations.sage " + " ".join(f"--screen {relative(path)}" for path in screen_paths),
}
output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"R17GENUS3RECONSTRUCT|primes={','.join(map(str, primes))}|traces={len(common_trace_indices)}"
    f"|crt_tuples={crt_tuple_count}|reconstructed={len(reconstructed_keys)}"
    f"|candidates={len(candidates)}|smooth_matches={payload['smooth_atlas_match_count']}"
    f"|collisions={len(collisions)}|output={relative(output_path)}|status={payload['status']}",
    flush=True,
)
