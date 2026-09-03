#!/usr/bin/env sage-python
"""Export exact common-cover systems for selected norm-four trace pairs."""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
DEFAULT_MODEL = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
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


def primitive_integer_polynomial(polynomial, integer_ring):
    denominator = ZZ(polynomial.denominator())
    result = integer_ring(denominator * polynomial)
    content = ZZ(result.content())
    if content:
        result //= content
    return result


def parse_pair(value):
    try:
        left, right = map(int, value.split(":"))
    except Exception as error:
        raise argparse.ArgumentTypeError("expected LEFT:RIGHT trace indices") from error
    if left == right or min(left, right) < 0:
        raise argparse.ArgumentTypeError("trace indices must be distinct and nonnegative")
    return left, right


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--max-l1", type=int, default=31)
parser.add_argument("--pair", type=parse_pair, action="append", required=True)
parser.add_argument("--output-dir", type=Path, required=True)
parser.add_argument("--summary", type=Path, required=True)
args = parser.parse_args()

chord = load_script("r17_direct_genus3_collision_msolve_chord", CHORD_SCRIPT)
model_path = args.model.resolve()
model = json.loads(model_path.read_text())
if model.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    raise ValueError("expected a certified direct model")

gram = matrix(ZZ, model["sections"]["height_gram"])
pari.allocatemem(2 * 1024**3)
short_vectors = matrix(ZZ, pari(gram).qfminim(4)[2])
best_norm_four = {}
for column in range(short_vectors.ncols()):
    value = short_vectors.column(column)
    if int(value * gram * value) != 4:
        continue
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
        best_norm_four[mask] = score, candidate
selected = [
    (mask, item[1])
    for mask, item in best_norm_four.items()
    if item[0][0] <= args.max_l1
]
selected.sort(key=lambda item: (sum(abs(entry) for entry in item[1]), tuple(item[1])))
requested = {index for pair in args.pair for index in pair}
if any(index >= len(selected) for index in requested):
    raise ValueError("trace index outside selected prefix")

u_ring = PolynomialRing(QQ, "u")
function_field = u_ring.fraction_field()
weierstrass = model["weierstrass_model"]
A = u_ring([QQ(value) for value in weierstrass["A_coefficients_low_to_high"]])
B = u_ring([QQ(value) for value in weierstrass["B_coefficients_low_to_high"]])
curve = EllipticCurve(function_field, [A, B])
basis = [
    curve(
        rational_function(record["X"], u_ring, function_field),
        rational_function(record["Y"], u_ring, function_field),
    )
    for record in model["sections"]["records"]
]

names = (
    "am0", "am1", "am2", "as0", "as1", "as2",
    "bm0", "bm1", "bm2", "bs0", "bs1", "bs2", "z",
)
coefficient_ring = PolynomialRing(QQ, names=names)
variables = coefficient_ring.gens()
integer_ring = PolynomialRing(ZZ, names=names)
symbolic_u_ring = PolynomialRing(coefficient_ring, "u")
u = symbolic_u_ring.gen()


def trace_system(trace_index, offset):
    mask, trace_vector = selected[trace_index]
    trace = sum(
        (
            coefficient * point
            for coefficient, point in zip(trace_vector, basis)
            if coefficient
        ),
        curve(0),
    )
    frame = chord.trace_chord_frame(trace[0], trace[1], u_ring)
    h, Nx, Ny = (frame[key] for key in ("h", "Nx", "Ny"))
    if h.degree() != 0:
        raise ArithmeticError("selected trace is not in the polynomial chart")
    m0, m1, m2, s0, s1, s2 = variables[offset : offset + 6]
    M = m0 + m1 * u + m2 * u**2
    q = (
        M**4
        - 6 * M**2 * symbolic_u_ring(Nx)
        - 8 * M * symbolic_u_ring(Ny)
        - 3 * symbolic_u_ring(Nx) ** 2
        - 4 * symbolic_u_ring(A)
    )
    q_coefficients = [coefficient_ring(q[index]) for index in range(9)]
    c2 = q_coefficients[8]
    c1 = q_coefficients[7] - 2 * s2 * c2
    c0 = q_coefficients[6] - 2 * s2 * c1 - (s2**2 + 2 * s1) * c2
    square_factor = u**3 + s2 * u**2 + s1 * u + s0
    residual = q - square_factor**2 * (c2 * u**2 + c1 * u + c0)
    equations = [coefficient_ring(residual[index]) for index in range(6)]
    return {
        "trace_index": trace_index,
        "mask": int(mask),
        "basis_coordinates": list(map(int, trace_vector)),
        "cover_coefficients": (c0, c1, c2),
        "equations": equations,
    }


output_dir = args.output_dir if args.output_dir.is_absolute() else ROOT / args.output_dir
output_dir.mkdir(parents=True, exist_ok=True)
systems = []
for left_index, right_index in args.pair:
    left = trace_system(left_index, 0)
    right = trace_system(right_index, 6)
    left_c0, left_c1, left_c2 = left["cover_coefficients"]
    right_c0, right_c1, right_c2 = right["cover_coefficients"]
    equations = left["equations"] + right["equations"]
    equations.extend(
        (
            left_c1 * right_c2 - right_c1 * left_c2,
            left_c0 * right_c2 - right_c0 * left_c2,
            variables[12] * left_c2 * right_c2 - 1,
        )
    )
    integer_equations = [
        primitive_integer_polynomial(equation, integer_ring)
        for equation in equations
    ]
    text = ",".join(names) + "\n0\n"
    text += ",\n".join(
        str(equation).replace("**", "^") for equation in integer_equations
    )
    text += "\n"
    path = output_dir / f"collision-{left_index:04d}-{right_index:04d}.ms"
    path.write_text(text)
    systems.append(
        {
            "trace_indices": [left_index, right_index],
            "translation_orbit_masks": [left["mask"], right["mask"]],
            "basis_coordinates": [
                left["basis_coordinates"], right["basis_coordinates"]
            ],
            "variable_names": list(names),
            "equation_count": len(integer_equations),
            "msolve_input": relative(path),
            "msolve_input_sha256": sha256(text.encode()).hexdigest(),
        }
    )

summary_path = args.summary if args.summary.is_absolute() else ROOT / args.summary
summary_path.parent.mkdir(parents=True, exist_ok=True)
payload = {
    "schema": "elkies-k3.r17-norm12-direct-genus3-collision-msolve-export.v1",
    "status": "PASS_EXACT_CHARACTERISTIC_ZERO_COLLISION_SYSTEMS_EXPORTED",
    "source_label": "norm12-orbit-11952",
    "trace_norm": 4,
    "coefficient_l1_bound": args.max_l1,
    "system_count": len(systems),
    "systems": systems,
    "proof_boundary": (
        "Each system imposes both monic-cubic normalization systems and equality "
        "of the two residual monic quadratic polynomial atoms, saturated by both "
        "degree-eight leading slope coefficients. A common atom is necessary but "
        "not sufficient for equality of the scalar squareclasses. Lower-degree "
        "branch charts and solution certification are separate gates."
    ),
    "inputs": {
        relative(model_path): digest(model_path),
        relative(CHORD_SCRIPT): digest(CHORD_SCRIPT),
    },
}
summary_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"R17GENUS3COLLISIONMSOLVEEXPORT|systems={len(systems)}"
    f"|output={relative(summary_path)}|status={payload['status']}",
    flush=True,
)
