#!/usr/bin/env sage-python
"""Export bounded mod-p searches for rational normalizations of genus-2 bisections.

For a norm-six trace section on a rootless ``24 I1`` K3, ``tau.O=1``.
Writing its finite-pole chord frame as ``tau=(Nx/h^2,Ny/h^3)``, every regular
residual-chord slope numerator is

    M = M0 + (l0 + l1*u)*h^2.

The resulting branch polynomial ``q(u)`` has degree at most six.  Its double
cover normalizes to a rational curve in the finite two-node chart precisely
when

    q(u) = (u^2 + a*u + b)^2 * (c2*u^2 + c1*u + c0).

The three ``c`` variables are eliminated from the top coefficients, leaving
four equations in ``l0,l1,a,b``.  This script exports those exact reductions
to msolve for the displayed coefficient-L1 prefix.  Modular solutions are a
discovery sieve only; characteristic-zero reconstruction and exact lift
verification are separate gates.
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
DEFAULT_OUTPUT_DIR = (
    ROOT / "artifacts/local/elkies-k3/r17-norm12-genus2-normalizations"
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
        raise ZeroDivisionError("prime divides a rational coefficient denominator")
    return field(value.numerator()) / field(value.denominator())


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--prime", type=int, required=True)
parser.add_argument("--max-l1", type=int, default=1)
parser.add_argument("--start", type=int, default=0)
parser.add_argument("--trace-limit", type=int)
parser.add_argument("--pari-stack-gb", type=int, default=2)
parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
parser.add_argument("--output", type=Path)
args = parser.parse_args()
if not ZZ(args.prime).is_prime():
    parser.error("--prime must be prime")
if args.max_l1 <= 0 or args.start < 0 or args.pari_stack_gb <= 0:
    parser.error("invalid nonpositive bound or negative start")
if args.trace_limit is not None and args.trace_limit <= 0:
    parser.error("--trace-limit must be positive")

chord = load_script("r17_direct_genus2_chord", CHORD_SCRIPT)
model_path = args.model.resolve()
model = json.loads(model_path.read_text())
if model.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    raise ValueError("expected a certified direct norm-12 model")

gram = matrix(ZZ, model["sections"]["height_gram"])
pari.allocatemem(args.pari_stack_gb * 1024**3)
minimum = pari(gram).qfminim(6)
short_vectors = matrix(ZZ, minimum[2])
minimum_norm_by_mask = {}
best_norm_six_by_mask = {}
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
    if mask not in best_norm_six_by_mask or score < best_norm_six_by_mask[mask][0]:
        best_norm_six_by_mask[mask] = (score, candidate)
best_norm_six_by_mask = {
    mask: item
    for mask, item in best_norm_six_by_mask.items()
    if minimum_norm_by_mask[mask] == 6
}
l1_histogram = Counter(item[0][0] for item in best_norm_six_by_mask.values())
selected = [
    (mask, item[1])
    for mask, item in best_norm_six_by_mask.items()
    if item[0][0] <= args.max_l1
]
selected.sort(key=lambda item: (sum(abs(entry) for entry in item[1]), tuple(item[1])))
available_in_prefix = len(selected)
selected = selected[args.start :]
if args.trace_limit is not None:
    selected = selected[: args.trace_limit]

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

parameter_ring = PolynomialRing(QQ, names=("l0", "l1"))
l0, l1 = parameter_ring.gens()
symbolic_u_ring = PolynomialRing(parameter_ring, "u")
symbolic_u = symbolic_u_ring.gen()

prime = args.prime
finite_field = GF(prime)
coefficient_ring = PolynomialRing(
    finite_field, names=("l0", "l1", "a", "b"), order="degrevlex"
)
l0p, l1p, ap, bp = coefficient_ring.gens()
modp_u_ring = PolynomialRing(coefficient_ring, "u")
modp_u = modp_u_ring.gen()


def parameter_polynomial_mod_p(polynomial):
    polynomial = parameter_ring(polynomial)
    return coefficient_ring(
        sum(
            mod_p(coefficient, finite_field) * l0p**exponents[0] * l1p**exponents[1]
            for exponents, coefficient in polynomial.dict().items()
        )
    )

tag = args.model.stem.replace("elkies-k3-r17-norm12-", "")
output_dir = (args.output_dir / tag / f"p{prime}").resolve()
output_dir.mkdir(parents=True, exist_ok=True)
systems = []
trace_records = []
for local_index, (mask, trace_vector) in enumerate(selected):
    global_index = args.start + local_index
    trace = sum(
        (
            coefficient * point
            for coefficient, point in zip(trace_vector, basis)
            if coefficient
        ),
        curve(0),
    )
    frame = chord.trace_chord_frame(trace[0], trace[1], ring)
    h, Nx, Ny, M0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
    trace_record = {
        "trace_index": global_index,
        "translation_orbit_mask": int(mask),
        "basis_coordinates": list(map(int, trace_vector)),
        "coefficient_l1": int(sum(abs(entry) for entry in trace_vector)),
        "finite_pole_degree": int(h.degree()),
    }
    if h.degree() != 1:
        trace_record["status"] = "SKIPPED_POLE_AT_INFINITY_CHART"
        trace_records.append(trace_record)
        continue

    h_symbolic = symbolic_u_ring(h)
    Nx_symbolic = symbolic_u_ring(Nx)
    Ny_symbolic = symbolic_u_ring(Ny)
    M = symbolic_u_ring(M0) + (l0 + l1 * symbolic_u) * h_symbolic**2
    numerator = (
        M**4
        - 6 * M**2 * Nx_symbolic
        - 8 * M * Ny_symbolic
        - 3 * Nx_symbolic**2
        - 4 * symbolic_u_ring(A) * h_symbolic**4
    )
    q, remainder = numerator.quo_rem(h_symbolic**6)
    if remainder or q.degree() > 6:
        raise ArithmeticError("symbolic genus-two branch division failed")

    try:
        q_modp = modp_u_ring(
            [parameter_polynomial_mod_p(coefficient) for coefficient in q]
        )
    except ZeroDivisionError:
        trace_record["status"] = "SKIPPED_BAD_REDUCTION_DENOMINATOR"
        trace_records.append(trace_record)
        continue
    q_coefficients = [q_modp[index] for index in range(7)]
    c2 = q_coefficients[6]
    c1 = q_coefficients[5] - 2 * ap * c2
    c0 = q_coefficients[4] - 2 * ap * c1 - (ap**2 + 2 * bp) * c2
    square_factor = modp_u**2 + ap * modp_u + bp
    residual_quadratic = c2 * modp_u**2 + c1 * modp_u + c0
    factor_residual = q_modp - square_factor**2 * residual_quadratic
    equations = [coefficient_ring(factor_residual[index]) for index in range(4)]
    if not all(equation != 0 for equation in equations):
        raise ArithmeticError("normalization system contains a zero equation")

    system_path = output_dir / f"trace-{global_index:05d}.ms"
    with system_path.open("w") as handle:
        handle.write("l0,l1,a,b\n")
        handle.write(f"{prime}\n")
        for equation_index, equation in enumerate(equations):
            handle.write(str(equation).replace("**", "^"))
            handle.write(",\n" if equation_index + 1 < len(equations) else "\n")
    systems.append(
        {
            "trace_index": global_index,
            "translation_orbit_mask": int(mask),
            "basis_coordinates": list(map(int, trace_vector)),
            "path": relative(system_path),
            "sha256": digest(system_path),
            "q_degree_over_QQ_l0_l1": int(q.degree()),
            "q_coefficients_in_QQ_l0_l1_low_to_high": [
                str(parameter_ring(coefficient)) for coefficient in q
            ],
        }
    )
    trace_record["status"] = "PASS_EXACT_MODP_TWO_NODE_SYSTEM_EXPORT"
    trace_record["system_sha256"] = systems[-1]["sha256"]
    trace_records.append(trace_record)

if args.output is None:
    output_path = output_dir / "export.json"
else:
    output_path = args.output if args.output.is_absolute() else ROOT / args.output
output_path.parent.mkdir(parents=True, exist_ok=True)
payload = {
    "schema": "elkies-k3.r17-norm12-direct-genus2-normalization-msolve-export.v1",
    "status": "PASS_EXACT_MODP_TWO_NODE_SYSTEM_EXPORT",
    "proof_boundary": (
        "Each exported system is the exact reduction modulo the displayed prime "
        "of the finite-chart factorization q=(u^2+a*u+b)^2*Q2. Modular solutions "
        "are only a discovery sieve; they do not establish a characteristic-zero "
        "normalization or section. Charts with a trace pole at infinity are skipped."
    ),
    "prime": prime,
    "search": {
        "norm": 6,
        "coefficient_l1_bound": args.max_l1,
        "start": args.start,
        "trace_limit": args.trace_limit,
        "exact_norm_six_representatives_up_to_sign": exact_norm_six_count,
        "minimum_norm_six_translation_classes": len(best_norm_six_by_mask),
        "available_in_l1_prefix": available_in_prefix,
        "selected_trace_count": len(selected),
        "exported_system_count": len(systems),
        "l1_histogram_best_representative_by_translation_class": {
            str(key): value for key, value in sorted(l1_histogram.items())
        },
    },
    "systems": systems,
    "trace_records": trace_records,
    "inputs": {
        relative(path): digest(path) for path in (model_path, CHORD_SCRIPT)
    },
    "reproducing_command": (
        "sage -python "
        "elkies-k3/scripts/export_r17_norm12_direct_genus2_normalizations_modp.sage "
        f"--prime {prime} --max-l1 {args.max_l1} --start {args.start}"
        + ("" if args.trace_limit is None else f" --trace-limit {args.trace_limit}")
    ),
}
output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"R17GENUS2EXPORT|p={prime}|selected={len(selected)}|systems={len(systems)}"
    f"|output={relative(output_path)}|status={payload['status']}",
    flush=True,
)
