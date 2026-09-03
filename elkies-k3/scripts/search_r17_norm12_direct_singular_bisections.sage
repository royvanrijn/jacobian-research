#!/usr/bin/env sage-python
"""Search low-complexity genus-one pencils for rational singular bisections.

On a rootless ``24 I1`` K3, a height-eight trace section ``tau`` has
``tau.O=2``.  In a finite-pole chart its residual-chord genus-one pencil is

    M = M0 + lambda*h^2.

The branch quartic ``q_lambda(u)`` is singular exactly when its discriminant
vanishes.  A rational root of odd multiplicity whose specialized squareclass
has degree two gives a rational singular bisection: after removing the square
factor, its normalization is an exact quadratic cover of the original base.

This script exhausts the displayed coefficient-L1 prefix of exact norm-eight
vectors, verifies every candidate lift over QQ, and compares its full rational
squareclass with the complete smooth rational-bisection atlas.  A no-hit result
is only a bounded equation-complexity statement.
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
DEFAULT_SMOOTH_COLLISIONS = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-alternate-bisection-collisions-full-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-singular-bisection-search-l1le3-v1.json"
)
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
    return (
        str(value.numerator())
        if value.denominator() == 1
        else f"{value.numerator()}/{value.denominator()}"
    )


def coefficients(polynomial):
    return [rational_text(polynomial[index]) for index in range(polynomial.degree() + 1)]


def rational_function(record, ring, field):
    numerator = ring(
        [QQ(value) for value in record["numerator_coefficients_low_to_high"]]
    )
    denominator = ring(
        [QQ(value) for value in record["denominator_coefficients_low_to_high"]]
    )
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


def primitive_integer_polynomial(polynomial):
    """Return the primitive ZZ polynomial defining the same projective roots."""

    denominator = ZZ(polynomial.denominator())
    integer_ring = PolynomialRing(ZZ, polynomial.parent().variable_name())
    result = integer_ring(denominator * polynomial)
    content = ZZ(result.content())
    if content:
        result //= content
    if result.leading_coefficient() < 0:
        result = -result
    return result


def projective_root_count_mod_prime(polynomial, prime):
    """Count distinct P1(F_p) roots of a primitive integer polynomial."""

    field = GF(prime)
    reduced_ring = PolynomialRing(field, polynomial.parent().variable_name())
    reduced = reduced_ring([field(coefficient) for coefficient in polynomial])
    if not reduced:
        return None
    finite_count = sum(reduced(value) == 0 for value in field)
    infinity_count = int(field(polynomial.leading_coefficient()) == 0)
    return int(finite_count + infinity_count)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--smooth-collisions", type=Path, default=DEFAULT_SMOOTH_COLLISIONS)
parser.add_argument("--max-l1", type=int, default=3)
parser.add_argument("--trace-limit", type=int, help="optional deterministic pilot limit")
parser.add_argument("--start", type=int, default=0, help="first row of the ordered prefix")
parser.add_argument("--pari-stack-gb", type=int, default=2)
parser.add_argument(
    "--rational-root-sieve-primes",
    default="17,19,23,29,31",
    help=(
        "comma-separated primes used as exact necessary local tests for a "
        "rational projective root of the squarefree pencil discriminant"
    ),
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if args.max_l1 <= 0 or args.pari_stack_gb <= 0:
    parser.error("--max-l1 and --pari-stack-gb must be positive")
if args.trace_limit is not None and args.trace_limit <= 0:
    parser.error("--trace-limit must be positive")
if args.start < 0:
    parser.error("--start must be nonnegative")
sieve_primes = tuple(ZZ(value) for value in args.rational_root_sieve_primes.split(","))
if not sieve_primes or any(not prime.is_prime() for prime in sieve_primes):
    parser.error("--rational-root-sieve-primes must be a nonempty prime list")

chord = load_script("r17_direct_singular_chord", CHORD_SCRIPT)
hasher = load_script("r17_direct_singular_hasher", HASH_SCRIPT)

model = json.loads(args.model.read_text())
if model.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    raise ValueError("expected a certified direct norm-12 model")
smooth = json.loads(args.smooth_collisions.read_text())
if (
    smooth.get("status") != "PASS_EXTENSION_CANONICALIZATION"
    or not smooth.get("compact_output")
    or int(smooth.get("collision_count", -1)) != 0
    or int(smooth.get("distinct_quadratic_extensions", -1))
    != int(smooth.get("input_bisection_count", -2))
):
    raise ValueError("expected the compact complete smooth-cover collision artifact")
smooth_by_digest = {
    item["extension_sha256"]: item for item in smooth["extension_manifest"]
}

gram = matrix(ZZ, model["sections"]["height_gram"])
pari.allocatemem(args.pari_stack_gb * 1024**3)
minimum = pari(gram).qfminim(8)
short_vectors = matrix(ZZ, minimum[2])
exact_norm_eight_count = 0
minimum_norm_by_mask = {}
best_norm_eight_by_mask = {}
for column in range(short_vectors.ncols()):
    value = short_vectors.column(column)
    norm = int(value * gram * value)
    mask = sum((int(entry) & 1) << index for index, entry in enumerate(value))
    minimum_norm_by_mask[mask] = min(norm, minimum_norm_by_mask.get(mask, norm))
    if norm != 8:
        continue
    exact_norm_eight_count += 1
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
best_norm_eight_by_mask = {
    mask: item
    for mask, item in best_norm_eight_by_mask.items()
    if minimum_norm_by_mask[mask] == 8
}
l1_histogram = Counter(item[0][0] for item in best_norm_eight_by_mask.values())
selected = [
    item[1]
    for item in best_norm_eight_by_mask.values()
    if item[0][0] <= args.max_l1
]
selected.sort(key=lambda value: (sum(abs(entry) for entry in value), tuple(value)))
available_in_prefix = len(selected)
selected = selected[args.start :]
if args.trace_limit is not None:
    selected = selected[: args.trace_limit]

ring = PolynomialRing(QQ, "u")
u = ring.gen()
field = ring.fraction_field()
weierstrass = model["weierstrass_model"]
A = ring([QQ(value) for value in weierstrass["A_coefficients_low_to_high"]])
B = ring([QQ(value) for value in weierstrass["B_coefficients_low_to_high"]])
surface_discriminant = ring(-16 * (4 * A**3 + 27 * B**2))
curve = EllipticCurve(field, [A, B])
basis = [
    curve(
        rational_function(record["X"], ring, field),
        rational_function(record["Y"], ring, field),
    )
    for record in model["sections"]["records"]
]

lambda_ring = PolynomialRing(QQ, "lambda")
lambda_variable = lambda_ring.gen()
bivariate_ring = PolynomialRing(lambda_ring, "u")

trace_records = []
candidates = []
factor_degree_histogram = Counter()
sieve_exclusion_histogram = Counter()
exact_factorization_count = 0
odd_discriminant_irreducibility_exclusion_count = 0
finite_pole_count = 0
for trace_index, trace_vector in enumerate(selected):
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
    record = {
        "trace_index": trace_index,
        "basis_coordinates": list(map(int, trace_vector)),
        "coefficient_l1": int(sum(abs(entry) for entry in trace_vector)),
        "finite_pole_degree": int(h.degree()),
    }
    if h.degree() != 2:
        record["status"] = "SKIPPED_POLE_AT_INFINITY_CHART"
        trace_records.append(record)
        continue
    finite_pole_count += 1

    h_symbolic = bivariate_ring(h)
    Nx_symbolic = bivariate_ring(Nx)
    Ny_symbolic = bivariate_ring(Ny)
    M_symbolic = bivariate_ring(M0) + lambda_variable * h_symbolic**2
    numerator = (
        M_symbolic**4
        - 6 * M_symbolic**2 * Nx_symbolic
        - 8 * M_symbolic * Ny_symbolic
        - 3 * Nx_symbolic**2
        - 4 * bivariate_ring(A) * h_symbolic**4
    )
    q_symbolic, remainder = numerator.quo_rem(h_symbolic**6)
    if remainder or q_symbolic.degree() != 4:
        raise ArithmeticError("symbolic genus-one branch division failed")
    pencil_discriminant = lambda_ring(q_symbolic.discriminant())
    # Sage's univariate ``squarefree_part`` is the representative modulo
    # polynomial squares: factors of even multiplicity disappear and factors
    # of odd multiplicity remain once.  This is exactly the locus relevant to
    # nonsplit normalizations; the ubiquitous double roots are split members.
    odd_multiplicity_discriminant = pencil_discriminant.squarefree_part()
    primitive_odd_discriminant = primitive_integer_polynomial(
        odd_multiplicity_discriminant
    )
    sieve_records = []
    excluded_prime = None
    for sieve_prime in sieve_primes:
        root_count = projective_root_count_mod_prime(
            primitive_odd_discriminant, sieve_prime
        )
        sieve_records.append(
            {
                "prime": int(sieve_prime),
                "projective_root_count": root_count,
            }
        )
        if root_count == 0:
            excluded_prime = int(sieve_prime)
            break
    record.update(
        {
            "pencil_discriminant_degree": int(pencil_discriminant.degree()),
            "odd_multiplicity_pencil_discriminant_degree": int(
                odd_multiplicity_discriminant.degree()
            ),
            "rational_root_sieve": sieve_records,
            "rational_singular_parameters": [],
        }
    )
    if excluded_prime is not None:
        record["status"] = "PASS_MODULAR_NO_RATIONAL_SINGULAR_PARAMETER"
        record["excluding_prime"] = excluded_prime
        sieve_exclusion_histogram[str(excluded_prime)] += 1
        trace_records.append(record)
        continue

    if (
        odd_multiplicity_discriminant.degree() > 1
        and odd_multiplicity_discriminant.is_irreducible()
    ):
        record["status"] = "PASS_EXACT_ODD_DISCRIMINANT_IRREDUCIBLE"
        odd_discriminant_irreducibility_exclusion_count += 1
        trace_records.append(record)
        continue

    exact_factorization_count += 1
    factorization = pencil_discriminant.factor()
    factor_profile = tuple(
        (int(factor.degree()), int(exponent)) for factor, exponent in factorization
    )
    factor_degree_histogram[str(factor_profile)] += 1
    record.update(
        {
            "status": "PASS_EXACT_PENCIL_DISCRIMINANT_FACTORIZATION",
            "factor_degree_multiplicities": [list(item) for item in factor_profile],
        }
    )

    for factor, exponent in factorization:
        if factor.degree() != 1:
            continue
        parameter = QQ(-factor[0] / factor[1])
        M = M0 + parameter * h**2
        specialized_numerator = (
            M**4 - 6 * M**2 * Nx - 8 * M * Ny - 3 * Nx**2 - 4 * A * h**4
        )
        q, specialized_remainder = specialized_numerator.quo_rem(h**6)
        if specialized_remainder:
            raise ArithmeticError("specialized branch division failed")
        square_part, q_reduced = squareclass_decomposition(q, ring)
        reduced_degree = 0 if q_reduced in QQ else int(q_reduced.degree())
        parameter_record = {
            "lambda": rational_text(parameter),
            "discriminant_root_multiplicity": int(exponent),
            "branch_squarefree_degree": reduced_degree,
        }
        record["rational_singular_parameters"].append(parameter_record)
        if int(exponent) % 2 == 0 or reduced_degree != 2:
            continue

        sum_x, sum_remainder = (M**2 - Nx).quo_rem(h**2)
        if sum_remainder:
            raise ArithmeticError("residual x-sum division failed")
        product_x = field(((M * Nx + Ny) ** 2 - B * h**6) / (h**4 * Nx))
        if product_x.denominator() != 1:
            raise ArithmeticError("residual x-product is not polynomial")
        product_x = ring(product_x)
        if sum_x**2 - 4 * product_x != h**2 * q:
            raise ArithmeticError("specialized residual discriminant identity failed")
        x0 = ring(sum_x / 2)
        x1 = ring(h * square_part / 2)
        intercept = field(-(Ny + M * Nx) / h**3)
        y0 = field(M / h) * x0 + intercept
        y1 = field(M / h) * ring(h / 2) * square_part
        if y0.denominator() != 1 or y1.denominator() != 1:
            raise ArithmeticError("normalized singular lift is not polynomial")
        y0, y1 = ring(y0), ring(y1)
        if (
            y0**2 + y1**2 * q_reduced
            != x0**3 + 3 * x0 * x1**2 * q_reduced + A * x0 + B
        ):
            raise ArithmeticError("constant normalized-lift identity failed")
        if (
            2 * y0 * y1
            != 3 * x0**2 * x1 + x1**3 * q_reduced + A * x1
        ):
            raise ArithmeticError("linear normalized-lift identity failed")
        if q_reduced.gcd(q_reduced.derivative()).degree() != 0:
            raise ArithmeticError("normalized branch quadratic is not squarefree")

        branch = {
            "numerator_coefficients": coefficients(q_reduced),
            "denominator_coefficients": ["1"],
        }
        sympy_variable = hasher.sp.Symbol("u")
        extension = hasher.extension_key(branch, sympy_variable)
        extension_digest = hasher.key_digest(extension)
        smooth_match = smooth_by_digest.get(extension_digest)
        candidate = {
            "label": f"singular-trace-{trace_index:04d}",
            "trace_index": trace_index,
            "trace_basis_coordinates": list(map(int, trace_vector)),
            "trace_translation_orbit_mask": sum(
                (int(entry) & 1) << index
                for index, entry in enumerate(trace_vector)
            ),
            "lambda": rational_text(parameter),
            "branch": branch,
            "removed_square_factor_coefficients_low_to_high": coefficients(square_part),
            "extension_squareclass": extension,
            "extension_sha256": extension_digest,
            "smooth_atlas_match": smooth_match,
            "branch_coprime_to_surface_discriminant": (
                q_reduced.gcd(surface_discriminant).degree() == 0
            ),
            "lifted_section": {
                "field": "QQ(u,r), r^2=q_reduced(u)",
                "x0_coefficients_low_to_high": coefficients(x0),
                "x1_coefficients_low_to_high": coefficients(x1),
                "y0_coefficients_low_to_high": coefficients(y0),
                "y1_coefficients_low_to_high": coefficients(y1),
                "constant_and_linear_identities_verified": True,
            },
        }
        candidates.append(candidate)
        parameter_record["candidate_label"] = candidate["label"]
        parameter_record["extension_sha256"] = extension_digest
    trace_records.append(record)

candidate_collisions = {}
for candidate in candidates:
    candidate_collisions.setdefault(candidate["extension_sha256"], []).append(
        candidate["label"]
    )
candidate_collisions = {
    key: labels for key, labels in candidate_collisions.items() if len(labels) >= 2
}

output = args.output if args.output.is_absolute() else ROOT / args.output
payload = {
    "schema": "elkies-k3.r17-norm12-direct-singular-bisection-search.v1",
    "status": (
        "PASS_EXACT_SINGULAR_QUADRATIC_COVER_CANDIDATES"
        if candidates
        else "PASS_BOUNDED_NO_NONSPLIT_RATIONAL_SINGULAR_MEMBER"
    ),
    "search": {
        "norm": 8,
        "coefficient_l1_bound": args.max_l1,
        "trace_limit": args.trace_limit,
        "start": args.start,
        "exact_norm_eight_representatives_up_to_sign": exact_norm_eight_count,
        "minimum_norm_eight_translation_classes": len(best_norm_eight_by_mask),
        "available_in_l1_prefix": available_in_prefix,
        "processed_trace_count": len(selected),
        "processed_half_open_range": [args.start, args.start + len(selected)],
        "finite_pole_trace_count": finite_pole_count,
        "rational_root_sieve_primes": list(map(int, sieve_primes)),
        "modular_rational_root_exclusion_count": sum(
            sieve_exclusion_histogram.values()
        ),
        "modular_rational_root_exclusion_prime_histogram": dict(
            sorted(sieve_exclusion_histogram.items(), key=lambda item: int(item[0]))
        ),
        "exact_discriminant_factorization_count": exact_factorization_count,
        "odd_discriminant_irreducibility_exclusion_count": (
            odd_discriminant_irreducibility_exclusion_count
        ),
        "l1_histogram_best_representative_by_translation_class": {
            str(key): value for key, value in sorted(l1_histogram.items())
        },
        "pencil_discriminant_factor_profile_histogram": dict(
            sorted(factor_degree_histogram.items())
        ),
    },
    "candidate_count": len(candidates),
    "smooth_atlas_match_count": sum(
        candidate["smooth_atlas_match"] is not None for candidate in candidates
    ),
    "candidate_collision_count": len(candidate_collisions),
    "candidate_collisions": candidate_collisions,
    "candidates": candidates,
    "trace_records": trace_records,
    "proof_boundary": (
        "The norm-eight vector enumeration is complete, but only representatives "
        "within the displayed coefficient-L1 prefix are processed. Even a complete "
        "prefix no-hit is not a global exclusion of singular rational bisections. "
        "A zero projective root count for the primitive squarefree discriminant "
        "modulo any displayed prime rigorously excludes rational discriminant roots. "
        "Survivors are factored over QQ. Only odd-multiplicity rational roots with "
        "quadratic normalized squareclass are accepted as nonsplit rational-normalization "
        "candidates."
    ),
    "inputs": {
        relative(path): digest(path)
        for path in (args.model, args.smooth_collisions, CHORD_SCRIPT, HASH_SCRIPT)
    },
    "reproducing_command": (
        "sage -python elkies-k3/scripts/search_r17_norm12_direct_singular_bisections.sage "
        f"--model {relative(args.model)} "
        f"--smooth-collisions {relative(args.smooth_collisions)} "
        f"--max-l1 {args.max_l1} --start {args.start}"
        + ("" if args.trace_limit is None else f" --trace-limit {args.trace_limit}")
        + f" --rational-root-sieve-primes {args.rational_root_sieve_primes}"
        + f" --output {relative(output)}"
    ),
}
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"R17DIRECTSINGULAR|l1={args.max_l1}|traces={len(selected)}"
    f"|finite={finite_pole_count}|candidates={len(candidates)}"
    f"|smooth_matches={payload['smooth_atlas_match_count']}"
    f"|collisions={len(candidate_collisions)}|output={relative(output)}"
    f"|status={payload['status']}",
    flush=True,
)
