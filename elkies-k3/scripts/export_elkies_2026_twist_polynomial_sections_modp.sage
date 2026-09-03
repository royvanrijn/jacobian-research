#!/usr/bin/env sage-python
"""Export reduced polynomial-section systems for one quadratic twist.

status: ACTIVE_COMPILER
claim: exact finite-field P.O=0 section-system export for one declared twist
inputs: a certified Weierstrass model and certified twist-character records
outputs: artifacts/local/elkies-k3/twist-polynomial-sections/<candidate>/p<prime>

For an elliptic surface of arithmetic genus ``chi`` in short Weierstrass
form, a section disjoint from the zero section has

    deg(X) <= 2*chi,   deg(Y) <= 3*chi.

After fixing the affine point met at infinity, the coefficients of ``Y`` are
recovered recursively from high degree.  The remaining coefficient equations
in the ``2*chi`` nonleading X-coefficients are exported to msolve, one system
per rational affine point on the infinity fibre.

This is an exhaustive finite-field polynomial-section scheme only when the
reported infinity fibre is smooth and has no rational 2-torsion.  It is not a
characteristic-zero rank certificate.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

from sage.all import GF, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

from screen_elkies_2026_quadratic_twist_ranks import (  # noqa: E402
    Candidate,
    DEFAULT_BISECTIONS,
    DEFAULT_MODEL,
    DEFAULT_PAIRS,
    load_candidates,
    multiply_integer_polynomials,
    square_equivalent_integer_polynomial,
    valuation,
)

DEFAULT_GENUS_ONE_CONSTRUCTIONS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
)


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def polynomial_coefficients_mod(values, prime):
    minimum = min(valuation(int(value), prime) for value in values)
    if minimum % 2:
        raise ValueError(f"twist has odd p-adic content at p={prime}")
    divisor = prime**minimum
    return [int(value // divisor) % prime for value in values]


parser = argparse.ArgumentParser(description=__doc__)
target = parser.add_mutually_exclusive_group(required=True)
target.add_argument("--singleton-mask", type=int)
target.add_argument("--product-key")
target.add_argument("--genus-one-label")
target.add_argument(
    "--direct-label",
    help="label in a direct norm-12 bisection-extension artifact",
)
target.add_argument(
    "--direct-product-key",
    help="two labels, separated by ':', in a direct bisection-extension artifact",
)
parser.add_argument("--prime", type=int, required=True)
parser.add_argument(
    "--allow-infinity-two-torsion",
    action="store_true",
    help=(
        "export only non-2-torsion leading blocks when every smooth chart fibre "
        "has rational 2-torsion; this is a discovery sieve, not an exhaustive scheme"
    ),
)
parser.add_argument("--bisections", type=Path, default=DEFAULT_BISECTIONS)
parser.add_argument("--pairs", type=Path, default=DEFAULT_PAIRS)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument(
    "--genus-one-constructions", type=Path, default=DEFAULT_GENUS_ONE_CONSTRUCTIONS
)
parser.add_argument(
    "--output-dir",
    type=Path,
    default=ROOT / "artifacts/local/elkies-k3/twist-polynomial-sections",
)
args = parser.parse_args()

prime = int(args.prime)
if prime < 5 or not GF(prime).is_prime_field():
    raise ValueError("the modular section exporter requires an odd prime at least five")
field = GF(prime)

if args.singleton_mask is not None:
    singletons, products, unused_schemas = load_candidates(args.bisections, args.pairs)
    key = str(args.singleton_mask)
    candidate = next((item for item in singletons if item.key == key), None)
    chi = 3
elif args.product_key is not None:
    singletons, products, unused_schemas = load_candidates(args.bisections, args.pairs)
    key = str(args.product_key)
    candidate = next((item for item in products if item.key == key), None)
    chi = 4
elif args.direct_label is not None or args.direct_product_key is not None:
    key = str(args.direct_label or args.direct_product_key)
    direct_covers = json.loads(args.bisections.read_text())
    if direct_covers.get("schema") != "elkies-k3.bisection-extension-input.v1":
        raise ValueError("direct cover input has the wrong generic schema")
    by_label = {item["label"]: item for item in direct_covers["bisections"]}
    if args.direct_label is not None:
        source_record = by_label.get(key)
        candidate = None if source_record is None else Candidate(
            kind="direct_singleton",
            key=key,
            masks=(int(source_record["lattice_orbit_mask"]),),
            coefficients=square_equivalent_integer_polynomial(
                source_record["branch"]["numerator_coefficients"]
            ),
            forced_twist_rank=1,
            metadata={"orbit_hex": f"0x{int(source_record['lattice_orbit_mask']):05x}"},
        )
        chi = 3
    else:
        labels = tuple(key.split(":"))
        if len(labels) != 2 or labels[0] == labels[1]:
            raise ValueError("--direct-product-key requires two distinct labels")
        source_records = tuple(by_label.get(label) for label in labels)
        if any(record is None for record in source_records):
            candidate = None
        else:
            factors = tuple(
                square_equivalent_integer_polynomial(
                    record["branch"]["numerator_coefficients"]
                )
                for record in source_records
            )
            candidate = Candidate(
                kind="direct_product",
                key=key,
                masks=tuple(int(record["lattice_orbit_mask"]) for record in source_records),
                coefficients=multiply_integer_polynomials(*factors),
                forced_twist_rank=0,
                metadata={
                    "labels": list(labels),
                    "orbit_hex": [
                        f"0x{int(record['lattice_orbit_mask']):05x}"
                        for record in source_records
                    ],
                },
            )
        chi = 4
else:
    key = str(args.genus_one_label)
    constructions = json.loads(args.genus_one_constructions.read_text())
    source_record = next(
        (
            item
            for item in constructions["construction"]["records"]
            if item["label"] == key
        ),
        None,
    )
    candidate = None if source_record is None else Candidate(
        kind="genus_one",
        key=key,
        masks=(int(source_record["lattice_orbit_mask"]),),
        coefficients=square_equivalent_integer_polynomial(
            source_record["branch_polynomial_q_coefficients_low_to_high"]
        ),
        forced_twist_rank=1,
        metadata={"orbit_hex": f"0x{int(source_record['lattice_orbit_mask']):05x}"},
    )
    chi = 4
if candidate is None:
    raise ValueError(f"unknown twist candidate {key}")

model = json.loads(args.model.read_text())
if model.get("status") == "PASS_TRANSCRIBED_PUBLISHED_R17_MODEL":
    model_coefficients = model
elif model.get("status") == "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    model_coefficients = model["weierstrass_model"]
else:
    raise ValueError("expected a certified published or direct norm-12 model")
base_a = [field(QQ(value)) for value in model_coefficients["A_coefficients_low_to_high"]]
base_b = [field(QQ(value)) for value in model_coefficients["B_coefficients_low_to_high"]]
twist_q = polynomial_coefficients_mod(candidate.coefficients, prime)

base_ring = PolynomialRing(field, "t")
t = base_ring.gen()
A0 = base_ring(base_a)
B0 = base_ring(base_b)
q = base_ring(twist_q)
if q == 0:
    raise ArithmeticError("twist polynomial reduces to zero")
if candidate.kind == "direct_product":
    base_discriminant = -field(16) * (field(4) * A0**3 + field(27) * B0**2)
    if (
        q.degree() != 4
        or not q.is_squarefree()
        or base_discriminant.degree() != 24
        or not base_discriminant.is_squarefree()
        or q.gcd(base_discriminant).degree()
    ):
        raise ArithmeticError(
            "direct product does not have squarefree good reduction away from the 24I1 fibres"
        )
A = A0 * q**2
B = B0 * q**3
if A.degree() > 4 * chi or B.degree() > 6 * chi:
    raise ArithmeticError("twist exceeds the expected Weierstrass degree bounds")

x_degree = 2 * chi
y_degree = 3 * chi


def infinity_fibre_data(coefficient_a, coefficient_b):
    discriminant = -field(16) * (
        field(4) * coefficient_a**3 + field(27) * coefficient_b**2
    )
    points = []
    two_torsion = []
    for leading_x in field:
        rhs = leading_x**3 + coefficient_a * leading_x + coefficient_b
        if rhs == 0:
            two_torsion.append(int(leading_x))
        for leading_y in field:
            if leading_y**2 == rhs:
                points.append((leading_x, leading_y))
    return discriminant, points, two_torsion


chart_parameter = None
complete_infinity_cover = True
a_infinity = field(A[4 * chi])
b_infinity = field(B[6 * chi])
infinity_discriminant, leading_points, rational_two_torsion_x = infinity_fibre_data(
    a_infinity, b_infinity
)
if not infinity_discriminant or rational_two_torsion_x:
    # Move a finite fibre t=c to infinity by t=c+1/s and the standard
    # Weierstrass scaling.  Global P.O=0 sections remain polynomial with the
    # same degree bounds in this chart.
    for candidate_parameter in field:
        candidate_a = field(A(candidate_parameter))
        candidate_b = field(B(candidate_parameter))
        candidate_discriminant, candidate_points, candidate_two_torsion = (
            infinity_fibre_data(candidate_a, candidate_b)
        )
        if candidate_discriminant and not candidate_two_torsion:
            chart_parameter = field(candidate_parameter)
            transformed_a = sum(
                A[index] * (chart_parameter * t + 1) ** index * t ** (4 * chi - index)
                for index in range(A.degree() + 1)
            )
            transformed_b = sum(
                B[index] * (chart_parameter * t + 1) ** index * t ** (6 * chi - index)
                for index in range(B.degree() + 1)
            )
            A = base_ring(transformed_a)
            B = base_ring(transformed_b)
            a_infinity = candidate_a
            b_infinity = candidate_b
            infinity_discriminant = candidate_discriminant
            leading_points = candidate_points
            rational_two_torsion_x = candidate_two_torsion
            break
if (not infinity_discriminant or rational_two_torsion_x) and args.allow_infinity_two_torsion:
    # Some small characteristics have rational 2-torsion on every smooth
    # rational chart fibre.  The high-to-low Y recursion still applies to
    # every leading point with y != 0.  Exporting those blocks is useful as a
    # discovery sieve, but sections meeting a rational 2-torsion point at
    # infinity are omitted and the union must not be called exhaustive.
    fallback = None
    for candidate_parameter in field:
        candidate_a = field(A(candidate_parameter))
        candidate_b = field(B(candidate_parameter))
        candidate_discriminant, candidate_points, candidate_two_torsion = (
            infinity_fibre_data(candidate_a, candidate_b)
        )
        non_two_torsion_points = [
            point for point in candidate_points if point[1] != 0
        ]
        if candidate_discriminant and non_two_torsion_points:
            fallback = (
                candidate_parameter,
                candidate_a,
                candidate_b,
                candidate_discriminant,
                non_two_torsion_points,
                candidate_two_torsion,
            )
            break
    if fallback is not None:
        (
            chart_parameter,
            a_infinity,
            b_infinity,
            infinity_discriminant,
            leading_points,
            rational_two_torsion_x,
        ) = fallback
        transformed_a = sum(
            A[index] * (chart_parameter * t + 1) ** index * t ** (4 * chi - index)
            for index in range(A.degree() + 1)
        )
        transformed_b = sum(
            B[index] * (chart_parameter * t + 1) ** index * t ** (6 * chi - index)
            for index in range(B.degree() + 1)
        )
        A = base_ring(transformed_a)
        B = base_ring(transformed_b)
        complete_infinity_cover = False
if not infinity_discriminant or (rational_two_torsion_x and complete_infinity_cover):
    raise ArithmeticError(
        "no smooth rational chart fibre without rational 2-torsion was found"
    )

tag = (
    f"singleton-{key}"
    if candidate.kind == "singleton"
    else f"product-{key.replace(':', '-')}"
    if candidate.kind == "product"
    else f"direct-singleton-{key}"
    if candidate.kind == "direct_singleton"
    else f"direct-product-{key.replace(':', '--')}"
    if candidate.kind == "direct_product"
    else f"genus-one-{key}"
)
output_dir = args.output_dir.resolve() / tag / f"p{prime}"
output_dir.mkdir(parents=True, exist_ok=True)

systems = []
names = tuple(f"x{index}" for index in range(x_degree - 1, -1, -1))
for block_index, (leading_x, leading_y) in enumerate(leading_points):
    coefficient_ring = PolynomialRing(field, names=names, order="degrevlex")
    variables = coefficient_ring.gens_dict()
    polynomial_ring = PolynomialRing(coefficient_ring, "T")
    T = polynomial_ring.gen()
    X = coefficient_ring(leading_x) * T**x_degree + sum(
        variables[f"x{index}"] * T**index for index in range(x_degree)
    )
    AS = polynomial_ring([coefficient_ring(value) for value in A])
    BS = polynomial_ring([coefficient_ring(value) for value in B])
    rhs = X**3 + AS * X + BS

    y_coefficients = {y_degree: coefficient_ring(leading_y)}
    for degree in range(2 * y_degree - 1, y_degree - 1, -1):
        index = degree - y_degree
        known = sum(
            y_coefficients[left] * y_coefficients[degree - left]
            for left in y_coefficients
            if degree - left in y_coefficients
        )
        y_coefficients[index] = (
            rhs[degree] - known
        ) / (2 * coefficient_ring(leading_y))
    Y = sum(y_coefficients[index] * T**index for index in range(y_degree + 1))
    residual = Y**2 - rhs
    equations = [coefficient_ring(residual[degree]) for degree in range(y_degree)]
    if not all(equation != 0 for equation in equations):
        raise ArithmeticError("recursive section system contains a zero equation")

    system_path = output_dir / f"block-{block_index:03d}.ms"
    with system_path.open("w") as handle:
        handle.write(",".join(names) + f"\n{prime}\n")
        for equation_index, equation in enumerate(equations):
            handle.write(str(equation).replace("**", "^"))
            handle.write(",\n" if equation_index + 1 < len(equations) else "\n")
    systems.append(
        {
            "block_index": block_index,
            "leading_x_y": [int(leading_x), int(leading_y)],
            "path": str(system_path.relative_to(ROOT)),
            "sha256": digest(system_path),
        }
    )

record = {
    "schema": "elkies-k3.elkies-2026-twist-polynomial-section-msolve-export.v1",
    "status": (
        "PASS_EXACT_MODP_REDUCED_POLYNOMIAL_SECTION_EXPORT"
        if complete_infinity_cover
        else "PASS_EXACT_MODP_NON_TWO_TORSION_BLOCK_EXPORT"
    ),
    "proof_boundary": (
        "The union of the exported systems is the complete polynomial P.O=0 section "
        "scheme over the displayed finite field because the infinity fibre is smooth and "
        "has no rational 2-torsion. Solver results and characteristic-zero lifting are separate."
        if complete_infinity_cover
        else "Only leading blocks with y != 0 are exported. Sections meeting a rational "
        "2-torsion point at infinity are omitted, so this is a discovery sieve and not a "
        "complete polynomial-section scheme."
    ),
    "candidate": {
        "kind": candidate.kind,
        "key": candidate.key,
        "masks": list(candidate.masks),
        "chi": chi,
        "x_degree_bound": x_degree,
        "y_degree_bound": y_degree,
    },
    "prime": prime,
    "reduced_twist_q_coefficients_low_to_high": [int(value) for value in q.list()],
    "twist_degrees_A_B": [int(A.degree()), int(B.degree())],
    "twist_A_coefficients_low_to_high": [
        int(A[index]) for index in range(4 * chi + 1)
    ],
    "twist_B_coefficients_low_to_high": [
        int(B[index]) for index in range(6 * chi + 1)
    ],
    "infinity_fibre": {
        "chart": (
            "original_infinity"
            if chart_parameter is None
            else f"original_t={int(chart_parameter)} via t={int(chart_parameter)}+1/s"
        ),
        "a_b": [int(a_infinity), int(b_infinity)],
        "discriminant": int(infinity_discriminant),
        "affine_point_count": len(leading_points),
        "rational_two_torsion_x": rational_two_torsion_x,
    },
    "systems": systems,
    "inputs": {
        str(path.resolve().relative_to(ROOT)): digest(path)
        for path in (
            args.bisections,
            args.model,
            *((args.pairs,) if candidate.kind in {"singleton", "product"} else ()),
            *((args.genus_one_constructions,) if candidate.kind == "genus_one" else ()),
        )
    },
}
record_path = output_dir / "export.json"
record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
print(
    f"TWISTPOLYEXPORT|kind={candidate.kind}|key={candidate.key}|p={prime}"
    f"|chi={chi}|blocks={len(systems)}|variables={len(names)}"
    f"|equations={y_degree}|output={record_path}|status=PASS_EXPORTED",
    flush=True,
)
