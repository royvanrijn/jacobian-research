#!/usr/bin/env sage -python
"""Compile the exact source-q8 smooth collision principal-part matrix.

At the four actual smooth P1.O collision fibres, use the transported frame

    q=(m-y(P1)/x(P1))/h,   X=h^2*x.

The q8 line-bundle condition in this frame is that every negative h-principal
part vanish.  The coefficient template is the exact QQ formula certified by
the principal-part derivation.  This program proves its rank by a full-column
good reduction, avoiding a needlessly enormous dense QQ h-adic expansion. It
does not declare a complete q8 cover: E8, the marked E7 chart, all generic E7
components, edge nodes, and overlaps are separate condition blocks.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, binomial, gcd, matrix


ROOT = Path(__file__).resolve().parents[2]
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
AMBIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-endpoint-rr-ambient.json"
FRAME = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-smooth-collision-frame.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-exact.json"

exec(compile(CORE.read_text(), str(CORE), "exec"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--extra-h-power", type=int, default=0,
    help="raise every endpoint denominator exponent by this nonnegative amount",
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if args.extra_h_power < 0:
    raise ValueError("extra-h-power must be nonnegative")

p1 = json.loads(P1.read_text())
endpoint = json.loads(AMBIENT.read_text())
frame = json.loads(FRAME.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1"
assert endpoint["status"] == "PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT"
assert frame["status"] == "PASS_EXACT_Q8_SMOOTH_COLLISION_FRAME"
assert endpoint["ambient_dimension"] == 54

# The endpoint construction is nested: the extra denominator power changes
# only k, leaving the actual E8 floor and marked-E7 upper inequality intact.
if args.extra_h_power:
    ambient_basis = []
    for family in endpoint["families"]:
        k = int(family["h_power"]) + args.extra_h_power
        lower = int(family["e8_minimal_u_power"])
        upper = 4*k + int(family["e7_allowed_t_denominator_power"])
        for i in range(lower, upper+1):
            ambient_basis.append({
                **family["generic_basis"], "u_power": i, "h_power": k,
            })
else:
    ambient_basis = endpoint["ambient_basis"]
assert len(ambient_basis) == 54+72*args.extra_h_power

ring = PolynomialRing(QQ, "u")
h = polynomial(ring, p1["structured_denominator"]["Z4_coefficients"])
assert h.degree() == 4 and gcd(h, h.derivative()) == 1
pole_bound = 15+args.extra_h_power

coordinates = tuple(
    (int(item["X_power"]), int(item["q_power"]))
    for item in frame["regular_degree_18_frame"]
)
assert coordinates == tuple([(0, power) for power in range(10)] + [(1, power) for power in range(8)])
coordinate_index = {entry: index for index, entry in enumerate(coordinates)}
# A full-column-rank good reduction is already an exact rank certificate. We
# therefore evaluate the exact template directly in GF(43), avoiding gigantic
# intermediate rational h-adic inverses in a dense QQ expansion.
rank_prime = 43
finite = GF(rank_prime)


def reduce_coefficient(value):
    value = QQ(value)
    denominator = finite(ZZ(value.denominator()))
    assert denominator
    return finite(ZZ(value.numerator()))/denominator


finite_ring = PolynomialRing(finite, "u")
u_finite = finite_ring.gen()
finite_field = finite_ring.fraction_field()
h_finite = finite_ring([reduce_coefficient(value) for value in h.list()])
finite_modulus = h_finite**pole_bound
assert gcd(h_finite, h_finite.derivative()) == 1
residue_dimension = finite_modulus.degree()
x_p_finite = finite_field(finite_ring([
    reduce_coefficient(value) for value in p1["x_entrance_base"]["numerator_coefficients"]
])) / finite_field(finite_ring([
    reduce_coefficient(value) for value in p1["x_entrance_base"]["denominator_coefficients"]
]))
y_p_finite = finite_field(finite_ring([
    reduce_coefficient(value) for value in p1["y_entrance_base"]["numerator_coefficients"]
])) / finite_field(finite_ring([
    reduce_coefficient(value) for value in p1["y_entrance_base"]["denominator_coefficients"]
]))
rho_finite = finite_field(h_finite)*y_p_finite/x_p_finite
assert gcd(h_finite, finite_ring(rho_finite.numerator())) == 1
assert gcd(h_finite, finite_ring(rho_finite.denominator())) == 1


def finite_residue(value):
    value = finite_field(value)
    numerator = finite_ring(value.numerator())
    denominator = finite_ring(value.denominator())
    assert gcd(denominator, finite_modulus) == 1
    return finite_ring((numerator*denominator.inverse_mod(finite_modulus)) % finite_modulus)


def finite_principal_part_coordinates(entry):
    a, b = int(entry["x_power"]), int(entry["m_power"])
    i, k = int(entry["u_power"]), int(entry["h_power"])
    result = [finite.zero()]*(len(coordinates)*residue_dimension)
    for j in range(b+1):
        exponent = 2*j-b-k-2*a
        if exponent >= 0:
            continue
        value = (finite(binomial(b, j))*u_finite**i*rho_finite**(b-j)
                 *h_finite**(pole_bound+exponent))
        remainder = finite_residue(value)
        offset = coordinate_index[(a, j)]*residue_dimension
        for degree, coefficient in enumerate(remainder.list()):
            result[offset+degree] += coefficient
    return result


def finite_principal_part_image(entry):
    """Sparse coordinates of the actual finite smooth quotient image."""

    values = finite_principal_part_coordinates(entry)
    return {
        (coordinates[index // residue_dimension][0],
         coordinates[index // residue_dimension][1],
         index % residue_dimension): value
        for index, value in enumerate(values)
        if value
    }


# The declared h-adic quotient has 18*degree(h^pole_bound) coordinates, but
# a finite ambient sees only a finite subset.  Compile that exact image rather
# than pretending that the quotient needs a hand-written dense matrix.
smooth_block = finite_ambient_image_condition(
    "q8 smooth collision principal parts modulo 43",
    tuple(range(len(ambient_basis))),
    lambda index: finite_principal_part_image(ambient_basis[index]),
    lambda coordinate: coordinate,
    finite,
    "actual q/X smooth-collision frame reduced at the rank-certifying prime",
)
reduced_matrix = smooth_block["matrix"]
reduced_rank = reduced_matrix.rank()
assert reduced_rank == len(ambient_basis)
condition_rows = len(coordinates)*residue_dimension
condition_columns = len(ambient_basis)
exact_rank = condition_columns
exact_kernel_dimension = 0

payload = {
    "schema": "elkies-k3.h92-q8-smooth-principal-parts-exact.v1",
    "status": "PASS_EXACT_Q8_SMOOTH_PRINCIPAL_PART_CONDITION_BLOCK",
    "inputs": {
        "p1": {"path": str(P1.relative_to(ROOT)), "sha256": digest(P1)},
        "compiler_core": {"path": str(CORE.relative_to(ROOT)), "sha256": digest(CORE)},
        "endpoint_ambient": {"path": str(AMBIENT.relative_to(ROOT)), "sha256": digest(AMBIENT)},
        "smooth_frame": {"path": str(FRAME.relative_to(ROOT)), "sha256": digest(FRAME)},
    },
    "parameters": {"extra_h_power": int(args.extra_h_power), "pole_bound": int(pole_bound)},
    "ambient_basis": ambient_basis,
    "condition": {
        "frame": "q=(m-y(P1)/x(P1))/h, X=h^2*x",
        "rule": "all negative h-principal parts vanish",
        "coordinate_count": int(len(coordinates)),
        "h_adic_residue_dimension": int(residue_dimension),
        "rows": int(condition_rows),
        "finite_ambient_image_rows": int(reduced_matrix.nrows()),
        "columns": int(condition_columns),
        "rank": int(exact_rank),
        "kernel_dimension": int(exact_kernel_dimension),
        "rank_certificate": {
            "prime": int(rank_prime),
            "reduced_rank": int(reduced_rank),
            "argument": (
                "The exact QQ matrix has at most its column count as rank; "
                "its reduction mod 43 has that full column rank, so its QQ "
                "rank is exactly the same."
            ),
            "compiler_block": {
                "name": smooth_block["name"],
                "provenance": smooth_block["provenance"],
                "finite_image_coordinate_count": int(reduced_matrix.nrows()),
            },
        },
    },
    "boundary": (
        "This is the exact smooth collision condition block only. It is not "
        "a complete q8 resolved-chart cover and does not certify h0, a pencil, "
        "a child equation, a bisection, or a rank statement."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8SMOOTHEXACT|extra_h={}|rows={}|columns={}|rank={}|kernel={}|"
    "status=PASS_EXACT_Q8_SMOOTH_PRINCIPAL_PART_CONDITION_BLOCK".format(
        args.extra_h_power, condition_rows, condition_columns,
        exact_rank, exact_kernel_dimension,
    ),
    flush=True,
)
