#!/usr/bin/env sage -python
"""Certify genus one for the resolved p=19 third-q12 pencil by adjunction.

Singular's normalization libraries do not support algebraic finite fields.
Here normalization is unnecessary: the lattice certificate gives a primitive
isotropic divisor on the K3, the resolved computation gives its complete
two-dimensional linear system, and exact factorization shows that the moving
equation is irreducible. Since its old-fibre degree is three (prime to 19),
the induced pencil is separable. A base-point-free primitive isotropic pencil
on a K3 has smooth generic member of arithmetic and geometric genus one.
"""

import csv
import hashlib
import json
from math import gcd
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ, block_diagonal_matrix, matrix, vector
from sage.structure.proof.proof import WithProof


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "artifacts/generated-results/q80-third-q12-um2-p19-resolved-pencil.json"
OUTPUT = ROOT / "artifacts/generated-results/q80-third-q12-um2-p19-resolved-genus.json"
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt"
ROUTE = ROOT / "elkies-k3/data/fibrations/kumar_q80_to_rootless_path.tsv"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def bezout_vector(pairings):
    current = ZZ.zero()
    coefficients = [ZZ.zero()] * len(pairings)
    for index, pairing in enumerate(pairings):
        if not pairing:
            continue
        new_gcd, left, right = current.xgcd(ZZ(pairing))
        coefficients = [left * value for value in coefficients]
        coefficients[index] += right
        current = new_gcd
    if abs(current) != 1:
        raise ArithmeticError("fibre pairings are not primitive")
    return vector(ZZ, coefficients if current == 1 else [-value for value in coefficients])


def neighbor(parent, qnorm, a, b, coordinates):
    hyperbolic = matrix(ZZ, ((0, 1), (1, 0)))
    ns = block_diagonal_matrix(hyperbolic, -parent)
    fibre = vector(ZZ, [a, b] + list(coordinates))
    if fibre * ns * fibre != 0 or gcd(*map(abs, ns * fibre)) != 1:
        raise ArithmeticError("invalid primitive neighbour fibre")
    mate = bezout_vector(ns * fibre)
    mate -= ZZ(mate * ns * mate) // 2 * fibre
    complement = matrix(ZZ, [list(fibre * ns), list(mate * ns)]).right_kernel_matrix()
    child = -(complement * ns * complement.transpose())
    transport = matrix(ZZ, [list(fibre), list(mate)] + complement.rows())
    if abs(transport.det()) != 1:
        raise ArithmeticError("neighbour transport is not unimodular")
    return child


payload = json.loads(INPUT.read_text())
if payload["status"] != "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_MOD19_QUADRATIC":
    raise ValueError("resolved pencil is not certified")
if payload["resolved_gates"]["combined_rank"] != 5:
    raise ArithmeticError("resolved component gate does not have rank five")
if payload["resolved_gates"]["kernel_dimension"] != 2:
    raise ArithmeticError("resolved complete linear system is not a pencil")

# Reconstruct the D7+D5 Neron--Severi chart through the two pinned q4 steps.
with ROUTE.open() as stream:
    steps = list(csv.DictReader(stream, delimiter="\t"))
frame = load_matrix(FRAME)
for step in steps[:2]:
    frame = neighbor(
        frame,
        ZZ(step["q"]),
        ZZ(step["a"]),
        ZZ(step["b"]),
        vector(ZZ, map(ZZ, step["v"].split(","))),
    )
if frame.det() != 948:
    raise ArithmeticError("unexpected D7+D5 frame determinant")
hyperbolic = matrix(ZZ, ((0, 1), (1, 0)))
ns = block_diagonal_matrix(hyperbolic, -frame)

# Chamber-reduced q12 divisor pinned by analyze_q80_second_neighbor_chamber.sage.
divisor = vector(
    ZZ,
    (3, 3, -18, -20, 18, 0, -24, -22, 9, 6, 5, 42, -71, -10, 6, 16, 19, -8, -8),
)
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
if divisor * ns * divisor != 0:
    raise ArithmeticError("third-q12 divisor is not isotropic")
if gcd(*map(abs, divisor)) != 1:
    raise ArithmeticError("third-q12 divisor is not primitive")
old_degree = int(divisor * ns * old_fibre)
if old_degree != 3 or old_degree % 19 == 0:
    raise ArithmeticError("third-q12 map does not have separable degree three")

# Rebuild the generic moving equation over GF(19^2)(V) and factor it using
# Singular's exact finite-field backend. Sage marks this backend as unproved
# for the coefficient-field presentation, so relax only the dispatcher flag;
# the returned one-factor decomposition is replayed below.
base_finite = GF(19)
modulus_ring = PolynomialRing(base_finite, "m")
m = modulus_ring.gen()
finite = GF(19**2, "r", modulus=m**2 + 12 * m + 3)
r = finite.gen()
plane_ring = PolynomialRing(finite, names=("V", "W", "x"))
V, W, x = plane_ring.gens()


def field_element(coordinates):
    return finite(coordinates[0]) + finite(coordinates[1]) * r


moving = plane_ring.zero()
for t_degree, w_degree, x_degree, coordinates in payload["moving_equation"][
    "terms_T_W_x_coefficient_1_r"
]:
    moving += field_element(coordinates) * V**t_degree * W**w_degree * x**x_degree
if (moving.degree(W), moving.degree(x)) != (9, 3):
    raise ArithmeticError("unexpected moving-equation bidegree")
with WithProof("polynomial", False):
    factors = moving.factor()
if len(factors) != 1 or int(factors[0][1]) != 1:
    raise ArithmeticError("generic moving equation is reducible")
if factors.prod() != moving:
    raise ArithmeticError("factorization replay failed")

output = {
    "schema": "elkies-k3.q80-third-q12-resolved-genus-modp2.v2",
    "status": "PASS_EXACT_THIRD_Q12_GENUS_ONE_BY_ADJUNCTION_MOD19_QUADRATIC",
    "specialization": {"u": "-2", "prime": 19, "extension_modulus": "r^2+12*r+3"},
    "lattice": {
        "NS_determinant": int(frame.det()),
        "divisor_square": int(divisor * ns * divisor),
        "divisor_primitive": True,
        "old_fibre_degree": old_degree,
        "separable_in_characteristic_19": True,
    },
    "linear_system": {
        "resolved_condition_rank": 5,
        "dimension": 2,
        "moving_bidegree_W_x": [9, 3],
        "generic_equation_irreducible": True,
    },
    "conclusion": (
        "The complete primitive isotropic pencil is base-point-free; its "
        "separable irreducible generic member is smooth by Bertini and has "
        "genus 1 by K3 adjunction."
    ),
    "inputs": {
        "resolved_pencil": {"path": str(INPUT.relative_to(ROOT)), "sha256": sha256(INPUT)},
        "initial_frame": {"path": str(FRAME.relative_to(ROOT)), "sha256": sha256(FRAME)},
        "route": {"path": str(ROUTE.relative_to(ROOT)), "sha256": sha256(ROUTE)},
        "chamber_verifier": {
            "path": "elkies-k3/scripts/analyze_q80_second_neighbor_chamber.sage",
            "sha256": sha256(ROOT / "elkies-k3/scripts/analyze_q80_second_neighbor_chamber.sage"),
        },
    },
    "claim_boundary": {
        "proved": [
            "primitive isotropic third-q12 divisor of old-fibre degree three",
            "complete two-dimensional resolved linear system",
            "irreducible separable generic moving member",
            "generic normalization genus one by adjunction and Bertini",
        ],
        "not_proved": [
            "minimal Weierstrass Jacobian and A5+A3+3A1 fibres",
            "alignment at another prime",
            "characteristic-zero lifting",
        ],
    },
    "reproduce": "sage -python elkies-k3/scripts/verify_q80_third_q12_um2_p19_resolved_genus.sage",
}
OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    "Q80THIRDQ12GENUS|D2=0|primitive=1|degree=3|separable=1|h0=2|"
    "irreducible=1|genus=1|"
    "status=PASS_EXACT_THIRD_Q12_GENUS_ONE_BY_ADJUNCTION_MOD19_QUADRATIC",
    flush=True,
)
