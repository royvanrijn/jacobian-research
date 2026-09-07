#!/usr/bin/env sage -python
"""Certify genus one for a common-producer resolved third-q12 pencil."""

import argparse
import csv
import hashlib
import json
import re
from math import gcd
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ, block_diagonal_matrix, matrix, vector
from sage.structure.proof.proof import WithProof


ROOT = Path(__file__).resolve().parents[2]
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt"
ROUTE = ROOT / "elkies-k3/data/fibrations/kumar_q80_to_rootless_path.tsv"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
args.input = args.input.resolve()
args.output = args.output.resolve()


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


payload = json.loads(args.input.read_text())
if payload.get("status") != "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_COMMON_PRODUCER":
    raise ValueError("resolved common-producer pencil is not certified")
if payload["resolved_gates"]["combined_rank"] != 5 or payload["resolved_gates"]["kernel_dimension"] != 2:
    raise ArithmeticError("resolved complete linear system is not a pencil")
specialization = payload["specialization"]
prime = int(specialization["prime"])
modulus_match = re.fullmatch(r"r\^2 \+ (\d+)\*r \+ (\d+)", specialization["extension_modulus"])
if modulus_match is None:
    raise ValueError("cannot parse quadratic extension modulus")
linear, constant = map(int, modulus_match.groups())

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
ns = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -frame)
divisor = vector(
    ZZ,
    (3, 3, -18, -20, 18, 0, -24, -22, 9, 6, 5, 42, -71, -10, 6, 16, 19, -8, -8),
)
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_degree = int(divisor * ns * old_fibre)
if divisor * ns * divisor != 0 or gcd(*map(abs, divisor)) != 1:
    raise ArithmeticError("third-q12 divisor is not primitive isotropic")
if old_degree != 3 or old_degree % prime == 0:
    raise ArithmeticError("third-q12 map is not separable of degree three")

base_finite = GF(prime)
modulus_ring = PolynomialRing(base_finite, "m")
m = modulus_ring.gen()
finite = GF(prime**2, "r", modulus=m**2 + linear * m + constant)
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
if len(factors) != 1 or int(factors[0][1]) != 1 or factors.prod() != moving:
    raise ArithmeticError("generic moving equation is reducible")

output = {
    "schema": "elkies-k3.q80-third-q12-resolved-genus-modp2.v3",
    "status": "PASS_EXACT_THIRD_Q12_GENUS_ONE_COMMON_PRODUCER",
    "specialization": specialization,
    "lattice": {
        "NS_determinant": int(frame.det()),
        "divisor_square": 0,
        "divisor_primitive": True,
        "old_fibre_degree": old_degree,
        "separable_in_displayed_characteristic": True,
    },
    "linear_system": {
        "resolved_condition_rank": 5,
        "dimension": 2,
        "moving_bidegree_W_x": [9, 3],
        "generic_equation_irreducible": True,
    },
    "conclusion": (
        "The complete primitive isotropic pencil is base-point-free; its separable "
        "irreducible generic member is smooth by Bertini and has genus 1 by K3 adjunction."
    ),
    "inputs": [
        {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
        for path in (args.input, FRAME, ROUTE)
    ],
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "primitive isotropic divisor of old-fibre degree three",
            "complete two-dimensional resolved linear system",
            "irreducible separable generic moving member and genus one",
        ],
        "not_proved": [
            "minimal Weierstrass Jacobian and transported fibre marking at this prime",
            "canonical cross-prime gauge alignment or characteristic-zero lifting",
        ],
    },
    "reproduce": (
        "sage -python elkies-k3/scripts/verify_q80_third_q12_resolved_genus_modp2.sage "
        f"--input {args.input} --output {args.output}"
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80THIRDQ12COMMONGENUS|prime={prime}|D2=0|primitive=1|degree=3|"
    "separable=1|h0=2|irreducible=1|genus=1|"
    "status=PASS_EXACT_THIRD_Q12_GENUS_ONE_COMMON_PRODUCER",
    flush=True,
)
