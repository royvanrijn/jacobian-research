#!/usr/bin/env sage-python
"""Prove geometric Picard rank 19 for the rational G720 ``3I6`` source.

The exact source certificate supplies 19 rational divisor classes.  Point
counts over F_p and F_(p^2) at p=17,19 show that each reduction has geometric
Picard rank 20.  Their Artin--Tate discriminant square classes are incompatible,
so van Luijk's two-prime argument forces characteristic-zero Picard rank 19.

The three split I6 fibres require a correction of ``15*q`` to the point count
of the singular short Weierstrass model over F_q.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = (
    ROOT / "artifacts/generated-results/elkies-k3-golay-det720-3a5-source-qq-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-3a5-picard19-v1.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mod_rational(value, p):
    value = QQ(value)
    numerator = int(value.numerator() % p)
    denominator = int(value.denominator() % p)
    if not denominator:
        raise ZeroDivisionError((value, p))
    return numerator * pow(denominator, -1, p) % p


def eval_fp(coefficients, value, p):
    answer = 0
    for coefficient in reversed(coefficients):
        answer = (answer * value + coefficient) % p
    return answer


def mul_fp2(left, right, p, nonsquare):
    a, b = left
    c, d = right
    return ((a * c + nonsquare * b * d) % p, (a * d + b * c) % p)


def eval_fp2(coefficients, value, p, nonsquare):
    answer = (0, 0)
    for coefficient in reversed(coefficients):
        answer = mul_fp2(answer, value, p, nonsquare)
        answer = ((answer[0] + coefficient) % p, answer[1])
    return answer


def first_nonsquare(p):
    return next(value for value in range(2, p) if pow(value, (p - 1) // 2, p) == p - 1)


def legendre_table(p):
    table = np.empty(p, dtype=np.int8)
    table[0] = 0
    for value in range(1, p):
        table[value] = 1 if pow(value, (p - 1) // 2, p) == 1 else -1
    return table


def count_raw_weierstrass_fp(A, B, p):
    chi = legendre_table(p)
    xs = np.arange(p, dtype=np.int64)
    x3 = xs * xs % p * xs % p
    character_sum = 0
    for t_value in range(p):
        rhs = (x3 + eval_fp(A, t_value, p) * xs + eval_fp(B, t_value, p)) % p
        character_sum += int(chi[rhs].sum())
    rhs = (x3 + A[8] * xs + B[12]) % p
    character_sum += int(chi[rhs].sum())
    return int((p + 1) ** 2 + character_sum)


def count_raw_weierstrass_fp2(A, B, p, nonsquare):
    chi = legendre_table(p)
    coordinates = np.arange(p, dtype=np.int64)
    xa = np.repeat(coordinates, p)
    xb = np.tile(coordinates, p)
    x2a = (xa * xa + nonsquare * xb * xb) % p
    x2b = (2 * xa * xb) % p
    x3a = (x2a * xa + nonsquare * x2b * xb) % p
    x3b = (x2a * xb + x2b * xa) % p

    def fibre_character_sum(av, bv):
        aa, ab = av
        ba, bb = bv
        rhs_a = (x3a + aa * xa + nonsquare * ab * xb + ba) % p
        rhs_b = (x3b + aa * xb + ab * xa + bb) % p
        norm = (rhs_a * rhs_a - nonsquare * rhs_b * rhs_b) % p
        return int(chi[norm].sum())

    character_sum = 0
    for ta in range(p):
        for tb in range(p):
            value = (ta, tb)
            character_sum += fibre_character_sum(
                eval_fp2(A, value, p, nonsquare),
                eval_fp2(B, value, p, nonsquare),
            )
    character_sum += fibre_character_sum((A[8], 0), (B[12], 0))
    q = p**2
    return int((q + 1) ** 2 + character_sum)


def reduction_record(A, B, residual, p):
    Rp = PolynomialRing(GF(p), "t")
    t = Rp.gen()
    A_mod = Rp([GF(p)(mod_rational(value, p)) for value in A.list()])
    B_mod = Rp([GF(p)(mod_rational(value, p)) for value in B.list()])
    residual_mod = Rp([GF(p)(mod_rational(value, p)) for value in residual.list()])
    assert A_mod.degree() == 8 and B_mod.degree() == 12
    assert residual_mod.degree() == 6 and residual_mod.is_squarefree()
    assert residual_mod.gcd(t * (t - 1) * A_mod) == 1
    assert A_mod(0) and A_mod(1) and A_mod[8]

    A_values = [int(value) for value in A_mod.list()]
    B_values = [int(value) for value in B_mod.list()]
    nonsquare = first_nonsquare(p)
    # A split I_n fibre has n*q points, while its nodal Weierstrass cubic has
    # q points.  Three I6 fibres therefore add 3*(6-1)*q = 15*q.
    n1_raw = count_raw_weierstrass_fp(A_values, B_values, p)
    n2_raw = count_raw_weierstrass_fp2(A_values, B_values, p, nonsquare)
    n1 = n1_raw + 15 * p
    n2 = n2_raw + 15 * p**2

    trace_h2_1 = ZZ(n1 - 1 - p**2)
    trace_h2_2 = ZZ(n2 - 1 - p**4)
    residual_trace_1 = trace_h2_1 - 19 * p
    residual_trace_2 = trace_h2_2 - 19 * p**2
    possible_signs = [
        epsilon
        for epsilon in (-1, 1)
        if residual_trace_1**2 - 2 * epsilon * p * residual_trace_1 == residual_trace_2
    ]
    assert len(possible_signs) == 1
    epsilon = possible_signs[0]
    pair_trace = residual_trace_1 - epsilon * p
    assert pair_trace not in (-2 * p, -p, 0, p, 2 * p)
    discriminant_class = ZZ(pair_trace**2 - 4 * p**2)
    return {
        "p": p,
        "fp2_nonsquare": nonsquare,
        "raw_weierstrass_point_count_fp": int(n1_raw),
        "raw_weierstrass_point_count_fp2": int(n2_raw),
        "split_I6_resolution_correction_fp": 15 * p,
        "split_I6_resolution_correction_fp2": 15 * p**2,
        "smooth_K3_point_count_fp": int(n1),
        "smooth_K3_point_count_fp2": int(n2),
        "trace_h2_fp": int(trace_h2_1),
        "trace_h2_frob_squared": int(trace_h2_2),
        "residual_trace_fp": int(residual_trace_1),
        "residual_trace_frob_squared": int(residual_trace_2),
        "linear_residual_eigenvalue_sign": epsilon,
        "remaining_pair_trace": int(pair_trace),
        "remaining_pair_normalized_root_of_unity_excluded": True,
        "geometric_picard_rank": 20,
        "artin_tate_discriminant_square_class_representative": int(discriminant_class),
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

source_path = arguments.source.resolve()
output_path = arguments.output.resolve()
source = json.loads(source_path.read_text())
assert source["status"] == "PASS_EXACT_QQ_3I6_MW2_RANK19_SUBLATTICE_DET720"
assert source["weierstrass_model"]["fibre_profile"] == "3I6+6I1"
assert all(row["split_over_Q"] for row in source["weierstrass_model"]["split_reducible_fibres"])
assert source["lattice"]["explicit_NS_sublattice_rank"] == 19

R = PolynomialRing(QQ, "t")
A = R(source["weierstrass_model"]["A_coefficients_low_to_high"])
B = R(source["weierstrass_model"]["B_coefficients_low_to_high"])
residual = R(source["weierstrass_model"]["residual_I1_polynomial_monic"])
reductions = [reduction_record(A, B, residual, p) for p in (17, 19)]
assert reductions[0]["smooth_K3_point_count_fp"] == 600
assert reductions[0]["smooth_K3_point_count_fp2"] == 89624
assert reductions[0]["artin_tate_discriminant_square_class_representative"] == -256
assert reductions[1]["smooth_K3_point_count_fp"] == 702
assert reductions[1]["smooth_K3_point_count_fp2"] == 136824
assert reductions[1]["artin_tate_discriminant_square_class_representative"] == -1440
ratio = QQ(-1440) / QQ(-256)
assert ratio == QQ(45) / 8 and not ratio.is_square()

payload = {
    "schema": "elkies-k3.golay-det720-3a5-picard19.v1",
    "status": "PASS_EXACT_TWO_PRIME_ARTIN_TATE_PICARD19",
    "reproduce": (
        "sage -python elkies-k3/scripts/certify_golay_det720_3a5_picard19.sage"
    ),
    "inputs": {relative(source_path): digest(source_path)},
    "known_characteristic_zero_divisor_rank": 19,
    "reductions": reductions,
    "artin_tate_discriminant_class_ratio": str(ratio),
    "ratio_is_square": False,
    "geometric_picard_rank_characteristic_zero": 19,
    "argument": (
        "The 19 displayed rational divisor classes give rho>=19.  At both good "
        "primes the two point counts leave one algebraic eigenvalue and one "
        "non-root-of-unity reciprocal pair, so rho=20 in reduction.  If rho over "
        "Qbar were 20, its NS discriminant square class would agree in both "
        "reductions; the Artin-Tate ratio 45/8 is nonsquare, so rho<=19."
    ),
    "proof_boundary": (
        "This proves Picard rank 19, but does not by itself prove saturation of "
        "the displayed Mordell-Weil basis or a physical neighbour corridor to "
        "the abstract MW17 target."
    ),
}

encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != encoded:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(encoded)

print(
    "GOLAY7203A5PICARD|p=17,19|rho_reductions=20,20|"
    "disc_ratio=45/8|rho_Qbar=19|status=PASS"
)
