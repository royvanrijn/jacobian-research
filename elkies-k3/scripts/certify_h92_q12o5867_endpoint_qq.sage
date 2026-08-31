#!/usr/bin/env sage -python
"""Certify the q12/o5867 endpoint identity, Picard rank, and full MW group.

This is deliberately a small endpoint calculation.  It does three things:

1. The q12 binary quartic is evaluated at the old finite I2 support v=0.
   The value is an exact square in QQ(u).  Pointing the quartic there gives a
   generalized Weierstrass equation whose short invariants differ from the
   stored endpoint only by x -> 9*x, y -> 27*y.  Hence the last Jacobian is
   the original genus-one fibration, not a nontrivial torsor.
2. Exact point counts at the good primes 131 and 137 give rank-20 reductions
   with Artin--Tate discriminant square classes -948 and -948*37.  Their
   incompatibility forces geometric Picard rank 19 in characteristic zero.
3. The determinant-948 section lattice has only one possible index-two
   integral overlattice.  Its new vector has odd norm 73, whereas the full
   rootless K3 Mordell--Weil height lattice is even.  Thus the 17 displayed
   sections are saturated and are the full Mordell--Weil group.

The finite-field counts use only vectorized arithmetic in F_p and the norm
character of F_{p^2}; no Groebner basis or surface elimination is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
Q8_MODEL = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"
Q12_MODEL = LOCAL / "q12o5867-smooth-rr-qq.json"
HEIGHTS = LOCAL / "q12o5867-rootless-height-basis-qq.json"
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=GENERATED / "elkies-k3-h3-q12o5867-endpoint-certificate.json",
    )
    return parser.parse_args()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hash_strings(values):
    raw = json.dumps(list(map(str, values)), separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def rational_bits(values):
    answer = 0
    for value in values:
        value = QQ(value)
        answer = max(
            answer,
            abs(ZZ(value.numerator())).nbits(),
            ZZ(value.denominator()).nbits(),
        )
    return int(answer)


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def mod_rational(value, p):
    value = QQ(value)
    numerator = int(value.numerator() % p)
    denominator = int(value.denominator() % p)
    if denominator == 0:
        raise ZeroDivisionError((value, p))
    return numerator * pow(denominator, -1, p) % p


def coefficients_mod_p(poly, p):
    return [mod_rational(value, p) for value in poly.list()]


def eval_fp(coefficients, value, p):
    answer = 0
    for coefficient in reversed(coefficients):
        answer = (answer * value + coefficient) % p
    return answer


def mul_fp2(left, right, p, nonsquare):
    a, b = left
    c, d = right
    return ((a*c+nonsquare*b*d) % p, (a*d+b*c) % p)


def eval_fp2(coefficients, value, p, nonsquare):
    answer = (0, 0)
    for coefficient in reversed(coefficients):
        answer = mul_fp2(answer, value, p, nonsquare)
        answer = ((answer[0]+coefficient) % p, answer[1])
    return answer


def first_nonsquare(p):
    for value in range(2, p):
        if pow(value, (p-1)//2, p) == p-1:
            return value
    raise ArithmeticError(p)


def legendre_table(p):
    table = np.empty(p, dtype=np.int8)
    table[0] = 0
    for value in range(1, p):
        symbol = pow(value, (p-1)//2, p)
        table[value] = 1 if symbol == 1 else -1
    return table


def count_surface_fp(A_coefficients, B_coefficients, p):
    """Count the smooth elliptic K3 over F_p, including the infinity fibre."""
    chi = legendre_table(p)
    xs = np.arange(p, dtype=np.int64)
    x3 = xs*xs % p * xs % p
    total_character = 0
    for t in range(p):
        av = eval_fp(A_coefficients, t, p)
        bv = eval_fp(B_coefficients, t, p)
        rhs = (x3+av*xs+bv) % p
        total_character += int(chi[rhs].sum())
    # The minimal model has weights (4,6); the infinity coefficients are the
    # leading coefficients of A and B because their degrees are exactly 8,12.
    rhs = (x3+A_coefficients[8]*xs+B_coefficients[12]) % p
    total_character += int(chi[rhs].sum())
    return int((p+1)**2+total_character)


def count_surface_fp2(A_coefficients, B_coefficients, p, nonsquare):
    """Count over F_{p^2}=F_p[w]/(w^2-nonsquare) by its norm character."""
    chi = legendre_table(p)
    coordinates = np.arange(p, dtype=np.int64)
    xa = np.repeat(coordinates, p)
    xb = np.tile(coordinates, p)
    x2a = (xa*xa+nonsquare*xb*xb) % p
    x2b = (2*xa*xb) % p
    x3a = (x2a*xa+nonsquare*x2b*xb) % p
    x3b = (x2a*xb+x2b*xa) % p

    def fibre_character_sum(av, bv):
        aa, ab = av
        ba, bb = bv
        rhs_a = (x3a+aa*xa+nonsquare*ab*xb+ba) % p
        rhs_b = (x3b+aa*xb+ab*xa+bb) % p
        norm = (rhs_a*rhs_a-nonsquare*rhs_b*rhs_b) % p
        return int(chi[norm].sum())

    total_character = 0
    for ta in range(p):
        for tb in range(p):
            t = (ta, tb)
            total_character += fibre_character_sum(
                eval_fp2(A_coefficients, t, p, nonsquare),
                eval_fp2(B_coefficients, t, p, nonsquare),
            )
    total_character += fibre_character_sum(
        (A_coefficients[8], 0), (B_coefficients[12], 0)
    )
    q = p*p
    return int((q+1)**2+total_character)


def good_reduction_and_counts(A, B, p):
    Rp = PolynomialRing(GF(p), "u")
    A_mod = Rp(coefficients_mod_p(A, p))
    B_mod = Rp(coefficients_mod_p(B, p))
    delta_mod = -16*(4*A_mod**3+27*B_mod**2)
    assert A_mod.degree() == 8 and B_mod.degree() == 12
    assert delta_mod.degree() == 24 and delta_mod.is_squarefree()
    assert delta_mod.gcd(A_mod) == 1  # all bad fibres remain nodal I1
    nonsquare = first_nonsquare(p)
    A_coefficients = [int(value) for value in A_mod.list()]
    B_coefficients = [int(value) for value in B_mod.list()]
    n1 = count_surface_fp(A_coefficients, B_coefficients, p)
    n2 = count_surface_fp2(A_coefficients, B_coefficients, p, nonsquare)

    trace_h2_1 = ZZ(n1-1-p**2)
    trace_h2_2 = ZZ(n2-1-p**4)
    residual_trace_1 = trace_h2_1-19*p
    residual_trace_2 = trace_h2_2-19*p**2
    possible_signs = [
        epsilon for epsilon in (-1, 1)
        if residual_trace_1**2-2*epsilon*p*residual_trace_1 == residual_trace_2
    ]
    assert len(possible_signs) == 1
    epsilon = possible_signs[0]
    pair_trace = residual_trace_1-epsilon*p
    assert pair_trace not in (-2*p, -p, 0, p, 2*p)
    discriminant_class = ZZ(pair_trace**2-4*p**2)
    return {
        "p": p,
        "fp2_nonsquare": nonsquare,
        "point_count_fp": int(n1),
        "point_count_fp2": int(n2),
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


args = parse_args()
started = time.monotonic()
q8 = json.loads(Q8_MODEL.read_text())
q12 = json.loads(Q12_MODEL.read_text())
heights = json.loads(HEIGHTS.read_text())
assert q8["status"] == "PASS_EXACT_QQ_Q8O376_4A1_RR_JACOBIAN_AND_P1229_ZERO"
assert q8["child"]["fibre_profile"] == "4I2+16I1"
assert any(item["factor"] == "u" for item in q8["child"]["finite_reducible_fibres"])
assert q12["status"] == "PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN"
assert q12["divisor"] == {
    "P_dot_O": 10,
    "class": "O+P-4F",
    "fibre_twist": -4,
    "vertical_layers": 0,
    "vertical_support": 0,
}
assert q12["child"]["degrees_A_B_Delta"] == [8, 12, 24]
assert q12["child"]["root_rank"] == 0
assert q12["child"]["infinity"]["kodaira"] == "smooth"
assert heights["status"] == "PASS_EXACT_QQ_Q12O5867_ROOTLESS_RANK17_HEIGHT_BASIS_PINNED"

R = PolynomialRing(QQ, "u")
K = R.fraction_field()
quartic_coefficients = [
    R([QQ(value) for value in row])
    for row in q12["binary_quartic"]["coefficients_in_old_v_low_to_high"]
]
assert len(quartic_coefficients) == 5
ee, dd, cc, bb, aa = quartic_coefficients
assert ee.is_square()
ordinate = R(ee.sqrt())
assert ordinate**2 == ee

# Standard nonbranch pointed-quartic model at (v,w)=(0,ordinate).
a1 = K(dd)/ordinate
a2 = K(cc)-K(dd)**2/(4*ordinate**2)
a3 = 2*ordinate*K(bb)
a4 = -4*ordinate**2*K(aa)
a6 = a2*a4
b2 = a1**2+4*a2
b4 = 2*a4+a1*a3
b6 = a3**2+4*a6
c4 = b2**2-24*b4
c6 = -b2**3+36*b2*b4-216*b6
A_pointed = -c4/48
B_pointed = -c6/864
A = R([QQ(value) for value in q12["child"]["minimal_A_coefficients_low_to_high"]])
B = R([QQ(value) for value in q12["child"]["minimal_B_coefficients_low_to_high"]])
assert 81*A_pointed == K(A)
assert 729*B_pointed == K(B)

reductions = [good_reduction_and_counts(A, B, p) for p in (131, 137)]
assert reductions[0]["point_count_fp"] == 19308
assert reductions[0]["point_count_fp2"] == 294853764
assert reductions[0]["artin_tate_discriminant_square_class_representative"] == -23700
assert reductions[1]["point_count_fp"] == 21036
assert reductions[1]["point_count_fp2"] == 352653204
assert reductions[1]["artin_tate_discriminant_square_class_representative"] == -35076
disc_ratio = QQ(
    reductions[1]["artin_tate_discriminant_square_class_representative"]
) / reductions[0]["artin_tate_discriminant_square_class_representative"]
assert disc_ratio == QQ(37)/25 and not disc_ratio.is_square()

gram = matrix(ZZ, heights["height_gram"])
pinned = load_matrix(PINNED)
assert gram.det() == pinned.det() == 948
assert gram.is_positive_definite() and pinned.is_positive_definite()
# Any proper finite-index integral overlattice must have index 2, since n^2
# divides 948.  The mod-2 radical gives its unique possible new coset.
possible_indices = [n for n in range(2, ZZ(948).isqrt()+1) if 948 % n**2 == 0]
assert possible_indices == [2]
kernel_mod2 = pinned.change_ring(GF(2)).right_kernel()
assert kernel_mod2.dimension() == 1
candidate = vector(ZZ, [ZZ(value) for value in kernel_mod2.basis()[0]])
assert all(value % 2 == 0 for value in pinned*candidate)
candidate_norm = QQ(candidate*pinned*candidate)/4
assert candidate_norm == 73 and ZZ(candidate_norm) % 2 == 1

input_paths = (Q8_MODEL, Q12_MODEL, HEIGHTS, PINNED)
payload = {
    "schema": "elkies-k3.h92-q12o5867-endpoint-certificate.v1",
    "status": "PASS_EXACT_Q12O5867_SOURCE_IDENTITY_RHO19_FULL_MW_R17",
    "reproducing_command": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/certify_h92_q12o5867_endpoint_qq.sage"
    ),
    "source_identity": {
        "old_q8_fibre_support": "v=0 finite I2",
        "quartic_point": {
            "old_base_coordinate": "0",
            "ordinate_coefficients_low_to_high": list(map(str, ordinate.list())),
            "ordinate_degree": int(ordinate.degree()),
            "exact_square_identity": True,
        },
        "pointed_generalized_weierstrass_a_invariants_sha256": hash_strings((a1, a2, a3, a4, a6)),
        "pointed_invariant_identities": [
            "81*A_pointed=A_endpoint",
            "729*B_pointed=B_endpoint",
        ],
        "short_isomorphism": "x_endpoint=9*x_pointed, y_endpoint=27*y_pointed",
        "conclusion": (
            "The q12 genus-one pencil has an exact QQ(u)-point, so its stored "
            "Jacobian is the same elliptic K3.  Composed with the already-pointed "
            "source-to-q8 forward chain, the endpoint is exactly the certified H3 source."
        ),
    },
    "picard_rank": {
        "known_characteristic_zero_classes": 19,
        "reductions": reductions,
        "tate_conjecture_use": (
            "For a K3 surface in odd characteristic, geometric Picard rank equals "
            "the number of Frobenius eigenvalues that are p times roots of unity."
        ),
        "artin_tate_use": (
            "Over F_(p^2), the rank-20 NS discriminant square class is represented "
            "by r^2-4p^2, where r is the trace of the remaining reciprocal pair."
        ),
        "discriminant_class_ratio": str(disc_ratio),
        "ratio_is_square": False,
        "specialization_argument": (
            "If rho over characteristic zero were 20, its NS discriminant square "
            "class would agree with both rank-20 reductions.  The ratio 37/25 is "
            "not a square, hence rho<=19; the 19 explicit classes give rho=19."
        ),
        "geometric_picard_rank_characteristic_zero": 19,
    },
    "mordell_weil": {
        "shioda_tate_rank": 17,
        "torsion": "trivial",
        "displayed_lattice_determinant": 948,
        "possible_proper_integral_overlattice_indices": possible_indices,
        "unique_index_two_coset_mod2": list(map(int, candidate)),
        "candidate_new_vector_norm": str(candidate_norm),
        "candidate_is_even": False,
        "evenness_reason": (
            "With only I1 fibres, section pairings are integral and every self-height "
            "is 4+2(P.O), so the full Mordell--Weil height lattice is even."
        ),
        "saturated": True,
        "full_geometric_mordell_weil_lattice": "pinned R17",
        "full_geometric_mordell_weil_rank": 17,
        "full_geometric_mordell_weil_determinant": 948,
    },
    "method": {
        "point_counts": "exact vectorized F_p and norm-character F_(p^2) sums",
        "large_groebner_required": False,
        "surface_elimination_required": False,
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "Exact source identity, geometric Picard rank 19, trivial torsion, and "
        "full saturated geometric Mordell--Weil lattice R17 of rank 17 and "
        "determinant 948. Arithmetic specialization is not attempted."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in input_paths],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in input_paths},
    },
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q12O5867ENDPOINT|source_identity=1|rho=19|MW_rank=17|torsion=0|"
    "det=948|saturated=1|counts=131:{},{};137:{},{}|status={}|seconds={:.3f}|output={}".format(
        reductions[0]["point_count_fp"], reductions[0]["point_count_fp2"],
        reductions[1]["point_count_fp"], reductions[1]["point_count_fp2"],
        payload["status"], payload["method"]["runtime_seconds"], args.output,
    ),
    flush=True,
)
