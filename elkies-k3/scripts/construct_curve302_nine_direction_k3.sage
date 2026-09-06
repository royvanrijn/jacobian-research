#!/usr/bin/env sage-python
"""Construct a pointed K3 family through 302 with nine certified directions.

This deliberately constructs a baseline, not provenance recovery. Eight points
define a cubic pencil; a line through its first basepoint and a ninth public
point defines a quadratic base change and a ninth section.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, gcd, lcm, matrix, vector
from sage.schemes.elliptic_curves.jacobian import Jacobian
from sage.env import SAGE_VERSION

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
PUBLIC = ROOT / "elliptic-curves/cas/icarm_curve302.py"
INDEPENDENCE = ROOT / "elliptic-curves/cas/verify_icarm_curve302_rank31.py"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-curve302-nine-direction-k3-v1.json"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def polynomial_record(polynomial):
    return [{"exponents": list(map(int, powers)), "coefficient": str(c)}
            for powers, c in sorted(polynomial.dict().items())]


def primitive(poly):
    result = poly*lcm([c.denominator() for c in poly.coefficients()])
    result /= gcd([ZZ(c) for c in result.coefficients()])
    if result.leading_coefficient() < 0:
        result = -result
    return poly.parent()(result)


def coefficients(poly):
    return list(map(str, poly.list()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    public = SourceFileLoader("curve302_nine_public", str(PUBLIC)).load_module()
    target = EllipticCurve(QQ, [QQ(str(c)) for c in public.GENERAL_WEIERSTRASS_COEFFICIENTS])
    points = [target(QQ(str(x)), QQ(str(y))) for x, y in public.POINTS]
    sys.path.insert(0, str(PUBLIC.parent))
    rank_helper = SourceFileLoader("curve302_nine_independence", str(INDEPENDENCE)).load_module()
    character_primes, character_rows, character_rank = rank_helper.quadratic_character_certificate(prime_bound=400)
    if character_rank != 31 or matrix(GF(2), character_rows).rank() != 31:
        raise ArithmeticError("public-point independence did not replay")
    c4, c6 = target.c_invariants()
    a, b = -27*c4, -54*c6
    short = EllipticCurve(QQ, [a, b])
    # Pinned public-to-short isomorphism, checked for all 31 input points.
    pp = [short(36*p[0]+15, 108*(2*p[1]+p[0]+1)) for p in points]
    ring = PolynomialRing(QQ, names=("X", "Y", "Z"))
    x, y, z = ring.gens()
    f0 = y*y*z-x**3-a*x*z*z-b*z**3
    monomials = [x**i*y**j*z**(3-i-j) for i in range(4) for j in range(4-i)]
    evaluation = matrix(QQ, [[m(*p) for m in monomials] for p in pp[:8]])
    if evaluation.rank() != 8:
        raise ArithmeticError("eight selected points do not impose eight conditions")
    kernel = evaluation.right_kernel().basis()
    f0_vector = vector(QQ, [f0.monomial_coefficient(m) for m in monomials])
    candidates = []
    for v in kernel:
        if matrix(QQ, [f0_vector, v]).rank() == 2:
            candidates.append(primitive(sum(c*m for c, m in zip(v, monomials))))
    f1 = min(candidates, key=lambda f: (max(abs(ZZ(c)).nbits() for c in f.coefficients()), str(f)))
    residual = -sum(pp[:8], short(0))
    basepoints = pp[:8]+[residual]
    if len(set(basepoints)) != 9:
        raise ArithmeticError("basepoints are not distinct")
    for p in basepoints:
        if f0(*p) or f1(*p):
            raise ArithmeticError("basepoint identity failed")
        gradients = matrix(QQ, [[f.derivative(v)(*p) for v in ring.gens()] for f in (f0, f1)])
        if gradients.rank() != 2:
            raise ArithmeticError("base locus is not transverse")
    print("NINEDIR302|base_locus=nine_distinct_transverse_rational_points", flush=True)

    ur = PolynomialRing(QQ, "u")
    u = ur.gen()
    cr = PolynomialRing(ur.fraction_field(), names=("X", "Y", "Z"))
    pencil = cr(f0)+u*cr(f1)
    jacobian = Jacobian(pencil)
    ja, jb = ur(jacobian.a4()), ur(jacobian.a6())
    if any(jacobian.a_invariants()[i] for i in (0, 1, 2)):
        raise ArithmeticError("Jacobian is not short")
    delta = -16*(4*ja**3+27*jb**2)
    if ja.degree() > 4 or jb.degree() > 6 or delta.degree() != 12:
        raise ArithmeticError("rational elliptic surface degree profile changed")
    if delta.gcd(delta.derivative()).degree() != 0 or ja.gcd(delta).degree() != 0:
        raise ArithmeticError("source pencil does not have twelve I1 fibres")
    specialized = EllipticCurve(QQ, [ja(0), jb(0)])
    if not specialized.is_isomorphic(short):
        raise ArithmeticError("source Jacobian does not specialize to the target over Q")
    print("NINEDIR302|source=12I1|source_geometric_rank=8", flush=True)

    tr = PolynomialRing(QQ, "t")
    t = tr.gen()
    line = [pp[0][0]+t*(pp[8][0]-pp[0][0]), pp[0][1]+t*(pp[8][1]-pp[0][1]), tr(1)]
    restriction0, restriction1 = tr(f0(*line)), tr(f1(*line))
    base_map = tr.fraction_field()(-restriction0/restriction1)
    numerator, denominator = base_map.numerator(), base_map.denominator()
    if max(numerator.degree(), denominator.degree()) != 2 or base_map(1) != 0:
        raise ArithmeticError("line does not give a degree-two cover through t=1")
    if denominator(1) == 0:
        raise ArithmeticError("target lift is not in the affine cover chart")
    # Homogenized pullback in the usual integral Weierstrass gauge.
    pull_a = sum(tr(ja[i])*numerator**i*denominator**(4-i) for i in range(5))
    pull_b = sum(tr(jb[i])*numerator**i*denominator**(6-i) for i in range(7))
    pull_delta = -16*(4*pull_a**3+27*pull_b**2)
    if pull_a.degree() > 8 or pull_b.degree() > 12 or pull_delta.degree() != 24:
        raise ArithmeticError("pullback does not have the expected K3 degree profile")
    if pull_delta.gcd(pull_delta.derivative()).degree() != 0 or pull_a.gcd(pull_delta).degree() != 0:
        raise ArithmeticError("pullback is not 24I1; branch or minimality gate needs analysis")
    # The affine moving point lies on the pulled-back cubic identically.
    if denominator*restriction0+numerator*restriction1:
        raise ArithmeticError("tautological ninth section identity failed")
    fibre = EllipticCurve(QQ, [pull_a(1), pull_b(1)])
    if not fibre.is_isomorphic(target):
        raise ArithmeticError("K3 fibre is not Q-isomorphic to curve 302")

    # These are specialization coordinates for the POINTED cubic with zero P1.
    rows = []
    for index in range(1, 8):
        row = [0]*31
        row[index], row[0] = 1, -1
        rows.append(row)
    rows.append([-2]+[-1]*7+[0]*23)  # residual basepoint minus P1
    row = [0]*31
    row[8], row[0] = 1, -1
    rows.append(row)
    image = matrix(ZZ, rows)
    specializations = [p-pp[0] for p in pp[1:8]]+[residual-pp[0], pp[8]-pp[0]]
    for q, row in zip(specializations, image.rows()):
        if q != sum((int(c)*p for c, p in zip(row, pp) if c), short(0)):
            raise ArithmeticError("specialized section group-law identity failed")
    if image.rank() != 9 or image[:8].rank() != 8:
        raise ArithmeticError("specialization overlap does not have ranks eight and nine")
    # In the rational surface, the strict transform of this line meets E1
    # once and the other eight exceptional sections zero times. Projection
    # formula gives S.O=1 and S.S_i=0 on the quadratic pullback. The constant
    # sections are pairwise disjoint and disjoint from O. With chi=2 and no
    # reducible fibres, Shioda's height formula now determines the full Gram.
    for p in basepoints[1:]:
        if matrix(QQ, [list(pp[0]), list(pp[8]), list(p)]).det() == 0:
            raise ArithmeticError("moving line meets an additional basepoint")
    height = matrix(ZZ, 9, 9, lambda i, j: (
        (4 if i == j else 2) if i < 8 and j < 8
        else (6 if i == j else 3)))
    if height.det() != 4608 or not height.is_positive_definite():
        raise ArithmeticError("nine-section height lattice failed")
    # A candidate point span of dimension nine maps injectively: any generic
    # relation would specialize to a relation among the independent P_i.
    result = {
        "schema": "elkies-k3.curve302-nine-direction-k3.v1",
        "status": "PASS_EXACT_CONSTRUCTED_K3_WITH_NINE_INDEPENDENT_SPECIALIZATIONS",
        "provenance_recovery": False,
        "target_curve": 302, "target_parameter": "1",
        "public_point_independence": {"method": "exact cubic-root quadratic characters at good primes",
                                      "prime_bound": 400, "primes": list(character_primes),
                                      "rows_mod_2": [list(r) for r in character_rows], "rank": 31},
        "source_point_indices_one_based": list(range(1, 10)),
        "pointed_cubic_family": {
            "equation": "D(t)*F0(X,Y,Z)+N(t)*F1(X,Y,Z)=0",
            "zero": list(map(str, pp[0])),
            "F0": polynomial_record(f0), "F1": polynomial_record(f1),
            "N_coefficients_low_to_high": coefficients(numerator),
            "D_coefficients_low_to_high": coefficients(denominator),
            "constant_sections": [list(map(str, p)) for p in basepoints[1:]],
            "moving_section_coordinates_low_to_high": [coefficients(p) for p in line],
            "specialization_map_to_public_model": "At t=1: (X,Y,Z) -> (X/Z,Y/Z)-P1 in E_short; then x=(X_short-15)/36, y=(Y_short/108-x-1)/2.",
        },
        "source_jacobian": {"A_coefficients_low_to_high": coefficients(ja),
                            "B_coefficients_low_to_high": coefficients(jb),
                            "singular_fibres": "12I1", "geometric_rank": 8},
        "weierstrass_family": {"equation": "y^2=x^3+A(t)*x+B(t)",
                               "A_coefficients_low_to_high": coefficients(pull_a),
                               "B_coefficients_low_to_high": coefficients(pull_b),
                               "singular_fibres": "24I1", "chi": 2,
                               "nonisotrivial": bool((pull_a**3).derivative()*pull_b**2 != pull_a**3*(pull_b**2).derivative()),
                               "fibre_to_target_isomorphism_u_r_s_t": list(map(str, fibre.isomorphism_to(target).tuple()))},
        "specialization_overlap": {"candidate": "span of the nine displayed pointed-cubic sections",
                                   "rank": 9, "matrix_9_by_31_rows": rows,
                                   "matrix_smith_factors": list(map(int, image.elementary_divisors())),
                                   "all_nine_group_law_identities": True,
                                   "full_generic_MW_rank_lower_bound": 9,
                                   "full_generic_MW_rank_upper_bound": 17,
                                   "full_generic_MW_rank": "UNKNOWN",
                                   "full_MW_overlap_with_M31": "UNKNOWN; at least 9",
                                   "specialized_remaining_quotient_rank": 22},
        "generic_section_lattice": {"height_gram": [list(map(int, r)) for r in height.rows()],
                                    "rank": 9, "determinant": int(height.det()),
                                    "saturation": "UNKNOWN",
                                    "NS_sublattice": "U plus negative of the displayed height Gram; determinant magnitude 4608",
                                    "intersection_input": {"constant_section_pairwise_intersections": 0,
                                                           "constant_section_zero_intersections": 0,
                                                           "moving_section_zero_intersection": 1,
                                                           "moving_section_constant_intersections": 0},
                                    "height_formula": "<P,Q>=chi+P.O+Q.O-P.Q for distinct sections; <P,P>=2chi+2P.O; all fibre corrections vanish"},
        "checks": {"public_points_on_curve": 31, "base_locus_transverse": True,
                   "source_discriminant_degree": int(delta.degree()),
                   "source_discriminant_squarefree": True,
                   "pullback_discriminant_degree": int(pull_delta.degree()),
                   "pullback_discriminant_squarefree": True,
                   "moving_section_polynomial_identity": True,
                   "target_Q_isomorphism": True},
        "proof_boundary": "This constructs a K3 through 302 from chosen known points. It is not the original parent, a 12-direction family, a saturated generic MW basis, an exact generic rank, or a rank upper bound for E302. The specialization matrix certifies rank-nine overlap for the displayed section sublattice; overlap for the entire generic group is only bounded below.",
        "inputs": {str(p.relative_to(ROOT)): digest(p) for p in (PUBLIC, INDEPENDENCE, Path(__file__))},
        "software": {"sage_version": SAGE_VERSION},
        "reproducing_command": "sage -python elkies-k3/scripts/construct_curve302_nine_direction_k3.sage",
    }
    rendered = json.dumps(result, sort_keys=True, indent=2)+"\n"
    if args.check:
        if args.output.read_text() != rendered:
            raise ArithmeticError("stored construction differs from exact replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(f"NINEDIR302|K3=24I1|overlap=9|full_generic_rank=UNKNOWN|output={args.output}", flush=True)


if __name__ == "__main__":
    main()
