#!/usr/bin/env sage-python
"""Rank pointed cubic-pencil quadratic pullbacks by known-point overlap.

Declared finite space: 31 cyclic eight-point anchors; lines through each of
their nine rational basepoints and each of the remaining 23 public points.
Equal quadratic function-field squareclasses, not numeric height resemblance,
are the gate for putting several directions on one common pullback.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, PolynomialRing, QQ, matrix, vector

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
PUBLIC = ROOT / "elliptic-curves/cas/icarm_curve302.py"
HELPER = ROOT / "elkies-k3/scripts/construct_curve302_nine_direction_k3.sage"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-curve302-cubic-pencil-overlap-v1.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    public = SourceFileLoader("pencil_overlap_public", str(PUBLIC)).load_module()
    helper = SourceFileLoader("pencil_overlap_helper", str(HELPER)).load_module()
    target = EllipticCurve(QQ, [QQ(str(c)) for c in public.GENERAL_WEIERSTRASS_COEFFICIENTS])
    c4, c6 = target.c_invariants()
    a, b = -27*c4, -54*c6
    curve = EllipticCurve(QQ, [a, b])
    points = [curve(36*QQ(str(px))+15, 108*(2*QQ(str(py))+QQ(str(px))+1)) for px, py in public.POINTS]
    ring = PolynomialRing(QQ, names=("X", "Y", "Z"))
    x, y, z = ring.gens()
    f0 = y*y*z-x**3-a*x*z*z-b*z**3
    monomials = [x**i*y**j*z**(3-i-j) for i in range(4) for j in range(4-i)]
    original = vector(QQ, [f0.monomial_coefficient(m) for m in monomials])
    sr = PolynomialRing(QQ, "s")
    s = sr.gen()
    ur = PolynomialRing(QQ, "u")
    u = ur.gen()
    all_pencils = []
    all_scores = Counter()
    unavailable = Counter()
    for start in range(31):
        indices = [(start+i) % 31 for i in range(8)]
        anchors = [points[i] for i in indices]
        evaluation = matrix(QQ, [[m(*p) for m in monomials] for p in anchors])
        if evaluation.rank() != 8:
            raise ArithmeticError("anchor evaluation rank is not eight")
        possibilities = [helper.primitive(sum(c*m for c, m in zip(v, monomials)))
                         for v in evaluation.right_kernel().basis()
                         if matrix(QQ, [original, v]).rank() == 2]
        f1 = min(possibilities, key=lambda f: (max(abs(c.numerator()).nbits() for c in f.coefficients()), str(f)))
        residual = -sum(anchors, curve(0))
        basepoints = anchors+[residual]
        if len(set(basepoints)) != 9:
            raise ArithmeticError("basepoint collision")
        for p in basepoints:
            if f0(*p) or f1(*p):
                raise ArithmeticError("basepoint not on both cubics")
            if matrix(QQ, [[f.derivative(v)(*p) for v in ring.gens()] for f in (f0, f1)]).rank() != 2:
                raise ArithmeticError("nontransverse base locus")
        buckets = {}
        for zero_index, zero in enumerate(basepoints):
            if not zero[2]:
                raise ArithmeticError("affine line chart unavailable")
            for target_index, q in enumerate(points):
                if target_index in indices:
                    continue
                line = [zero[0]+s*(q[0]-zero[0]), zero[1]+s*(q[1]-zero[1]), sr(1)]
                g0, r0 = sr(f0(*line)).quo_rem(s)
                g1, r1 = sr(f1(*line)).quo_rem(s)
                if r0 or r1 or max(g0.degree(), g1.degree()) > 2:
                    raise ArithmeticError("line restriction is not residual quadratic")
                if g0.gcd(g1).degree() != 0:
                    unavailable["line_contains_another_basepoint"] += 1
                    continue
                aa, bb, cc = (ur(g0[i])+u*ur(g1[i]) for i in (2, 1, 0))
                disc = bb**2-4*aa*cc
                if disc.degree() != 2 or disc.gcd(disc.derivative()).degree() != 0 or not disc(0):
                    unavailable["degenerate_or_ramified_anchor_cover"] += 1
                    continue
                # Primitive monic branch polynomial groups identical branch
                # divisors. The remaining scalar must be a rational square.
                normalized = disc/disc.leading_coefficient()
                key = tuple(map(str, normalized.list()))
                groups = buckets.setdefault(key, [])
                record = {"basepoint_index": zero_index, "public_point_index_one_based": target_index+1,
                          "discriminant_scalar": str(disc.leading_coefficient()),
                          "residual_F0_coefficients_low_to_high": helper.coefficients(g0),
                          "residual_F1_coefficients_low_to_high": helper.coefficients(g1)}
                for group in groups:
                    ratio = disc.leading_coefficient()/QQ(group[0]["discriminant_scalar"])
                    if ratio.is_square():
                        group.append(record)
                        break
                else:
                    groups.append([record])
        groups_out = []
        for branch, groups in sorted(buckets.items()):
            for group in groups:
                extra = sorted({r["public_point_index_one_based"] for r in group})
                # Modulo the eight anchor coordinates, each distinct new
                # public index is an independent direction in M31.
                rank = 8+len(extra)
                if rank > 17:
                    raise ArithmeticError("apparent K3 rank exceeds arithmetic ceiling; do not promote")
                all_scores[rank] += 1
                groups_out.append({"certified_point_span_rank": rank,
                                   "branch_polynomial_monic_low_to_high": list(branch),
                                   "extra_public_indices_one_based": extra, "covers": group})
        groups_out.sort(key=lambda g: (-g["certified_point_span_rank"], g["extra_public_indices_one_based"],
                                      g["branch_polynomial_monic_low_to_high"]))
        best = groups_out[0]
        # Preserve all grouping data in a compact checkpoint; the summary
        # certificate retains every branch fingerprint and member label.
        checkpoint = ROOT / f"artifacts/local/elkies-k3/curve302-cubic-pencil-overlap-v1/pencil-{start:02d}.json"
        checkpoint.parent.mkdir(parents=True, exist_ok=True)
        checkpoint.write_text(json.dumps({"anchors_one_based": [i+1 for i in indices],
                                         "F1": helper.polynomial_record(f1), "groups": groups_out}, sort_keys=True)+"\n")
        all_pencils.append({"anchors_one_based": [i+1 for i in indices],
                            "cover_class_count": len(groups_out),
                            "best_certified_point_span_rank": best["certified_point_span_rank"],
                            "best_extra_public_indices_one_based": best["extra_public_indices_one_based"],
                            "classes": [{"branch_sha256": sha256(json.dumps(g["branch_polynomial_monic_low_to_high"]).encode()).hexdigest(),
                                         "extra_public_indices_one_based": g["extra_public_indices_one_based"],
                                         "basepoint_indices": [r["basepoint_index"] for r in g["covers"]],
                                         "certified_point_span_rank": g["certified_point_span_rank"]} for g in groups_out],
                            "checkpoint": str(checkpoint.relative_to(ROOT)), "checkpoint_sha256": helper.digest(checkpoint)})
        print(f"PENCILOVERLAP302|pencil={start+1}/31|covers={len(groups_out)}|best={best['certified_point_span_rank']}", flush=True)
    result = {"schema": "elkies-k3.curve302-cubic-pencil-overlap.v1",
              "status": "EXACT_BOUNDED_POINT_DIRECTED_COVER_GROUPING",
              "target_curve": 302, "anchor_pencil_count": 31,
              "attempted_line_covers": 31*9*23,
              "known_direction_count": 31, "generic_rank_filter": None,
              "overlap_score_histogram": dict(sorted(all_scores.items())),
              "maximum_displayed_point_span_rank": max(all_scores),
              "unavailable": dict(unavailable), "pencils": all_pencils,
              "ranking_scope": "Within each anchored rational cubic pencil, quadratic covers are deduplicated by exact QQ(u) squareclass. Each class lifts the eight-dimensional anchor span and every listed extra public direction. Full generic MW groups and cross-pencil fibration equivalence remain unknown.",
              "proof_boundary": "The grouping is exact for the declared 6417 line covers, not for all multisections or all point combinations. Nine-direction K3 geometry and an actual t-family are separately certified in elkies-k3-curve302-nine-direction-k3-v1.json. Other groups are cover-incidence candidates until their geometry and maps are materialized; their count is not a count of distinct K3 fibrations.",
              "inputs": {str(p.relative_to(ROOT)): helper.digest(p) for p in (Path(__file__), PUBLIC, HELPER)},
              "reproducing_command": "sage -python elkies-k3/scripts/search_curve302_cubic_pencil_overlap.sage"}
    rendered = json.dumps(result, indent=2, sort_keys=True)+"\n"
    if args.check:
        if args.output.read_text() != rendered:
            raise ArithmeticError("cover grouping replay changed")
    else:
        args.output.write_text(rendered)
    print(f"PENCILOVERLAP302|maximum={max(all_scores)}|classes={sum(all_scores.values())}|output={args.output}", flush=True)


if __name__ == "__main__":
    main()
