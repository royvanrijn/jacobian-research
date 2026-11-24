#!/usr/bin/env python3
"""Bounded exact controls for the DZ-coordinate fibre rigidity theorem.

The unrestricted stable-field and fibre-exhaustion arguments are written
proofs. This checker does not infer them from the displayed examples.
"""

from __future__ import annotations

if not __debug__:
    raise RuntimeError("Run this exact checker without -O or -OO.")

import argparse
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import comb
from pathlib import Path
import json
import platform
import sys


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def load_arithmetic(repo):
    path = repo / "research/foundations/extended-geometry/verify_odd_start_fixed_coordinate_fibres.py"
    spec = spec_from_file_location("dz_fibre_exact_arithmetic", path)
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def dyck_words(edges):
    """Each word is one actual ordered rooted plane tree, by edge traversal."""
    def visit(prefix, opened, closed):
        if closed == edges:
            yield prefix
        else:
            if opened < edges:
                yield from visit(prefix + "(", opened + 1, closed)
            if closed < opened:
                yield from visit(prefix + ")", opened, closed + 1)
    yield from visit("", 0, 0)


def black_degrees(word):
    # The white root is what remains beside the removed marked black leaf.
    children = [0]
    result = []
    for token in word:
        if token == "(":
            children[-1] += 1
            children.append(0)
        else:
            depth = len(children) - 1
            count = children.pop()
            if depth % 2:
                result.append(count + 1)
    assert len(children) == 1
    return result


def coefficient_truncated_power(a, r, m):
    p = [1] + [0] * m
    for _ in range(r):
        q = [0] * (m + 1)
        for i, c in enumerate(p):
            for j in range(min(a, m - i) + 1):
                q[i + j] += c
        p = q
    return p[m]


def tree_controls():
    rows = []
    for a in (1, 2, 3, 4, 5):
        weight_choices = {}
        for degree in range(1, a + 1):
            # Enumerate positive edge weights, rather than substituting a
            # binomial count into the rooted-tree recursion.
            weight_choices[degree] = sum(
                sum(weights) <= a
                for weights in product(range(1, a + 1), repeat=degree))
        for r in range(1, 8):
            ordinary, weighted, tree_count = 0, 0, 0
            for word in dyck_words(r - 1):
                tree_count += 1
                degrees = black_degrees(word)
                if any(d > a for d in degrees):
                    continue
                ordinary += 1
                ways = 1
                for d in degrees:
                    ways *= weight_choices[d]
                weighted += ways
            assert tree_count == comb(2 * r - 2, r - 1) // r
            ordinary_formula = coefficient_truncated_power(a, r, r - 1) // r
            assert coefficient_truncated_power(a, r, r - 1) % r == 0
            inclusion = sum((-1) ** j * comb(r, j) *
                            comb(2 * r - 2 - (a + 1) * j, r - 1)
                            for j in range((r - 1) // (a + 1) + 1))
            assert inclusion % r == 0 and inclusion // r == ordinary
            assert ordinary == ordinary_formula
            assert weighted == comb(a * r, r - 1) // r
            assert weighted >= ordinary
            rows.append({"a": a, "r": r, "ordinary_plane_trees": ordinary,
                         "weighted_plane_trees": weighted,
                         "singleton_coordinate_classes_at_odd_height":
                         weighted - ordinary if a % 2 else None})
    a3 = [row for row in rows if row["a"] == 3 and row["r"] <= 6]
    assert [row["ordinary_plane_trees"] for row in a3] == [1, 1, 2, 5, 13, 36]
    assert [row["weighted_plane_trees"] for row in a3] == [1, 3, 12, 55, 273, 1428]
    degree_counts = []
    for a in (1, 3, 5):
        lookup = {row["r"]: row for row in rows if row["a"] == a}
        for s in range(1, 8):
            rigid = lookup[s]["weighted_plane_trees"] - lookup[s]["ordinary_plane_trees"]
            explicit = rigid + sum(lookup[r]["ordinary_plane_trees"]
                                   for r in range(1, s + 1)
                                   for k in range(s)
                                   if r * (k + 1) == s)
            formula = lookup[s]["weighted_plane_trees"] + sum(
                lookup[r]["ordinary_plane_trees"]
                for r in range(1, s) if s % r == 0)
            assert explicit == formula
            if a == 1:
                assert formula == sum(s % r == 0 for r in range(1, s + 1))
            degree_counts.append({"a": a, "s": s,
                                  "degree": (a + 2) * (s + 1),
                                  "phase_source_classes": formula})
    return rows, degree_counts


def algebra_controls(t):
    x, one = t.RX, t.RONE
    counts = {k: 0 for k in (
        "pure_pairs", "raw_moments", "dz_identities", "height_rejections",
        "even_height_nonconstant_endpoint_rejections",
        "weighted_deformation_rejections", "outside_class_pairs")}
    negatives = []
    samples = []

    def dz_data(u, r, N):
        p = u / (x * u.diff())
        assert p.polynomial() and len(p.as_poly()) - 1 == r and p.at(0) == 1
        assert len(t.gcd_poly(p.as_poly(), t.derivative(p.as_poly()))) == 1
        assert len(u.n) - 1 == len(u.d) - 1 == N
        aa = t.scale(u.n, 1 / u.n[-1])
        bb = t.scale(u.d, 1 / u.d[-1])
        assert aa[0] == 0 and aa[1] != 0 and bb[0] != 0
        assert len(t.sub(aa, bb)) - 1 == N - r
        assert u.diff() / u == 1 / (x * p)
        counts["dz_identities"] += 1
        return p

    def check_pair(u, eta, a, degree, outside=False):
        b, xx, c, d = t.construct(u, eta, a)
        assert c.polynomial() and d.polynomial()
        cp, dp = c.as_poly(), d.as_poly()
        assert cp[0] == 1 and dp[0] == -1
        assert x * xx ** a / c == u and -d / c == xx ** 2
        assert max(len(cp) + a - 2, len(dp) + a) == degree
        for m in range(1, 7):
            assert t.raw_moment(cp, dp, a, m) == 0
            counts["raw_moments"] += 1
        if outside:
            assert not b.polynomial()
            counts["outside_class_pairs"] += 1
        else:
            assert b.polynomial()
        counts["pure_pairs"] += 1
        return cp, dp

    # Safe coordinate includes the a=0 boundary.
    for a in range(6):
        check_pair(x, one, a, a + 1)

    # Entire Shabat star fibres, over both signs of source scale.
    for lam in (Q(1), Q(-2)):
        for r in range(1, 5):
            rr = 2 * (1 + lam * x) ** r - 1
            u = (2 / rr.diff().at(0)) * (rr - 1) / (rr + 1)
            p = dz_data(u, r, r)
            assert p == (rr ** 2 - 1) / (2 * x * rr.diff())
            for a in (1, 3, 5):
                for k in range(4):
                    check_pair(u, rr ** k, a, (a + 2) * (r * (k + 1) + 1))
            for a in (2, 4):
                check_pair(u, one, a, (a + 2) * (r + 1))
                _, _, c, d = t.construct(u, rr, a)
                assert c.polynomial() and d.polynomial()
                # The coefficient of s in Phi is -1/c, so even-height
                # purity already fails at the FIRST actual moment.
                c_infty = u.n[-1] / u.d[-1]
                assert t.raw_moment(c.as_poly(), d.as_poly(), a, 1) == \
                    -t.kappa(a, 1) / c_infty != 0
                counts["even_height_nonconstant_endpoint_rejections"] += 1
            if r == 1:
                check_pair(u, one, 0, 4)
            else:
                _, _, c, d = t.construct(u, one, 0)
                assert not c.polynomial() or not d.polynomial()
                counts["height_rejections"] += 1

    # The non-star Chebyshev tree: an unmarked double zero requires a>=2.
    rr = 4 * (1 + x) ** 3 - 3 * (1 + x)
    u = (2 / rr.diff().at(0)) * (rr - 1) / (rr + 1)
    dz_data(u, 3, 3)
    assert (rr - 1).at(Q(-3, 2)) == 0
    assert rr.diff().at(Q(-3, 2)) == 0
    assert rr.diff().diff().at(Q(-3, 2)) != 0
    for a in (0, 1):
        for k in range(4):
            _, _, c, d = t.construct(u, rr ** k, a)
            assert not c.polynomial() or not d.polynomial()
            counts["height_rejections"] += 1
    for a in (2, 4):
        check_pair(u, one, a, (a + 2) * 4)
    for a in (3, 5):
        for k in range(4):
            check_pair(u, rr ** k, a, (a + 2) * (3 * (k + 1) + 1))
    negatives.append("nonconstant eta cannot repair the unmarked double-zero height cap")

    # A weighted DZ coordinate with N>r: the proposed nonconstant
    # H(u) remains reciprocal but introduces an actual forbidden pole.
    uu, vv = 1 + x, 1 + 2 * x / 3
    u = x * vv ** 2 / uu ** 3
    p = dz_data(u, 2, 3)
    assert p == uu * vv
    c_infty = Q(4, 9)
    rr_bad = (c_infty + u) / (c_infty - u)
    assert not rr_bad.polynomial()
    assert t.evaluate(rr_bad.d, Q(-4, 3)) == 0
    for a in range(6):
        _, _, c, d = t.construct(u, one, a)
        if a < 2:
            assert not c.polynomial() or not d.polynomial()
            counts["height_rejections"] += 1
        else:
            assert c == vv ** (a - 2) * uu ** (a + 3)
            assert d == -(vv ** a * uu ** (a + 5))
            check_pair(u, one, a, 3 * (a + 2))
    for a in (3, 5):
        for k in (1, 2, 3):
            _, _, c, d = t.construct(u, rr_bad ** k, a)
            assert not c.polynomial() and not d.polynomial()
            counts["weighted_deformation_rejections"] += 1
    negatives.append("reciprocal deformation of weighted DZ coordinate has a finite u=c pole")
    negatives.append("even height requires constant Phi, not only even-part cancellation")

    # A deliberately non-Shabat polynomial does not pass condition (S1).
    rr = 1 + x + x ** 2
    u = 2 * (rr - 1) / (rr + 1)
    assert not (u / (x * u.diff())).polynomial()
    assert rr.diff().at(Q(-1, 2)) == 0
    assert rr.at(Q(-1, 2)) not in (1, -1)
    negatives.append("critical value outside {minus1,plus1} fails polynomial logarithmic quotient")

    # Root's degree-three reciprocal quotient: finite critical points +/-i.
    z = 1 + x
    w = (z - 1) / (z + 1)
    u = Q(2, 3) * (w ** 3 + 3 * w)
    assert u == Q(8, 3) * (z ** 3 - 1) / (z + 1) ** 3
    assert u.diff().at(0) == 1
    b = u / (x * u.diff())
    assert b == (z ** 2 + z + 1) * (z + 1) / (3 * (z ** 2 + 1))
    assert not b.polynomial()
    for a in (1, 3, 5):
        for k in (2, 6, 10):
            check_pair(u, z ** k, a, (a + 2) * (k + 2), outside=True)
    for k in (0, 1, 3, 4, 5):
        _, _, c, d = t.construct(u, z ** k, 3)
        assert not c.polynomial() or not d.polynomial()
    negatives.append("degree-three non-Phi1 coordinate invalidates unrestricted Shabat rigidity")

    # The separate degree-two quotient has primitive sixth-root critical points.
    u = -(z ** 2 - 1) / (z ** 2 - 4 * z + 1)
    assert u.diff().at(0) == 1
    b = u / (x * u.diff())
    assert not b.polynomial()
    modulus = (z ** 2 - z + 1).as_poly()
    assert not t.divmod_poly((z ** 3 + 1).n, modulus)[1]
    for a in (1, 3, 5):
        for k in (3, 9):
            check_pair(u, z ** k, a, (a + 2) * (k + 2), outside=True)
    negatives.append("degree-two non-Phi1 coordinate also has non-Shabat infinite deformations")

    # A corrupted member remains normalized and polynomial, but is not pure.
    u = x / (1 + x)
    eta = (1 + 2 * x) ** 2 + x
    _, _, c, d = t.construct(u, eta, 3)
    assert c.polynomial() and d.polynomial()
    values = [t.raw_moment(c.as_poly(), d.as_poly(), 3, m) for m in (2, 4, 6)]
    assert any(values)
    samples.append({"corrupted_pair_even_moments": list(map(str, values))})
    negatives.append("corrupted normalized polynomial pair fails raw moments")
    return counts, negatives, samples


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    directory = Path("research/foundations/extended-geometry")
    proof = directory / "DAVENPORT_ZANNIER_COORDINATE_FIBRES_AND_SHABAT_RIGIDITY.md"
    checker = directory / Path(__file__).name
    ledger_path = directory / "davenport_zannier_coordinate_fibres_source_ledger.json"
    ledger = json.loads((repo / ledger_path).read_text())
    for entry in ledger["sources"]:
        assert digest(repo / entry["path"]) == entry["sha256"], entry["path"]
    arithmetic = load_arithmetic(repo)
    counts, negatives, samples = algebra_controls(arithmetic)
    trees, degree_counts = tree_controls()
    result = {
        "schema": "davenport_zannier_coordinate_fibres_verification_v1",
        "status": "PASS", "python": platform.python_version(),
        "arithmetic": "exact Fraction and integer arithmetic; standard library only",
        "counts": counts, "negative_controls": negatives,
        "tree_enumeration_rows": trees, "degree_count_rows": degree_counts,
        "samples": samples, "source_pins_checked": len(ledger["sources"]),
        "packet_pins": {str(p): digest(repo / p) for p in (proof, checker, ledger_path)},
        "scope": "Bounded exact algebraic, moment and plane-tree controls. "
                 "The written proof supplies stable-field rigidity and exhaustive "
                 "all-height empty/singleton/infinite classification. "
                 "No classification outside polynomial u/(xu') is claimed.",
        "assurance": "Internal author replay; no external human review or formal verification."
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "counts": counts,
                      "tree_count_cases": len(trees),
                      "degree_count_cases": len(degree_counts),
                      "negative_controls": len(negatives),
                      "source_pins_checked": len(ledger["sources"]),
                      "receipt": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
