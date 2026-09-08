#!/usr/bin/env sage-python
"""Replay all saved model maps, sections, finite certificates and quartic identities.

No minimization, coordinate optimization, CVP enumeration or point search is
needed for this narrow independent certificate replay.
"""
import argparse
import gzip
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))
from mod2_reduction_independence import (
    mod2_reduction_signature, gf2_rank,
    short_curve_has_no_rational_2_torsion_modular_certificate,
)


def verify_benchmark():
    art = ROOT / "artifacts/generated-results/elliptic-curves"
    path = art / "mw16_short_models_chart_benchmark_v1.json"
    benchmark = json.loads(path.read_text())
    compressed = (art / "mw16_short_models_h300_v1.json.gz").read_bytes()
    assert sha256(compressed).hexdigest() == benchmark["input_sha256"]
    data = json.loads(gzip.decompress(compressed))
    source_path = art / "a1_mw16_target_free_parameter_candidates_h300_v1.json"
    source = json.loads(source_path.read_text())
    lookup = {r["candidate_id"]: r for r in data["candidates"]}
    raw_lookup = {r["candidate_id"]: r for r in source["candidates"]}
    witnesses, unknown = [], []
    assert benchmark["height_bound"] == 10000 and benchmark["timeout_seconds"] == 2
    assert len(benchmark["results"]) == len(benchmark["frozen_charts"]) == 9
    for frozen, result in zip(benchmark["frozen_charts"], benchmark["results"]):
        assert (frozen["candidate_id"], frozen["mask"]) == (result["candidate_id"], result["mask"])
        assert len(result["trials"]) == 4
        raw = raw_lookup[result["candidate_id"]]
        row = lookup[result["candidate_id"]]
        chart = next(c for c in row["charts"] if c["mask"] == result["mask"])
        E = EllipticCurve(QQ, raw["raw_short_model"])
        points = [E(QQ(p["x"]), QQ(p["y"])) for p in raw["raw_generic_points"]]
        rep = chart["representative"]
        Q = sum((n*p for n, p in zip(rep, points)), E(0))
        known = {E(0): [0]*16, Q: rep}
        for i, p in enumerate(points):
            for sign in (-1, 1):
                vector = [sign if j == i else 0 for j in range(16)]
                known[sign*p] = vector
                known[Q-sign*p] = [a-b for a, b in zip(rep, vector)]
        for trial in result["trials"]:
            assert trial["status"] == "bounded_search_complete"
            for i, point in enumerate(trial["points_transported_to_raw"]):
                p = E(QQ(point["x"]), QQ(point["y"]))
                witness = {"candidate_id": result["candidate_id"], "mask": result["mask"],
                           "mode": trial["mode"], "point_index": i}
                if p not in known:
                    unknown.append(witness)
                    continue
                vector = known[p]
                assert p == sum((n*q for n, q in zip(vector, points)), E(0))
                witnesses.append({**witness, "generic_coefficients": vector})
    output = {"status": "PASS_ALL_BENCHMARK_POINTS_IN_MW16" if not unknown else "UNCLASSIFIED_BENCHMARK_POINTS",
        "benchmark_sha256": sha256(path.read_bytes()).hexdigest(),
        "models_sha256": sha256(compressed).hexdigest(),
        "verifier_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "exact_generic_group_witnesses": witnesses, "unclassified_points": unknown,
        "claim_boundary": "Exact containment of reported benchmark points only; no rank upper bound or statement about unsearched boxes."}
    (art / "mw16_short_models_benchmark_verification_v1.json").write_text(json.dumps(output, indent=2, sort_keys=True)+"\n")
    assert not unknown, "returned benchmark points require further exact classification"
    print(f"PASS: all {len(witnesses)} reported benchmark points have exact MW16 group-law witnesses.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-only", action="store_true")
    args = parser.parse_args()
    if args.benchmark_only:
        verify_benchmark()
        return
    art = ROOT / "artifacts/generated-results/elliptic-curves"
    summary = json.loads((art / "mw16_short_models_h300_summary_v1.json").read_text())
    compressed = (ROOT / summary["certificate"]).read_bytes()
    assert sha256(compressed).hexdigest() == summary["certificate_sha256"]
    data = json.loads(gzip.decompress(compressed))
    for path, expected in data["inputs"].items():
        assert sha256((ROOT / path).read_bytes()).hexdigest() == expected, path
    source = json.loads((art / "a1_mw16_target_free_parameter_candidates_h300_v1.json").read_text())
    previous = json.loads((art / "a1_mw16_target_free_parameter_search_h300_v1.json").read_text())
    masks = {r["candidate_id"]: r["deepest_masks"] for r in previous["results"]}
    assert [c["candidate_id"] for c in data["candidates"]] == [c["candidate_id"] for c in source["candidates"]]
    assert len(data["candidates"]) == 104
    ring = PolynomialRing(QQ, "w")
    w = ring.gen()
    chart_count = 0
    for raw, row in zip(source["candidates"], data["candidates"]):
        E = EllipticCurve(QQ, raw["raw_short_model"])
        original = [E(QQ(p["x"]), QQ(p["y"])) for p in raw["raw_generic_points"]]
        assert len(original) == 16
        for model_key, points_key, map_key in (
            ("global_minimal_model", "minimal_sections", "raw_to_minimal"),
            ("selected_short_model", "selected_short_sections", "raw_to_selected_short")):
            F = EllipticCurve(QQ, row[model_key])
            transported = [F(QQ(p["x"]), QQ(p["y"])) for p in row[points_key]]
            u, r, s, t = map(QQ, row[map_key]["u_r_s_t"])
            assert u and len(transported) == 16
            assert E.discriminant() == u**12*F.discriminant()
            assert E.j_invariant() == F.j_invariant()
            for old, new in zip(original, transported):
                assert old[0] == u*u*new[0]+r
                assert old[1] == u**3*new[1]+s*u*u*new[0]+t
        # F and transported now denote the selected short model.
        model = tuple(Fraction(str(a)) for a in F.ainvs())
        points = tuple(tuple(Fraction(str(v)) for v in p.xy()) for p in transported)
        all_rows = []
        for sig in row["independence"]["signatures"]:
            actual = mod2_reduction_signature(model, points, sig["prime"])
            assert [list(v) for v in actual.rows] == sig["rows"]
            assert actual.group_order == sig["group_order"]
            assert actual.doubled_subgroup_order == sig["doubled_subgroup_order"]
            all_rows += sig["rows"]
        assert gf2_rank(all_rows, 16) == 16
        assert short_curve_has_no_rational_2_torsion_modular_certificate(
            model, row["independence"]["no_rational_two_torsion_prime"])
        assert [c["mask"] for c in row["charts"]] == masks[row["candidate_id"]]
        gram = [[QQ(v) for v in line] for line in raw["generic_height_gram"]]
        for chart in row["charts"]:
            rep = chart["representative"]
            assert len(rep) == 16
            assert sum((v % 2) << i for i, v in enumerate(rep)) == chart["mask"]
            assert 2*sum(rep[i]*gram[i][j]*rep[j] for i in range(16) for j in range(16)) == 23
            p = sum((n*q for n, q in zip(rep, transported)), F(0))
            assert (str(p[0]), str(p[1])) == (chart["base_point"]["x"], chart["base_point"]["y"])
            u, r, s, t = map(QQ, row["raw_to_selected_short"]["u_r_s_t"])
            rawp = E(u*u*p[0]+r, u**3*p[1]+s*u*u*p[0]+t)
            assert rawp == sum((n*q for n, q in zip(rep, original)), E(0))
            selected = chart["selected"]
            g = ring(selected["integral_coefficients_ascending"])
            assert all(v.denominator() == 1 for v in g)
            a, b, c, d = map(QQ, selected["matrix_a_b_c_d"])
            assert a*d-b*c
            N, D = a*w+b, c*w+d
            ftrans = N**4-6*p[0]*N*N*D*D-8*p[1]*N*D**3-(3*p[0]**2+4*F.a4())*D**4
            assert g == QQ(selected["ordinate_scale"])**2*ftrans
            heights = [max(abs(g[i].numerator()).nbits(), g[i].denominator().nbits()) for i in range(5)]
            assert max(heights) == selected["maximum_bits"]
            assert sum(heights) == selected["total_bits"]
            assert (max(heights), sum(heights), selected["name"]) == min(
                (trial["maximum_bits"], trial["total_bits"], trial["name"]) for trial in selected["trials"])
            chart_count += 1
        j = E.j_invariant()
        lower = max(1, (abs(j.numerator()).nbits()-20)//6, (j.denominator().nbits()-11)//6)
        assert lower == row["any_integral_binary_quartic_coefficient_bit_lower_bound"]
    assert chart_count == summary["chart_count"] == 856
    print("PASS: 104 model maps; 1664 sections on both models; 104 fresh rank-16 certificates; 856 exact quartic maps.")
    verify_benchmark()


if __name__ == "__main__":
    main()
