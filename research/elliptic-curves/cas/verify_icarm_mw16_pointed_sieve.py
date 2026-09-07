#!/usr/bin/env python3
"""Package/check the bounded MW16 sieve run; optionally replay every QQ map.

No local checkpoints are needed for --check: the full raw evidence is retained
in a deterministic gzip certificate. --replay-charts reconstructs all chart
base points and verifies the five-coefficient substitution identities and hits.
This checker does not turn a bounded miss into a statement about E(Q).
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as Q
import gzip
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import statistics
import subprocess

import half_lattice_pointed_sieve as sieve
from alternate_quartic_covers import alternate_cover
from mod2_reduction_independence import (
    combined_mod2_rank, find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime, mod2_reduction_signature,
)

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
INPUT = ART / "icarm_mw16_nagao_finalist_specializations_h300_v1.json"
CONTROLS = ART / "icarm_mw16_pointed_sieve_controls_v1.json"
CERTIFICATE = ART / "icarm_mw16_pointed_sieve_h10000_v1.json.gz"
SUMMARY = ART / "icarm_mw16_pointed_sieve_h10000_summary_v1.json"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def check_sources(payload):
    for name, expected in payload["inputs"].items():
        if digest(ROOT / name) != expected:
            raise ArithmeticError(f"changed certificate source: {name}")


def verify_chart(model, points, cover, replay):
    record = cover["search"]
    if record["backend"] != sieve.BACKEND_NAME:
        raise ArithmeticError("wrong quartic backend")
    if record["hyperellminimalmodel_called"] or record["hyperellred_called"]:
        raise ArithmeticError("generic hyperelliptic preprocessing was used")
    if record["status"] != "bounded_search_complete":
        raise ArithmeticError("an incomplete box cannot enter the complete certificate")
    h = record["height_bound"]
    if (record["denominator_start"],record["denominator_end"],record["completed_denominator"]) != (1,h,h):
        raise ArithmeticError("denominator coverage mismatch")
    if record["integer_pairs_covered"] != h*(2*h+1) or not record["infinity_checked"]:
        raise ArithmeticError("projective box coverage mismatch")
    if not 0 <= record["exact_square_tests"] <= record["all_prime_survivors"] <= record["word_sieve_survivors"] <= h*(2*h+1):
        raise ArithmeticError("sieve stage census mismatch")
    if not replay:
        return
    base = sieve.linear_combination(model,points,cover["specialized_representative"])
    if record["base_point"] != {"x":str(base[0]),"y":str(base[1])}:
        raise ArithmeticError("chart base point is not its recorded MW16 combination")
    chart = sieve.make_chart(model,base)
    if chart.record() != record["chart"]:
        raise ArithmeticError("regenerated pointed chart changed")
    # Check the entire rational polynomial, not just j or I,J.
    a,b,c,d = chart.matrix
    den,k,u = chart.denominator,chart.shift,chart.curve_scale
    raw_matrix = (Q(den*den*a+k*c,den)/u,Q(den*den*b+k*d,den)/u,c,d)
    transported = tuple(v*u**4/den**2 for v in sieve.binary_transform(
        alternate_cover(model,base).coefficients,raw_matrix))
    if transported != chart.coefficients:
        raise ArithmeticError("exact five-coefficient horizontal substitution failed")
    mapped = set()
    for n,den,r in (tuple(map(int,p)) for p in record["primitive_square_hits"]):
        if r*r != sum(f*n**i*den**(4-i) for i,f in enumerate(chart.coefficients)):
            raise ArithmeticError("recorded square hit left its binary quartic")
        for root in {r,-r}:
            point = chart.map_point(n,den,root)
            if point is not None:
                mapped.add(point)
    expected = {(Q(p["x"]),Q(p["y"])) for p in record["finite_curve_points"]}
    if expected != mapped:
        raise ArithmeticError("recorded hit transport census changed")


def verify_control_group(model, points, result):
    """Independent group-law and finite-reduction replay of the initial gains."""
    no_two_torsion_prime = find_two_torsion_certificate_prime(model,prime_bound=1000)
    classification = result["discovered_group_saturation"]
    basis = list(points)
    for event in classification["events"]:
        if event["type"] != "NEW_Q_INDEPENDENT_DIRECTION":
            raise ArithmeticError("this retained control checker requires direction-only basis growth")
        basis.append((Q(event["point"]["x"]),Q(event["point"]["y"])))
        if event["basis_rank_after"] != len(basis):
            raise ArithmeticError("control basis rank changed")
    if len(basis) != 16+result["exact_quotient_rank_recovered"]:
        raise ArithmeticError("control gain differs from its explicit basis")
    if classification["events"]:
        certificate = classification["events"][-1]["finite_reduction_certificate"]
        signatures = []
        for saved in certificate["signatures"]:
            signature = mod2_reduction_signature(model,basis,saved["prime"])
            actual = {"prime":signature.prime,"group_order":signature.group_order,
                      "doubled_subgroup_order":signature.doubled_subgroup_order,
                      "quotient_dimension":signature.quotient_dimension,
                      "rows":[list(r) for r in signature.rows]}
            if actual != saved:
                raise ArithmeticError("control finite-reduction signature failed replay")
            signatures.append(signature)
    else:
        signatures = find_mod2_reduction_certificate(model,basis,prime_bound=1000)
    if combined_mod2_rank(signatures,len(basis)) != len(basis):
        raise ArithmeticError("control directions are not certified independent")
    for relation in classification["exact_integral_relations"]:
        point = (Q(relation["point"]["x"]),Q(relation["point"]["y"]))
        if sieve.linear_combination(model,basis,relation["coordinates"]) != point:
            raise ArithmeticError("control integral relation failed exact group law")
    return no_two_torsion_prime


def summarize(raw, controls, replay=False):
    check_sources(raw)
    check_sources(controls)
    inputs = json.loads(INPUT.read_text())
    candidates = inputs["candidates"]
    if raw["status"] != "PASS_COMPLETE_FROZEN_NAGAO_FINALIST_HALF_LATTICE_GATE":
        raise ArithmeticError("full prospective gate is not complete")
    if [c["candidate_id"] for c in candidates] != [r["candidate_id"] for r in raw["results"]]:
        raise ArithmeticError("frozen prospective candidate order changed")
    if controls["status"] != "PASS_COMPLETE_INITIAL_POINTED_SIEVE_CONTROLS":
        raise ArithmeticError("initial control run did not complete")
    records, rows = [], []
    for candidate,result in zip(candidates,raw["results"]):
        model = tuple(map(Q,candidate["raw_short_model"]))
        points = tuple((Q(p["x"]),Q(p["y"])) for p in candidate["raw_generic_points"])
        if result["generic_mod2_independence_rank"] != 16 or result["generic_half_lattice"]["complete_class_count"] != 65536:
            raise ArithmeticError("generic independence/census prerequisite changed")
        masks = result["generic_half_lattice"]["deepest_masks"]
        if sorted(masks) != sorted(c["mask"] for c in result["cover_records"]) or len(set(masks)) != len(masks):
            raise ArithmeticError("maximum-depth stratum coverage changed")
        for cover in result["cover_records"]:
            verify_chart(model,points,cover,replay)
            if cover["search"]["height_bound"] != raw["declared_budget"]["height_bound_each_quartic"]:
                raise ArithmeticError("chart height differs from the declared campaign bound")
            records.append(cover["search"])
        rows.append({"candidate_id":candidate["candidate_id"],"chart_count":len(masks),
                     "completed_chart_count":len(masks),
                     "exact_quotient_rank_recovered":result["exact_quotient_rank_recovered"]})
    if len(candidates) != 104 or len(records) != 856:
        raise ArithmeticError("frozen 104-fibre/856-chart census changed")
    control_fixture = json.loads((ROOT / "elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json").read_text())
    parents = {p["parent_id"]:p for p in control_fixture["parents"]}
    control_rows = []
    for result in controls["results"]:
        parent = parents[result["parent_id"]]
        model = tuple(map(Q,parent["target_short_model"]))
        points = tuple((Q(p["x"]),Q(p["y"])) for p in parent["specialized_generic_points"])
        for cover in result["cover_records"]:
            verify_chart(model,points,cover,replay)
        if result["discovered_group_saturation"]["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
            raise ArithmeticError("control points are not classified exactly")
        no_two_torsion_prime = verify_control_group(model,points,result)
        control_rows.append({"curve_id":result["curve_id"],"chart_count":len(result["cover_records"]),
                             "exact_quotient_rank_recovered":result["exact_quotient_rank_recovered"],
                             "no_rational_two_torsion_prime":no_two_torsion_prime})
    bits = [r["chart"]["maximum_coefficient_bits"] for r in records]
    times = [r["worker_seconds"] for r in records]
    walls = [r["wall_seconds"] for r in records]
    return {
        "schema":"elliptic-curves.icarm-mw16-pointed-sieve-summary.v1",
        "status":"PASS_COMPLETE_BOUNDED_POINTED_SIEVE",
        "backend":sieve.BACKEND_NAME,
        "candidate_count":len(candidates),"chart_count":len(records),
        "completed_chart_count":len(records),"timeout_count":0,
        "height_bound":raw["declared_budget"]["height_bound_each_quartic"],
        "search_budget_seconds_each_chart":raw["declared_budget"]["timeout_seconds_each_quartic"],
        "integer_pairs_covered":sum(r["integer_pairs_covered"] for r in records),
        "exact_square_tests":sum(r["exact_square_tests"] for r in records),
        "finite_curve_points_reported":sum(len(r["finite_curve_points"]) for r in records),
        "positive_candidate_count":raw["positive_candidate_count"],
        "coefficient_bits":{"minimum":min(bits),"median":statistics.median(bits),"maximum":max(bits)},
        "worker_seconds":{"sum":sum(times),"median":statistics.median(times),"maximum":max(times)},
        "chart_wall_seconds_including_group_and_transform":{
            "sum":sum(walls),"median":statistics.median(walls),"maximum":max(walls)},
        "software":raw["software"],
        "initial_controls":control_rows,
        "initial_control_directions":sum(r["exact_quotient_rank_recovered"] for r in control_rows),
        "prospective_candidates":rows,
        "inputs":{**sieve.provenance(),**{str(p.relative_to(ROOT)):digest(p) for p in
                  (Path(__file__).resolve(),Path(__file__).with_name("mod2_reduction_independence.py"),
                   INPUT,CONTROLS,CERTIFICATE)}},
        "claim_boundary":[
            "All 856 boxes complete in the recorded exact slope coordinates; these are different boxes from PARI-reduced or raw-coordinate searches.",
            "The initial control run is not a replay of the historical adaptive 54/55 recovery.",
            "The first prospective timeout campaign is retained as historical, wholly censored evidence.",
            "A bounded miss gives no rank upper bound, point absence, saturation, or Selmer information.",
            "The complete search still stops before adaptive lifts, unrestricted searches, and residual Selmer calculations."],
        "reproducing_command":"python3 elliptic-curves/cas/verify_icarm_mw16_pointed_sieve.py --check --replay-charts",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw",type=Path)
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--replay-charts",action="store_true")
    args = parser.parse_args()
    if args.check:
        expected = json.loads(SUMMARY.read_text())
        check_sources(expected)
        raw = json.loads(gzip.decompress(CERTIFICATE.read_bytes()))
    else:
        if args.raw is None:
            parser.error("--raw is required when packaging the certificate")
        raw = json.loads(args.raw.read_text())
        raw["software"]["sage_version_at_packaging"] = subprocess.check_output(
            ["sage","--version"],text=True).strip()
        version = subprocess.run(["gp","--version"],text=True,capture_output=True,check=True)
        raw["software"]["pari_version_at_packaging"] = (version.stdout+version.stderr).strip()
        CERTIFICATE.write_bytes(gzip.compress((json.dumps(raw,indent=2,sort_keys=True)+'\n').encode(),mtime=0))
    summary = summarize(raw,json.loads(CONTROLS.read_text()),args.replay_charts)
    if args.check:
        if summary != expected:
            raise ArithmeticError("summary does not match the retained raw evidence")
    else:
        SUMMARY.write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(f"POINTEDSIEVE|charts={summary['chart_count']}|complete={summary['completed_chart_count']}|"
          f"positive={summary['positive_candidate_count']}|QQ_replay={args.replay_charts}")


if __name__ == "__main__":
    main()
