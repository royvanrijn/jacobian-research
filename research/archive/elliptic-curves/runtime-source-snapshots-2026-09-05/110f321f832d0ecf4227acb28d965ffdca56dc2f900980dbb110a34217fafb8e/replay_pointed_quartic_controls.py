#!/usr/bin/env python3
"""Re-sieve the 1,034 frozen calibrated MW16 control boxes through the API.

Compare exact coefficients, square hits and mapped points chart by chart;
replay the saved integral relations and finite-reduction independence. This
is a backend regression on blind, already frozen centres, not a new search
policy calibration or a prospective campaign. No public complement is read.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
from hashlib import sha256
import json
from pathlib import Path

import pointed_quartic_search as pq
import verify_mw16_sensitivity as exact

ROOT = pq.ROOT
BUNDLE = ROOT/"artifacts/generated-results/elliptic-curves/mw16_sensitivity_controls_v1.json.gz"
SUMMARY = ROOT/"artifacts/generated-results/elliptic-curves/mw16_sensitivity_controls_summary_v1.json"
PREFIX = "artifacts/local/elliptic-curves/mw16-sensitivity/"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--checkpoint-dir", type=Path, default=ROOT/"artifacts/local/elliptic-curves/pointed-quartic-search/control-regression")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        raise ValueError("use one to eight workers")
    digest = sha256(BUNDLE.read_bytes()).hexdigest()
    if digest != json.loads(SUMMARY.read_text())["bundle_sha256"]:
        raise ArithmeticError("control evidence bundle differs from its pinned summary")
    bundle = json.loads(gzip.decompress(BUNDLE.read_bytes()))
    documents = bundle["files"]
    # The saved bundle carries executed sources. They need not equal the
    # amended working tree; their bytes remain independently hash checked.
    import base64
    for name, record in documents.items():
        data = base64.b64decode(record["base64"]) if "base64" in record else record["text"].encode()
        if sha256(data).hexdigest() != record["sha256"]:
            raise ArithmeticError("source/evidence snapshot corrupt: "+name)
    fixture = json.loads(documents[exact.FIXTURE]["text"])
    parents = {r["parent_id"]: r for r in fixture["parents"]}
    initial = json.loads(documents[PREFIX+"frozen-initial.json"]["text"])
    adaptive = json.loads(documents[PREFIX+"frozen-adaptive.json"]["text"])
    jobs = []
    for payload in (initial, adaptive):
        for result in payload["results"]:
            parent = parents[result["parent_id"]]
            for setting in result["settings"]:
                if setting["specification"] != "metric:16":
                    raise ArithmeticError("frozen control regression left calibrated coordinates")
                for record in setting["charts"]:
                    jobs.append((parent["target_short_model"], result["search_basis"], record,
                                 result["curve_id"], setting["centre"]))
    if len(jobs) != 1034:
        raise ArithmeticError("calibrated control box census changed")

    def replay(job):
        model, subgroup, saved, curve_id, centre_kind = job
        search = pq.PointedQuarticSearch(model, subgroup, {"coefficients": saved["representative"]}, "metric:16")
        outcome = search.search(saved["height_bound"], 20, checkpoint_dir=args.checkpoint_dir)
        if outcome.record["status"] != "bounded_search_complete":
            raise ArithmeticError("control box timed out; regression is incomplete")
        for field in ("pointed_chart", "horizontal_matrix", "ordinate_scale", "coefficients",
                      "primitive_square_hits", "finite_curve_points", "integer_pairs_covered"):
            if outcome.record[field] != saved[field]:
                raise ArithmeticError(f"control {curve_id} differs in {field}")
        return {"curve_id": curve_id, "centre": centre_kind, "mask": saved["mask"],
                "height": saved["height_bound"], "coefficient_bits": outcome.record["maximum_coefficient_bits"],
                "square_hits": len(outcome.record["primitive_square_hits"]),
                "finite_points": len(outcome.curve_points),
                "integer_pairs_covered": outcome.record["integer_pairs_covered"],
                "exact_result_sha256": sha256(pq.canonical({k: outcome.record[k] for k in
                    ("coefficients", "primitive_square_hits", "finite_curve_points")}).encode()).hexdigest()}

    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for row in pool.map(replay, jobs):
            rows.append(row)
            if len(rows) % 50 == 0 or len(rows) == len(jobs):
                print(f"POINTED_CONTROL|exact_matches={len(rows)}/{len(jobs)}", flush=True)
    # Exact quotient rank is certified only here, outside the search API.
    summaries = [exact.verify(payload, documents) for payload in (initial, adaptive)]
    gains = {r["curve_id"]: r["exact_quotient_rank_recovered"] for r in summaries[0]["results"]}
    gains.update({r["curve_id"]: r["exact_quotient_rank_recovered"] for r in summaries[1]["results"]})
    if gains != {398: 14, 400: 12, 401: 11, 542: 10, 548: 8}:
        raise ArithmeticError("exact control quotient ladder changed")
    report = {"schema": "elliptic-curves.universal-pointed-control-regression.v1",
              "status": "PASS_EXACT_BOX_AND_QUOTIENT_REPLAY", "backend": pq.BACKEND_NAME,
              "sources": {**pq.sources(), **{str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest()
                  for p in (Path(__file__).resolve(), Path(exact.__file__).resolve(),
                            Path(exact.__file__).with_name('mod2_reduction_independence.py'))}},
              "input_bundle": str(BUNDLE.relative_to(ROOT)), "input_bundle_sha256": digest,
              "chart_count": len(rows), "exact_quotient_gains": gains,
              "total_control_directions": sum(gains.values()), "charts": rows,
              "claim_boundary": "Same 55 certified control directions on the same 1034 frozen boxes. No new prospective campaign, gain, rank upper bound, or calibration claim for MW17/MW18."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    print("POINTED_CONTROL|PASS|directions=55", flush=True)


if __name__ == "__main__":
    main()
