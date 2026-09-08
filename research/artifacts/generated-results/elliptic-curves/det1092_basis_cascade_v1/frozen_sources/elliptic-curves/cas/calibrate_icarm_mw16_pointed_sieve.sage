#!/usr/bin/env sage-python
"""Complement-blind initial-wave controls for the specialized slope sieve.

One representative of each of the five supplied parent fibrations. Recompute
the full generic census and specialized representatives, then classify every
discovery with the existing exact relation/finite-reduction engine. No public
exceptional points, target ranks, or old search points enter the search.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0,str(CAS))
import half_lattice_pointed_sieve as sieve

INPUT = ROOT / "elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_pointed_sieve_controls_v1.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=OUTPUT)
    parser.add_argument("--height-bound",type=int,default=10000)
    parser.add_argument("--timeout-seconds",type=float,default=2)
    args = parser.parse_args()
    if not 1 <= args.height_bound <= 1000000 or not 0 < args.timeout_seconds <= 60:
        raise SystemExit("invalid bounded search budget")
    args.stack_bytes = 1000000000
    args.relation_chunk_size = 64
    args.relation_timeout_seconds = 180
    ladder = SourceFileLoader("pointed_control_ladder",str(CAS / "run_icarm_mw16_parent_ladder_blind.sage")).load_module()
    legacy = SourceFileLoader("pointed_control_legacy",str(CAS / "run_curve385_iterated_half_lattice_search.sage")).load_module()
    legacy.GENERIC_DIMENSION = 16
    legacy.engine = sieve.CheckpointedBackend(ROOT / "artifacts/local/elliptic-curves/mw16-pointed-sieve-control-charts")
    fixture = json.loads(INPUT.read_text())
    if fixture["status"] != "PASS_EXACT_COMPLEMENT_BLIND_NINE_PARENT_INPUTS":
        raise ArithmeticError("control fixture is not passing")
    selected = {}
    for row in sorted(fixture["parents"],key=lambda p:int(p["priority_rank"])):
        selected.setdefault(int(row["curve_id"]),row)
    paths = [Path(__file__),INPUT,CAS / "run_icarm_mw16_parent_ladder_blind.sage",
             CAS / "run_curve385_iterated_half_lattice_search.sage"]
    payload = {
        "schema":"elliptic-curves.icarm-mw16-pointed-sieve-controls.v1",
        "status":"SEARCHING", "backend":sieve.BACKEND_NAME,
        "inputs":{**sieve.provenance(),**{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths}},
        "declared_budget":{"height_bound":args.height_bound,"search_seconds_each_chart":args.timeout_seconds,
                           "adaptive_lifts":0,"height_coordinate":"recorded pointed/Gauss slope chart"},
        "software":subprocess.check_output([str(sieve.compiled_worker()),"--version"],text=True).strip(),
        "blindness_boundary":fixture["blindness_boundary"],
        "results":[],
        "claim_boundary":["This validates the initial search path, not the historical adaptive 54/55 recovery.",
                          "Changing slope coordinates changes finite height boxes.",
                          "No rank upper bound or saturation in the full Mordell-Weil group is asserted."],
    }
    if args.output.exists():
        cached = json.loads(args.output.read_text())
        for key in ("schema","inputs","declared_budget"):
            if payload[key] != cached[key]:
                raise ArithmeticError("control checkpoint configuration changed")
        payload = cached
    for curve_id,row in sorted(selected.items()):
        if any(r["curve_id"] == curve_id for r in payload["results"]):
            continue
        result = ladder.run_parent(row,legacy,args)
        payload["results"].append(result)
        ladder.write_payload(args.output,payload)
        print(f"POINTEDCONTROL|curve={curve_id}|gain={result['exact_quotient_rank_recovered']}",flush=True)
    complete = all(c["search"]["status"] == "bounded_search_complete"
                   for r in payload["results"] for c in r["cover_records"])
    payload["status"] = "PASS_COMPLETE_INITIAL_POINTED_SIEVE_CONTROLS" if complete else "INCOMPLETE_POINTED_SIEVE_CONTROLS"
    payload["reproducing_command"] = f"sage -python {Path(__file__).relative_to(ROOT)} --height-bound {args.height_bound} --timeout-seconds {args.timeout_seconds}"
    ladder.write_payload(args.output,payload)


if __name__ == "__main__":
    main()
