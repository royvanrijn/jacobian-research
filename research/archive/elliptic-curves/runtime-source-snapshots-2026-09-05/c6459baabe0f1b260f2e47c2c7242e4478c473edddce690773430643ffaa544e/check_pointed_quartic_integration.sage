#!/usr/bin/env sage-python
"""Cheap executable checks of active family adapters on one MW17 chart."""
import argparse
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import sys
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT/"elliptic-curves/cas"
sys.path.insert(0, str(CAS))
import pointed_quartic_search as pq
from pointed_quartic_migration import runtime_search, require_runtime


def load(name, path):
    return SourceFileLoader(name, str(path)).load_module()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    legacy = load("pointed_integration_legacy", CAS/"run_curve385_iterated_half_lattice_search.sage")
    mw17 = load("pointed_integration_mw17", CAS/"run_mw17_jump_v2.sage")
    rescue = load("pointed_integration_rescue", CAS/"run_mw17_jump_v2_zero_gain_rescue.sage")
    sensitivity = load("pointed_integration_mw16", CAS/"run_mw16_sensitivity.sage")
    crt = load("pointed_integration_crt", ROOT/"elkies-k3/scripts/run_r17_prospective_crt_half_lattice_search.sage")
    mw17.load_campaign()
    rescue.load_protocol()
    crt.load_inputs()
    for version in ("", "_v2"):
        sparse = load("pointed_integration_sparse"+version,
                      CAS/("run_curve385_sparse_quotient_rank32_search"+version+".sage"))
        sparse.load_protocol()
    assert legacy.engine.run_quartic_search is pq.run_quartic_search
    assert crt.engine.run_quartic_search is pq.run_quartic_search
    assert sensitivity.backend is pq
    fixture = ROOT/"elliptic-curves/data/r17_refresh_jump_ladder_blind_inputs_v1.json"
    case = json.loads(fixture.read_text())["cases"][0]
    subgroup = tuple(pq.point(p) for p in case["generic_points"])
    assert len(subgroup) == 17
    kwargs = dict(mask=1, representative=[1]+[0]*16, short_model=case["short_model"],
                  generic_points=subgroup, height_bound=17, timeout_seconds=2, stack_bytes=0)
    results = []
    original_run = pq.base.subprocess.run

    def no_pari(command, *positional, **keywords):
        if Path(command[0]).name == "gp":
            raise AssertionError("active pointed search called PARI")
        return original_run(command, *positional, **keywords)

    with patch.object(pq.base.subprocess, "run", side_effect=no_pari):
        for label, search in (
            ("curve-specific/shared-lattice", lambda: legacy.engine.run_quartic_search(**kwargs)),
            ("MW17", lambda: mw17.run_quartic_search_raw(None, **kwargs)),
            ("zero-gain-rescue", lambda: rescue.load_base_module().run_quartic_search_raw(None, **kwargs)),
            ("R17-CRT", lambda: crt.engine.run_quartic_search(**kwargs)),
        ):
            outcome = search()
            assert outcome.record["backend"] == pq.BACKEND_NAME
            assert outcome.record["status"] == "bounded_search_complete"
            compact = crt.compact_cover_record("canary", 1, kwargs["representative"], outcome, 0)
            assert compact["pointed_search"] == outcome.record
            results.append({"adapter": label, "finite_points": len(outcome.curve_points),
                            "coefficient_bits": outcome.record["maximum_coefficient_bits"]})
    try:
        require_runtime({"status": "COMPLETE_CHUNK"})
    except ArithmeticError:
        pass
    else:
        raise AssertionError("historical checkpoint accepted by active campaign")
    result = {"status": "PASS_UNIVERSAL_POINTED_ADAPTERS", "runtime_search": runtime_search(),
              "MW16_sensitivity_uses_shared_module": True, "fixture_curve_id": case["curve_id"],
              "height": 17, "seconds_per_chart": 2, "adapters": results,
              "claim_boundary": "Bounded integration canary only; no MW17/MW18 sensitivity or prospective gain is claimed."}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print("POINTED_INTEGRATION|PASS", flush=True)


if __name__ == "__main__":
    main()
