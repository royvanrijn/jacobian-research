import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as F

os.environ["CURVE302_TEST_SKIP_FSYNC"] = "1"  # sandbox CI mounts may reject fsync; production stays strict

HERE = Path(__file__).resolve()
CAS = HERE.parents[1] / "cas"
sys.path.insert(0, str(CAS))

import curve302_chart_exposure as core
import run_curve302_chart_exposure as runner
from curve302_short_core_controls import Vocabulary


NAMES = runner.NAMES


def dump(path, payload):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n")


def sha(path):
    return runner.sha(path)


def candidate_vectors():
    n = 14
    rows = []
    # All coordinate axes first (norm 1), sign-canonical and lexicographically sorted.
    for i in range(n):
        v = [0]*n; v[i] = 1; rows.append((1, tuple(v)))
    # Add simple primitive pairs at norm 2.
    for i in range(10):
        v = [0]*n; v[i] = 1; v[(i+1) % n] = 1; rows.append((2, tuple(v)))
    rows.sort(key=lambda x: (F(x[0]), x[1]))
    return rows


def make_fixture(root, *, point_only=False):
    root = Path(root)
    short = root/"short"; structure = root/"structure"; ledger = root/"ledger.json"
    short.mkdir(parents=True); structure.mkdir(parents=True)
    directions = candidate_vectors()
    tsv = short/"primitive-directions.tsv"
    with tsv.open("w") as out:
        out.write("norm\tvector\n")
        for norm, v in directions:
            out.write(f"{norm}\t{','.join(map(str,v))}\n")
    enum = {"status":"PASS_COMPLETE_EXACT_ENUMERATION", "enumeration_sha256":sha(tsv),
            "direction_count":len(directions), "bound":"2"}
    dump(short/"enumeration.json", enum)
    cores = []
    axes = []
    for d in range(1,15):
        v = [0]*14; v[d-1] = 1; axes.append(v)
        cores.append({"quotient_dimension":d, "rank":d, "run_count":14 if d <= 12 else 13,
                      "basis":[list(x) for x in axes]})
    filtration = {"status":"PASS_INTRINSIC_AND_OBSERVED_FILTRATIONS", "observed_common_integral_cores":cores}
    dump(short/"filtration.json", filtration)
    dump(short/"ranks-basins.json", {"status":"PASS_COMPLETE_RANK_AND_BASIN_CENSUS"})

    form = [[int(i==j) for j in range(14)] for i in range(14)]
    rel = {"status":"PASS_QUOTIENT_RELATION_ANALYSIS", "direction_ids":list(NAMES), "schur_quotient":form}
    dump(structure/"quotient-relations.json", rel)
    runs = []; ledger_runs = []
    for seed in range(14):
        others = [j for j in range(14) if j != seed]
        if seed == 0: others = others[:-2]
        stages = []; lstages = []
        for epoch, j in enumerate(others):
            v = [0]*14; v[j] = 1
            stages.append({"epoch":epoch, "new":[{"integral":True, "denominator":1,
                                                   "quotient_word":v, "primitive_quotient_word":v}]})
            chart = {"chart_id":f"{seed}-{epoch}", "order":0, "score_band":"A", "complete":True,
                     "quartic_coefficients":[1, seed+1, epoch+1, 2, 3], "search_bound":100}
            if point_only:
                chart["points"] = [{"x":"1/2", "y":"3/4", "parameter_height":5}]
            else:
                # actual hit plus one alternative short axis, if unknown
                alt = (j+1) % 14
                if alt == seed: alt = (alt+1) % 14
                av = [0]*14; av[alt] = 1
                chart["exposures"] = [
                    {"quotient_word":v, "parameter_height":3, "within_bound":True},
                    {"quotient_word":av, "parameter_height":9, "within_bound":True},
                ]
            lstages.append({"epoch":epoch, "charts":[chart]})
        runs.append({"seed":NAMES[seed], "final_rank":18+len(others), "charts":len(others), "stages":stages})
        ledger_runs.append({"seed":NAMES[seed], "stages":lstages})
    traj = {"status":"PASS_ALL_14_SEEDED_TRAJECTORIES_RECONCILED", "direction_ids":list(NAMES), "runs":runs}
    dump(structure/"trajectories.json", traj)
    dump(structure/"REPORT.json", {"status":"PASS_THREE_CLOSURE_EXPERIMENTS",
                                    "outputs":{"quotient-relations":sha(structure/"quotient-relations.json"),
                                               "trajectories":sha(structure/"trajectories.json")}})
    dump(short/"plan.json", {"source":str(structure), "source_hashes":{}})
    dump(short/"REPORT.json", {"status":"PASS_THREE_SHORT_VECTOR_CORE_EXPERIMENTS",
                                "enumeration_sha256":sha(tsv),
                                "outputs":{"enumeration":sha(short/"enumeration.json"),
                                           "filtration":sha(short/"filtration.json"),
                                           "ranks-basins":sha(short/"ranks-basins.json")}})
    payload = {"schema":"curve302-chart-exposure-ledger.v1", "direction_ids":list(NAMES), "runs":ledger_runs}
    dump(ledger, payload)
    return short, structure, ledger


class CoreTests(unittest.TestCase):
    def test_nominal_word_is_not_exposure(self):
        v = [1]+[0]*13
        chart = core.normalize_chart({"chart_id":"x", "quotient_word":v}, 0, 14)
        self.assertEqual(chart["exposures"], [])

    def test_explicit_result_word_is_exposure(self):
        v = [1]+[0]*13
        chart = core.normalize_chart({"chart_id":"x", "exposures":[{"quotient_word":v, "parameter":"3/7"}]}, 0, 14)
        self.assertEqual(len(chart["exposures"]), 1)
        self.assertEqual(chart["exposures"][0]["parameter_height"], 7)

    def test_false_within_bound_is_excluded(self):
        v = [1]+[0]*13
        chart = core.normalize_chart({"exposures":[{"quotient_word":v, "parameter_height":11, "within_bound":False}]}, 0, 14)
        self.assertFalse(chart["exposures"])

    def test_point_only_is_diagnostic_not_exposure(self):
        chart = core.normalize_chart({"points":[{"x":"1/2","y":"2/3"}]}, 0, 14)
        self.assertEqual(chart["unsupported_point_only_records"], 1)
        self.assertEqual(chart["exposures"], [])

    def test_contradictory_word_aliases_fail(self):
        with self.assertRaises(ValueError):
            core.explicit_word({"quotient_word":[1]+[0]*13, "candidate_word":[0,1]+[0]*12}, 14)

    def test_integer_containment(self):
        self.assertTrue(core.integer_contains([(1,0),(0,1)], (2,3)))
        self.assertFalse(core.integer_contains([(2,0),(0,1)], (1,0)))

    def test_candidate_population_adds_actual_outside_limit(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"v.tsv"; rows=candidate_vectors()
            with p.open("w") as out:
                out.write("norm\tvector\n")
                for norm,v in rows: out.write(f"{norm}\t{','.join(map(str,v))}\n")
            actual=rows[-1][1]
            vocab=Vocabulary.load(p,14,len(rows),F(2),[actual])
            runs=[{"events":[{"primitive":actual}]}]
            vectors,pos=core.build_candidate_population(vocab,runs,static_limit=5)
            self.assertIn(actual,vectors); self.assertEqual(pos[actual],len(rows)-1)

    def test_counterfactual_reverse_can_change_gain(self):
        axes=[tuple(int(i==j) for i in range(14)) for j in range(14)]
        stage={"seed":"s","epoch":0,"actual":{"word":axes[1]}, "coverage_complete":True,
               "chart_rows":[{"chart_id":"a","original_order":0,"score_band":"x","quartic_coefficient_bits":10,"complete_flag":True,"exposed":[axes[1]]},
                             {"chart_id":"b","original_order":1,"score_band":"x","quartic_coefficient_bits":5,"complete_flag":True,"exposed":[axes[2]]}]}
        result=core.counterfactual_stage(stage,[axes[0]],{2:(axes[1],)},NAMES,{axes[1]:0,axes[2]:1},random_orders=4)
        by={x["policy"]:x for x in result["policies"]}
        self.assertEqual(by["original"]["matches_actual"],1)
        self.assertEqual(by["reverse"]["matches_actual"],0)


    def test_incomplete_stage_counterfactual_is_unknown(self):
        axes=[tuple(int(i==j) for i in range(14)) for j in range(14)]
        stage={"seed":"s","epoch":0,"actual":{"word":axes[1]}, "coverage_complete":False,
               "chart_rows":[{"chart_id":"a","original_order":0,"score_band":"x","quartic_coefficient_bits":10,"complete_flag":None,"exposed":[axes[1]]}]}
        result=core.counterfactual_stage(stage,[axes[0]],{},NAMES,{axes[1]:0},random_orders=2)
        self.assertTrue(all(p["status"] == "UNKNOWN_INCOMPLETE_CHART_COVERAGE" for p in result["policies"]))

    def test_tied_percentile(self):
        self.assertEqual(core.tied_percentile(2,[1,2,3],higher_better=True), str(F(1,2)))


class EndToEndTests(unittest.TestCase):
    def test_run_and_check(self):
        with tempfile.TemporaryDirectory() as td:
            td=Path(td); short, structure, ledger=make_fixture(td)
            out=td/"out"
            runner.prepare(out, short, structure, ledger_path=ledger, static_limit=20, random_orders=8,
                           stage_seconds=120, memory_bytes=2*1024**3)
            runner.resume(out)
            self.assertEqual(runner.read(out/"REPORT.json")["status"],"PASS_THREE_CHART_EXPOSURE_EXPERIMENTS")
            census=runner.read(out/"exposure-census.json")
            self.assertEqual(census["summary"]["stages"],180)
            self.assertEqual(census["summary"]["actual_exposed"],180)
            runner.check(out)

    def test_tampered_ledger_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            td=Path(td); short, structure, ledger=make_fixture(td); out=td/"out"
            runner.prepare(out, short, structure, ledger_path=ledger, static_limit=20, random_orders=2,
                           stage_seconds=120, memory_bytes=2*1024**3)
            p=out/"inputs/chart-exposure-ledger.json"; obj=runner.read(p); obj["runs"][0]["stages"][0]["charts"][0]["chart_id"]="tampered"; dump(p,obj)
            with self.assertRaises(ValueError): runner.guard(out)

    def test_point_only_prepare_fails_closed_with_probe(self):
        with tempfile.TemporaryDirectory() as td:
            td=Path(td); short, structure, ledger=make_fixture(td, point_only=True); out=td/"out"
            with self.assertRaises(ValueError):
                runner.prepare(out, short, structure, ledger_path=ledger, static_limit=20, random_orders=2,
                               stage_seconds=120, memory_bytes=2*1024**3)
            self.assertTrue((out/"PREPARE_FAILED.json").is_file())
            self.assertTrue((out/"SCHEMA_PROBE.json").is_file())


if __name__ == "__main__": unittest.main()
