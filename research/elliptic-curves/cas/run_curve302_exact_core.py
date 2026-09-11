#!/usr/bin/env python3
"""Exact parity, rational and integral cores of the sealed Curve302 trajectories.

Use Sage Python. One bounded stage, frozen source/data, exact witnesses, and
deterministic recomputation. No heights, point searches or new EC identities.
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from sage.all import QQ, ZZ, gcd, lcm, matrix, vector
import sage.version

HERE = Path(__file__).resolve().parent
SUPPORT = HERE / "support.py"
if not SUPPORT.exists():
    SUPPORT = HERE / "run_curve302_closure_followup.py"
spec = importlib.util.spec_from_file_location("exact_core_support", SUPPORT)
lab = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lab)

SCHEMA = "curve302-exact-core.v1"
BOUNDARY = (
    "Exact linear algebra on previously EC-verified integer quotient words in displayed D/M17. "
    "Parity, rational span and integral generation are distinct. Certificates express identities "
    "modulo M17, not new EC proofs or statements about saturation in the full Mordell-Weil group. "
    "At each dimension only runs attaining that dimension are compared; no missing state is imputed."
)


def rows_json(m):
    return [[int(x) for x in row] for row in m.rows()]


def word_witness(generators, target):
    """Unique rational coefficients and the least positive integral multiple."""
    g = generators.change_ring(QQ)
    target = vector(QQ, target)
    if target not in g.row_space():
        return None
    word = g.solve_left(target)
    lab.require(word * g == target, "rational word identity failed")
    multiple = ZZ(lcm([x.denominator() for x in word]))
    integers = vector(ZZ, [multiple*x for x in word])
    lab.require(integers * generators == multiple*target, "integral multiple identity failed")
    lab.require(gcd([multiple, *integers]) == 1, "multiple witness was not minimal")
    return {"rational_coefficients": [str(x) for x in word],
            "minimum_positive_multiple": int(multiple),
            "integer_coefficients_for_multiple": list(map(int, integers)),
            "integral": multiple == 1}


def primitive_rational_basis(space):
    answer = []
    for row in space.basis_matrix().rows():
        denominator = lcm([x.denominator() for x in row])
        answer.append(list(lab.primitive([int(denominator*x) for x in row])))
    return answer


def common_core(entries, n):
    """Full Z-lattice intersection plus Q intersection and per-run witnesses."""
    rational = entries[0][1].change_ring(QQ).row_space()
    integral = entries[0][1].row_module()
    parity_rows = lab.rref([lab.parity(list(map(int, row))) for row in entries[0][1]], n)
    for _, g in entries[1:]:
        rational = rational.intersection(g.change_ring(QQ).row_space())
        integral = integral.intersection(g.row_module())
        parity_rows = lab.intersection(parity_rows, [lab.parity(list(map(int, row))) for row in g], n)
    z_basis = integral.basis_matrix()
    lab.require(z_basis.change_ring(QQ).row_space() == rational,
                "integral intersection has the wrong rational span")
    rational_rows = primitive_rational_basis(rational)
    integral_rows = rows_json(z_basis)
    witnesses = {}
    for seed, g in entries:
        q_words = [word_witness(g, row) for row in rational_rows]
        z_words = [word_witness(g, row) for row in integral_rows]
        lab.require(all(w is not None for w in q_words), "common Q basis missing in run")
        lab.require(all(w is not None and w["integral"] for w in z_words), "common Z basis missing in run")
        witnesses[seed] = {"rational_basis": q_words, "integral_basis": z_words}
    return {"runs_present": len(entries), "seeds": [s for s, _ in entries],
            "common_mod2_dimension": len(parity_rows), "common_mod2_basis": parity_rows,
            "common_rational_dimension": int(rational.dimension()),
            "common_rational_basis_primitive_integer_rows": rational_rows,
            "common_integral_rank": int(integral.rank()), "common_integral_basis_hnf": integral_rows,
            "witnesses_in_each_run": witnesses}


def analyze(data):
    names = data["names"]
    n = len(names)
    target_indices = [names.index(f"recovered-local-{i:02d}") for i in (2, 1, 4)]
    by_dimension, terminals, runs = {}, [], []
    for run in data["runs"]:
        rows = [[int(i == run["seed_index"]) for i in range(n)]]
        states = []
        for step in range(len(run["events"])+1):
            if step:
                rows.append(list(run["events"][step-1]["word"]))
            g = matrix(ZZ, rows)
            lab.require(g.rank() == len(rows), "dependent exact quotient generators")
            binary = [lab.parity(row) for row in rows]
            lab.require(len(lab.rref(binary, n)) == len(rows), "dependent parity generators")
            smith = g.smith_form(transformation=False)
            invariants = [abs(int(smith[i,i])) for i in range(len(rows))]
            index = 1
            for value in invariants:
                index *= value
            lab.require(index % 2 == 1, "even saturation index contradicts parity rank")
            axes = {}
            for j in target_indices:
                axis = [int(i == j) for i in range(n)]
                w = word_witness(g, axis)
                in_mod2 = lab.in_span(1 << j, binary, n)
                # Odd lattice index guarantees Q-contained integral axes have
                # odd clearing multiples, hence also parity containment.
                if w is not None:
                    lab.require(w["minimum_positive_multiple"] % 2 == 1 and in_mod2,
                                "axis containment disagrees with odd index")
                axes[names[j]] = {"mod2": in_mod2, "rational": w is not None,
                                  "integral": w is not None and w["integral"], "witness": w}
            state = {"step": step, "quotient_dimension": len(rows), "generators": rows_json(g),
                     "smith_nonzero_invariants": invariants,
                     "index_in_saturation_inside_Z14": index, "named_axes": axes}
            states.append(state)
            by_dimension.setdefault(len(rows), []).append((run["seed"], g))
        terminals.append((run["seed"], g))
        first = {name: {level: next((s["quotient_dimension"] for s in states if s["named_axes"][name][level]), None)
                        for level in ("mod2", "rational", "integral")}
                 for name in (names[j] for j in target_indices)}
        runs.append({"seed": run["seed"], "final_rank": run["final_rank"],
                     "first_containment_dimension": first, "states": states})
    cores = [{"quotient_dimension": r, **common_core(entries, n)}
             for r, entries in sorted(by_dimension.items())]
    return {"schema": SCHEMA, "status": "COMPLETE_EXACT_CORE", "direction_ids": names,
            "source_acquisitions": sum(len(r["events"]) for r in data["runs"]),
            "runs": runs, "common_cores_by_dimension": cores,
            "common_terminal_core_all_runs": common_core(terminals, n), "boundary": BOUNDARY}


def software():
    return {**lab.runtime(), "sage": sage.version.version}


def guard(folder):
    plan = lab.read(folder/"plan.json")
    lab.require(plan["schema"] == SCHEMA, "wrong schema")
    lab.require(lab.digest(folder/"plan.json") == lab.read(folder/"manifest.json")["plan_sha256"], "plan changed")
    lab.require(plan["software"] == software(), "runtime changed")
    for name, digest in plan["frozen_sha256"].items():
        lab.require(lab.digest(folder/name) == digest, "frozen input/source changed: "+name)
    return plan


def stage(folder, phase, output):
    plan = guard(folder)
    phase.mkdir()
    lab.atomic(phase/"started.json", {"plan_sha256": lab.digest(folder/"plan.json")})
    receipt = lab.supervise([sys.executable, str(folder/"source/runner.py"), "_worker",
                            "--folder", str(folder), "--output", str(output)],
                           folder, phase, plan["seconds"], plan["memory_bytes"])
    guard(folder)
    lab.require(receipt["outcome"] == "COMPLETED", f"stage failed/censored; preserve {phase}")
    lab.require(lab.read(output)["status"] == "COMPLETE_EXACT_CORE", "missing complete output")
    lab.atomic(phase/"seal.json", {"output_sha256": lab.digest(output),
                                   "receipt_sha256": lab.digest(phase/"receipt.json")})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("run", "check", "status", "_worker"))
    root = HERE.parents[2]
    parser.add_argument("--folder", type=Path, default=root/"research/artifacts/local/elliptic-curves/curve302-exact-core-v1")
    parser.add_argument("--inputs", type=Path, default=root/"research/artifacts/local/elliptic-curves/curve302-closure-followup-v1-rref-compat/inputs")
    parser.add_argument("--output", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    folder = args.folder.resolve()
    if args.command == "status":
        print(lab.read(folder/"REPORT.json") if (folder/"REPORT.json").exists() else
              {"status": "NOT_PREPARED" if not folder.exists() else "INCOMPLETE_INSPECT_RECEIPTS"})
        return
    if args.command == "_worker":
        guard(folder)
        output = args.output.resolve()
        lab.require(output == folder/"exact-core.json" or output.is_relative_to(folder/"checks"), "invalid output location")
        lab.require(not output.exists(), "refuse to overwrite evidence")
        result = analyze(lab.validate_results(folder/"inputs"))
        guard(folder)
        lab.atomic(output, result)
        return
    if args.command == "run":
        lab.require(not folder.exists(), "folder exists; inspect, never rerun in place")
        # Read-only preflight before creating any output.
        lab.validate_results(args.inputs)
        folder.mkdir(parents=True)
        with lab.lock(folder):
            (folder/"inputs").mkdir()
            (folder/"source").mkdir()
            sources = {f"inputs/{n}": args.inputs/n for n in
                       ("landscape.json", "closure-laws.json", "quotient-relations.json", "trajectories.json", "REPORT.json")}
            sources.update({"source/runner.py": Path(__file__).resolve(), "source/support.py": SUPPORT})
            frozen = {}
            for name, source in sources.items():
                before = lab.digest(source)
                shutil.copyfile(source, folder/name)
                lab.require(before == lab.digest(source) == lab.digest(folder/name), "source changed during snapshot")
                frozen[name] = before
            plan = {"schema": SCHEMA, "software": software(), "frozen_sha256": frozen,
                    "origins": {n: str(p.resolve()) for n, p in sources.items()},
                    "seconds": 1800, "memory_bytes": 3*1024**3,
                    "source_commit_observed": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
                    "boundary": BOUNDARY}
            lab.atomic(folder/"plan.json", plan)
            lab.atomic(folder/"manifest.json", {"plan_sha256": lab.digest(folder/"plan.json")})
            try:
                stage(folder, folder/"phase", folder/"exact-core.json")
            except BaseException:
                lab.atomic(folder/"REPORT.json", {"status": "UNKNOWN_FAILED_OR_CENSORED"})
                raise
            lab.atomic(folder/"REPORT.json", {"status": "COMPLETE_EXACT_CORE",
                                               "output_sha256": lab.digest(folder/"exact-core.json"), "boundary": BOUNDARY})
    else:
        with lab.lock(folder):
            guard(folder)
            report = lab.read(folder/"REPORT.json")
            seal = lab.read(folder/"phase/seal.json")
            receipt = lab.read(folder/"phase/receipt.json")
            lab.require(report["status"] == "COMPLETE_EXACT_CORE" and receipt["outcome"] == "COMPLETED", "incomplete run")
            lab.require(seal["receipt_sha256"] == lab.digest(folder/"phase/receipt.json"), "receipt changed")
            lab.require(seal["output_sha256"] == report["output_sha256"] == lab.digest(folder/"exact-core.json"), "output changed")
            (folder/"checks").mkdir(exist_ok=True)
            check = Path(tempfile.mkdtemp(prefix="attempt-", dir=folder/"checks"))
            stage(folder, check/"phase", check/"exact-core.json")
            lab.require(lab.read(check/"exact-core.json") == lab.read(folder/"exact-core.json"), "recomputation mismatch")
            lab.atomic(check/"CHECK.json", {"status": "PASS_DETERMINISTIC_RECOMPUTATION", "report_sha256": lab.digest(folder/"REPORT.json")})
    print("EXACT_CORE_"+args.command.upper()+"_PASS", flush=True)


if __name__ == "__main__":
    main()
