#!/usr/bin/env python3
"""Replay exported exact-core certificates using stacked annihilators.

This checks quotient words against the prior sealed trajectory export, then
checks exact witnesses and completeness of the common Z and Q cores. It does
not repeat the original elliptic-curve coordinate recognition.
"""
import hashlib
import json
from pathlib import Path
import sys

from sage.all import GF, QQ, ZZ, gcd, matrix, vector


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def witness(g, target, w):
    target = vector(ZZ, target)
    c = vector(QQ, w["rational_coefficients"])
    a = vector(ZZ, w["integer_coefficients_for_multiple"])
    d = ZZ(w["minimum_positive_multiple"])
    require(c*g == target and a*g == d*target and a == d*c, "word identity failed")
    require(d > 0 and gcd([d, *a]) == 1, "nonminimal denominator")
    require(w["integral"] == (d == 1), "integrality flag differs")


def annihilator_intersection(generators, field, n):
    constraints = [list(row) for g in generators for row in g.change_ring(field).right_kernel().basis_matrix()]
    return matrix(field, constraints, ncols=n).right_kernel()


def verify(folder):
    manifest = read(folder/"manifest.json")
    for name, digest in manifest["files_sha256"].items():
        require(sha(folder/name) == digest, "export hash mismatch: "+name)
    report = read(folder/"source-REPORT.json")
    require(report["status"] == "PASS_THREE_CLOSURE_EXPERIMENTS", "source report not passed")
    require(report["outputs"]["trajectories"] == sha(folder/"source-trajectories.json"), "source trajectory binding differs")
    source = read(folder/"source-trajectories.json")
    result = read(folder/"exact-core.json")
    names = source["direction_ids"]
    n = len(names)
    require(n == 14 and result["status"] == "COMPLETE_EXACT_CORE" and result["direction_ids"] == names, "wrong result")
    source_runs = {r["seed"]: r for r in source["runs"]}
    require(len(result["runs"]) == len(source_runs) == 14, "incomplete roster")
    require({r["seed"] for r in result["runs"]} == set(source_runs), "duplicate or missing seed")
    targets = {f"recovered-local-{i:02d}" for i in (1, 2, 4)}
    per_dim, terminal = {}, {}
    for run in result["runs"]:
        seed = run["seed"]
        events = [v for s in source_runs[seed]["stages"] for v in s["new"]]
        rows = [[int(i == names.index(seed)) for i in range(n)]]
        require(len(run["states"]) == len(events)+1, "missing prefix")
        for k, state in enumerate(run["states"]):
            if k:
                require(events[k-1]["denominator"] == 1, "nonintegral source")
                rows.append(events[k-1]["quotient_word"])
            g = matrix(ZZ, rows)
            require(state["generators"] == rows and state["quotient_dimension"] == k+1, "source/prefix mismatch")
            require(g.rank() == g.change_ring(GF(2)).rank() == k+1, "dependent prefix")
            smith = g.smith_form(transformation=False)
            require(all(abs(smith[i,i]) == 1 for i in range(k+1)), "prefix not saturated in displayed Z14")
            require(state["smith_nonzero_invariants"] == [1]*(k+1) and state["index_in_saturation_inside_Z14"] == 1, "index report differs")
            require(set(state["named_axes"]) == targets, "missing named axis")
            for name, info in state["named_axes"].items():
                axis = [int(i == names.index(name)) for i in range(n)]
                q_in = vector(QQ, axis) in g.change_ring(QQ).row_space()
                f_in = vector(GF(2), axis) in g.change_ring(GF(2)).row_space()
                require(info["rational"] == info["integral"] == q_in and info["mod2"] == f_in, "axis containment differs")
                if q_in:
                    witness(g, axis, info["witness"])
                else:
                    require(info["witness"] is None, "spurious axis witness")
            per_dim.setdefault(k+1, {})[seed] = g
        require(run["final_rank"] == source_runs[seed]["final_rank"] == 17+len(rows), "endpoint differs")
        expected_first = {name: {level: next((s["quotient_dimension"] for s in run["states"] if s["named_axes"][name][level]), None)
                                 for level in ("mod2", "rational", "integral")} for name in targets}
        require(run["first_containment_dimension"] == expected_first, "first-containment times differ")
        terminal[seed] = g
    require(len(result["common_cores_by_dimension"]) == len(per_dim), "missing dimension")
    for block, entries in [*( (b, per_dim[b["quotient_dimension"]]) for b in result["common_cores_by_dimension"]),
                           (result["common_terminal_core_all_runs"], terminal)]:
        require(block["runs_present"] == len(entries) and set(block["seeds"]) == set(entries), "core roster differs")
        qspace = annihilator_intersection(list(entries.values()), QQ, n)
        fspace = annihilator_intersection(list(entries.values()), GF(2), n)
        q = matrix(QQ, block["common_rational_basis_primitive_integer_rows"], ncols=n)
        z = matrix(ZZ, block["common_integral_basis_hnf"], ncols=n)
        require(q.row_space() == z.change_ring(QQ).row_space() == qspace, "common rational space differs")
        require(z.row_module().saturation() == z.row_module(), "common Z lattice not saturated")
        require(q.nrows() == z.nrows() == qspace.dimension() == block["common_rational_dimension"] == block["common_integral_rank"], "core dimension differs")
        parity = matrix(GF(2), [[(v >> j) & 1 for j in range(n)] for v in block["common_mod2_basis"]], ncols=n)
        require(parity.row_space() == fspace and parity.nrows() == block["common_mod2_dimension"] == fspace.dimension(), "common parity space differs")
        require(z.change_ring(GF(2)).row_space() == fspace, "common parity core does not fully lift")
        for seed, g in entries.items():
            for key, basis in (("rational_basis", q), ("integral_basis", z)):
                words = block["witnesses_in_each_run"][seed][key]
                require(len(words) == basis.nrows(), "missing core witness")
                for target, w in zip(basis, words):
                    witness(g, target, w)
                    if key == "integral_basis":
                        require(w["integral"], "nonintegral core witness")
    print("PASS_EXACT_CORE_ANNIHILATOR_AND_INTEGER_WITNESSES|runs=14|acquisitions=180|prefixes=194")


if __name__ == "__main__":
    default = Path(__file__).resolve().parents[2]/"artifacts/generated-results/elliptic-curves/curve302_exact_core_v1"
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else default)
