#!/usr/bin/env python3
"""Exact selected-relation replay for the cubic cube-zero one-generator algebra."""
if not __debug__:
    raise SystemExit("Exact certificate checks require ordinary Python execution.")

from argparse import ArgumentParser
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import gzip
import json
import sys
import time

P = 1000000009
DIMENSIONS = [1, 1, 1, 1, 2, 3, 4, 6, 9, 13, 19, 26, 34, 40, 35, 16, 13, 9, 2]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    folder = args.repo_root / "research/jacobian/subtree-shuffling"
    prefix_path = folder / "cubic_cube_zero_shuffle_certificates.json.gz"
    grade29_path = folder / "cubic_cube_zero_grade29_certificate.json.gz"
    ledger_path = folder / "cubic_cube_zero_shuffle_source_ledger.json"
    ledger = json.loads(ledger_path.read_text())
    for entry in ledger["inputs"]:
        base = args.repo_root if entry["scope"] == "repository" else folder
        path = base / entry["path"]
        require(sha256(path.read_bytes()).hexdigest() == entry["sha256"],
                "Changed pinned prefix source: " + entry["path"])
    prefix = json.loads(gzip.decompress(prefix_path.read_bytes()))
    grade29 = json.loads(gzip.decompress(grade29_path.read_bytes()))
    tail_stream = gzip.open(args.input, "rt", encoding="utf-8")
    tail = json.loads(tail_stream.readline())
    require(prefix["arity"] == prefix["path_length"] == 3, "Wrong prefix problem")
    require([c["grade"] for c in prefix["exact_grades"]] == list(range(1, 19)),
            "Incomplete exact prefix")
    require([c["grade"] for c in prefix["modular_zero_grades"]] == list(range(19, 29)),
            "Incomplete original zero prefix")
    require(grade29["schema_version"] == 1 and grade29["status"] == "FULL",
            "Incomplete grade29 certificate")
    require(grade29["grade"] == 29 and grade29["prime"] == P, "Wrong grade29 problem")
    require(grade29["columns"] == grade29["rank"] == len(grade29["relations"]) == 46928,
            "Wrong grade29 certificate size")
    prefix_hash = sha256(prefix_path.read_bytes()).hexdigest()
    grade29_hash = sha256(grade29_path.read_bytes()).hexdigest()
    require(grade29["prefix_certificate_sha256"] == prefix_hash, "Wrong grade29 prefix")
    require(tail["schema_version"] == 1 and tail["arity"] == tail["path_length"] == 3,
            "Wrong tail problem")
    require(tail["prime"] == P, "Wrong tail field")
    require(tail["prefix_certificate_sha256"] == prefix_hash and
            tail["grade29_certificate_sha256"] == grade29_hash, "Wrong tail prefix")
    last_grade = tail["last_grade"]
    require(isinstance(last_grade, int) and 30 <= last_grade <= 55, "Wrong tail endpoint")
    require(tail["first_grade"] == 30, "Wrong first tail grade")
    require(tail["format"] == "header_then_one_grade_per_line", "Wrong streaming format")

    s = load("whole_cubezero_kernel", folder / "verify_subtree_shuffle.py")
    v = load("whole_cubezero_prefix", folder / "verify_cubic_cube_zero_shuffle.py")
    q = s.ShuffleQuotient(3, 3)
    exact, modular_pivots = [], [{}]
    for item in prefix["exact_grades"]:
        row, pivots = v.verify_exact(q, item, s)
        exact.append(row)
        modular_pivots.append(pivots)
        print("EXACT", json.dumps(row), flush=True)
    require([len(basis) for basis in q.basis] == DIMENSIONS, "Wrong exact dimensions")
    require(sum(DIMENSIONS) == 235 and len(q.basis[18]) == 2, "Wrong sharp dimension")

    # A literal 37-leaf monomial, not an inverse-coefficient combination.
    generator = 0
    a = q.graft((generator, generator, generator))
    b = q.graft((a, a, a))
    witness = b
    for _ in range(4):
        witness = q.graft((generator, generator, witness))
        witness = q.graft((generator, a, witness))
    for _ in range(2):
        witness = q.graft((generator, generator, witness))
    require(witness == 3569 and q.grade[witness] == 18, "Wrong 37-leaf monomial")
    require(v.reduce_rref({witness: Fraction(1)}, q.pivots[18]) == {witness: 1},
            "The sharp monomial is not a nonzero free basis element")
    witness_record = {"tree_id": witness, "grade": 18, "leaf_count": 37,
                      "expression": "L^2 (M L)^4 b",
                      "definitions": "a=mu(x,x,x), b=mu(a,a,a), L(y)=mu(x,x,y), M(y)=mu(x,a,y)",
                      "rational_remainder": {str(witness): "1"}}

    class ModularQuotient(s.ShuffleQuotient):
        def lower_shuffle(self, side, bottom):
            if not side:
                return {bottom: 1}
            grade = self.grade[bottom] + sum(self.grade[i] for i in side) + len(side) // 2
            require(grade < len(self.basis), "Using an unproved lower grade")
            if grade >= 19:
                require(not self.basis[grade], "Zero shortcut on a nonzero piece")
                return {}
            return self.nonzero_lower_shuffle(side, bottom, grade)

        @lru_cache(maxsize=250000)
        def nonzero_lower_shuffle(self, side, bottom, grade):
            return v.nf_mod(self.root_shuffles(side, bottom), self.pivots[grade])

    class TriangularPrefix:
        """Gaussian row space with exact propagation of proved unit rows."""

        def __init__(self, columns):
            self.rows = {}
            self.reverse = {}
            self.units = set()
            self.first, self.last = min(columns), max(columns)
            require(self.last - self.first + 1 == len(columns), "Nonconsecutive columns")
            self.frontier = self.first

        def insert(self, row):
            require(all(self.first <= i <= self.last for i in row),
                    "Relation outside the current prequotient")
            reduced = v.nf_mod(row, self.rows)
            if not reduced:
                return False
            pivot = max(reduced)
            scale = pow(reduced[pivot], -1, P)
            normalized = {j: coefficient * scale % P for j, coefficient in reduced.items()}
            require(normalized[pivot] == 1 and all(self.first <= j <= pivot for j in normalized),
                    "Invalid normalized pivot row")
            self.rows[pivot] = normalized
            pending_units = []
            if len(normalized) == 1:
                pending_units.append(pivot)
            else:
                for j in normalized:
                    if j != pivot:
                        self.reverse.setdefault(j, set()).add(pivot)
            while pending_units:
                unit = pending_units.pop()
                if unit in self.units:
                    continue
                require(self.rows[unit] == {unit: 1}, "An unproved unit row was used")
                self.units.add(unit)
                for parent in self.reverse.pop(unit, ()):
                    require(parent > unit and self.rows[parent].get(unit, 0) != 0,
                            "Invalid unit-row dependency")
                    # Subtract this coefficient times the existing unit row.
                    # The row span and the parent's leading pivot are unchanged.
                    del self.rows[parent][unit]
                    if len(self.rows[parent]) == 1:
                        pending_units.append(parent)
            while self.frontier in self.units:
                require(self.rows[self.frontier] == {self.frontier: 1}, "Invalid unit prefix")
                self.frontier += 1
            return True

    control = TriangularPrefix([4, 5])
    require(control.insert({5: 1, 4: 3}) and control.frontier == 4,
            "A gap was incorrectly treated as a unit prefix")
    require(control.insert({4: 2}) and control.rows == {5: {5: 1}, 4: {4: 1}},
            "Valid triangular cleanup failed")
    require(not control.insert({5: 9, 4: 7}), "A dependent relation was accepted")
    detached = TriangularPrefix([4, 5, 6])
    require(detached.insert({6: 1, 5: 2}) and detached.insert({5: 3}),
            "Independent detached rows were rejected")
    require(detached.frontier == 4 and detached.units == {5, 6} and
            detached.rows == {5: {5: 1}, 6: {6: 1}},
            "Unit propagation across a missing smaller column failed")
    invalid_column_rejected = False
    try:
        control.insert({3: 1})
    except ValueError:
        invalid_column_rejected = True
    require(invalid_column_rejected, "An out-of-prequotient relation was accepted")

    mq = ModularQuotient(3, 3)
    mq.grade = list(q.grade)
    mq.tree_ids = dict(q.tree_ids)
    mq.basis = list(q.basis)
    mq.pivots = modular_pivots
    zeros = []
    def zero_certificates():
        yield from prefix["modular_zero_grades"]
        yield grade29
        for expected_grade in range(30, last_grade + 1):
            line = tail_stream.readline()
            require(bool(line), "Missing streamed tail grade")
            item = json.loads(line)
            require(item["grade"] == expected_grade, "Reordered streamed tail grade")
            yield item
        require(not tail_stream.read().strip(), "Unexpected extra tail records")
        tail_stream.close()

    for item in zero_certificates():
        k = item["grade"]
        require(item["prime"] == P and k == len(mq.basis), "Nonconsecutive zero piece")
        children = v.make_children(mq, k)
        require(len(item["relations"]) == len(children), "Wrong relation count")
        if "columns" in item:
            require(item["columns"] == len(children), "Wrong prequotient dimension")
        pivots = TriangularPrefix(children)
        grade_start = last_log = time.monotonic()
        for side, bottom in item["relations"]:
            side = tuple(side)
            v.validate_substitution(mq, k, side, bottom)
            require(pivots.insert(mq.root_shuffles(side, bottom)),
                    "Dependent selected zero-piece relation")
            if k >= 29 and time.monotonic() - last_log >= 15:
                print("PROGRESS", json.dumps({"grade": k, "rank": len(pivots.rows),
                      "columns": len(children), "grade_seconds": time.monotonic() - grade_start,
                      "total_seconds": time.monotonic() - start}), flush=True)
                last_log = time.monotonic()
        require(len(pivots.rows) == len(children) and pivots.frontier == pivots.last + 1,
                "The whole prequotient did not reduce to zero")
        mq.pivots.append({})
        mq.basis.append(())
        # Every use of a zero lower grade is intercepted above. Its old
        # tree dictionary entries are therefore unnecessary, not relations.
        for cc in children.values():
            del mq.tree_ids[cc]
        row = {"grade": k, "columns": len(children),
               "independent_valid_relations": len(pivots.rows), "dimension": 0,
               "seconds": time.monotonic() - grade_start}
        zeros.append(row)
        print("ZERO", json.dumps(row), flush=True)

    complete = last_grade == 55
    require([r["grade"] for r in zeros] == list(range(19, last_grade + 1)),
            "Incomplete zero interval")
    if complete:
        require(1 + 3 * 18 == 55 and all(not mq.basis[k] for k in range(19, 56)),
                "Root-grafting propagation hypotheses failed")
    result = {"status": "PASS", "arithmetic": "QQ and F_1000000009; standard library only",
              "exact_grades": exact, "modular_zero_grades": zeros,
              "last_verified_zero_grade": last_grade, "grade18_witness": witness_record,
              "full_algebra_dimension": 235 if complete else None,
              "optimal_repeated_element_leaf_threshold": 39 if complete else None,
              "mixed_input_nilpotence_for_arbitrary_algebras_claimed": False,
              "input_sha256": sha256(args.input.read_bytes()).hexdigest(),
              "prefix_certificate_sha256": prefix_hash,
              "grade29_certificate_sha256": grade29_hash,
              "prefix_ledger_sha256": sha256(ledger_path.read_bytes()).hexdigest(),
              "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "controls": {"unit_row_propagation": True, "detached_unit_block": True, "gap_not_skipped": True,
                           "dependent_relation_rejected": True,
                           "out_of_prequotient_column_rejected": True,
                           "inherited": v.controls()},
              "seconds": time.monotonic() - start}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS", json.dumps({"last_zero_grade": last_grade,
          "full_algebra_dimension": result["full_algebra_dimension"],
          "seconds": result["seconds"]}), flush=True)


if __name__ == "__main__":
    main()
