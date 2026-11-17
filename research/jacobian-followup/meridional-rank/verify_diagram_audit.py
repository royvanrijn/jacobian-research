#!/usr/bin/env python3
"""Standard-library named-diagram transfer and independent tangle construction.

No SnapPy/Spherogram or producer matching code is imported. Closed planar
maps are compared by a dart bijection preserving edge pairing, cyclic
rotation (globally possibly reversed), and crossing data (globally
possibly mirrored). Rational tangles are constructed from elementary
paired-endpoint graphs, with continued fractions checked exactly.
"""

from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def edge_pairs(pd):
    ends = {}
    for c, crossing in enumerate(pd):
        for h, edge in enumerate(crossing):
            ends.setdefault(edge, []).append(4 * c + h)
    assert len(ends) == 2 * len(pd) and all(len(v) == 2 for v in ends.values())
    return {a: b for pair in ends.values() for a, b in [pair, pair[::-1]]}


def walk(pd, start=0):
    pairs = edge_pairs(pd)
    current, visits = start, []
    for _ in range(2 * len(pd)):
        assert current not in visits
        visits.append(current)
        opposite = 4 * (current // 4) + (current + 2) % 4
        current = pairs[opposite]
    assert current == start and {d // 4 for d in visits} == set(range(len(pd)))
    return visits


def plane_map(pd):
    pairs = edge_pairs(pd)
    unseen = set(pairs)
    faces = 0
    while unseen:
        current = start = min(unseen)
        while True:
            unseen.remove(current)
            d = pairs[current]
            current = 4 * (d // 4) + (d + 1) % 4
            if current == start:
                break
        faces += 1
    assert len(pd) - 2 * len(pd) + faces == 2


def map_witness(left, right):
    assert len(left) == len(right)
    a, b = edge_pairs(left), edge_pairs(right)
    size = len(a)
    for orientation in [1, -1]:
        for first in range(size):
            swap = first % 2
            mapping, todo, failed = {}, [(0, first)], False
            while todo and not failed:
                x, y = todo.pop()
                if x in mapping:
                    failed = mapping[x] != y
                    continue
                if (x % 2) ^ (y % 2) != swap:
                    failed = True
                    break
                mapping[x] = y
                todo.append((a[x], b[y]))
                todo.append((4 * (x // 4) + (x + 1) % 4,
                             4 * (y // 4) + (y + orientation) % 4))
            if not failed and len(mapping) == len(set(mapping.values())) == size:
                assert all(mapping[a[x]] == b[mapping[x]] for x in a)
                return {"orientation": orientation, "over_under_swap": swap,
                        "dart_images": [mapping[i] for i in range(size)]}
    raise AssertionError("no exact crossing-preserving planar-map bijection")


class Tangle:
    # Endpoints 4*c+h are crossing darts. Boundary ports -1,-2,-3,-4
    # are bottom-left, bottom-right, top-left, top-right respectively.
    def __init__(self, count, pairs):
        self.n = count
        self.pairs = dict(pairs)
        assert set(self.pairs) == set(range(4 * count)) | {-1, -2, -3, -4}
        assert all(self.pairs[self.pairs[x]] == x and self.pairs[x] != x for x in self.pairs)

    @classmethod
    def simple(cls, kind):
        if kind == "zero":
            pairs, count = [(-1, -2), (-3, -4)], 0
        elif kind == "infinity":
            pairs, count = [(-1, -3), (-2, -4)], 0
        else:
            assert kind == "one"
            pairs, count = [(-1, 0), (-2, 1), (-3, 3), (-4, 2)], 1
        return cls(count, {a: b for pair in pairs for a, b in [pair, pair[::-1]]})

    def relabel(self, function):
        return Tangle(self.n, {function(x): function(y) for x, y in self.pairs.items()})

    def rotate(self, turns):
        order = [0, 1, 3, 2]
        # New boundary p receives old boundary order[index(p)+turns].
        old_to_new = {order[(i + turns) % 4]: order[i] for i in range(4)}
        return self.relabel(lambda x: -1 - old_to_new[-1 - x] if x < 0 else x)

    def inverse(self):
        rotated = self.rotate(1)
        return rotated.relabel(lambda x: x if x < 0 else 4 * (x // 4) + (x - 1) % 4)

    def add(self, other):
        def bshift(x):
            return x + 4 * self.n if x >= 0 else x - 4
        pairs = {**self.pairs, **{bshift(x): bshift(y) for x, y in other.pairs.items()}}
        for x, y in [(-2, -5), (-4, -7)]:
            p, q = pairs.pop(x), pairs.pop(y)
            assert p != y and q != x, "closed component created in tangle sum"
            pairs[p], pairs[q] = q, p
        boundary = {-1: -1, -6: -2, -3: -3, -8: -4}
        return Tangle(self.n + other.n, {boundary.get(x, x): boundary.get(y, y) for x, y in pairs.items()})

    @classmethod
    def integer(cls, value):
        assert value >= 0
        result = cls.simple("zero")
        for _ in range(value):
            result = result.add(cls.simple("one"))
        return result

    @classmethod
    def rational(cls, a, b):
        wanted = Fraction(a, b)
        assert wanted > 0
        quotients = []
        while b:
            q, rest = divmod(a, b)
            quotients.append(q)
            a, b = b, rest
        result = cls.simple("infinity")
        fraction = None  # infinity
        for q in reversed(quotients):
            result = cls.integer(q).add(result.inverse())
            fraction = Fraction(q) if fraction is None else Fraction(q) + 1 / fraction
        assert fraction == wanted
        return result

    def closure_pd(self):
        pairs = self.pairs.copy()
        for x, y in [(-1, -2), (-3, -4)]:
            p, q = pairs.pop(x), pairs.pop(y)
            assert p != y and q != x
            pairs[p], pairs[q] = q, p
        labels = {}
        for x in sorted(pairs):
            if x not in labels:
                label = len(labels) // 2
                labels[x] = labels[pairs[x]] = label
        return [[labels[4 * c + h] for h in range(4)] for c in range(self.n)]


def source_arc_assignment(pd, code):
    wanted, ids = [], {}
    for c in code:
        ids.setdefault(abs(c), len(ids))
        wanted.append((ids[abs(c)], int(c > 0)))
    for start in range(4 * len(pd)):
        visits = walk(pd, start)
        for mirror in [0, 1]:
            ids, got = {}, []
            for d in visits:
                ids.setdefault(d // 4, len(ids))
                got.append((ids[d // 4], (d % 2) ^ mirror))
            if got != wanted:
                continue
            # The outgoing edge after each event belongs to its last-under arc.
            previous, edges = len(pd) - 1, {}
            for value, dart in zip(code, visits):
                if value < 0:
                    previous = (previous + 1) % len(pd)
                edge = pd[dart // 4][(dart + 2) % 4]
                assert edge not in edges
                edges[edge] = previous
            return edges, visits, mirror
    raise AssertionError("census over/under walk does not match source")


def rank(rows, prime):
    a = [[x % prime for x in row] for row in rows]
    pivot = 0
    for c in range(len(a[0])):
        found = next((i for i in range(pivot, len(a)) if a[i][c]), None)
        if found is None:
            continue
        a[pivot], a[found] = a[found], a[pivot]
        f = pow(a[pivot][c], -1, prime)
        a[pivot] = [x * f % prime for x in a[pivot]]
        for i in range(pivot + 1, len(a)):
            f = a[i][c]
            a[i] = [(x - f * y) % prime for x, y in zip(a[i], a[pivot])]
        pivot += 1
    return pivot


def check_transferred_certificate(pd, certificate):
    arc_of, visits, mirror = source_arc_assignment(pd, certificate["gauss_code"])
    entrances = {}
    for d in visits:
        entrances.setdefault(d // 4, []).append(d % 4)
    crossings = []
    for c, slots in entrances.items():
        under = next(h for h in slots if h % 2 == mirror)
        over = next(h for h in slots if h % 2 != mirror)
        assert arc_of[pd[c][over]] == arc_of[pd[c][(over + 2) % 4]]
        crossings.append((arc_of[pd[c][over]], arc_of[pd[c][under]], arc_of[pd[c][(under + 2) % 4]]))
    for basis in certificate["lift_basis"]:
        for o, a, b in crossings:
            permutation = list(range(4))
            translation = [0] * 4
            for generator in (o, a, o, b):
                p, q = certificate["arc_transpositions"][generator]
                translation[permutation[p]] += basis[generator]
                translation[permutation[q]] -= basis[generator]
                permutation[p], permutation[q] = permutation[q], permutation[p]
            assert permutation == list(range(4)) and all(x % 5 == 0 for x in translation)
    assert rank(certificate["lift_basis"], 5) == 5
    seen = set(certificate["seeds"])
    assert len(seen) == 5
    for o, old, new in certificate["seed_trace"]:
        assert o in seen and old in seen and new not in seen
        assert (o, old, new) in crossings or (o, new, old) in crossings
        seen.add(new)
    assert len(seen) == len(pd)
    return {"start_dart": visits[0], "mirror": mirror,
            "edge_to_source_arc": sorted(arc_of.items())}


def main():
    audit = json.loads((HERE / "diagram_identity_audit.json").read_text())
    rotations = json.loads((HERE / "montesinos_rotation_results.json").read_text())
    certificates = json.loads((HERE / "sharp_five_certificates.json").read_text())["certificates"]
    by_name = {c["name"]: c for c in certificates}
    assert audit["status"] == rotations["status"] == "COMPLETE"
    assert set(by_name) == set(rotations["matches"]) and len(by_name) == 16
    for receipt in [audit, rotations]:
        for path, digest in receipt["input_sha256"].items():
            assert hashlib.sha256((HERE / path).read_bytes()).hexdigest() == digest
    results = []
    for census in audit["records"]:
        name = census["name"]
        pd = census["census_pd"]
        plane_map(pd)
        transfer = check_transferred_certificate(pd, by_name[name])
        candidate = rotations["matches"][name]
        assert candidate["integer_twists"] == 0 and len(candidate["fractions"]) == 5
        tangle = Tangle.integer(0)
        for i, (a, b) in enumerate(candidate["fractions"]):
            assert 0 < a < b and Fraction(a, b).denominator > 1
            term = Tangle.rational(a, b)
            if candidate["half_turn_mask"] >> i & 1:
                term = term.rotate(2)
            tangle = tangle.add(term)
        independently_built_pd = tangle.closure_pd()
        plane_map(independently_built_pd)
        generated_map = map_witness(independently_built_pd, candidate["pd"])
        census_map = map_witness(candidate["pd"], pd)
        results.append({"name": name, "source_certificate_transfer": transfer,
                        "fractions": candidate["fractions"], "rotation_mask": candidate["half_turn_mask"],
                        "independent_tangle_pd": independently_built_pd,
                        "independent_construction_to_producer": generated_map,
                        "producer_to_named_census": census_map,
                        "bounds": {"mr": 5, "bridge": 5},
                        "prior_family": "Five-noninteger-tangle Montesinos, covered by Boileau-Zieschang1985"})
    # A local crossing switch must destroy these exact diagram witnesses.
    bad = [row[:] for row in audit["records"][0]["census_pd"]]
    bad[0] = bad[0][1:] + bad[0][:1]
    try:
        map_witness(bad, audit["records"][0]["census_pd"])
    except AssertionError:
        pass
    else:
        raise AssertionError("single crossing-switch control accepted")
    inputs = ["verify_diagram_audit.py", "sharp_five_certificates.json", "diagram_identity_audit.json",
              "montesinos_rotation_results.json", "montesinos_audit_plan.json", "montesinos_rotation_plan.json"]
    report = {"status": "PASS", "named_certificate_transfers": 16,
              "independent_Montesinos_constructions": 16, "rejected_single_crossing_switch": True,
              "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in inputs},
              "results": results,
              "claim_boundary": "All16 named sharp-five certificates transfer to independent HT census diagrams and all16 have explicit five-tangle Montesinos presentations, up to global reflection. The MRC equalities are therefore already covered by the classical family theorem. No new knot case or full MRC solution. Finite checks are not a formal verification or full original-paper proof audit."}
    (HERE / "diagram_audit_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ["status", "named_certificate_transfers",
          "independent_Montesinos_constructions", "rejected_single_crossing_switch"]}, indent=2))


if __name__ == "__main__":
    main()
