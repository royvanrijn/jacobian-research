#!/usr/bin/env python3
"""Independent Sage replay of the input freeze; no field or point search."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts/generated-results/elliptic-curves/fresh_constructor_transfer_v1"


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    from sage.all import QQ, ZZ, PolynomialRing, EllipticCurve
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU, (60, 65))
    start = time.process_time()
    selection = read(OUT/"selection.json")
    assert selection["status"] == "INPUTS_FROZEN_ARITHMETIC_NOT_STARTED"
    assert selection["no_replacements_after_arithmetic"]
    for name, h in selection["sources"].items():
        assert digest(ROOT/name) == h, name
    for name, h in selection["bindings"].items():
        assert digest(OUT/name) == h, name
    family = read(OUT/"family.json")
    aliases = read(OUT/"equation-aliases.json")
    source = read(ROOT/"elliptic-curves/data/research_curves/database.json")
    assert digest(ROOT/"elliptic-curves/data/research_curves/database.json") == aliases["source_sha256"]
    assert [e["ainvs"] for e in aliases["equations"]] == [c["ainvs"] for c in source["curves"]]
    assert all(set(e) == {"ainvs", "c4", "c6"} for e in aliases["equations"])
    for equation in aliases["equations"]:
        E = EllipticCurve(list(map(QQ, equation["ainvs"])))
        assert list(map(str, E.c_invariants())) == [equation["c4"], equation["c6"]]
    lo, hi = selection["height_band_inclusive"]
    candidates = {QQ(m)/n for n in range(1, hi+1) for m in range(-hi, hi+1)
                  if ZZ(m).gcd(n) == 1 and lo <= max(abs(m), n) <= hi}
    assert len(candidates) == selection["population_size"]
    def order(t):
        return hashlib.sha256((selection["seed"]+"\n"+str(t)).encode()).hexdigest()
    candidates = sorted(candidates, key=lambda t: (order(t), t.numerator(), t.denominator()))
    R = PolynomialRing(QQ, "t")
    A, B = (R(family[k+"_coefficients_low_to_high"]) for k in ("A", "B"))
    control = EllipticCurve([A(QQ(3)/17), B(QQ(3)/17)])
    known = [EllipticCurve(list(map(QQ, e["ainvs"]))) for e in aliases["equations"]]
    # j merely screens comparisons; Sage Q-isomorphism decides aliases.
    known_j = {}
    for i, E in enumerate(known):
        known_j.setdefault(E.j_invariant(), []).append(i)
    retained, accepted, equations = [], [], []
    for number, entry in enumerate(selection["ledger"]):
        t = candidates[number]
        assert entry["ordinal"] == number and entry["parameter"] == str(t)
        assert entry["ordering_sha256"] == order(t)
        a, b = A(t), B(t)
        disposition, witness = "SELECTED", None
        if 4*a**3+27*b**2 == 0:
            disposition = "SINGULAR"
        else:
            E = EllipticCurve([a, b])
            if E.is_isomorphic(control):
                disposition = "COMMISSIONING_CONTROL_ALIAS"
            else:
                for i in known_j.get(E.j_invariant(), []):
                    if E.is_isomorphic(known[i]):
                        disposition, witness = "FROZEN_EQUATION_ALIAS", i
                        break
                if disposition == "SELECTED":
                    for i, previous in enumerate(equations):
                        if E.is_isomorphic(previous):
                            disposition, witness = "EARLIER_SELECTED_EQUATION_ALIAS", i
                            break
        assert (disposition, witness) == (entry["disposition"], entry["alias_row"])
        if disposition == "SELECTED":
            packet = read(OUT/"inputs"/(entry["id"]+".json"))
            assert list(map(QQ, packet["short_model"])) == [a, b]
            assert packet["parameter"] == str(t)
            assert packet["generic_height_gram"] == family["generic_height_gram"]
            assert packet["generic_independence"] == "PENDING_FRESH_CERTIFICATE"
            assert len(packet["generic_points"]) == 16
            for section, point in zip(family["sections"], packet["generic_points"]):
                coordinates = []
                for key in ("X", "Y"):
                    row = section[key]
                    coordinates.append(R(row["numerator_coefficients_low_to_high"])(t) /
                                       R(row["denominator_coefficients_low_to_high"])(t))
                assert list(map(QQ, point)) == coordinates
                assert E.is_on_curve(*coordinates)
            retained.append(str(t))
            accepted.append(entry["id"])
            equations.append(E)
    assert len(retained) == 8 and retained == selection["fresh_parameters"]
    assert selection["ledger"][-1]["disposition"] == "SELECTED"
    cost = read(OUT/"historical-cost-audit.json")
    total = 0.0
    keys = set()
    for row in cost["receipts"]:
        path = ROOT/row["path"]
        assert digest(path) == row["sha256"]
        value = read(path)
        for key in row["field"].split("."):
            value = value[key]
        assert float(value) == row["recorded_wall_seconds"]
        assert (row["path"], row["field"]) not in keys
        keys.add((row["path"], row["field"]))
        total += float(value)
    assert total == cost["recorded_nonoverlapping_wall_seconds_subtotal"]
    assert cost["status"] == "INCOMPLETE_HISTORICAL_METERING" and cost["unmetered_components"]
    result = {"schema": "fresh-constructor-transfer.input-replay.v1", "status": "PASS",
              "fresh_fibres": 8, "generic_point_equations_checked": 128,
              "frozen_alias_equations": len(known), "new_field_computations": 0,
              "new_relation_searches": 0, "point_searches": 0,
              "historical_cost_is_complete": False, "historical_receipts_checked": len(keys),
              "source_sha256": digest(Path(__file__)),
              "bindings": {name: digest(OUT/name) for name in (
                  "selection.json", "family.json", "equation-aliases.json", "historical-cost-audit.json")},
              "cpu_seconds": time.process_time()-start,
              "boundary": "Independent input-selection and bookkeeping replay only. No new class, rational lift, rank gain or full reference-cost measurement."}
    if args.write:
        with (OUT/"independent-input-replay.json").open("x") as stream:
            json.dump(result, stream, sort_keys=True, indent=2)
            stream.write("\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
