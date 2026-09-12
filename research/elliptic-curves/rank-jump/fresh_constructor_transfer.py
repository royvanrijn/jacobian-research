#!/usr/bin/env python3
"""Rank-blind input freeze and cost audit for the eight-fibre transfer gate.

This entry point performs rational arithmetic and reads retained receipts only.
It does not initialize a number field, collect a relation, or search for points.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from math import gcd, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
OUT = ART / "fresh_constructor_transfer_v1"
ATLAS = ART / "compact_five_mw16_atlas_v1.json"
DATABASE = ROOT / "elliptic-curves/data/research_curves/database.json"
SEED = "MW16-05 fresh dependency transfer v1; 2026-09-12"
BAND = (8, 32)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(value, stream, sort_keys=True, indent=2)
        stream.write("\n")


def read(path):
    return json.loads(Path(path).read_text())


def evaluate(coefficients, t):
    result = Q(0)
    for coefficient in reversed(coefficients):
        result = result * t + Q(coefficient)
    return result


def section_value(coordinate, t):
    return evaluate(coordinate["numerator_coefficients_low_to_high"], t) / evaluate(
        coordinate["denominator_coefficients_low_to_high"], t)


def invariants(ainvs):
    a1, a2, a3, a4, a6 = map(Q, ainvs)
    b2, b4, b6 = a1*a1+4*a2, a1*a3+2*a4, a3*a3+4*a6
    return b2*b2-24*b4, -b2**3+36*b2*b4-216*b6


def rational_square(q):
    q = Q(q)
    return q >= 0 and isqrt(q.numerator)**2 == q.numerator and isqrt(q.denominator)**2 == q.denominator


def rational_power(q, degree):
    q = Q(q)
    if q <= 0:
        return False
    def power(n):
        lo, hi = 0, 1 << ((n.bit_length()+degree-1)//degree)
        while lo+1 < hi:
            m = (lo+hi)//2
            if m**degree < n:
                lo = m
            else:
                hi = m
        return hi**degree == n
    return power(q.numerator) and power(q.denominator)


def q_isomorphic(left, right):
    """Exact Q-isomorphism for smooth curves via c4,c6, including j=0,1728."""
    a, b = map(Q, left)
    c, d = map(Q, right)
    if not a or not c:
        return a == c == 0 and bool(b and d) and rational_power(b/d, 6)
    if not b or not d:
        return b == d == 0 and rational_power(a/c, 4)
    r4, r6 = a/c, b/d
    u2 = r6/r4
    return rational_square(u2) and u2*u2 == r4 and u2**3 == r6


def fibre(family, t):
    a = evaluate(family["A_coefficients_low_to_high"], t)
    b = evaluate(family["B_coefficients_low_to_high"], t)
    points = []
    for section in family["sections"]:
        x, y = (section_value(section[k], t) for k in ("X", "Y"))
        assert y*y == x*x*x+a*x+b
        points.append([str(x), str(y)])
    assert len(points) == 16
    return {"parameter": str(t), "short_model": [str(a), str(b)],
            "generic_points": points, "generic_height_gram": family["generic_height_gram"],
            "generic_independence": "PENDING_FRESH_CERTIFICATE"}


def freeze_inputs():
    if OUT.exists():
        raise FileExistsError("Retain the frozen input panel; do not reselect it")
    raw_family = next(f for f in read(ATLAS)["families"] if f["fibration_id"] == "a1-fibration-05")
    family = {key: raw_family[key] for key in (
        "fibration_id", "presentation_id", "A_coefficients_low_to_high",
        "B_coefficients_low_to_high", "sections", "generic_height_gram")}
    # The selection stage alone reads the historical database. Export only
    # equations, never ranks, point coordinates, scores or fibre labels.
    equations = []
    for row in read(DATABASE)["curves"]:
        model = list(map(str, row["ainvs"]))
        c4, c6 = invariants(model)
        assert c4**3 != c6*c6
        equations.append({"ainvs": model, "c4": str(c4), "c6": str(c6)})
    assert len(equations) == len(read(DATABASE)["curves"])
    aliases = {"source_sha256": sha(DATABASE), "equations": equations,
               "boundary": "Frozen equation projection only; no historical ranks, parameters, scores or point coordinates."}
    population = [Q(m, n) for n in range(1, BAND[1]+1) for m in range(-BAND[1], BAND[1]+1)
                  if gcd(m, n) == 1 and BAND[0] <= max(abs(m), n) <= BAND[1]]
    assert len(population) == len(set(population))
    def key(t):
        return hashlib.sha256((SEED+"\n"+str(t)).encode()).hexdigest()
    population.sort(key=lambda t: (key(t), t.numerator, t.denominator))
    selected, ledger = [], []
    c = fibre(family, Q(3, 17))
    control_inv = (-48*Q(c["short_model"][0]), -864*Q(c["short_model"][1]))
    known = [(Q(e["c4"]), Q(e["c6"])) for e in equations]
    accepted = []
    for ordinal, t in enumerate(population):
        a = evaluate(family["A_coefficients_low_to_high"], t)
        b = evaluate(family["B_coefficients_low_to_high"], t)
        inv = (-48*a, -864*b)
        reason, witness = None, None
        if 4*a**3+27*b*b == 0:
            reason = "SINGULAR"
        elif t == Q(3, 17) or q_isomorphic(inv, control_inv):
            reason = "COMMISSIONING_CONTROL_ALIAS"
        else:
            for i, old in enumerate(known):
                if q_isomorphic(inv, old):
                    reason, witness = "FROZEN_EQUATION_ALIAS", i
                    break
            if reason is None:
                for i, old in enumerate(accepted):
                    if q_isomorphic(inv, old):
                        reason, witness = "EARLIER_SELECTED_EQUATION_ALIAS", i
                        break
        record = {"ordinal": ordinal, "parameter": str(t), "ordering_sha256": key(t),
                  "disposition": reason or "SELECTED", "alias_row": witness}
        if reason is None:
            packet = fibre(family, t)
            packet["id"] = "fresh-%02d" % (len(selected)+1)
            selected.append(packet)
            accepted.append(inv)
            record["id"] = packet["id"]
        ledger.append(record)
        if len(selected) == 8:
            break
    assert len(selected) == 8
    write(OUT/"family.json", family)
    write(OUT/"equation-aliases.json", aliases)
    write(OUT/"commissioning-input.json", {"id": "control", **c})
    for packet in selected:
        write(OUT/"inputs"/(packet["id"]+".json"), packet)
    write(OUT/"selection.json", {
        "schema": "fresh-constructor-transfer.selection.v1",
        "status": "INPUTS_FROZEN_ARITHMETIC_NOT_STARTED",
        "family": "MW16-05", "height": "max(abs(m),n), gcd(m,n)=1, n>0",
        "height_band_inclusive": list(BAND), "population_size": len(population),
        "ordering": "SHA256(seed + newline + canonical Fraction string), then m,n",
        "seed": SEED, "ledger": ledger,
        "fresh_parameters": [p["parameter"] for p in selected],
        "no_replacements_after_arithmetic": True,
        "sources": {str(p.relative_to(ROOT)): sha(p) for p in (ATLAS, DATABASE, Path(__file__))},
        "bindings": {str(p.relative_to(OUT)): sha(p) for p in sorted(OUT.rglob("*.json"))},
        "boundary": "Rank-blind finite panel selection and exact generic point substitution only. No fresh number-field, local, relation, class or point search has run. Known-alias exclusion is complete relative to this frozen equation list, not every historical experiment."})
    print(json.dumps({"selected": [p["parameter"] for p in selected], "considered": len(ledger)}))


def historical_cost():
    """An honest subtotal: never turn absent metering into zero cost."""
    rows = []
    def add(path, field, stage):
        p = ROOT/path
        value = read(p)
        for part in field.split("."):
            value = value[part]
        rows.append({"stage": stage, "path": path, "sha256": sha(p),
                     "field": field, "recorded_wall_seconds": float(value)})
    for stem in ("norm-batch", "special-primes", "small-base-targets"):
        for phase in ("worker", "check", "audit"):
            path = f"artifacts/local/elliptic-curves/small-conductor-{stem}-v1/{phase}.supervisor.json"
            if (ROOT/path).exists():
                add(path, "wall_seconds", f"early {stem}: {phase}")
    for kind, count in (("wave", 4), ("strip_wave", 5)):
        for i in range(1, count+1):
            path = ART/f"small_conductor_class_target_{kind}_{i:03d}_v1.json"
            d = read(path)
            # Each chunk time is disjoint within this worker. No enclosing
            # worker time is added again. Setup/audit time remains missing.
            for chunk in d["chunks"]:
                add(chunk["path"], "wall_seconds", f"early {kind} {i}: target")
    add("artifacts/local/rank-jump-reference-strict-class-construction-v1/execution.json", "elapsed_seconds", "failed unseeded BNF")
    for name, field in (
        ("rank_jump_seeded_reference_class_v1.json", "reference.elapsed_seconds"),
        ("rank_jump_seeded_reference_relations_v1.json", "terminals.reference.elapsed_seconds"),
        ("rank_jump_seeded_reference_relations_memory_v1.json", "terminals.reference.elapsed_seconds"),
        ("rank_jump_reference_class_targeted_relations_v1.json", "elapsed_seconds"),
        ("rank_jump_reference_adaptive_class_completion_v1.json", "elapsed_seconds"),
        ("rank_jump_reference_small_representative_class_v1.json", "elapsed_seconds")):
        add(str((ART/name).relative_to(ROOT)), field, name.removesuffix("_v1.json"))
    missing = [
        "Cold equation/factorization/prime proofs and maximal-order setup before retained hints",
        "Initial factor-base construction and generic-only local/anchor preparation",
        "Target-wave selection/setup and full audit overhead not inside timed chunks",
        "Standalone strict extraction, large-support re-extraction and independent Hilbert/ideal verification",
        "Half-ideal compaction, cover compilation, reduction, failed lifting preparations and all final verifications",
        "Implementation/debugging time; never presented as measured arithmetic"]
    result = {"schema": "fresh-constructor-transfer.historical-cost.v1",
              "status": "INCOMPLETE_HISTORICAL_METERING",
              "recorded_nonoverlapping_wall_seconds_subtotal": sum(r["recorded_wall_seconds"] for r in rows),
              "receipts": rows, "unmetered_components": missing,
              "budget_gate": "Do not size the transfer from this subtotal or the 18.401-second post-construction recovery. Meter commissioning end to end before freezing numeric fresh-fibre allowances.",
              "source_sha256": sha(Path(__file__))}
    write(OUT/"historical-cost-audit.json", result)
    print(json.dumps({k:v for k,v in result.items() if k != "receipts"}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze-inputs", "historical-cost"))
    args = parser.parse_args()
    {"freeze-inputs": freeze_inputs, "historical-cost": historical_cost}[args.mode]()
