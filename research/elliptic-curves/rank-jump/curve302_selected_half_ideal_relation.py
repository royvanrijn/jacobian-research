#!/usr/bin/env python3
"""Algebraic principal-square solve for the frozen first E302 twenty-ideal set."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import time

import retrospective as r


HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "CURVE302_SELECTED_HALF_IDEAL_RELATION_PROTOCOL.json"
SELECTED = r.OUT / "rank_jump_curve302_square_half_ideal_constructor_v1.json"
BASE_PATH = HERE / "curve302_square_half_ideal_constructor.py"
OUTPUT = r.OUT / "rank_jump_curve302_selected_half_ideal_relation_v1.json"
WORK = r.ROOT / "artifacts/local/rank-jump-curve302-selected-half-ideal-relation-v1"


def base_module():
    spec = importlib.util.spec_from_file_location("curve302_selected_half_base", BASE_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def bindings(base):
    paths = [Path(__file__), PROTOCOL, SELECTED, BASE_PATH, base.ARITH,
             r.ROOT / "elliptic-curves/cas/research_runtime/local_kummer.py"]
    return {str(path.relative_to(r.ROOT)): r.digest(path.read_bytes()) for path in paths}


def setup(base):
    frozen = r.read(SELECTED)
    assert frozen["status"] == "PASS" and frozen["tested_subset_count"] == str(2**20 - 1)
    assert frozen["candidate_count"] == 0 and len(frozen["targets"]) == 20
    context = base.setup()
    policy = r.read(PROTOCOL)
    pari, nf = context["pari"], context["nf"]
    targets = []
    for record in frozen["targets"]:
        p = record["p"]
        prime = next(prime for prime in pari.idealprimedec(nf, p)
                     if str(pari.idealhnf(nf, prime)) == record["hnf"])
        targets.append({**record, "ideal": prime})
    context["targets"] = targets
    context["policy"] = policy
    return context


def class_condition(mask, coordinates, cyclic):
    sums = [0] * len(cyclic)
    for index, row in enumerate(coordinates):
        if mask >> index & 1:
            for column, value in enumerate(row):
                sums[column] = (sums[column] + value) % cyclic[column]
    return all((2 * value) % modulus == 0 for value, modulus in zip(sums, cyclic)), sums


def half_ideal(context, mask):
    pari, nf = context["pari"], context["nf"]
    answer = pari.idealhnf(nf, 1)
    for index, target in enumerate(context["targets"]):
        if mask >> index & 1:
            answer = pari.idealmul(nf, answer, target["ideal"])
    return answer


def worker():
    base = base_module()
    WORK.mkdir(parents=True, exist_ok=True)
    context = setup(base)
    pari, nf = context["pari"], context["nf"]
    policy = context["policy"]
    pari.allocatemem(64000000, policy["bounds"]["pari_stack_bytes"], silent=True)
    print("BNF_START", flush=True)
    try:
        bnf = pari.bnfinit(nf, 1, pari(policy["bounds"]["compact_relation_tech"]))
    except Exception as error:
        r.write_new(WORK / f"bnf_stack_{policy['bounds']['pari_stack_bytes']}.json", {
            "schema": "rank-jump.curve302-selected-half-ideal-relation-bnf-attempt.v1",
            "selected_half_ideals": len(context["targets"]), "pari_stack_bytes": policy["bounds"]["pari_stack_bytes"],
            "status": "FAILED", "exception_type": type(error).__name__, "reason": str(error),
        })
        raise
    cyclic = list(map(int, bnf.bnf_get_cyc()))
    coordinates = [list(map(int, pari.bnfisprincipal(bnf, target["ideal"], 0))) for target in context["targets"]]
    assert all(len(row) == len(cyclic) for row in coordinates)
    print("BNF_RETURNED", cyclic, flush=True)
    masks = []
    for mask in range(1, 1 << len(context["targets"])):
        passed, coordinates_sum = class_condition(mask, coordinates, cyclic)
        if passed:
            masks.append((mask, coordinates_sum))
    print("SQUARE_CLASS_MASKS", len(masks), flush=True)
    candidates = []
    exact_principal_masks = []
    for mask, coordinates_sum in masks:
        half = half_ideal(context, mask)
        square = pari.idealpow(nf, half, 2)
        relation = pari.bnfisprincipal(bnf, square, 1)
        if any(relation[0]):
            raise ArithmeticError("class-coordinate square condition did not principalize")
        alpha = pari.nfbasistoalg(nf, relation[1])
        assert pari.idealhnf(nf, alpha) == square
        exact_principal_masks.append(str(mask))
        record = base.candidate_record(context, mask, half, square, alpha)
        if record is not None:
            record["class_coordinate_sum"] = list(map(int, coordinates_sum))
            candidates.append(record)
            print("STRICT_CLASS_CANDIDATE", mask, flush=True)
            break
    certified = bool(pari.bnfcertify(bnf) == 1)
    result = {
        "schema": "rank-jump.curve302-selected-half-ideal-relation.v1",
        "status": "PASS" if certified else "UNCERTIFIED_RELATION_SPACE",
        "relation_space_status": "CERTIFIED" if certified else "UNVERIFIED",
        "positive_endpoint": "STRICT_CLASS_CANDIDATE_PENDING_POSTHOC_M24_AND_COVER" if candidates else "NOT_REACHED",
        "bindings": bindings(base),
        "class_group_cyclic": cyclic,
        "compact_relation_tech": policy["bounds"]["compact_relation_tech"],
        "selected_prime_class_coordinates": coordinates,
        "principal_square_mask_count": len(masks),
        "exact_principal_square_masks": exact_principal_masks,
        "strict_class_candidates": candidates,
        "boundary": policy["failure_semantics"],
    }
    r.write_new(OUTPUT, result)
    print(result["status"], result["positive_endpoint"], flush=True)


def check():
    base = base_module()
    stored = r.read(OUTPUT)
    assert stored["bindings"] == bindings(base)
    assert stored["status"] in ("PASS", "UNCERTIFIED_RELATION_SPACE")
    print("PASS", stored["relation_space_status"], stored["positive_endpoint"], flush=True)


def record_compact_timeout():
    """Retain the deliberately terminated compact-BNF attempt as UNKNOWN."""
    base = base_module()
    context = setup(base)
    policy = context["policy"]
    WORK.mkdir(parents=True, exist_ok=True)
    path = WORK / "compact_bnf_timeout.json"
    r.write_new(path, {
        "schema": "rank-jump.curve302-selected-half-ideal-relation-bnf-attempt.v1",
        "selected_half_ideals": len(context["targets"]),
        "pari_stack_bytes": policy["bounds"]["pari_stack_bytes"],
        "compact_relation_tech": policy["bounds"]["compact_relation_tech"],
        "wall_seconds_lower_bound": 1080,
        "status": "TIMEOUT_UNKNOWN",
        "reason": "No class coordinates were emitted before the bounded compact-BNF attempt was terminated.",
        "recorded_at_unix": int(time.time()),
    })
    print(path, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("worker", "check", "record-compact-timeout"))
    args = parser.parse_args()
    if args.mode == "worker":
        worker()
    elif args.mode == "record-compact-timeout":
        record_compact_timeout()
    else:
        check()
