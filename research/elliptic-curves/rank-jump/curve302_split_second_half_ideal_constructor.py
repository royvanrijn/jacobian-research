#!/usr/bin/env python3
"""Second-factor, equation/generic-only principal-square constructor for E302."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import retrospective as r


HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "CURVE302_SPLIT_SECOND_HALF_IDEAL_CONSTRUCTOR_PROTOCOL.json"
BASE_PATH = HERE / "curve302_square_half_ideal_constructor.py"
OUTPUT = r.OUT / "rank_jump_curve302_split_second_half_ideal_constructor_v1.json"
WORK = r.ROOT / "artifacts/local/rank-jump-curve302-split-second-half-ideal-constructor-v1"


def base_module():
    spec = importlib.util.spec_from_file_location("curve302_square_half_base", BASE_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def bindings(base):
    paths = [Path(__file__), PROTOCOL, BASE_PATH, base.ARITH, r.ROOT / "elliptic-curves/cas/research_runtime/local_kummer.py"]
    return {str(path.relative_to(r.ROOT)): r.digest(path.read_bytes()) for path in paths}


def setup(base):
    context = base.setup()
    policy = r.read(PROTOCOL)
    from sage.all import prime_range

    targets = []
    for p0 in prime_range(3, policy["bounds"]["prime_search_bound"] + 1):
        p = int(p0)
        if p in set(context["data"]["S_finite"]) or context["cubic"].discriminant() % p == 0:
            continue
        primes = list(context["pari"].idealprimedec(context["nf"], p))
        degree_one = [(index, prime) for index, prime in enumerate(primes) if int(prime[3]) == 1]
        if len(degree_one) < 2:
            continue
        index, prime = degree_one[1]
        targets.append({"p": p, "prime_index": index, "hnf": str(context["pari"].idealhnf(context["nf"], prime)), "ideal": prime})
        if len(targets) == policy["bounds"]["half_ideal_primes"]:
            break
    assert len(targets) == policy["bounds"]["half_ideal_primes"]
    context["policy"] = policy
    context["targets"] = targets
    return context


def worker():
    base = base_module()
    WORK.mkdir(parents=True, exist_ok=True)
    context = setup(base)
    targets = context["targets"]
    pari, nf = context["pari"], context["nf"]
    total = 1 << len(targets)
    low_count = len(targets) // 2
    unit = pari.idealhnf(nf, 1)

    def tables(primes):
        halves, squares = [unit], [unit]
        for mask in range(1, 1 << len(primes)):
            bit = mask & -mask
            index = bit.bit_length() - 1
            prior = mask ^ bit
            halves.append(pari.idealmul(nf, halves[prior], primes[index]["ideal"]))
            squares.append(pari.idealmul(nf, squares[prior], pari.idealpow(nf, primes[index]["ideal"], 2)))
        return halves, squares

    low_half, low_square = tables(targets[:low_count])
    high_half, high_square = tables(targets[low_count:])
    candidates = []
    for high_mask in range(1 << (len(targets) - low_count)):
        for low_mask in range(1 << low_count):
            mask = (high_mask << low_count) | low_mask
            if not mask:
                continue
            square = pari.idealmul(nf, high_square[high_mask], low_square[low_mask])
            reduced, alpha = pari.idealred(nf, [square, 1])
            if reduced == unit:
                half = pari.idealmul(nf, high_half[high_mask], low_half[low_mask])
                record = base.candidate_record(context, mask, half, square, alpha)
                if record is not None:
                    candidates.append(record)
                    print("STRICT_CLASS_CANDIDATE", mask, flush=True)
                    break
            if mask % context["policy"]["bounds"]["page_masks"] == 0:
                page = mask // context["policy"]["bounds"]["page_masks"] - 1
                r.write_new(WORK / f"page_{page:04d}.json", {
                    "schema": "rank-jump.curve302-split-second-half-ideal-page.v1",
                    "mask_start": str(mask - context["policy"]["bounds"]["page_masks"] + 1),
                    "mask_end": str(mask), "strict_class_candidates": candidates,
                })
                print("PAGE", page, mask, flush=True)
        if candidates:
            break
    result = {
        "schema": "rank-jump.curve302-split-second-half-ideal-constructor.v1", "status": "PASS",
        "positive_endpoint": "STRICT_CLASS_CANDIDATE_PENDING_POSTHOC_M24_AND_COVER" if candidates else "NOT_REACHED",
        "bindings": bindings(base),
        "targets": [{key: value for key, value in target.items() if key != "ideal"} for target in targets],
        "tested_subset_count": str(mask), "strict_class_candidates": candidates, "candidate_count": len(candidates),
        "boundary": context["policy"]["failure_semantics"],
    }
    r.write_new(OUTPUT, result)
    print(result["positive_endpoint"], result["tested_subset_count"], flush=True)


def check():
    base = base_module()
    stored = r.read(OUTPUT)
    assert stored["status"] == "PASS" and stored["bindings"] == bindings(base)
    print("PASS", stored["positive_endpoint"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("worker", "check"))
    args = parser.parse_args()
    if args.mode == "worker":
        worker()
    else:
        check()
