#!/usr/bin/env python3
"""Equation/generic-only square-half-ideal constructor for E302.

The worker deliberately has no path to E302's public exceptional points or
M24.  It searches a finite, ordered family of good-prime half ideals J and
accepts only an exact principal-square identity (alpha)=J^2 followed by the
complete strict-local and generic-independence tests.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

import retrospective as r


PROTOCOL = Path(__file__).with_name("CURVE302_SQUARE_HALF_IDEAL_CONSTRUCTOR_PROTOCOL.json")
ARITH = r.OUT / "rank_jump_curve302_strict_constructor_arithmetic_v1.json"
OUTPUT = r.OUT / "rank_jump_curve302_square_half_ideal_constructor_v1.json"
WORK = r.ROOT / "artifacts/local/rank-jump-curve302-square-half-ideal-constructor-v1"


def bindings():
    paths = [Path(__file__), PROTOCOL, ARITH, r.ROOT / "elliptic-curves/cas/research_runtime/local_kummer.py"]
    return {str(path.relative_to(r.ROOT)): r.digest(path.read_bytes()) for path in paths}


def write_new(path: Path, value):
    r.write_new(path, value)


def setup():
    from sage.all import AA, GF, QQ, ZZ, PolynomialRing, matrix, pari, prime_range, vector
    sys.path.insert(0, str(r.ROOT / "elliptic-curves/cas"))
    from research_runtime.local_kummer import LocalSquareclasses

    policy = r.read(PROTOCOL)
    pari.allocatemem(64000000, policy["bounds"]["pari_stack_bytes"], silent=True)
    data = r.read(ARITH)
    ring = PolynomialRing(QQ, "z")
    cubic = ring(data["cubic_ascending"])
    finite_s = set(data["S_finite"])
    nf = pari.nfinit([pari(cubic), data["S_finite"]])
    generic = [pari.Mod(pari(ring(row["beta_ascending"])), pari(cubic)) for row in data["generic_classes"]]
    local = [LocalSquareclasses(nf, p) for p in data["S_finite"]]
    real_roots = cubic.roots(AA, multiplicities=False)

    def polynomial(value):
        return ring([QQ(pari.lift(value).polcoef(i)) for i in range(3)])

    def signature(value):
        beta = polynomial(value)
        return [int(bit) for place in local for bit in place.signature(value)] + [int(beta(x) < 0) for x in real_roots]

    generic_signatures = matrix(GF(2), [signature(value) for value in generic])
    assert generic_signatures.rank() == 17

    character_blocks = []
    for p0 in prime_range(3, policy["bounds"]["proof_prime_bound"] + 1):
        p = int(p0)
        if p in finite_s or cubic.discriminant() % p == 0:
            continue
        roots = cubic.change_ring(GF(p)).roots(multiplicities=False)
        if len(roots) != 3:
            continue
        if all(polynomial(value).change_ring(GF(p))(root) != 0 for value in generic for root in roots):
            character_blocks.append((p, roots))

    def character(value):
        beta = polynomial(value)
        bits = []
        for p, roots in character_blocks:
            reduced = beta.change_ring(GF(p))
            for root in roots:
                residue = reduced(root)
                if not residue:
                    return None
                bits.append(int(not residue.is_square()))
        return vector(GF(2), bits)

    generic_characters = matrix(GF(2), [character(value) for value in generic])
    assert generic_characters.rank() == 17

    targets = []
    for p0 in prime_range(3, policy["bounds"]["prime_search_bound"] + 1):
        p = int(p0)
        if p in finite_s or cubic.discriminant() % p == 0:
            continue
        primes = list(pari.idealprimedec(nf, p))
        for index, prime in enumerate(primes):
            if int(prime[3]) != 1:
                continue
            targets.append({"p": p, "prime_index": index, "hnf": str(pari.idealhnf(nf, prime)), "ideal": prime})
            break
        if len(targets) == policy["bounds"]["half_ideal_primes"]:
            break
    assert len(targets) == policy["bounds"]["half_ideal_primes"]
    return {
        "policy": policy,
        "data": data,
        "ring": ring,
        "cubic": cubic,
        "nf": nf,
        "pari": pari,
        "generic": generic,
        "generic_signatures": generic_signatures,
        "generic_characters": generic_characters,
        "signature": signature,
        "character": character,
        "polynomial": polynomial,
        "targets": targets,
    }


def candidate_record(context, mask, half, square, alpha):
    from sage.all import GF, QQ, vector, pari

    nf = context["nf"]
    assert pari.idealpow(nf, half, 2) == square
    assert pari.idealhnf(nf, alpha) == square
    norm = QQ(pari.nfeltnorm(nf, alpha))
    if norm < 0:
        alpha = -alpha
        norm = -norm
    assert norm > 0 and norm.is_square()
    coordinates = [QQ(pari.lift(alpha).polcoef(i)) for i in range(3)]
    # A rational square carries no new Kummer direction, regardless of its
    # half-ideal presentation.
    if coordinates[1:] == [0, 0] and coordinates[0].is_square():
        return None
    local = vector(GF(2), context["signature"](alpha))
    generic_locals = context["generic_signatures"]
    if local not in generic_locals.row_space():
        return None
    correction = generic_locals.transpose().solve_right(local)
    strict = alpha
    for index, bit in enumerate(correction):
        if bit:
            strict *= context["generic"][index]
    assert vector(GF(2), context["signature"](strict)).is_zero()
    char = context["character"](strict)
    if char is None or char in context["generic_characters"].row_space():
        return None
    return {
        "subset_mask": str(mask),
        "half_ideal_target_indices": [index for index in range(len(context["targets"])) if mask >> index & 1],
        "half_ideal_targets": [{key: value for key, value in context["targets"][index].items() if key != "ideal"}
                               for index in range(len(context["targets"])) if mask >> index & 1],
        "half_ideal_hnf": str(half),
        "square_ideal_hnf": str(square),
        "principal_generator_ascending": [str(value) for value in coordinates],
        "norm": str(norm),
        "norm_square_root": str(norm.sqrt()),
        "generic_correction_indices": [index for index, bit in enumerate(correction) if bit],
        "strict_representative_ascending": [str(value) for value in context["polynomial"](strict).list()],
        "strict_local_signature": [0] * len(local),
        "proof_character": [int(bit) for bit in char],
        "status": "STRICT_CLASS_CANDIDATE",
    }


def worker():
    WORK.mkdir(parents=True, exist_ok=True)
    context = setup()
    policy = context["policy"]
    targets = context["targets"]
    nf = context["nf"]
    count = len(targets)
    total = 1 << count
    page_width = policy["bounds"]["page_masks"]
    pari = context["pari"]
    unit = pari.idealhnf(nf, 1)
    # Split the subset product into two ten-prime tables.  This preserves
    # binary-mask order while avoiding a million large PARI ideal objects in
    # memory; it is a finite relation enumeration, not a class-group call.
    low_count = count // 2
    low_primes, high_primes = targets[:low_count], targets[low_count:]

    def subset_tables(primes):
        halves = [unit]
        squares = [unit]
        for mask in range(1, 1 << len(primes)):
            bit = mask & -mask
            index = bit.bit_length() - 1
            prior = mask ^ bit
            halves.append(pari.idealmul(nf, halves[prior], primes[index]["ideal"]))
            squares.append(pari.idealmul(nf, squares[prior], pari.idealpow(nf, primes[index]["ideal"], 2)))
        return halves, squares

    low_halves, low_squares = subset_tables(low_primes)
    high_halves, high_squares = subset_tables(high_primes)
    candidates = []
    mask = 0
    for high_mask in range(1 << len(high_primes)):
        for low_mask in range(1 << len(low_primes)):
            mask = (high_mask << low_count) | low_mask
            if not mask:
                continue
            square = pari.idealmul(nf, high_squares[high_mask], low_squares[low_mask])
            reduced, alpha = pari.idealred(nf, [square, 1])
            if reduced == unit:
                half = pari.idealmul(nf, high_halves[high_mask], low_halves[low_mask])
                record = candidate_record(context, mask, half, square, alpha)
                if record is not None:
                    candidates.append(record)
                    print("STRICT_CLASS_CANDIDATE", mask, flush=True)
                    break
            if mask % page_width == 0:
                page = mask // page_width - 1
                write_new(WORK / f"page_{page:04d}.json", {
                    "schema": "rank-jump.curve302-square-half-ideal-constructor-page.v1",
                    "mask_start": str(mask - page_width + 1), "mask_end": str(mask),
                    "tested": page_width, "strict_class_candidates": candidates,
                })
                print("PAGE", page, mask, flush=True)
        if candidates:
            break
    result = {
        "schema": "rank-jump.curve302-square-half-ideal-constructor.v1",
        "status": "PASS",
        "positive_endpoint": "STRICT_CLASS_CANDIDATE_PENDING_POSTHOC_M24_AND_COVER" if candidates else "NOT_REACHED",
        "bindings": bindings(),
        "targets": [{key: value for key, value in target.items() if key != "ideal"} for target in targets],
        "tested_subset_count": str(mask),
        "strict_class_candidates": candidates,
        "candidate_count": len(candidates),
        "boundary": policy["failure_semantics"],
    }
    write_new(OUTPUT, result)
    print(result["positive_endpoint"], result["tested_subset_count"], flush=True)


def capture():
    WORK.mkdir(parents=True, exist_ok=True)
    with (WORK / "worker.log").open("x") as log:
        try:
            process = subprocess.run([sys.executable, __file__, "worker"], stdout=log, stderr=log,
                                     timeout=r.read(PROTOCOL)["bounds"]["worker_seconds"])
            error = None if process.returncode == 0 else "worker failure"
        except subprocess.TimeoutExpired:
            error = "bounded timeout"
    if error and not OUTPUT.exists():
        write_new(OUTPUT, {"schema": "rank-jump.curve302-square-half-ideal-constructor.v1", "status": "UNKNOWN",
                           "reason": error, "bindings": bindings()})
    print(r.read(OUTPUT)["status"], flush=True)


def check():
    stored = r.read(OUTPUT)
    assert stored["bindings"] == bindings()
    assert stored["status"] == "PASS"
    if stored["strict_class_candidates"]:
        context = setup()
        assert stored["positive_endpoint"] == "STRICT_CLASS_CANDIDATE_PENDING_POSTHOC_M24_AND_COVER"
        for item in stored["strict_class_candidates"]:
            mask = int(item["subset_mask"])
            assert mask and len(item["proof_character"]) > 0
    print("PASS", stored["positive_endpoint"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("capture", "worker", "check"))
    args = parser.parse_args()
    if args.mode == "worker":
        worker()
    elif args.mode == "capture":
        capture()
    else:
        check()
