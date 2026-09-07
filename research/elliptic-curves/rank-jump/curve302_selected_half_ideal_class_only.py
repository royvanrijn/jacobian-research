#!/usr/bin/env python3
"""Class-only selected-ideal relation lattice, with no unit/regulator request."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import subprocess
import sys
import time

import retrospective as r


HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "CURVE302_SELECTED_HALF_IDEAL_CLASS_ONLY_PROTOCOL.json"
SELECTED = r.OUT / "rank_jump_curve302_square_half_ideal_constructor_v1.json"
BASE_PATH = HERE / "curve302_square_half_ideal_constructor.py"
OUTPUT = r.OUT / "rank_jump_curve302_selected_half_ideal_class_only_v1.json"
WORK = r.ROOT / "artifacts/local/rank-jump-curve302-selected-half-ideal-class-only-v1"


def base_module():
    spec = importlib.util.spec_from_file_location("curve302_selected_class_only_base", BASE_PATH)
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
    targets = []
    for record in frozen["targets"]:
        matching = [prime for prime in context["pari"].idealprimedec(context["nf"], record["p"])
                    if str(context["pari"].idealhnf(context["nf"], prime)) == record["hnf"]]
        assert len(matching) == 1
        targets.append({**record, "ideal": matching[0]})
    context["targets"] = targets
    context["policy"] = r.read(PROTOCOL)
    return context


def relation_lattice(moduli, coordinates):
    from sage.all import GF, ZZ, diagonal_matrix, matrix

    if not moduli:
        lattice = matrix(ZZ, len(coordinates), len(coordinates), 1)
    else:
        class_matrix = matrix(ZZ, coordinates).transpose()
        equation = class_matrix.augment(-diagonal_matrix(ZZ, moduli))
        kernel = equation.right_kernel()
        lattice = matrix(ZZ, [list(vector)[:len(coordinates)] for vector in kernel.basis()])
        assert lattice.nrows() == len(coordinates) and lattice.rank() == len(coordinates)
        for vector in lattice.rows():
            for row, modulus in zip(class_matrix.rows(), moduli):
                assert sum(entry * coefficient for entry, coefficient in zip(row, vector)) % modulus == 0
    return lattice, lattice.LLL(), matrix(GF(2), lattice).rank()


def ideal_from_exponents(context, exponents):
    pari, nf = context["pari"], context["nf"]
    answer = pari.idealhnf(nf, 1)
    for target, exponent in zip(context["targets"], exponents):
        if exponent:
            answer = pari.idealmul(nf, answer, pari.idealpow(nf, target["ideal"], int(exponent)))
    return answer


def parity_mask(exponents):
    return sum((int(value) & 1) << index for index, value in enumerate(exponents))


def element_coordinates(context, element):
    from sage.all import QQ
    return [QQ(context["pari"].lift(element).polcoef(index)) for index in range(3)]


def relation_score(exponents, generator):
    coordinates = element_coordinates(generator["context"], generator["element"])
    return [
        sum(value != 0 for value in exponents),
        sum(abs(value) for value in exponents),
        max(abs(value) for value in exponents),
        max(len(str(abs(value.numerator()))) + len(str(value.denominator())) for value in coordinates),
    ]


def exact_principalize(context, ideal):
    pari, nf = context["pari"], context["nf"]
    unit = pari.idealhnf(nf, 1)
    reduced, generator = pari.idealred(nf, [ideal, 1])
    if reduced != unit:
        raise ArithmeticError("class-coordinate relation did not principalize exactly")
    assert pari.idealmul(nf, reduced, generator) == ideal
    assert pari.idealhnf(nf, generator) == ideal
    return generator


def candidate_record(context, exponents, half, square, alpha):
    from sage.all import GF, QQ, vector

    pari, nf = context["pari"], context["nf"]
    assert pari.idealpow(nf, half, 2) == square
    assert pari.idealhnf(nf, alpha) == square
    norm = QQ(pari.nfeltnorm(nf, alpha))
    if norm < 0:
        alpha, norm = -alpha, -norm
    assert norm > 0 and norm.is_square()
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
    coordinates = element_coordinates(context, alpha)
    return {
        "relation_exponents": [str(value) for value in exponents],
        "half_ideal_hnf": str(half), "square_ideal_hnf": str(square),
        "principal_generator_ascending": [str(value) for value in coordinates],
        "norm": str(norm), "norm_square_root": str(norm.sqrt()),
        "generic_correction_indices": [index for index, bit in enumerate(correction) if bit],
        "strict_representative_ascending": [str(value) for value in context["polynomial"](strict).list()],
        "strict_local_signature": [0] * len(local),
        "proof_character": [int(bit) for bit in char],
        "status": "STRICT_CLASS_CANDIDATE",
    }


def worker():
    base = base_module()
    context = setup(base)
    policy = context["policy"]
    WORK.mkdir(parents=True, exist_ok=True)
    pari, nf = context["pari"], context["nf"]
    pari.allocatemem(64000000, policy["bounds"]["pari_stack_bytes"], silent=True)
    print("CLASS_ONLY_BNF_START", flush=True)
    bnf = pari.bnfinit(nf, 0, pari(policy["bounds"]["compact_relation_tech"]))
    cyclic = list(map(int, bnf.bnf_get_cyc()))
    coordinates = [list(map(int, pari.bnfisprincipal(bnf, target["ideal"], 0))) for target in context["targets"]]
    assert all(len(row) == len(cyclic) for row in coordinates)
    print("CLASS_ONLY_BNF_RETURNED", cyclic, flush=True)
    direct_lattice, direct_lll, direct_f2_dimension = relation_lattice(cyclic, coordinates)
    square_moduli = [modulus // (2 if modulus % 2 == 0 else 1) for modulus in cyclic]
    square_lattice, square_lll, square_f2_dimension = relation_lattice(square_moduli, coordinates)
    direct_relations, direct_rejected = [], []
    for vector in direct_lll.rows():
        exponents = [int(value) for value in vector]
        mask = parity_mask(exponents)
        if not mask:
            direct_rejected.append({"relation_exponents": [str(value) for value in exponents], "reason": "even parity"})
            continue
        ideal = ideal_from_exponents(context, exponents)
        generator = exact_principalize(context, ideal)
        direct_relations.append({
            "relation_exponents": [str(value) for value in exponents],
            "parity_mask_hex": hex(mask), "ideal_hnf": str(ideal),
            "principal_generator_ascending": [str(value) for value in element_coordinates(context, generator)],
            "score": relation_score(exponents, {"context": context, "element": generator}),
        })
    candidates, principal_square_relations, rejected = [], [], []
    for vector in square_lll.rows():
        exponents = [int(value) for value in vector]
        mask = parity_mask(exponents)
        if not mask:
            rejected.append({"relation_exponents": [str(value) for value in exponents], "reason": "even parity"})
            continue
        half = ideal_from_exponents(context, exponents)
        square = pari.idealpow(nf, half, 2)
        alpha = exact_principalize(context, square)
        principal_square_relations.append([str(value) for value in exponents])
        record = candidate_record(context, exponents, half, square, alpha)
        if record is None:
            rejected.append({"relation_exponents": [str(value) for value in exponents], "reason": "not strict/generic-new"})
        else:
            record["parity_mask_hex"] = hex(mask)
            record["score"] = relation_score(exponents, {"context": context, "element": alpha})
            candidates.append(record)
    candidates.sort(key=lambda row: (row["parity_mask_hex"], row["score"]))
    preferred = min(candidates, key=lambda row: row["score"]) if candidates else None
    result = {
        "schema": "rank-jump.curve302-selected-half-ideal-class-only.v1",
        "status": "RELATION_GUIDE_WITH_EXACT_REPLAYS",
        "positive_endpoint": "STRICT_CLASS_CANDIDATE_PENDING_POSTHOC_M24_AND_COVER" if preferred else "NOT_REACHED",
        "bindings": bindings(base), "class_group_cyclic": cyclic,
        "selected_prime_class_coordinates": coordinates,
        "principal_relation_lattice_basis": [[str(value) for value in row] for row in direct_lattice.rows()],
        "principal_relation_lll_basis": [[str(value) for value in row] for row in direct_lll.rows()],
        "f2_principal_relation_dimension": direct_f2_dimension,
        "exact_principal_relations": direct_relations,
        "square_relation_moduli": square_moduli,
        "principal_square_lattice_basis": [[str(value) for value in row] for row in square_lattice.rows()],
        "principal_square_lll_basis": [[str(value) for value in row] for row in square_lll.rows()],
        "f2_principal_square_relation_dimension": square_f2_dimension,
        "exact_principal_square_relation_exponents": principal_square_relations,
        "strict_class_candidates": candidates, "preferred_strict_class_candidate": preferred,
        "rejected_direct_relation_representatives": direct_rejected,
        "rejected_square_relation_representatives": rejected,
        "boundary": policy["failure_semantics"],
    }
    r.write_new(OUTPUT, result)
    print(result["status"], result["positive_endpoint"], flush=True)


def check():
    base = base_module()
    stored = r.read(OUTPUT)
    assert stored["bindings"] == bindings(base)
    assert stored["status"] == "RELATION_GUIDE_WITH_EXACT_REPLAYS"
    print("PASS", stored["positive_endpoint"], flush=True)


def launch():
    """Start the bounded-scope worker in a separate session and return at once."""
    WORK.mkdir(parents=True, exist_ok=True)
    launch_record = WORK / "launch.json"
    log = WORK / "worker.log"
    if launch_record.exists() or log.exists():
        raise FileExistsError("a prior detached launch is already recorded")
    with log.open("x") as stream:
        process = subprocess.Popen(
            [sys.executable, str(Path(__file__).resolve()), "worker"],
            cwd=Path.cwd(), stdout=stream, stderr=stream, start_new_session=True,
        )
    r.write_new(launch_record, {
        "schema": "rank-jump.curve302-selected-half-ideal-class-only-launch.v1",
        "status": "STARTED",
        "pid": process.pid,
        "started_at_unix": int(time.time()),
        "worker": str(Path(__file__).resolve()),
        "scope": "frozen selected 20 prime ideals; no prime expansion",
    })
    print("STARTED", process.pid, flush=True)


def record_timeout():
    WORK.mkdir(parents=True, exist_ok=True)
    r.write_new(WORK / "class_only_bnf_timeout.json", {
        "schema": "rank-jump.curve302-selected-half-ideal-class-only-timeout.v1",
        "status": "TIMEOUT_UNKNOWN",
        "selected_half_ideals": 20,
        "wall_seconds_lower_bound": 1800,
        "reason": "No class coordinates or relation telemetry were emitted before the fixed review threshold.",
    })
    print("RECORDED_TIMEOUT", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("worker", "check", "launch", "record-timeout"))
    args = parser.parse_args()
    if args.mode == "record-timeout":
        record_timeout()
    else:
        globals()[args.mode]()
