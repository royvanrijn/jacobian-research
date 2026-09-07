#!/usr/bin/env python3
"""Short reduced-ideal relation lattice for the frozen E302 20-prime family."""

from __future__ import annotations

import argparse
from itertools import combinations, product
import importlib.util
from pathlib import Path
import subprocess
import sys
import time

import retrospective as r


HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "CURVE302_SELECTED_HALF_IDEAL_REDUCED_RELATION_PROTOCOL.json"
SELECTED = r.OUT / "rank_jump_curve302_square_half_ideal_constructor_v1.json"
BASE_PATH = HERE / "curve302_square_half_ideal_constructor.py"
OUTPUT = r.OUT / "rank_jump_curve302_selected_half_ideal_reduced_relation_v1.json"
WORK = r.ROOT / "artifacts/local/rank-jump-curve302-selected-half-ideal-reduced-relation-v1"


def base_module():
    spec = importlib.util.spec_from_file_location("curve302_reduced_relation_base", BASE_PATH)
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
        matches = [prime for prime in context["pari"].idealprimedec(context["nf"], record["p"])
                   if str(context["pari"].idealhnf(context["nf"], prime)) == record["hnf"]]
        assert len(matches) == 1
        targets.append({**record, "ideal": matches[0]})
    context["targets"] = targets
    context["policy"] = r.read(PROTOCOL)
    return context


def compositions(total, width):
    if width == 1:
        yield (total,)
    else:
        for first in range(1, total - width + 2):
            for rest in compositions(total - first, width - 1):
                yield (first,) + rest


def short_states(width, radius):
    for norm in range(1, radius + 1):
        for support_size in range(1, min(norm, width) + 1):
            for support in combinations(range(width), support_size):
                for magnitudes in compositions(norm, support_size):
                    for signs in product((-1, 1), repeat=support_size):
                        row = [0] * width
                        for index, magnitude, sign in zip(support, magnitudes, signs):
                            row[index] = sign * magnitude
                        yield tuple(row)


def canonical_sign(row):
    row = tuple(int(value) for value in row)
    for value in row:
        if value:
            return row if value > 0 else tuple(-entry for entry in row)
    raise ValueError("zero relation")


def parity_mask(row):
    return sum((int(value) & 1) << index for index, value in enumerate(row))


def ideal_from_exponents(context, exponents):
    pari, nf = context["pari"], context["nf"]
    answer = pari.idealhnf(nf, 1)
    for target, exponent in zip(context["targets"], exponents):
        if exponent:
            answer = pari.idealmul(nf, answer, pari.idealpow(nf, target["ideal"], int(exponent)))
    return answer


def reduced_key(context, ideal):
    pari, nf = context["pari"], context["nf"]
    reduced, multiplier = pari.idealred(nf, [ideal, 1])
    assert pari.idealmul(nf, reduced, multiplier) == ideal
    return str(reduced)


def exact_principalize(context, ideal):
    pari, nf = context["pari"], context["nf"]
    unit = pari.idealhnf(nf, 1)
    reduced, generator = pari.idealred(nf, [ideal, 1])
    if reduced != unit:
        raise ArithmeticError("reduced-ideal collision did not principalize exactly")
    assert pari.idealmul(nf, reduced, generator) == ideal
    assert pari.idealhnf(nf, generator) == ideal
    return generator


def element_coordinates(context, element):
    from sage.all import QQ
    return [QQ(context["pari"].lift(element).polcoef(index)) for index in range(3)]


def score(context, exponents, generator):
    coordinates = element_coordinates(context, generator)
    return [
        sum(value != 0 for value in exponents),
        sum(abs(value) for value in exponents),
        max(abs(value) for value in exponents),
        max(len(str(abs(value.numerator()))) + len(str(value.denominator())) for value in coordinates),
    ]


def relation_basis(rows):
    from sage.all import ZZ, matrix
    if not rows:
        return matrix(ZZ, 0, 20)
    return matrix(ZZ, rows).row_module().basis_matrix()


def half_relation_lattice(basis):
    """Return {e : 2e lies in the direct relation lattice}, reduced by LLL."""
    from sage.all import ZZ, identity_matrix, matrix
    if basis.nrows() == 0:
        return matrix(ZZ, 0, 20), matrix(ZZ, 0, 20)
    rank = basis.nrows()
    equation = basis.transpose().augment(-2 * identity_matrix(ZZ, 20))
    kernel = equation.right_kernel()
    coefficients = matrix(ZZ, [list(vector)[:rank] for vector in kernel.basis()])
    halves = matrix(ZZ, [[sum(coefficients[row, index] * basis[index, column]
                                    for index in range(rank)) // 2
                                for column in range(20)]
                               for row in range(coefficients.nrows())])
    assert all(all(2 * value == sum(coefficients[row, index] * basis[index, column]
                                    for index in range(rank))
                   for column, value in enumerate(halves.row(row)))
               for row in range(halves.nrows()))
    return halves, halves.LLL()


def short_lattice_candidates(lll_basis):
    rows = [tuple(map(int, row)) for row in lll_basis.rows()]
    answers = set(rows)
    for left, right in combinations(rows, 2):
        answers.add(tuple(a + b for a, b in zip(left, right)))
        answers.add(tuple(a - b for a, b in zip(left, right)))
    answers.discard((0,) * 20)
    return sorted((row for row in answers if parity_mask(row)), key=lambda row: (sum(value != 0 for value in row),
                                                                                  sum(abs(value) for value in row),
                                                                                  max(abs(value) for value in row), row))


def strict_candidate(context, exponents, half, square, alpha):
    from sage.all import GF, QQ, vector
    pari, nf = context["pari"], context["nf"]
    assert pari.idealpow(nf, half, 2) == square and pari.idealhnf(nf, alpha) == square
    norm = QQ(pari.nfeltnorm(nf, alpha))
    if norm < 0:
        alpha, norm = -alpha, -norm
    if not (norm > 0 and norm.is_square()):
        return None
    local = vector(GF(2), context["signature"](alpha))
    generic = context["generic_signatures"]
    if local not in generic.row_space():
        return None
    correction = generic.transpose().solve_right(local)
    strict = alpha
    for index, bit in enumerate(correction):
        if bit:
            strict *= context["generic"][index]
    assert vector(GF(2), context["signature"](strict)).is_zero()
    character = context["character"](strict)
    if character is None or character in context["generic_characters"].row_space():
        return None
    return {
        "relation_exponents": [str(value) for value in exponents], "parity_mask_hex": hex(parity_mask(exponents)),
        "half_ideal_hnf": str(half), "square_ideal_hnf": str(square),
        "principal_generator_ascending": [str(value) for value in element_coordinates(context, alpha)],
        "norm": str(norm), "norm_square_root": str(norm.sqrt()),
        "generic_correction_indices": [index for index, bit in enumerate(correction) if bit],
        "strict_representative_ascending": [str(value) for value in context["polynomial"](strict).list()],
        "strict_local_signature": [0] * len(local), "proof_character": [int(bit) for bit in character],
        "score": score(context, exponents, alpha), "status": "STRICT_CLASS_CANDIDATE",
    }


def worker():
    base = base_module()
    context = setup(base)
    policy = context["policy"]
    WORK.mkdir(parents=True, exist_ok=True)
    states = list(short_states(20, policy["bounds"]["signed_l1_norm"]))
    assert len(states) == 11520
    unit = context["pari"].idealhnf(context["nf"], 1)
    seen, relation_rows, relation_records = {str(unit): (0,) * 20}, set(), []
    for count, state in enumerate(states, 1):
        ideal = ideal_from_exponents(context, state)
        key = reduced_key(context, ideal)
        earlier = seen.get(key)
        if earlier is None:
            seen[key] = state
        else:
            relation = canonical_sign(tuple(left - right for left, right in zip(state, earlier)))
            if relation not in relation_rows:
                generator = exact_principalize(context, ideal_from_exponents(context, relation))
                relation_rows.add(relation)
                relation_records.append({"relation_exponents": [str(value) for value in relation],
                    "parity_mask_hex": hex(parity_mask(relation)),
                    "principal_generator_ascending": [str(value) for value in element_coordinates(context, generator)],
                    "score": score(context, relation, generator)})
        if count % 512 == 0:
            r.write_new(WORK / f"page_{count // 512:03d}.json", {"schema": "rank-jump.curve302-reduced-relation-page.v1",
                "states_completed": count, "distinct_reduced_ideals": len(seen), "exact_direct_relations": len(relation_rows)})
            print("PAGE", count, len(relation_rows), flush=True)
    basis = relation_basis(sorted(relation_rows))
    half_lattice, half_lll = half_relation_lattice(basis)
    square_records, rejected, candidates = [], [], []
    for exponents in short_lattice_candidates(half_lll):
        half = ideal_from_exponents(context, exponents)
        square = context["pari"].idealpow(context["nf"], half, 2)
        alpha = exact_principalize(context, square)
        square_records.append({"relation_exponents": [str(value) for value in exponents],
            "parity_mask_hex": hex(parity_mask(exponents)), "score": score(context, exponents, alpha)})
        candidate = strict_candidate(context, exponents, half, square, alpha)
        if candidate is None:
            rejected.append({"relation_exponents": [str(value) for value in exponents], "reason": "not strict/generic-new"})
        else:
            candidates.append(candidate)
    candidates.sort(key=lambda row: row["score"])
    result = {
        "schema": "rank-jump.curve302-selected-half-ideal-reduced-relation.v1", "status": "PASS",
        "positive_endpoint": "STRICT_CLASS_CANDIDATE_PENDING_POSTHOC_M24_AND_COVER" if candidates else "NOT_REACHED",
        "bindings": bindings(base), "state_count": len(states), "exact_direct_relation_count": len(relation_rows),
        "exact_direct_relations": sorted(relation_records, key=lambda row: row["score"]),
        "direct_relation_hnf_basis": [[str(value) for value in row] for row in basis.rows()],
        "direct_relation_lll_basis": [[str(value) for value in row] for row in basis.LLL().rows()],
        "half_relation_lattice_basis": [[str(value) for value in row] for row in half_lattice.rows()],
        "half_relation_lll_basis": [[str(value) for value in row] for row in half_lll.rows()],
        "exact_principal_square_candidates": square_records,
        "strict_class_candidates": candidates, "preferred_strict_class_candidate": candidates[0] if candidates else None,
        "rejected_square_candidates": rejected, "boundary": policy["failure_semantics"],
    }
    r.write_new(OUTPUT, result)
    print(result["status"], result["positive_endpoint"], flush=True)


def launch():
    WORK.mkdir(parents=True, exist_ok=True)
    log, record = WORK / "worker.log", WORK / "launch.json"
    if log.exists() or record.exists():
        raise FileExistsError("prior detached launch already recorded")
    with log.open("x") as stream:
        process = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), "worker"], cwd=Path.cwd(),
                                   stdout=stream, stderr=stream, start_new_session=True)
    r.write_new(record, {"schema": "rank-jump.curve302-reduced-relation-launch.v1", "status": "STARTED",
        "pid": process.pid, "started_at_unix": int(time.time()), "scope": "same frozen 20 ideals; no prime expansion"})
    print("STARTED", process.pid, flush=True)


def check():
    base = base_module()
    stored = r.read(OUTPUT)
    assert stored["bindings"] == bindings(base) and stored["status"] == "PASS"
    print("PASS", stored["positive_endpoint"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("worker", "launch", "check"))
    args = parser.parse_args()
    globals()[args.mode]()
