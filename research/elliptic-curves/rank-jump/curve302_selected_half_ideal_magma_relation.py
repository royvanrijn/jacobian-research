#!/usr/bin/env python3
"""Selected-prime algebraic relation lattice for the frozen E302 20-ideal set.

This is deliberately not another subset sweep.  It obtains class coordinates,
forms the integer lattice of all exponent vectors e with J(e)^2 principal,
and tries a deterministic LLL basis of that lattice.  The lattice permits the
even-exponent corrections which a binary subset calibration necessarily misses.
"""

from __future__ import annotations

import argparse
import importlib.util
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from hashlib import sha256
from pathlib import Path

import retrospective as r


HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "CURVE302_SELECTED_HALF_IDEAL_MAGMA_RELATION_PROTOCOL.json"
SELECTED = r.OUT / "rank_jump_curve302_square_half_ideal_constructor_v1.json"
BASE_PATH = HERE / "curve302_square_half_ideal_constructor.py"
OUTPUT = r.OUT / "rank_jump_curve302_selected_half_ideal_magma_relation_v1.json"
WORK = r.ROOT / "artifacts/local/rank-jump-curve302-selected-half-ideal-magma-relation-v1"


def base_module():
    spec = importlib.util.spec_from_file_location("curve302_selected_magma_base", BASE_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def bindings(base):
    paths = [Path(__file__), PROTOCOL, SELECTED, BASE_PATH, base.ARITH,
             r.ROOT / "elliptic-curves/cas/research_runtime/local_kummer.py"]
    return {str(path.relative_to(r.ROOT)): r.digest(path.read_bytes()) for path in paths}


def target_root(record):
    """Recover theta's residue from PARI's column-HNF of the frozen ideal."""
    values = list(map(int, re.findall(r"-?\d+", record["hnf"])))
    assert values[0] == record["p"]
    # The maximal-order basis starts 1, theta+2, ... .  In the HNF columns,
    # theta-r=(-2-r,1,0) lies in the ideal precisely for r=-2-a mod p.
    return (-2 - values[1]) % record["p"]


def setup(base):
    frozen = r.read(SELECTED)
    assert frozen["status"] == "PASS" and frozen["tested_subset_count"] == str(2**20 - 1)
    assert frozen["candidate_count"] == 0 and len(frozen["targets"]) == 20
    context = base.setup()
    assert len(context["targets"]) == 20
    targets = []
    for record in frozen["targets"]:
        root = target_root(record)
        assert context["cubic"](root) % record["p"] == 0
        matching = [prime for prime in context["pari"].idealprimedec(context["nf"], record["p"])
                    if str(context["pari"].idealhnf(context["nf"], prime)) == record["hnf"]]
        assert len(matching) == 1
        targets.append({**record, "root": root, "ideal": matching[0]})
    context["targets"] = targets
    context["policy"] = r.read(PROTOCOL)
    return context


def magma_program(context):
    terms = []
    for power, coefficient in enumerate(context["data"]["cubic_ascending"]):
        if coefficient:
            terms.append(f"({coefficient})*x^{power}")
    target_rows = ",".join(f"[{target['p']},{target['root']}]" for target in context["targets"])
    lines = [
        "SetColumns(0); SetSeed(1); Q:=Rationals(); R<x>:=PolynomialRing(Q);",
        "f:=" + " + ".join(terms) + "; K<a>:=NumberField(f); O:=MaximalOrder(K);",
        "A,mA:=ClassGroup(O); print \"CLASS_INVARIANTS\",Invariants(A);",
        "targets:=[" + target_rows + "];",
        "for data in targets do",
        "  p:=data[1]; rr:=data[2]; P:=ideal<O | p,a-rr>;",
        "  assert Norm(P) eq p; c:=P @@ mA;",
        "  print \"PRIME_COORD\",p,rr,Eltseq(c);",
        "end for;",
        "print \"DONE\";",
    ]
    return "\n".join(lines)


def remote(program, timeout):
    request = urllib.request.Request(
        "https://magma.maths.usyd.edu.au/xml/calculator.xml",
        data=urllib.parse.urlencode({"input": program}).encode(),
        headers={"Referer": "https://magma.maths.usyd.edu.au/calc/", "Content-Type": "application/x-www-form-urlencoded"},
    )
    started = time.monotonic()
    raw = urllib.request.urlopen(request, timeout=timeout).read().decode()
    root = ET.fromstring(raw)
    text = "\n".join("".join(node.itertext()) for node in root.findall(".//results/line"))
    return raw, text, time.monotonic() - started, len(root.findall(".//warning"))


def parse_list(source):
    assert re.fullmatch(r"\s*\[?\s*(?:-?\d+(?:\s*,\s*-?\d+)*)?\s*\]?\s*", source)
    return list(map(int, re.findall(r"-?\d+", source)))


def parse_coordinates(text, targets):
    match = re.search(r"^CLASS_INVARIANTS\s*(\[[^\n]*\])\s*$", text, re.MULTILINE)
    if match is None:
        raise ValueError("missing finite-abelian invariants")
    cyclic = parse_list(match.group(1))
    coordinates = []
    for target in targets:
        pattern = rf"^PRIME_COORD\s+{target['p']}\s+{target['root']}\s*(\[[^\n]*\])\s*$"
        match = re.search(pattern, text, re.MULTILINE)
        if match is None:
            raise ValueError(f"missing class coordinate for p={target['p']}")
        row = parse_list(match.group(1))
        if len(row) != len(cyclic):
            raise ValueError("class-coordinate width mismatch")
        coordinates.append([value % modulus for value, modulus in zip(row, cyclic)])
    if "DONE" not in text:
        raise ValueError("remote class-group run did not complete")
    return cyclic, coordinates


def relation_lattice(cyclic, coordinates):
    from sage.all import GF, ZZ, diagonal_matrix, matrix

    # 2c=0 in Z/dZ iff c=0 in Z/(d/gcd(d,2))Z.
    moduli = [modulus // (2 if modulus % 2 == 0 else 1) for modulus in cyclic]
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
    lll = lattice.LLL()
    return moduli, lattice, lll, matrix(GF(2), lattice).rank()


def ideal_from_exponents(context, exponents):
    pari, nf = context["pari"], context["nf"]
    answer = pari.idealhnf(nf, 1)
    for target, exponent in zip(context["targets"], exponents):
        if exponent:
            answer = pari.idealmul(nf, answer, pari.idealpow(nf, target["ideal"], int(exponent)))
    return answer


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
    coordinates = [QQ(pari.lift(alpha).polcoef(i)) for i in range(3)]
    return {
        "relation_exponents": [str(value) for value in exponents],
        "half_ideal_hnf": str(half),
        "square_ideal_hnf": str(square),
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
    program = magma_program(context)
    remote_record = {"schema": "rank-jump.curve302-selected-half-ideal-magma-transport.v1",
                     "program_sha256": sha256(program.encode()).hexdigest(), "status": "UNKNOWN"}
    try:
        raw, text, elapsed, warnings = remote(program, policy["bounds"]["remote_request_timeout_seconds"])
        remote_record.update({"status": "COMPLETED", "wall_seconds": elapsed, "warnings": warnings,
                              "xml_sha256": sha256(raw.encode()).hexdigest(), "text": text})
        cyclic, coordinates = parse_coordinates(text, context["targets"])
    except Exception as error:
        remote_record.update({"status": "TRANSPORT_OR_PARSE_FAILURE", "error": repr(error)})
        r.write_new(WORK / "remote.json", remote_record)
        print(remote_record["status"], flush=True)
        return
    r.write_new(WORK / "remote.json", remote_record)
    moduli, lattice, lll, f2_dimension = relation_lattice(cyclic, coordinates)
    candidates, principal_relations, rejected = [], [], []
    pari, nf = context["pari"], context["nf"]
    unit = pari.idealhnf(nf, 1)
    for vector in lll.rows():
        exponents = [int(value) for value in vector]
        if not any(value % 2 for value in exponents):
            rejected.append({"relation_exponents": [str(value) for value in exponents], "reason": "even parity"})
            continue
        half = ideal_from_exponents(context, exponents)
        square = pari.idealpow(nf, half, 2)
        reduced, alpha = pari.idealred(nf, [square, 1])
        if reduced != unit:
            raise ArithmeticError("class-coordinate relation did not principalize exactly")
        assert pari.idealmul(nf, reduced, alpha) == square
        principal_relations.append([str(value) for value in exponents])
        record = candidate_record(context, exponents, half, square, alpha)
        if record is None:
            rejected.append({"relation_exponents": [str(value) for value in exponents], "reason": "not strict/generic-new"})
        else:
            candidates.append(record)
            print("STRICT_CLASS_CANDIDATE", flush=True)
            break
    result = {
        "schema": "rank-jump.curve302-selected-half-ideal-magma-relation.v1",
        "status": "PASS",
        "positive_endpoint": "STRICT_CLASS_CANDIDATE_PENDING_POSTHOC_M24_AND_COVER" if candidates else "NOT_REACHED",
        "bindings": bindings(base),
        "external_class_coordinate_status": "UNVERIFIED_EXTERNAL_OUTPUT_EXACTLY_REPLAYED_FOR_TESTED_RELATIONS",
        "class_group_cyclic": cyclic,
        "selected_prime_class_coordinates": coordinates,
        "principal_square_moduli": moduli,
        "relation_lattice_basis": [[str(value) for value in row] for row in lattice.rows()],
        "lll_relation_basis": [[str(value) for value in row] for row in lll.rows()],
        "f2_principal_square_relation_dimension": f2_dimension,
        "tested_relation_representative_count": len(principal_relations) + len(rejected),
        "exact_principal_relation_exponents": principal_relations,
        "strict_class_candidates": candidates,
        "rejected_relation_representatives": rejected,
        "boundary": policy["failure_semantics"],
    }
    r.write_new(OUTPUT, result)
    print(result["status"], result["positive_endpoint"], flush=True)


def check():
    base = base_module()
    stored = r.read(OUTPUT)
    assert stored["bindings"] == bindings(base)
    assert stored["status"] == "PASS"
    assert stored["f2_principal_square_relation_dimension"] >= 0
    print("PASS", stored["positive_endpoint"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("worker", "check"))
    args = parser.parse_args()
    globals()[args.mode]()
