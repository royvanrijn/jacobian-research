#!/usr/bin/env python3
"""Check two integer A8 ideal identities by exact conjugation orbit sums."""
from collections import Counter
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path

from verify_uniform_fox_ideal_models import convolution, reconstruct

HERE = Path(__file__).resolve().parent
CERTIFICATE = HERE / "integral_a8_fox_ideal_certificates.json"
ONE = bytes(range(8))
BYTE_TAIL = bytes(range(8, 256))
SUBGROUP_ORDERS = {0: 20160, 3: 64, 5: 1344, 7: 720}


def product(left, right):
    # Apply left to each of the values of right.
    return right.translate(left + BYTE_TAIL)


def inverse(p):
    return bytes(p.index(i) for i in range(8))


def signed_fox_rows(word, top):
    """Evaluate both signs, also checking integral Fox's fundamental identity."""
    result = {}
    for sign in (1, -1):
        row = [Counter() for _ in top]
        prefix = ONE
        scalar = 1
        for x in word:
            i = abs(x) - 1
            if x > 0:
                row[i][prefix] += scalar
            prefix = product(prefix, top[i])
            scalar *= sign
            if x < 0:
                row[i][prefix] -= scalar
        assert prefix == ONE and scalar == 1
        fundamental = Counter()
        for p, derivative in zip(top, row):
            for g, c in derivative.items():
                fundamental[product(g, p)] += sign * c
                fundamental[g] -= c
        assert not any(fundamental.values())
        result[sign] = [{g: c for g, c in derivative.items() if c} for derivative in row]
    return result


def verify_parity_surjection(top, group):
    seen = {(ONE, 0)}
    queue = [(ONE, 0)]
    for g, parity in queue:
        for p in top:
            value = (product(g, p), parity ^ 1)
            if value not in seen:
                seen.add(value)
                queue.append(value)
    assert len(seen) == 2 * len(group) == 40320


def orbit_map(generators, group):
    assert generators and all(g in group for g in generators)
    subgroup = {ONE}
    queue = [ONE]
    for g in queue:
        for h in generators:
            value = product(g, h)
            if value not in subgroup:
                subgroup.add(value)
                queue.append(value)
    moves = [(g, inverse(g)) for g in generators]
    orbit = {}
    sizes = []
    for g in sorted(group):
        if g in orbit:
            continue
        i = len(sizes)
        orbit[g] = i
        queue = [g]
        for h in queue:
            for k, inv in moves:
                value = product(product(k, h), inv)
                if value not in orbit:
                    orbit[value] = i
                    queue.append(value)
                assert orbit[value] == i
        sizes.append(len(queue))
    assert orbit[ONE] == 0 and sizes[0] == 1
    assert sum(sizes) == len(group) and all(len(subgroup) % s == 0 for s in sizes)
    return len(subgroup), orbit, sizes


def verify_average(record, row, group):
    p = record["modulus"]
    assert p in SUBGROUP_ORDERS
    order, orbit, sizes = orbit_map(list(map(bytes, record["subgroup_generators"])), group)
    assert order == record["subgroup_order"] == SUBGROUP_ORDERS[p]
    assert len(sizes) == record["orbit_count"]
    if p:
        assert order % p
    total = [0] * len(sizes)
    for multiplier, index, left, right in record["terms"]:
        a, b = bytes(left), bytes(right)
        assert type(multiplier) is int and multiplier != 0
        assert type(index) is int and 0 <= index < 4 and a in group and b in group
        for value, c in row[index].items():
            g = product(product(a, value), b)
            total[orbit[g]] += multiplier * c
    # Coefficients of the conjugation sum on an orbit O are
    # the orbit total multiplied by |H|/|O|.
    assert ([x % p for x in total] if p else total) == [1] + [0] * (len(sizes) - 1)
    return {"meridian_sign": record["meridian_sign"], "modulus": p, "subgroup_order": order,
            "orbit_count": len(sizes), "translate_terms": len(record["terms"]),
            "averaged_identity_coefficient": order % p if p else order}


def verify():
    data = json.loads(CERTIFICATE.read_text())
    for name, digest in data["sources"].items():
        assert sha256((HERE / name).read_bytes()).hexdigest() == digest, name
    literal, original, group_tuples, _ = reconstruct()
    group = set(map(bytes, group_tuples))
    top = list(map(bytes, literal["ordered_A8_top"]))
    verify_parity_surjection(top, group)
    word = literal["all_four_closure_relators"][0]
    rows = signed_fox_rows(word, top)
    assert rows[1] == [{bytes(g): c for g, c in row.items()} for row in original[:4]]
    assert all((rows[1][i].get(g, 0) - rows[-1][i].get(g, 0)) % 2 == 0
               for i in range(4) for g in rows[1][i].keys() | rows[-1][i].keys())
    char2_name = "uniform_fox_a8_char2_certificate.json"
    cert2 = json.loads((HERE / char2_name).read_text())
    assert cert2["derivative_pairs"] == [[0, 0], [0, 1]]
    binary_supports = convolution(cert2, original, group_tuples)
    records = data["certificates"]
    assert [(r["meridian_sign"], r["modulus"]) for r in records] == [
        (sign, p) for sign in (1, -1) for p in SUBGROUP_ORDERS]
    checks = [verify_average(r, rows[r["meridian_sign"]], group) for r in records]
    rejected = []
    for label, row in (("changed_integral_coefficient", records[0]),
                       ("changed_modular_coefficient", records[1])):
        corrupt = deepcopy(row)
        corrupt["terms"][0][0] += 1
        try:
            verify_average(corrupt, rows[1], group)
        except AssertionError:
            rejected.append(label)
        else:
            raise AssertionError(f"Accepted {label}")
    try:
        verify_average(records[4], rows[1], group)
    except AssertionError:
        rejected.append("wrong_meridian_sign")
    else:
        raise AssertionError("Accepted a certificate with the wrong sign")
    return {"status": "PASS integral fixed A8 and A8 x C2 first-row Fox ideal certificates",
            "certificate_sha256": sha256(CERTIFICATE.read_bytes()).hexdigest(),
            "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "literal_source_sha256": data["sources"],
            "parity_quotient_order": 40320, "binary_supports": binary_supports,
            "orbit_sum_identities": checks, "rejected_controls": rejected,
            "integer_descent": "20160 annihilates each sign quotient. Multiplication by each of 2, 3, 5, 7 is surjective, so multiplication by 20160 is surjective and both quotients vanish.",
            "parity_descent": "Lift 1 from each sign ideal; the product of the two complementary lifts lies in (z-1)(z+1)=0.",
            "assurance": {"whole_proof_independent_replay": False, "formal_verification": False,
                          "external_specialist_review": False}}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
