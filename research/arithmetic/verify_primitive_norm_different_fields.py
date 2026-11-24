#!/usr/bin/env python3
"""Exact bounded controls for different-field primitive norm class errors.

Finite residue spaces, finite mark probabilities and finite filter bijections
only. Analytic bounds and universal quantifiers are in the written proof.
All cases are deterministic below; no generated input file is required.
"""
from argparse import ArgumentParser
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
from math import gcd, prod
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
PRIMES = (3, 5, 7, 11, 13)
PROBABILITIES = (Q(0), Q(1, 2), Q(9, 10), Q(999, 1000), Q(1))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def local_survivors(p, target, inert1, inert2, marked_slot=None):
    result = set()
    for a in range(p * p):
        b = (target - a) % (p * p)
        if (inert1 and a % p == 0) or (inert2 and b % p == 0):
            continue
        if marked_slot is not None:
            value = a if marked_slot == 0 else b
            if value % p == 0 and value != 0:
                continue
        result.add(a)
    return result


def expected_increment(p, target, opposite_inert):
    if not opposite_inert:
        return Q(p - 1, p * p)
    return Q(0) if target % p == 0 else Q(1, p)


def local_controls():
    baseline_cases = marked_cases = 0
    maximal_root_counts = {}
    for p in PRIMES:
        max_roots = 0
        for target in range(p * p):
            for inert1, inert2 in product((False, True), repeat=2):
                base = local_survivors(p, target, inert1, inert2)
                rho = int(inert1) + int(inert2) - int(inert1 and inert2 and target % p == 0)
                require(Q(len(base), p * p) == 1 - Q(rho, p),
                        "Baseline root count mismatch.")
                generic = 1 - Q(int(inert1) + int(inert2), p)
                target_factor = Q(p - 1, p - 2) if inert1 and inert2 and target % p == 0 else Q(1)
                require(Q(len(base), p * p) == generic * target_factor,
                        "Common-inert target correction mismatch.")
                baseline_cases += 1
                for slot, own_inert, other_inert in ((0, inert1, inert2), (1, inert2, inert1)):
                    if own_inert:
                        continue
                    marked = local_survivors(p, target, inert1, inert2, slot)
                    increment = Q(len(base) - len(marked), len(base))
                    require(increment == expected_increment(p, target, other_inert),
                            "Incremental marking overlap mismatch.")
                    roots = p * p - len(marked)
                    require(roots <= 2 * p - 1 and Q(roots, p * p) <= Q(2, p),
                            "Combined prime-square root budget fails.")
                    max_roots = max(max_roots, roots)
                    marked_cases += 1
        require(max_roots == 2 * p - 1, "The larger root count was not attained.")
        maximal_root_counts[str(p)] = max_roots
    return dict(baseline_cases=baseline_cases, marked_cases=marked_cases,
                complete_target_classes_mod_p_squared=True,
                maximum_prime_square_forbidden_counts=maximal_root_counts)


def mark_product_controls():
    primes = PRIMES[:4]
    subsets = tuple(product((False, True), repeat=len(primes)))
    cases = 0
    for target in (1, 3, 15, 105):
        for opposite_inert in subsets:
            baseline = [Q(len(local_survivors(p, target, False, inert)), p * p)
                        for p, inert in zip(primes, opposite_inert)]
            increments = [expected_increment(p, target, inert)
                          for p, inert in zip(primes, opposite_inert)]
            target_weight = prod(Q(p, p - 1) for p, inert in zip(primes, opposite_inert)
                                 if inert and target % p == 0)
            for t in PROBABILITIES:
                expectation = Q(0)
                for selected in subsets:
                    weight = t ** sum(selected) * (1 - t) ** (len(primes) - sum(selected))
                    survived = prod(base * (1 - inc if mark else 1)
                                    for base, inc, mark in zip(baseline, increments, selected))
                    expectation += weight * survived
                expected = prod(baseline) * prod(1 - t * b for b in increments)
                require(expectation == expected, "Finite marking expectation fails.")
                reference = prod(1 - t * Q(p - 1, p * p) for p in primes)
                require(prod(1 - t * b for b in increments) <= reference * target_weight,
                        "Exact target-factor upper bound fails.")
                cases += 1
    return dict(primes=list(primes), targets=[1, 3, 15, 105],
                probabilities=[str(t) for t in PROBABILITIES], cases=cases,
                mark_sets_each=len(subsets), opposite_inert_patterns=len(subsets))


def crt_control():
    primes, target = (3, 5), 1
    forbidden = [{a for a in range(p * p)
                  if a not in local_survivors(p, target, False, True, 0)} for p in primes]
    modulus = prod(p * p for p in primes)
    roots = [a for a in range(modulus)
             if all(a % (p * p) in local for p, local in zip(primes, forbidden))]
    require(len(roots) == prod(len(s) for s in forbidden) == 45, "CRT root count fails.")
    require(len(roots) <= 2 ** len(primes) * prod(primes), "Corrected CRT bound fails.")
    return dict(primes=list(primes), modulus=modulus, forbidden_root_count=len(roots),
                old_bound=prod(primes), corrected_bound=2 ** len(primes) * prod(primes))


def subsets(values):
    values = sorted(values)
    return [frozenset(x for x, flag in zip(values, flags) if flag)
            for flags in product((False, True), repeat=len(values))]


def filter_controls():
    # Two distinct imaginary quadratic characters: -4 and -8 modulo 8.
    L, H1, H2 = 8, frozenset((1, 5)), frozenset((1, 3))
    normalization_cases = transport_cases = 0
    for A, B in ((1, 1), (1, 2), (4, 1), (1, 8)):
        require(gcd(A, B) == 1, "Bad fixed coefficient pair.")
        modulus = B * L
        lift1 = [x for x in range(modulus) if x % L in H1]
        lift2 = [y for y in range(modulus) if y % L in H2]
        require(len(lift1) == B * len(H1) and len(lift2) == B * len(H2),
                "Lifted kernel sizes fail.")
        for U1, U2 in product(subsets(H1), subsets(H2)):
            for target in range(modulus):
                roots = {(x, y) for x in U1 for y in U2 if (A * x + B * y - target) % L == 0}
                lifted = {(x, y) for x in lift1 for y in lift2
                          if x % L in U1 and y % L in U2
                          and (A * x + B * y - target) % modulus == 0}
                require(len(lifted) == B * len(roots), "Finite root lift count fails.")
                gamma = Q(L * len(roots), len(H1) * len(H2))
                require(gamma == Q(modulus * len(lifted), len(lift1) * len(lift2)),
                        "Period normalization creates a spurious B factor.")
                first_filter = [x for x in lift1 if x % L in U1
                                and (target - A * x) % B == 0
                                and ((target - A * x) // B) % L in U2]
                require(len(first_filter) == len(roots), "First-variable filter loses roots.")
                normalization_cases += 1
                for u, v in product(H1, H2):
                    residual1 = {pow(u, -1, L) * x % L for x in U1}
                    residual2 = {pow(v, -1, L) * y % L for y in U2}
                    residual = {(x, y) for x in residual1 for y in residual2
                                if (A * u * x + B * v * y - target) % L == 0}
                    restored = {(u * x % L, v * y % L) for x, y in residual}
                    require(restored == roots, "Transported masks do not give a bijection.")
                    transport_cases += 1
    # A nontrivial common residue tests simultaneous target division by d.
    larger_modulus, d, u, v, A, B = 40, 9, 13, 3, 1, 2
    U1, U2 = {1, 9}, {3, 11}
    def in_first(z):
        return gcd(z, larger_modulus) == 1 and z % 4 == 1
    def in_second(z):
        return gcd(z, larger_modulus) == 1 and z % 8 in (1, 3)
    require(in_first(d) and in_second(d) and in_first(u) and in_second(v)
            and all(in_first(a) for a in U1) and all(in_second(b) for b in U2),
            "The common-factor example violates its core/mask hypotheses.")
    common_cases = 0
    for m in range(larger_modulus):
        target = d * m
        roots = {(a, b) for a in U1 for b in U2
                 if (A * a + B * b - target) % larger_modulus == 0}
        V1 = {pow(d * u, -1, larger_modulus) * a % larger_modulus for a in U1}
        V2 = {pow(d * v, -1, larger_modulus) * b % larger_modulus for b in U2}
        residual = {(x, y) for x in V1 for y in V2
                    if (A * u * x + B * v * y - m) % larger_modulus == 0}
        require({(d * u * x % larger_modulus, d * v * y % larger_modulus)
                 for x, y in residual} == roots,
                "Common-factor target division and mask transport fail.")
        common_cases += 1
    return dict(modulus=L, first_kernel=sorted(H1), second_kernel=sorted(H2),
                coefficient_pairs=[[1, 1], [1, 2], [4, 1], [1, 8]],
                normalization_cases=normalization_cases, mask_transport_cases=transport_cases,
                common_factor_transport_cases=common_cases,
                common_factor_example=dict(modulus=larger_modulus, d=d, u=u, v=v))


def concrete_pair_controls():
    modulus = 92
    joint = {(a, b): 0 for a, b in product((-1, 1), repeat=2)}
    for r in range(modulus):
        if gcd(r, modulus) != 1:
            continue
        chi23 = 1 if pow(r % 23, 11, 23) == 1 else -1
        chi4 = 1 if r % 4 == 1 else -1
        joint[(chi23, chi4)] += 1
    require(set(joint.values()) == {11}, "Distinct quadratic masks not equinumerous.")
    first_mod8 = {(x * x + x * y + 6 * y * y) % 8
                  for x, y in product(range(8), repeat=2) if x % 2 or y % 2}
    require(first_mod8 == set(range(8)), "Minus-23 primitive mod8 image incomplete.")
    gaussian23 = {(x * x + y * y) % 23 for x, y in product(range(23), repeat=2)
                  if x or y}
    require(gaussian23 == set(range(1, 23)), "Gaussian unit norm image at23 fails.")
    require(all(any((n - a) % 23 in gaussian23 for a in (1, 4)) for n in range(23)),
            "Concrete local target cover at23 fails.")
    return dict(joint_quadratic_modulus=modulus,
                joint_sign_counts={f"{a},{b}": n for (a, b), n in joint.items()},
                minus23_primitive_mod8_image=sorted(first_mod8),
                gaussian_mod23_nonzero_image=sorted(gaussian23))


def negative_controls():
    rejected = []
    def reject(name, proposed):
        require(not proposed, "Corrupted control accepted: " + name)
        rejected.append(name)
    reject("reuse_old_nu_d_at_most_d", 45 <= 15)
    base = local_survivors(3, 3, False, True)
    marked = local_survivors(3, 3, False, True, 0)
    reject("claim_marker_gain_at_opposite_inert_target_prime",
           Q(len(base) - len(marked), len(base)) == Q(1, 3))
    reject("put_exclusive_target_prime_in_common_inert_J",
           Q(len(base), 9) == (1 - Q(1, 3)) * Q(2, 1))
    # U1=U2={1}, original target2. After multiplying side1 by5,
    # the residual first mask must be {5}; keeping it {1} loses the root.
    reject("hold_masks_fixed_during_cofactor_restoration", (5 + 1 - 2) % 8 == 0)
    base_roots = int((1 + 2 * 3 - 3) % 8 == 0)
    shifted_roots = int((5 * 1 + 2 * 3 - 3) % 8 == 0)
    require(base_roots == 0 and shifted_roots == 1, "Shifted-Gamma countermodel fails.")
    reject("freeze_base_Gamma_for_an_untransported_shift", base_roots == shifted_roots)
    return rejected


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "primitive_norm_different_fields_checks.json")
    args = parser.parse_args()
    paths = [Path(__file__).resolve(), HERE / "PRIMITIVE_NORM_DIFFERENT_FIELDS_CLASS_ERROR.md",
             HERE / "PRIMITIVE_NORM_ENDPOINT_CLASS_ERROR.md"]
    packet = dict(schema_version=1, passed=True,
                  arithmetic="Python standard library; exact integers and fractions",
                  local_controls=local_controls(), mark_controls=mark_product_controls(),
                  crt_control=crt_control(), finite_filter_controls=filter_controls(),
                  concrete_pair_controls=concrete_pair_controls(),
                  negative_controls_rejected=negative_controls(),
                  source_sha256={str(p.relative_to(REPO)): sha256(p.read_bytes()).hexdigest() for p in paths},
                  scope="Complete displayed finite spaces only. No numerical sampling proves the sieve, Mertens estimate, lower-count bridge or asymptotic theorem. No external review or formal verification.")
    args.output.write_text(json.dumps(packet, indent=2, sort_keys=True) + "\n")
    print(json.dumps(dict(passed=True, local=packet['local_controls'],
        mark_cases=packet['mark_controls']['cases'],
        filters=packet['finite_filter_controls'], negative_controls=len(packet['negative_controls_rejected']))))


if __name__ == "__main__":
    main()
