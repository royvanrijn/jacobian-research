#!/usr/bin/env python3
"""Exact finite controls for the endpoint class-error proof.

Checks finite probabilities, valuation-one residues, a derivative identity,
proper-subgroup profiles and bounded subset sums. It does not sample an
infinite asymptotic or certify imported analytic theorems.
"""
from argparse import ArgumentParser
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import isqrt
from pathlib import Path
import json

from verify_class_group_error_cover import (
    FiniteAbelianGroup, VerificationFailure, doubled_subgroup,
    enumerate_subgroups, inverse_families, maximal_subgroups, require,
)

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def product_fraction(values):
    result = Fraction(1)
    for value in values:
        result *= value
    return result


def local_density(p):
    require(is_prime(p) and p >= 3, "Local domain requires odd primes.")
    residues = [a for a in range(p * p) if a % p == 0 and a != 0]
    density = Fraction(len(residues), p * p)
    require(density == Fraction(p - 1, p * p), "Wrong valuation-one density.")
    return density, residues


def exact_mark_probabilities(primes, probabilities):
    require(len(set(primes)) == len(primes), "Repeated mark prime.")
    densities = [local_density(p)[0] for p in primes]
    subsets = list(product((False, True), repeat=len(primes)))
    rows, survival_checks, cap_checks = [], 0, 0
    for text in probabilities:
        t = Fraction(text)
        require(0 <= t <= 1, "Invalid probability.")
        weights = [t ** sum(marked) * (1 - t) ** (len(primes) - sum(marked))
                   for marked in subsets]
        require(sum(weights) == 1, "Probabilities do not sum to one.")
        average = sum(weight * product_fraction(
            1 - f for f, selected in zip(densities, marked) if selected)
            for marked, weight in zip(subsets, weights))
        require(average == product_fraction(1 - t * f for f in densities),
                "Mark-product expectation fails.")
        for active in subsets:
            survival = sum(weight for marked, weight in zip(subsets, weights)
                           if not any(a and m for a, m in zip(active, marked)))
            r = sum(active)
            require(survival == (1 - t) ** r, "Active-set survival fails.")
            survival_checks += 1
            if t < 1:
                for cap in range(r, len(primes) + 1):
                    require(survival >= (1 - t) ** cap, "Deficiency inequality fails.")
                    require(survival / (1 - t) ** cap >= 1,
                            "Survival normalization insufficient.")
                    cap_checks += 1
            else:
                require(survival == int(r == 0), "Probability-one boundary fails.")
        rows.append({"t": str(t), "mark_sets": len(subsets),
                     "product_expectation": str(average)})
    return rows, survival_checks, cap_checks


def poly_add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            result[exponent] = result.get(exponent, Fraction(0)) + coefficient
    return {e: c for e, c in result.items() if c}


def poly_scale(polynomial, scalar):
    return {e: c * scalar for e, c in polynomial.items() if c * scalar}


def poly_multiply(left, right):
    result = {}
    for (t1, u1), a in left.items():
        for (t2, u2), b in right.items():
            e = (t1 + t2, u1 + u2)
            result[e] = result.get(e, Fraction(0)) + a * b
    return {e: c for e, c in result.items() if c}


def derivative_identity():
    one = {(0, 0): Fraction(1)}
    t, u = {(1, 0): Fraction(1)}, {(0, 1): Fraction(1)}
    u2 = poly_multiply(u, u)
    f = poly_add(u, poly_scale(u2, -1))
    denominator = poly_add(one, poly_scale(poly_multiply(t, f), -1))
    numerator = poly_add(poly_multiply(u, denominator), poly_scale(f, -1))
    expected = poly_multiply(u2, poly_add(one, poly_scale(t, -1),
                                         poly_multiply(t, u)))
    require(numerator == expected, "Derivative numerator identity fails.")
    return {"ring": "Q[t,u]", "cleared_numerator": "u^2*(1-t+t*u)",
            "coefficients_t_degree_u_degree":
            [[a, b, str(c)] for (a, b), c in sorted(numerator.items())],
            "sign_scope": "Written proof uses 0<=t<=1 and 0<u<=1/3, not a grid sign test."}


def group_profiles(specification):
    group = FiniteAbelianGroup(specification["moduli"])
    B, families = doubled_subgroup(group), inverse_families(group)
    maxima, subgroups = maximal_subgroups(group), enumerate_subgroups(group, B)
    direct_maxima = {H for H in subgroups if H != B
                     and not any(H < K < B for K in subgroups)}
    require({M for _, M in maxima} == direct_maxima,
            "Homomorphism kernels disagree with subgroup enumeration.")
    tally = Counter()
    for index, M in maxima:
        outside = [g for g, _ in families if group.scale(2, g) not in M]
        density = Fraction(len(outside), group.size)
        require(density == Fraction(index - 1, 2 * index), "Wrong class density.")
        cap = sum(group.order(group.scale(2, g)) - 2 for g in outside)
        require(cap >= 0, "Negative familywise cap.")
        tally[(index, str(density), cap)] += 1
    profiles = [dict(index_prime=i, density=d, cap=k, multiplicity=n)
                for (i, d, k), n in sorted(tally.items())]
    require(profiles == specification["expected_profiles"], "Familywise profile mismatch.")
    require(bool(profiles) == (len(B) > 1), "Wrong trivial doubled-group case.")
    tally = Counter()
    for H in subgroups:
        if H == B:
            continue
        index = len(B) // len(H)
        outside = [g for g, _ in families if group.scale(2, g) not in H]
        density = Fraction(len(outside), group.size)
        require(density == Fraction(index - 1, 2 * index),
                "Wrong proper-subgroup class density.")
        tally[(index, str(density), index - 2)] += 1
    stronger = [dict(index=i, density=d, cap=k, multiplicity=n)
                for (i, d, k), n in sorted(tally.items())]
    require(stronger == specification["expected_stabilizer_profiles"],
            "Stabilizer profile mismatch.")
    compressed = None
    if stronger:
        q = min(row["index"] for row in stronger)
        require(is_prime(q), "Minimum proper-subgroup index is not prime.")
        compressed = dict(smallest_prime=q, density=str(Fraction(q - 1, 2 * q)),
                          loglog_cap=q - 2)
    return dict(name=specification["name"], moduli=specification["moduli"],
                group_order=group.size, doubled_group_order=len(B),
                maximal_subgroups=len(maxima), profiles=profiles,
                proper_subgroups=len(subgroups) - 1, stabilizer_profiles=stronger,
                compressed_stabilizer_profile=compressed)


def subset_sums(group, sequence):
    result = set()
    for selected in product((False, True), repeat=len(sequence)):
        total = group.zero
        for value, active in zip(sequence, selected):
            if active:
                total = group.add(total, value)
        result.add(total)
    return frozenset(result)


def stabilizer(group, B, subset):
    return frozenset(h for h in B if {group.add(s, h) for s in subset} == subset)


def subset_sum_controls(specification, max_length):
    group = FiniteAbelianGroup(specification["moduli"])
    B = doubled_subgroup(group)
    sequences = deficient = prefixes = translations = 0
    for length in range(max_length + 1):
        for sequence in product(sorted(B), repeat=length):
            sequences += 1
            sums = subset_sums(group, sequence)
            H = stabilizer(group, B, sums)
            require(H <= B and group.zero in H, "Invalid stabilizer.")
            if sums == B:
                require(H == B, "Full sumset stabilizer is not full.")
                continue
            deficient += 1
            require(H < B, "Proper sumset has full stabilizer.")
            index = len(B) // len(H)
            outside = sum(v not in H for v in sequence)
            require(outside <= index - 2, "Deficiency cap fails.")
            require(1 + outside <= len(sums) // len(H) <= index - 1,
                    "Quotient cardinality bound fails.")
            previous_size = 1
            for i in range(length + 1):
                prefix = subset_sums(group, sequence[:i])
                saturated = {group.add(a, h) for a in prefix for h in H}
                require(stabilizer(group, B, saturated) == H,
                        "A quotient prefix has a nonzero final period.")
                size = len(saturated) // len(H)
                if i:
                    require(size >= previous_size + int(sequence[i - 1] not in H),
                            "A nonzero quotient increment did not grow its prefix.")
                previous_size = size
                prefixes += 1
            for residual in group.elements:
                target = {group.add(residual, b) for b in B}
                attained = {group.add(residual, s) for s in sums}
                require(attained < target, "Translation changed target deficiency.")
                require(stabilizer(group, B, attained) == H,
                        "Translation changed the stabilizer.")
                translations += 1
    return dict(name=specification["name"], max_sequence_length=max_length,
                sequences=sequences, proper_sumsets=deficient,
                quotient_prefix_checks=prefixes, translated_target_checks=translations)


def negative_controls(inputs):
    rejected = []
    def must_reject(name, check):
        try:
            check()
        except VerificationFailure:
            rejected.append(name)
            return
        raise VerificationFailure("Invalid control accepted: " + name)
    must_reject("replace_valuation_one_density_by_one_over_p", lambda:
                require(local_density(3)[0] == Fraction(1, 3), "Zero residue is excluded."))
    def singleton(t, active_count, normalization):
        expectation = Fraction(0)
        for marked in product((False, True), repeat=active_count):
            weight = t ** sum(marked) * (1 - t) ** (active_count - sum(marked))
            if not any(marked):
                expectation += weight
        require(normalization * expectation >= 1, "Normalization undercounts singleton.")
    singleton(Fraction(9, 10), 1, Fraction(10))
    must_reject("omit_inverse_survival_normalization", lambda:
                singleton(Fraction(9, 10), 1, Fraction(1)))
    must_reject("use_probability_one_with_positive_cap", lambda:
                singleton(Fraction(1), 1, Fraction(1)))
    wrong = deepcopy(next(g for g in inputs["groups"] if g["name"] == "C3"))
    wrong["expected_profiles"][0]["cap"] = 0
    must_reject("erase_C3_loglog_cap", lambda: group_profiles(wrong))
    group = FiniteAbelianGroup([8])
    B, sequence = doubled_subgroup(group), ((2,), (2,))
    sums = subset_sums(group, sequence)
    require(sums < B and stabilizer(group, B, sums) == frozenset(((0,),)),
            "Invalid maximal-only countermodel.")
    must_reject("restrict_index_minus_two_cover_to_maximal_subgroups", lambda:
                require(any(sum(v not in M for v in sequence) <= index - 2
                            for index, M in maximal_subgroups(group)),
                        "Composite-index stabilizers cannot be discarded."))
    return rejected


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path,
                        default=HERE / "primitive_norm_endpoint_class_error_inputs.json")
    parser.add_argument("--output", type=Path,
                        default=HERE / "primitive_norm_endpoint_class_error_checks.json")
    args = parser.parse_args()
    inputs = json.loads(args.inputs.read_text())
    require(inputs["schema_version"] == 1, "Unknown input schema.")
    local = []
    for p in inputs["local_primes"]:
        density, residues = local_density(p)
        local.append(dict(p=p, modulus=p * p, valuation_one_residues=residues,
                          density=str(density)))
    rows, survival, normalizations = exact_mark_probabilities(
        inputs["mark_primes"], inputs["mark_probabilities"])
    groups = [group_profiles(g) for g in inputs["groups"]]
    sums = [subset_sum_controls(g, inputs["max_sequence_length"]) for g in inputs["groups"]]
    paths = [Path(__file__).resolve(), args.inputs.resolve(),
             HERE / "verify_class_group_error_cover.py",
             HERE / "PRIMITIVE_NORM_ENDPOINT_CLASS_ERROR.md"]
    packet = dict(schema_version=1, passed=True,
        arithmetic="Python standard library; exact integers and fractions",
        local_controls=local, mark_product_controls=rows,
        active_set_survival_checks=survival, cap_normalization_checks=normalizations,
        uniform_derivative_identity=derivative_identity(), group_controls=groups,
        subset_sum_controls=sums, negative_controls_rejected=negative_controls(inputs),
        source_sha256={str(p.relative_to(REPO)): digest(p) for p in paths},
        scope="Finite algebra and probability only. Unbounded class error, Mertens, sieve and form transfer are the written proof. No formal verification or priority claim.")
    args.output.write_text(json.dumps(packet, indent=2, sort_keys=True) + "\n")
    print(json.dumps(dict(passed=True, groups=len(groups), local_primes=len(local),
        survival_checks=survival, cap_normalization_checks=normalizations,
        subset_sum_sequences=sum(r["sequences"] for r in sums),
        deficient_subset_sums=sum(r["proper_sumsets"] for r in sums),
        negative_controls=len(packet["negative_controls_rejected"]), output=str(args.output))))


if __name__ == "__main__":
    main()
