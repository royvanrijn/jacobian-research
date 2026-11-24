#!/usr/bin/env python3
"""Exact finite controls for class-specific mixed-field target factors.

This checks finite groups, complete stated residue spaces, and rational local
probabilities. Prime distribution, maximal orders, and all-target estimates
remain written theorems with the primary premises in the source ledger.
"""
from argparse import ArgumentParser
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
from math import gcd
from pathlib import Path
import json
import platform

HERE = Path(__file__).resolve().parent
PREFIX = Path("research/arithmetic")
NOTE = "PRIMITIVE_NORM_TARGET_FACTOR_REFINEMENT.md"
CHECKER = "verify_primitive_norm_target_factor_refinement.py"
LEDGER = "primitive_norm_target_factor_refinement_source_ledger.json"
RECEIPT = "primitive_norm_target_factor_refinement_checks.json"
SOURCES = {
    "proof": (NOTE, "T.1--T.12: subgroup target products, enlarged-order overlap densities, compressed endpoints, exact maximal order, and normalization limits."),
    "checker": (CHECKER, "Deterministic finite input domains are contained in this exact standard-library checker."),
    "mixed_field_counting": ("PRIMITIVE_NORM_DIFFERENT_FIELDS_CLASS_ERROR.md", "D.7 subgroup cover, D.10 local probabilities, D.11 uniform marked sieve, D.15--D.16 complete remainder, D.19--D.21 classical counts and concrete admissibility."),
    "endpoint_cover": ("PRIMITIVE_NORM_ENDPOINT_CLASS_ERROR.md", "All-proper-subgroup stabilizer cap and fixed-order prime-class Mertens; endpoint marking and subgroup compression."),
    "class_number_three": ("PRIMITIVE_NORM_CLASS_NUMBER_THREE.md", "Discriminant minus23 has class group C3; concrete example only."),
    "class_number_four": ("PRIMITIVE_NORM_CLASS_NUMBER_FOUR.md", "Discriminant minus56 has class group C4 and chi8 detects class parity; concrete examples only."),
    "class_prime_primary": ("sources/bourgain-fuchs-1105.3698v2.pdf", "Section2 part(5), Lemma2, equation(5.2), page15: fixed form-class prime count, inverse-class normalization and exponential error; actual retained statement reread."),
    "genus_primary": ("sources/ibukiyama-genus-characters-2303.14983v1.pdf", "Proposition3.4 and Theorem3.5, printed/PDF pages12--14; prime-to-conductor ideal correspondence on page13. Actual retained statements and proof inspected; norm formula used only at good primes."),
}
ASSURANCE = {
    "formal_verification": False,
    "external_review": False,
    "independent_replay": False,
    "historical_priority": False,
    "finite_controls_prove_prime_distribution": False,
    "optimality_of_actual_class_error": False,
}
PRIMES = (3, 5, 7)
PROBABILITIES = (Q(0), Q(1, 3), Q(2, 3), Q(9, 10), Q(1))
GROUPS = ((2,), (3,), (4,), (5,), (6,), (8,), (9,), (10,), (12,), (16,),
          (2, 2), (2, 4), (2, 6), (2, 8), (3, 3), (4, 4),
          (2, 2, 2), (2, 2, 4))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def check_sources(repo):
    ledger = json.loads((HERE / LEDGER).read_text())
    require(ledger["schema_version"] == 1 and ledger["assurance"] == ASSURANCE,
            "Changed ledger schema or assurance scope")
    rows = ledger["sources"]
    require(len(rows) == len(SOURCES) and {r["role"] for r in rows} == set(SOURCES),
            "Incorrect source roles")
    for row in rows:
        name, scope = SOURCES[row["role"]]
        require(row["path"] == str(PREFIX / name) and row["scope"] == scope,
                "Changed source path or scope")
        require(digest(repo / row["path"]) == row["sha256"],
                "Stale source: " + row["role"])
    return rows


def survivors(p, target, inert1, inert2, marked_slot=None):
    kept = []
    for a in range(p * p):
        b = (target - a) % (p * p)
        if (inert1 and a % p == 0) or (inert2 and b % p == 0):
            continue
        if marked_slot is not None:
            value = (a, b)[marked_slot]
            if value % p == 0 and value != 0:
                continue
        kept.append(a)
    return kept


def local_controls():
    normalized_cases = mark_cases = 0
    for p in PRIMES:
        for target in range(p * p):
            for inert1, inert2 in product((False, True), repeat=2):
                base = survivors(p, target, inert1, inert2)
                density = Q(len(base), p * p)
                sigma = density / ((1 - Q(int(inert1), p)) *
                                   (1 - Q(int(inert2), p)))
                if inert1 and inert2:
                    expected = (Q(p, p - 1) if target % p == 0 else
                                Q(p * (p - 2), (p - 1) ** 2))
                else:
                    expected = Q(1)
                require(sigma == expected, "Normalized shadow factor mismatch")
                normalized_cases += 1
                for slot, own, other in ((0, inert1, inert2), (1, inert2, inert1)):
                    if own:
                        continue
                    marked = survivors(p, target, inert1, inert2, slot)
                    b = Q(len(base) - len(marked), len(base))
                    reference = Q(1, p) if other else Q(p - 1, p * p)
                    blocked = other and target % p == 0
                    require(b == (0 if blocked else reference), "Marker overlap mismatch")
                    for t in PROBABILITIES:
                        ratio = (1 - t * b) / (1 - t * reference)
                        expected_ratio = 1 / (1 - t / p) if blocked else Q(1)
                        require(ratio == expected_ratio, "Exact T.4 factor mismatch")
                        require(ratio <= (Q(p, p - 1) if blocked else 1),
                                "Uniform full-marking upper bound mismatch")
                        mark_cases += 1
    return {"primes": list(PRIMES), "targets": "Every residue modulo p^2",
            "probabilities": [str(t) for t in PROBABILITIES],
            "normalized_shadow_cases": normalized_cases,
            "rational_mark_ratio_cases": mark_cases}


class FiniteGroup:
    def __init__(self, moduli):
        self.moduli = tuple(moduli)
        self.elements = tuple(product(*(range(n) for n in self.moduli)))
        self.zero = (0,) * len(self.moduli)

    def add(self, a, b):
        return tuple((x + y) % n for x, y, n in zip(a, b, self.moduli))

    def double(self, a):
        return self.add(a, a)

    def neg(self, a):
        return tuple(-x % n for x, n in zip(a, self.moduli))

    def characters(self):
        domains = [(0, 1) if n % 2 == 0 else (0,) for n in self.moduli]
        return [flags for flags in product(*domains) if any(flags)]

    @staticmethod
    def character(flags, element):
        return -1 if sum(a * b for a, b in zip(flags, element)) % 2 else 1

    def subgroups(self, universe=None):
        universe = frozenset(self.elements if universe is None else universe)
        require(self.zero in universe, "Subgroup universe lacks zero")
        seen = {frozenset((self.zero,))}
        pending = list(seen)
        while pending:
            old = pending.pop()
            for generator in sorted(universe - old):
                cyclic, element = [], self.zero
                while not cyclic or element != self.zero:
                    cyclic.append(element)
                    element = self.add(element, generator)
                new = frozenset(self.add(a, b) for a in old for b in cyclic)
                require(new <= universe, "Specified universe is not a subgroup")
                if new not in seen:
                    seen.add(new)
                    pending.append(new)
        return sorted(seen, key=lambda h: (len(h), tuple(sorted(h))))


def overlap_density(group, subgroup, flags):
    selected = [g for g in group.elements if g not in subgroup and
                group.character(flags, g) == -1]
    return Q(len(selected), 2 * len(group.elements)), selected


def fibre_controls():
    cases = aligned = odd_index = inverse_pairs = 0
    group_records = []
    for moduli in GROUPS:
        group = FiniteGroup(moduli)
        torsion = {g for g in group.elements if group.double(g) == group.zero}
        subgroups = [h for h in group.subgroups() if len(h) < len(group.elements)
                     and torsion <= h]
        local_cases = 0
        for subgroup in subgroups:
            r = len(group.elements) // len(subgroup)
            for flags in group.characters():
                density, selected = overlap_density(group, subgroup, flags)
                alignment = all(group.character(flags, a) == 1 for a in subgroup)
                expected = Q(1, 4) if alignment else Q(r - 1, 4 * r)
                require(density == expected, "Character-fibre density mismatch")
                require(all(group.neg(g) in selected and group.neg(g) != g
                            for g in selected), "Inverse-pair class count failed")
                require(not alignment or r % 2 == 0, "An odd quotient has a sign character")
                cases += 1
                local_cases += 1
                aligned += int(alignment)
                odd_index += int(r % 2 == 1)
                inverse_pairs += len(selected) // 2
        group_records.append({"cyclic_factors": list(moduli),
                              "relevant_subgroups": len(subgroups),
                              "nontrivial_sign_characters": len(group.characters()),
                              "cases": local_cases})
    require(aligned and odd_index, "Missing one of the two density mechanisms")
    return {"groups": group_records, "character_subgroup_cases": cases,
            "aligned_cases": aligned, "odd_index_cases": odd_index,
            "inverse_pairs_counted": inverse_pairs}


def least_prime_divisor(n):
    require(n > 1, "Trivial doubled group has no least prime")
    p = 2
    while n % p:
        p += 1
    return p


def compression_controls():
    cases = subgroup_cases = 0
    profiles = set()
    for moduli in GROUPS:
        target = FiniteGroup(moduli)
        doubled = frozenset(target.double(g) for g in target.elements)
        if len(doubled) == 1:
            continue
        hs = [h for h in target.subgroups(doubled) if h != doubled]
        q = least_prime_divisor(len(doubled))
        delta = Q(q - 1, 2 * q)
        models = [(target, flags, False) for flags in target.characters()]
        enlarged = FiniteGroup((*moduli, 2))
        models.append((enlarged, (*((0,) * len(moduli)), 1), True))
        for source, flags, projection in models:
            def pi(g):
                return g[:-1] if projection else g
            kernel = {g for g in source.elements if target.double(pi(g)) == target.zero}
            factors = all(source.character(flags, g) == 1 for g in kernel)
            minimal_densities = []
            for h in hs:
                preimage = frozenset(g for g in source.elements if target.double(pi(g)) in h)
                r = len(doubled) // len(h)
                density, _ = overlap_density(source, preimage, flags)
                aligned = all(source.character(flags, g) == 1 for g in preimage)
                require(density == (Q(1, 4) if aligned else Q(r - 1, 4 * r)),
                        "Actual doubling-preimage density failed")
                require(not aligned or factors, "Alignment did not kill the doubling kernel")
                if r == q:
                    minimal_densities.append(density)
                else:
                    require(Q(r - 1, 2 * r) > delta, "Larger index failed strict log saving")
                subgroup_cases += 1
            expected = delta / 2 if q % 2 else (Q(1, 4) if factors else Q(1, 8))
            require(max(minimal_densities) == expected, "Compressed beta classification failed")
            profiles.add((q, str(delta), str(expected)))
            cases += 1
    require((2, "1/4", "1/4") in profiles and (2, "1/4", "1/8") in profiles
            and (3, "1/3", "1/6") in profiles, "Missing compressed profile type")
    return {"surjection_character_models": cases, "proper_subgroup_cases": subgroup_cases,
            "profiles": [list(p) for p in sorted(profiles)]}


def chi8(n):
    require(n % 2 == 1, "chi8 input not a unit")
    return 1 if n % 8 in (1, 7) else -1


def chi_minus7(n):
    require(n % 7 != 0, "chi-7 input not a unit")
    return 1 if n % 7 in (1, 2, 4) else -1


def chi_minus4(n):
    require(n % 2 == 1, "chi-4 input not a unit")
    return 1 if n % 4 == 1 else -1


def concrete_controls():
    units = [a for a in range(56) if gcd(a, 56) == 1]
    markers = [a for a in units if chi8(a) == -1 and chi8(a) * chi_minus7(a) == 1]
    aligned = [a for a in markers if chi_minus7(a) == -1]
    independent = [a for a in markers if chi_minus4(a) == -1]
    require(len(units) == 24 and markers == aligned and len(aligned) == 6
            and len(independent) == 3, "Minus56 overlap residue densities failed")
    form_cases = 0
    for A, B, C, parity in ((1, 0, 14, 1), (2, 0, 7, 1),
                            (3, 2, 5, -1), (3, -2, 5, -1)):
        for x, y in product(range(8), repeat=2):
            n = A * x * x + B * x * y + C * y * y
            if n % 2:
                require(chi8(n) == parity, "Class-parity residue formula failed")
                form_cases += 1
    require(3 * 3 + 14 == 23 and chi_minus4(23) == -1 and
            chi8(23) * chi_minus7(23) == 1, "Principal-prime nonalignment witness failed")
    require(4 ** 2 * (-23) == -368 and (-368) // (-4) == 92 and 92 % 4 == 0,
            "Enlarged-order discriminant interface failed")
    table = [("minus23_Gaussian", Q(1, 3), Q(1, 6), 1),
             ("minus56_minus7", Q(1, 4), Q(1, 4), 0),
             ("minus56_Gaussian", Q(1, 4), Q(1, 8), 0)]
    require([cap + beta for _, _, beta, cap in table] == [Q(7, 6), Q(1, 4), Q(1, 8)],
            "Concrete endpoint exponents failed")
    return {"unit_classes_mod56": units, "minus56_markers": markers,
            "minus7_overlap": aligned, "Gaussian_overlap": independent,
            "form_parity_cases": form_cases,
            "endpoint_examples": [{"name": name, "delta": str(delta),
                                   "beta": str(beta), "cap": cap,
                                   "relative_loglog_exponent": str(cap + beta)}
                                  for name, delta, beta, cap in table]}


def negative_controls():
    rejected = []
    def reject(wrong_assertion, label):
        require(not wrong_assertion, "Corrupted control accepted: " + label)
        rejected.append(label)
    cyclic4 = FiniteGroup((4,))
    density4, _ = overlap_density(cyclic4, frozenset(((0,), (2,))), (1,))
    reject(density4 == Q(1, 8), "Halve an aligned index-two obstruction density")
    cyclic6 = FiniteGroup((6,))
    density6, _ = overlap_density(cyclic6, frozenset(((0,), (3,))), (1,))
    reject(density6 == Q(1, 4), "Assign exact quarter density to an odd-index obstruction")
    p = 3
    base = survivors(p, 0, False, True)
    marked = survivors(p, 0, False, True, 0)
    sigma = Q(len(base), p * p) / (1 - Q(1, p))
    reject(sigma == Q(p, p - 1), "Insert an exclusive-prime target factor in the shadow density")
    reject(Q(len(base) - len(marked), len(base)) == Q(1, p),
           "Retain a marker gain at a target-dividing opposite-inert prime")
    product_group = FiniteGroup((2, 4))
    doubled_kernel = frozenset(g for g in product_group.elements
                              if product_group.double(g) == product_group.zero)
    density, _ = overlap_density(product_group, doubled_kernel, (1, 0))
    reject(density == Q(1, 4), "Replace descent to the doubled group by descent to the original group")
    reject(density6 == Q(1, 3), "Count oriented classes without dividing into inverse pairs")
    return rejected


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=HERE.parent.parent)
    parser.add_argument("--output", type=Path, default=HERE / RECEIPT)
    args = parser.parse_args()
    sources = check_sources(args.repo_root.resolve())
    checks = {"local": local_controls(), "finite_group_fibres": fibre_controls(),
              "compression": compression_controls(), "concrete_examples": concrete_controls(),
              "negative_controls": negative_controls()}
    receipt = {"schema_version": 1, "status": "PASS", "python": platform.python_version(),
               "arithmetic": "Exact integers and rational probabilities; standard library only.",
               "source_ledger_sha256": digest(HERE / LEDGER), "sources": sources,
               "checks": checks, "assurance": ASSURANCE,
               "scope": "Complete stated finite domains. Distribution, all-target bounds and maximal-order limits are the written proof with the pinned primary statements."}
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print("PASS:", checks["local"]["normalized_shadow_cases"], "local normalizations;",
          checks["local"]["rational_mark_ratio_cases"], "mark ratios;",
          checks["finite_group_fibres"]["character_subgroup_cases"], "class-fibre cases;",
          checks["compression"]["proper_subgroup_cases"], "compression cases;",
          len(checks["negative_controls"]), "negative controls.")
    print(args.output)


if __name__ == "__main__":
    main()
