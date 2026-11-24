#!/usr/bin/env python3
"""Separate stdlib replay of every terminal-11/12 Bezout certificate.

This program does not import verify_stepanov_terminal_extension.py, FLINT, or SymPy. Moment
coefficients are reconstructed by the binomial recurrence and consecutive
rising-factorial blocks, without factorial tables or factorial quotients.
The analytic/complement theorems remain written mathematical inputs; their
current bytes and the exact numerical hypotheses used here are checked.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys


CERTIFICATE_SHA256 = "401fe2b233dd45dff604aa955b0b6134662ccfd894675b2155375f082d699ae1"
PROOF_DIRECTORY = Path("research/factorial/image-rigidity")
PROOF_SHA256 = {
    "lower_exponent_initial_and_boundary.tex": "e66e02d0354de5cd92d33eef54efea15c387b9f6a0c1b62ed379295e13f07d45",
    "GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md": "22700eab6cdb655dc676d75158f17ca5350fe0b7684fc6fe933c5a1195f2db9c",
    "GAMMA_PRODUCT_ROOT_MODULUS.tex": "b6cb212dd2b461a12c5f687375089da7839495031737ef04d0de2cc8151950b9",
    "GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md": "5ddbe130e436416808dda833a06b6512a0d2b62eda682ccd5029e8414f32f473",
    "ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex": "604e1687ea759785fe5d97a2cf214c84732be10ac0df0cae208d7bd01435dcd7",
}
# terminal: (previous radius, upper radius, first automatic positive-a gap,
#            first automatic a=0 gap, complete finite case count)
REGIONS = {11: (14, 16, 57, 9, 14316), 12: (16, 17, 65, 10, 20745)}
ALLOWED_PRIMES = (65537, 1000003, 1000033)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def without_duplicate_keys(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "Duplicate JSON key")
        result[key] = value
    return result


def reject_json_constant(value):
    raise ValueError("Non-finite JSON constant: " + value)


def parse_record(line):
    return json.loads(line, object_pairs_hook=without_duplicate_keys,
                      parse_constant=reject_json_constant)


def is_prime(value):
    return (type(value) is int and value >= 2
            and all(value % q for q in range(2, math.isqrt(value) + 1)))


def proof_sources(repo):
    actual = {}
    for name, expected in PROOF_SHA256.items():
        path = repo / PROOF_DIRECTORY / name
        require(path.is_file(), "Missing proof input: " + name)
        actual[name] = sha256(path)
        require(actual[name] == expected, "Changed proof input: " + name)
    return actual


def independent_domain():
    # Enumerate a first, then d, unlike the certificate generator.  The
    # inequality is tested directly, rather than sharing its enumerator.
    expected = set()
    counts = {}
    for n, (_, _, gap, constant_gap, wanted) in REGIONS.items():
        positive = {
            (n, a, d)
            for a in range(1, (gap - 1) * (n - 2))
            for d in range(1, gap)
            if a < d * (n - 2)
        }
        boundary = {(n, 0, d) for d in range(1, constant_gap)}
        require(not positive & boundary, "Disjoint support boundary")
        current = positive | boundary
        direct_sum = sum(d * (n - 2) - 1 for d in range(1, gap))
        require(len(positive) == direct_sum, "Positive support count")
        require(len(current) == wanted, "Complete terminal count")
        expected.update(current)
        counts[str(n)] = {"positive_lower_cases": len(positive),
                          "constant_boundary_cases": len(boundary),
                          "total": len(current)}
    require(len(expected) == 35061, "Complete two-terminal count")
    maximum = max(n * (a + d) for n, a, d in expected)
    require(maximum == 8436, "Maximum factorial argument")
    return expected, counts, maximum


def numerical_complements():
    result = []
    for n, (c0, c1, gap, constant_gap, _) in REGIONS.items():
        radius_margins = [2 * c0 ** (n - 1) - (c0 + 1) ** (n - 1),
                          2 * c1 ** n - (c1 + 1) ** n]
        require(all(value > 0 for value in radius_margins),
                "Strict rational root-annulus radius")
        # B_j=((gap+1)j)!/j!: reconstruct as a direct integer product.
        b_previous = math.prod(range(n, (gap + 1) * (n - 1) + 1))
        b_upper = math.prod(range(n + 1, (gap + 1) * n + 1))
        difference = (b_upper ** (n - 1)
                      - (c0 * c1) ** (n * (n - 1)) * b_previous ** n)
        require(difference > 0, "Strict endpoint-center separation")
        require(2 ** constant_gap >= 4 * n * n,
                "Constant-plus-power large-gap criterion")
        # Upper-order inner bound is strict. At a=d(N-2), its Gamma
        # start equals the previous order's outer-bound Gamma start.
        for d in range(1, gap):
            a = d * (n - 2)
            require(a * n + 1 == a * (n - 1) + 1 + d * (n - 2),
                    "Sharp-soft boundary equality uses actual shapes")
        encoded = difference.to_bytes((difference.bit_length() + 7) // 8, "big")
        result.append({"terminal": n, "radii": [c0, c1],
                       "strict_radius_integer_margins": radius_margins,
                       "positive_lower_gap_cutoff": gap,
                       "constant_gap_cutoff": constant_gap,
                       "annulus_difference_bits": difference.bit_length(),
                       "annulus_difference_sha256": hashlib.sha256(encoded).hexdigest()})
    return result


def coefficients_from_shape(m, gamma, d, p):
    """Ascending coefficients of sum binom(m,k)(gamma)_(dk) z^(m-k)."""
    require(type(m) is int and m >= 1 and gamma >= 1 and d >= 1,
            "Coefficient parameter domain")
    require(gamma + d * m - 1 < p, "Every rising-factorial factor is a unit")
    coefficient = [0] * (m + 1)
    binomial = 1
    rising = 1
    for k in range(m + 1):
        coefficient[m - k] = binomial * rising % p
        if k != m:
            numerator = binomial * (m - k)
            require(numerator % (k + 1) == 0, "Exact binomial recurrence")
            binomial = numerator // (k + 1)
            for factor in range(gamma + d * k, gamma + d * (k + 1)):
                rising = rising * factor % p
    require(coefficient[-1] == 1 and coefficient[0] != 0,
            "Reversed moment polynomial has unit endpoints")
    return coefficient


def coefficients(m, a, d, p):
    return coefficients_from_shape(m, a * m + 1, d, p)


def check_record(row, expected, seen):
    require(type(row) is dict and set(row) == {"N", "a", "d", "p", "u", "v"},
            "Exact record schema")
    require(all(type(row[key]) is int for key in ("N", "a", "d", "p")),
            "Integer parameters; booleans are not integers")
    case = row["N"], row["a"], row["d"]
    require(case in expected, "Record outside the complete residual domain")
    require(case not in seen, "Repeated certificate case")
    n, a, d = case
    p = row["p"]
    require(p in ALLOWED_PRIMES, "Allowed verified prime")
    require(p > n * (a + d), "Full good-reduction domain")
    for key, length_limit in (("u", n - 1), ("v", n)):
        values = row[key]
        require(type(values) is list and 1 <= len(values) <= length_limit,
                "Witness vector and degree bound")
        require(all(type(value) is int and 0 <= value < p for value in values),
                "Canonical integer witness coefficients")
        require(values[-1] != 0, "Canonical nonzero leading witness coefficient")
    return case


def check_identity(row, upper, lower):
    n, p, u, v = row["N"], row["p"], row["u"], row["v"]
    require(len(upper) == n + 1 and len(lower) == n,
            "Actual upper and lower polynomial degrees")
    require(upper[-1] == lower[-1] == 1 and upper[0] and lower[0],
            "Full-degree monic moment polynomials")
    # Check by coefficientwise dot products, not the generator's
    # polynomial multiplication, gcd, or resultant implementation.
    for degree in range(2 * n - 1):
        value = sum(u[i] * upper[degree - i]
                    for i in range(max(0, degree - n), min(len(u), degree + 1)))
        value += sum(v[i] * lower[degree - i]
                     for i in range(max(0, degree - n + 1), min(len(v), degree + 1)))
        require(value % p == (1 if degree == 0 else 0),
                "Full Bezout coefficient identity")


def check_coverage(expected, seen):
    require(seen == expected, "Complete residual domain, with no omitted cases")


def negative_controls(sample, expected):
    rejected = []

    def reject(name, action):
        try:
            action()
        except (ValueError, TypeError, KeyError):
            rejected.append(name)
        else:
            raise ValueError("Negative control was accepted: " + name)

    case = sample["N"], sample["a"], sample["d"]
    n, a, d = case
    p = sample["p"]
    upper = coefficients(n, a, d, p)
    lower = coefficients(n - 1, a, d, p)
    check_record(sample, expected, set())
    check_identity(sample, upper, lower)

    altered = copy.deepcopy(sample)
    altered["N"] = True
    reject("boolean-terminal", lambda: check_record(altered, expected, set()))
    altered = copy.deepcopy(sample)
    del altered["u"]
    reject("missing-witness", lambda: check_record(altered, expected, set()))
    altered = copy.deepcopy(sample)
    altered["u"][0] = float(altered["u"][0])
    reject("noninteger-witness-coefficient", lambda: check_record(altered, expected, set()))
    altered = copy.deepcopy(sample)
    altered["u"][0] = p
    reject("noncanonical-witness-residue", lambda: check_record(altered, expected, set()))
    altered = copy.deepcopy(sample)
    altered["u"] = [1] * n
    reject("excess-witness-degree", lambda: check_record(altered, expected, set()))
    reject("duplicate-case", lambda: check_record(sample, expected, {case}))
    altered = copy.deepcopy(sample)
    altered["d"] = REGIONS[n][2]
    reject("out-of-domain-case", lambda: check_record(altered, expected, set()))
    reject("missing-case", lambda: check_coverage(expected, expected - {case}))
    reject("duplicate-json-key", lambda: parse_record('{"N":11,"N":12}'))
    reject("invalid-json-number", lambda: parse_record('{"N":NaN}'))

    # The following controls reach the polynomial-identity layer; they
    # do not merely fail the fixed compressed-file SHA256 check.
    altered = copy.deepcopy(sample)
    altered["u"][0] = (altered["u"][0] + 1) % p
    reject("altered-bezout-witness", lambda: check_identity(altered, upper, lower))
    wrong_upper = upper.copy()
    wrong_upper[0] = (wrong_upper[0] + 1) % p
    reject("altered-upper-moment-coefficient",
           lambda: check_identity(sample, wrong_upper, lower))
    wrong_lower = coefficients_from_shape(n - 1, a * n + 1, d, p)
    require(wrong_lower != lower and a > 0, "Real two-shape negative-control difference")
    reject("frozen-upper-gamma-shape-at-lower-order",
           lambda: check_identity(sample, upper, wrong_lower))
    require(not is_prime(65535), "Composite good-reduction negative")
    rejected.append("composite-prime")
    require(2 ** 9 < 4 * 12 ** 2, "N12 constant boundary d9 cannot be omitted")
    rejected.append("N12-constant-d9-is-not-in-large-gap-complement")
    return rejected


def replay(path, expected):
    seen = set()
    per_terminal = {str(n): 0 for n in REGIONS}
    per_prime = {str(p): 0 for p in ALLOWED_PRIMES}
    case_digest = hashlib.sha256()
    coefficient_digest = hashlib.sha256()
    identity_coefficients = 0
    boundary_appell = 0
    sample = None
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            row = parse_record(line)
            case = check_record(row, expected, seen)
            n, a, d = case
            p = row["p"]
            upper = coefficients(n, a, d, p)
            lower = coefficients(n - 1, a, d, p)
            check_identity(row, upper, lower)
            if a == 0:
                require(all((j + 1) * upper[j + 1] % p == n * lower[j] % p
                            for j in range(n)), "Exact constant-family Appell identity")
                boundary_appell += 1
            if sample is None and a > 0:
                sample = row
            seen.add(case)
            per_terminal[str(n)] += 1
            per_prime[str(p)] += 1
            identity_coefficients += 2 * n - 1
            case_digest.update((f"{n},{a},{d},{p}\n").encode("ascii"))
            coefficient_digest.update(json.dumps([n, a, d, p, upper, lower],
                                                  separators=(",", ":")).encode("ascii") + b"\n")
    check_coverage(expected, seen)
    require(sample is not None, "Actual positive-lower case for negative controls")
    require(len(seen) == 35061 and boundary_appell == 17, "Complete finite and boundary counts")
    return {
        "full_bezout_identities": len(seen),
        "independently_reconstructed_moment_vectors": 2 * len(seen),
        "bezout_coefficient_equalities": identity_coefficients,
        "constant_boundary_appell_identities": boundary_appell,
        "per_terminal": per_terminal,
        "per_prime": per_prime,
        "case_stream_sha256": case_digest.hexdigest(),
        "rising_coefficient_stream_sha256": coefficient_digest.hexdigest(),
        "negative_controls": negative_controls(sample, expected),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--certificate", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    certificate = args.certificate.resolve()
    sources = proof_sources(repo)
    require(sha256(certificate) == CERTIFICATE_SHA256, "Exact audited compressed certificate")
    for p in ALLOWED_PRIMES:
        require(is_prime(p) and p > 8436, "Trial-division primality and uniform good reduction")
    expected, counts, maximum = independent_domain()
    complements = numerical_complements()
    finite = replay(certificate, expected)
    require(sources == proof_sources(repo), "Proof inputs unchanged during replay")
    require(sha256(certificate) == CERTIFICATE_SHA256, "Certificate unchanged during replay")
    result = {
        "schema": 1,
        "status": "PASS",
        "terminals": [11, 12],
        "conclusion": "Complete finite certificate for every residual support at terminals 11 and 12; written complement theorems cover all other ordinary univariate supports at those two terminals",
        "algorithmically_separate_internal_replay": True,
        "independent_external_review": False,
        "formal_verification": False,
        "standard_library_only": True,
        "generator_imported": False,
        "coefficient_method": "Exact binomial recurrence and modular rising-factorial blocks at each own-order Gamma shape; no factorial tables or quotients",
        "identity_method": "Every coefficient of U H_N + V H_(N-1) checked directly as a dot product modulo the certified prime",
        "all_order_SFC": False,
        "multivariate_extension": False,
        "uniform_eventual_theorem_used": False,
        "proof_source_sha256": sources,
        "certificate_sha256": CERTIFICATE_SHA256,
        "certificate_bytes": certificate.stat().st_size,
        "checker_sha256": sha256(Path(__file__)),
        "python_version": sys.version.split()[0],
        "domain": counts,
        "maximum_factorial_argument": maximum,
        "complement_hypotheses": complements,
        "finite_replay": finite,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "PASS", "cases": finite["full_bezout_identities"],
                      "moment_vectors": finite["independently_reconstructed_moment_vectors"],
                      "coefficient_equalities": finite["bezout_coefficient_equalities"],
                      "negative_controls": len(finite["negative_controls"]),
                      "output": str(args.output)}))


if __name__ == "__main__":
    main()
