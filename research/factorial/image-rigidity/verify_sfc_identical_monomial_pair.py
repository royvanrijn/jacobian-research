#!/usr/bin/env python3
"""Exact controls for SFC_IDENTICAL_MONOMIAL_PAIR.md (standard library only).

The finite checks support the identities, not the uniform external interlacing
theorem. Every reported coprimality is independently checked by a full rational
polynomial Bezout identity. No floating-point root finding is used.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json
import math


NOTE_SHA = "7157a05550267289e1cbfc2c01f08659eb2d85e4b9d62661d3e4e57c2fb0e04e"
PDF_SHA = "bf4731acc6af6348c4c5b306eaffff7d0be660f67945b2234458e742ac816e13"
TEXT_SHA = "9c4f1c16648a3e193c6a0024669e72ec85a8e3c543adcad0598c7b4ded00b5b3"
PROOF_SHA = "de91eb83ce4ec7a2e5bd7ff5be59a3c7023c19f7f42da9eb2dcd78f6d405a445"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(a, b):
    return trim([(a[k] if k < len(a) else 0) +
                 (b[k] if k < len(b) else 0)
                 for k in range(max(len(a), len(b)))])


def scale(a, scalar):
    return trim([scalar * x for x in a])


def mul(a, b):
    result = [Q(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return trim(result)


def divide(a, b):
    a, b = trim(a), trim(b)
    require(b != [0], "Polynomial divisor is nonzero")
    quotient = [Q(0)] * max(1, len(a)-len(b)+1)
    while a != [0] and len(a) >= len(b):
        degree = len(a)-len(b)
        factor = a[-1]/b[-1]
        quotient[degree] += factor
        for j, value in enumerate(b):
            a[degree+j] -= factor*value
        a = trim(a)
    return trim(quotient), a


def bezout(a, b):
    r0, r1 = a, b
    s0, s1 = [Q(1)], [Q(0)]
    t0, t1 = [Q(0)], [Q(1)]
    while r1 != [0]:
        quotient, remainder = divide(r0, r1)
        r0, r1 = r1, remainder
        s0, s1 = s1, add(s0, scale(mul(quotient, s1), -1))
        t0, t1 = t1, add(t0, scale(mul(quotient, t1), -1))
    require(len(r0) == 1 and r0[0] != 0, "Consecutive moment polynomials are coprime")
    s0, t0 = scale(s0, 1/r0[0]), scale(t0, 1/r0[0])
    require(add(mul(a, s0), mul(b, t0)) == [1], "Full exact Bezout identity")
    return s0, t0


def evaluate(p, z):
    answer = Q(0)
    for value in reversed(p):
        answer = answer*z+value
    return answer


def factorial_monomial(exponents, k):
    return math.prod(math.factorial(r*k) for r in exponents)


def coefficients(exponents, maximum):
    require(bool(exponents), "A template must have variables")
    require(all(type(r) is int and r >= 0 for r in exponents), "Nonnegative integer exponents")
    require(sum(exponents) >= 1, "Nonconstant template")
    a = [Q(factorial_monomial(exponents, k), math.factorial(k))
         for k in range(maximum+1)]
    delta = [Q(math.factorial(k), a[k]) for k in range(maximum+1)]
    return a, delta


def moment_polynomial(a, n):
    return [a[k]*a[n-k] for k in range(n+1)]


def direct_moment(exponents, c, d, n):
    return sum(math.comb(n, k)*c**(n-k)*d**k *
               factorial_monomial(exponents, n-k)*factorial_monomial(exponents, k)
               for k in range(n+1))


def is_logconcave(values):
    return all(values[k]**2 >= values[k-1]*values[k+1]
               for k in range(1, len(values)-1))


def validate_reciprocal_coefficients(a, n, values):
    require(len(values) == n+1, "Correct reciprocal degree")
    require(all(values[k]*a[k]*a[n-k] == 1 for k in range(n+1)),
            "Exact coefficientwise reciprocal identity")


def rejected(operation, label):
    try:
        operation()
    except ValueError:
        return label
    raise ValueError("A required rejection did not happen: " + label)


def exact_checks():
    templates = ((1,), (2,), (3,), (5,), (1, 1), (1, 2), (2, 2),
                 (1, 1, 1), (1, 3), (2, 3), (0, 2), (2, 0, 1))
    maximum = 10
    counts = dict(ratio_identities=0, logarithmic_derivative_bounds=0,
                  normalized_coefficient_identities=0, normalized_logconcavity_checks=0,
                  moment_identities=0, bezout_identities=0, degree_one_identities=0,
                  zero_coefficient_cases=0)
    rows = []
    for exponents in templates:
        degree = sum(exponents)
        a, delta = coefficients(exponents, maximum)
        for k in range(maximum):
            ratio = Q((k+1)**2, math.prod(r*k+j for r in exponents for j in range(1, r+1)))
            require(delta[k+1]/delta[k] == ratio, "Factorial ratio identity")
            counts["ratio_identities"] += 1
        if degree >= 2:
            require(is_logconcave(delta), "Uniform delta log-concavity control")
            for x in (Q(0), Q(1, 2), Q(7, 3)):
                log_derivative = Q(2, x+1) - sum(Q(r, r*x+j)
                    for r in exponents for j in range(1, r+1))
                require(log_derivative <= Q(2-degree, x+1) <= 0,
                        "Exact logarithmic derivative inequality")
                counts["logarithmic_derivative_bounds"] += 1
        fs = [moment_polynomial(a, n) for n in range(maximum+1)]
        for n in range(1, maximum+1):
            reciprocal = [1/value for value in fs[n]]
            validate_reciprocal_coefficients(a, n, reciprocal)
            normalized = [reciprocal[k]/math.comb(n, k) for k in range(n+1)]
            for k in range(n+1):
                require(normalized[k] == delta[k]*delta[n-k]/math.factorial(n),
                        "Normalized coefficient identity")
                counts["normalized_coefficient_identities"] += 1
            require(reciprocal == reciprocal[::-1] and all(v > 0 for v in reciprocal),
                    "Positive palindromic reciprocal coefficients")
            if degree >= 2:
                require(is_logconcave(normalized), "Imported theorem hypotheses")
                counts["normalized_logconcavity_checks"] += 1
            for c, d in ((Q(1), Q(-1)), (Q(2), Q(3)), (Q(-3, 2), Q(2, 5))):
                direct = direct_moment(exponents, c, d, n)
                reconstructed = math.factorial(n)*c**n*evaluate(fs[n], d/c)
                require(direct == reconstructed, "Moment polynomial identity")
                counts["moment_identities"] += 1
            for c, d in ((Q(0), Q(3)), (Q(-2), Q(0))):
                require(direct_moment(exponents, c, d, n) ==
                        (c+d)**n*factorial_monomial(exponents, n) != 0,
                        "One active monomial has nonzero moments")
                counts["zero_coefficient_cases"] += 1
        for n in range(1, 9):
            bezout(fs[n], fs[n+1])
            counts["bezout_identities"] += 1
            if degree == 1:
                require(add(fs[n+1], scale([Q(0)]+fs[n], -1)) == [1],
                        "Degree-one geometric identity")
                counts["degree_one_identities"] += 1
        rows.append(dict(exponents=exponents, total_degree=degree,
                         checked_moment_orders=[1, maximum], checked_adjacent_pairs=[1, 8]))

    controls = []
    geometric = [Q(1)]*5
    require(is_logconcave(geometric), "Constant sequence is also log-concave")
    require(all(geometric[k]**2 <= geometric[k-1]*geometric[k+1] for k in range(1, 4)),
            "The generalized-kernel control is log-convex")
    false_normalized = [Q(1), Q(1, 2), Q(1)]
    require(not is_logconcave(false_normalized), "Arbitrary log-convexity does not verify the theorem hypothesis")
    controls.append("positive log-convex a_k=1 fails the normalized-logconcavity hypothesis at n=2; not a coprimality counterexample")
    require(not is_logconcave([Q(math.factorial(k)) for k in range(5)]),
            "Degree one must not be fed into the D>=2 argument")
    controls.append("degree-one delta=k! rejects log-concavity, while its separate Bezout proof passes")
    a, _ = coefficients((2,), 4)
    good = [1/value for value in moment_polynomial(a, 4)]
    bad = list(good)
    bad[2] += 1
    controls.append(rejected(lambda: validate_reciprocal_coefficients(a, 4, bad), "altered middle coefficient"))
    controls.append(rejected(lambda: validate_reciprocal_coefficients(a, 4, good[:-1]), "malformed coefficient array"))
    controls.append(rejected(lambda: coefficients((-1, 2), 4), "negative template exponent"))
    controls.append(rejected(lambda: coefficients((0, 0), 4), "constant template"))
    return dict(templates=rows, counts=counts, rejecting_controls=controls)


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=here/"sfc_identical_monomial_pair_verification.json")
    parser.add_argument("--source-pdf", type=Path)
    args = parser.parse_args()
    note = here/"SFC_IDENTICAL_MONOMIAL_PAIR.md"
    require(hashlib.sha256(note.read_bytes()).hexdigest() == NOTE_SHA, "Reviewed note SHA-256")
    source = {"url": "https://arxiv.org/abs/2609.34551v1", "statement": "Theorem 1.1 with the explicit weaker hypothesis of Remark 3.4",
              "pdf_sha256": PDF_SHA, "extracted_text_sha256": TEXT_SHA,
              "reviewed_sections_2_3_sha256": PROOF_SHA, "pdf_checked_this_run": False}
    source_pdf = args.source_pdf
    if source_pdf is None:
        candidate = here.parents[2]/"research/artifacts/local/stepanov_closure/observable_mao_wang_zheng_2609.34551.pdf"
        if candidate.exists():
            source_pdf = candidate
    if source_pdf is not None:
        require(hashlib.sha256(source_pdf.read_bytes()).hexdigest() == PDF_SHA, "Inspected primary-source PDF SHA-256")
        source["pdf_checked_this_run"] = True
    result = {"status": "PASS", "scope": "finite exact controls; all-order proof is the written factorial-ratio argument and attributed interlacing input",
              "note": note.name, "note_sha256": NOTE_SHA,
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "primary_source": source, "verification": exact_checks()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"status": result["status"], "counts": result["verification"]["counts"],
                      "rejecting_controls": len(result["verification"]["rejecting_controls"]),
                      "source_pdf_checked": source["pdf_checked_this_run"], "receipt": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
