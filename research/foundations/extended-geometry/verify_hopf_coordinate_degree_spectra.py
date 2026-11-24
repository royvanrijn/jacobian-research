#!/usr/bin/env python3
"""Exact bounded controls for the written quadric degree/spectrum theorem.

The primary generator theorem and the all-word degree induction are not
inferred from these examples. Python standard library only.
"""
import argparse
import hashlib
import importlib.util
import json
import platform
from fractions import Fraction as F
from math import comb, gcd, isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
PREFIX = Path("research/foundations/extended-geometry")
NOTE = "HOPF_COORDINATE_DEGREES_AND_FULL_HYPERBOLIC_SPECTRA.md"
CHECKER = "verify_hopf_coordinate_degree_spectra.py"
LEDGER = "hopf_coordinate_degree_spectra_source_ledger.json"
RECEIPT = "hopf_coordinate_degree_spectra_verification.json"
SOURCES = {
    "written_proof": (NOTE, "All-word degree induction, explicit realizations, full displayed-family spectra and divisor criterion; uses the separately pinned classical generator theorem."),
    "exact_checker": (CHECKER, "Deterministic bounded exact controls, with complete input-generation parameters in this file."),
    "arithmetic_helper": ("verify_hopf_lnd_orbit.py", "Sparse rational polynomial arithmetic and the independent xy-reduction normal form only."),
    "primary_generation_theorem": ("sources/blanc-dubouloz-a1-fibered-surfaces-2011.pdf", "Theorem 5.4.5(7),(8), printed pages 5915-5919, with field convention page 5889: generation for the split quadratic surface. The printed missing-v shear term is corrected directly in C.4."),
    "retained_family_context": ("HOPF_AUTOMORPHISM_ORBIT_SPECTRA.md", "Section 8 identifies the retained family and its separate classification/Image interpretation; no such result is needed for the new algebraic degree proof."),
}
PRIMARY = {
    "authors": "Jeremy Blanc and Adrien Dubouloz",
    "title": "Automorphisms of A1-fibered affine surfaces",
    "publication": "Transactions of the AMS 363 (2011), 5887-5924",
    "url": "https://algebra.dmi.unibas.ch/blanc/articles/autoa1fibrations.pdf",
    "locator": "Theorem 5.4.5(7),(8), printed pages 5915-5919; field convention page 5889",
    "role": "mathematical_input",
    "inspected": "Actual primary statement, factorization proof, explicit generator proof, and split-root field hypotheses.",
    "printed_formula_correction": "The second coordinate in T and T0 omits the original v summand. C.4 is derived from uv=P(w), independently verified, and agrees with the plane-shear lift in the proof on page 5919.",
    "novelty_claim": False,
}
ASSURANCE = {
    "bounded_exact_controls_only": True,
    "finite_controls_prove_universal_quantifiers": False,
    "primary_generator_theorem_is_imported": True,
    "independent_replay": False,
    "independent_mathematical_review": False,
    "formal_verification": False,
    "external_review": False,
    "worldwide_priority": False,
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(repo):
    obj = json.loads((HERE / LEDGER).read_text())
    require(set(obj) == {"schema", "sources", "primary_input", "assurance"}, "Ledger keys")
    require(type(obj["schema"]) is int and obj["schema"] == 1, "Ledger schema")
    require(obj["primary_input"] == PRIMARY and obj["assurance"] == ASSURANCE, "Scope metadata")
    rows = obj["sources"]
    require(len(rows) == len(SOURCES) and {r["role"] for r in rows} == set(SOURCES), "Source roles")
    for row in rows:
        require(set(row) == {"role", "path", "scope", "sha256"}, "Source row keys")
        name, scope = SOURCES[row["role"]]
        require(row["path"] == str(PREFIX / name) and row["scope"] == scope, "Source path/scope")
        require(sha(repo / row["path"]) == row["sha256"], "Stale source: " + row["role"])
    return rows


spec = importlib.util.spec_from_file_location("hopf_sparse_arithmetic", HERE / "verify_hopf_lnd_orbit.py")
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
x, y, t, one = helper.X, helper.Y, helper.T, helper.ONE
IDENTITY = [x, y, t]
checks = 0


def nf(p):
    """Independent reduction t^2=1-2xy; helper instead removes xy."""
    out = {}
    for (a, b, e), coefficient in p.items():
        k, parity = divmod(e, 2)
        for j in range(k + 1):
            key = (a + j, b + j, parity)
            out[key] = out.get(key, F(0)) + coefficient * comb(k, j) * (-2) ** j
    return {key: F(v) for key, v in out.items() if v}


def add(*ps):
    return nf(helper.add(*ps))


def scale(p, coefficient):
    return helper.scale(p, coefficient)


def mul(a, b):
    return nf(helper.mul(a, b))


def power(p, n):
    out = one
    while n:
        if n & 1:
            out = mul(out, p)
        n //= 2
        if n:
            p = mul(p, p)
    return out


def evaluate(p, images):
    cache = [{0: one} for _ in images]
    out = {}
    for exponents, coefficient in p.items():
        term = scale(one, coefficient)
        for i, exponent in enumerate(exponents):
            if exponent not in cache[i]:
                cache[i][exponent] = power(images[i], exponent)
            term = mul(term, cache[i][exponent])
        out = add(out, term)
    return out


def equal(a, b, message):
    global checks
    delta = helper.add(a, helper.scale(b, -1))
    require(not nf(delta) and not helper.nf(delta), message)
    checks += 1


def degree(p):
    result = max((sum(m) for m in nf(p)), default=-1)
    require(result == helper.degree(p), "Independent quotient degrees disagree")
    return result


def degrees(triple):
    return tuple(degree(p) for p in triple)


def univariate(coefficients, value):
    out = {}
    for coefficient in reversed(coefficients):
        out = add(mul(out, value), scale(one, coefficient))
    return out


def shear(triple, coefficients):
    X, Y, T = triple
    Q = univariate(coefficients, X)
    return [X, add(Y, scale(mul(T, Q), -1), scale(mul(X, mul(Q, Q)), F(-1, 2))),
            add(T, mul(X, Q))]


def mono_shear(triple, q):
    return shear(triple, [0] * (q - 1) + [1])


def swap(triple):
    return [triple[1], triple[0], triple[2]]


def linear(triple, name):
    if name == "identity":
        return triple
    if name == "swap":
        return swap(triple)
    if name == "generic":
        return shear(swap(triple), [1])
    if name == "diagonal":
        return [scale(triple[0], 2), scale(triple[1], F(1, 2)), scale(triple[2], -1)]
    raise ValueError("Unknown linear map")


def theta(r, q, swapped=False):
    result = mono_shear(swap(mono_shear(IDENTITY, (r + 1) // 2)), q)
    return swap(result) if swapped else result


def permitted(a, c, b):
    r = min(a, c)
    return a >= 1 and c >= 1 and b >= 1 and a + c == 2 * b and r % 2 == 1 and b % r == 0


def preserve(triple, message):
    equal(add(mul(triple[2], triple[2]), scale(mul(triple[0], triple[1]), 2)), one, message)


def formal_shear_controls():
    # Free parameters U,V are invariant along x-shears. No quotient reduction.
    X, Y, T, U, V = [helper.var(i, 5) for i in range(5)]
    variables = [X, Y, T, U, V]
    quadratic = helper.add(helper.mul(T, T), helper.scale(helper.mul(X, Y), 2))

    def raw(parameter):
        return [X, helper.add(Y, helper.scale(helper.mul(T, parameter), -1),
                              helper.scale(helper.mul(X, helper.mul(parameter, parameter)), F(-1, 2))),
                helper.add(T, helper.mul(X, parameter)), U, V]

    require(helper.sub(quadratic, raw(U)) == quadratic, "Formal ambient quadratic preservation")
    for i in range(3):
        require(helper.sub(raw(U)[i], raw(helper.scale(U, -1))) == variables[i], "Formal inverse shear")
        require(helper.sub(raw(U)[i], raw(V)) == raw(helper.add(U, V))[i], "Formal shear addition")
    return {"free_variables": 5, "ambient_polynomial_identities": 7}


def coordinate_controls():
    rows = []
    for r in (1, 3, 5):
        for q in (1, 2, 3):
            for flipped in (False, True):
                triple = theta(r, q, flipped)
                wanted = (r, (2 * q - 1) * r, q * r)
                if flipped:
                    wanted = (wanted[1], wanted[0], wanted[2])
                require(degrees(triple) == wanted and permitted(*wanted), "Explicit realization degree")
                preserve(triple, "Explicit realization relation")
                rows.append({"r": r, "q": q, "swapped": flipped, "degrees_x_y_t": list(wanted)})
    words = []
    for q1, q2 in ((2, 2), (2, 3), (3, 2)):
        for middle in ("swap", "generic"):
            triple = mono_shear(IDENTITY, q1)
            triple = linear(triple, middle)
            # A lower-degree term tests more than monomial shears.
            coefficients = [F(1, 2)] + [0] * (q2 - 2) + [F(-2)]
            triple = shear(triple, coefficients)
            r = 2 * q1 - 1
            for final in ("identity", "swap", "generic"):
                result = linear(triple, final)
                wanted = (r, (2 * q2 - 1) * r, q2 * r)
                if final == "swap":
                    wanted = (wanted[1], wanted[0], wanted[2])
                if final == "generic":
                    wanted = ((2 * q2 - 1) * r,) * 3
                require(degrees(result) == wanted and permitted(*wanted), "Reduced-word/final-linear degrees")
                preserve(result, "Reduced word relation")
                words.append({"q1": q1, "q2": q2, "middle": middle, "final": final,
                              "degrees_x_y_t": list(wanted)})
    # Three nonlinear blocks: the written multiplicative induction has two steps.
    triple = mono_shear(swap(mono_shear(swap(mono_shear(IDENTITY, 2)), 2)), 2)
    require(degrees(triple) == (9, 27, 18), "Three nonlinear blocks")
    preserve(triple, "Three-block relation")
    # Build an inverse word independently, then reverse the actual forward steps.
    first = [0, 1]
    second = [F(1, 2), 0, -2]
    forward = shear(swap(shear(IDENTITY, first)), second)
    inverse = shear(swap(shear(IDENTITY, [-c for c in second])), [-c for c in first])
    require(max(degrees(forward)) == max(degrees(inverse)) == 15, "Inverse degree symmetry")
    restored = shear(swap(shear(forward, [-c for c in second])), [-c for c in first])
    for i in range(3):
        equal(restored[i], IDENTITY[i], "Inverse-word composition")
    cancellation = shear(shear(IDENTITY, [0, 1]), [0, -1])
    require(degrees(cancellation) == (1, 1, 1), "Unreduced inverse cancellation")
    for i in range(3):
        equal(cancellation[i], IDENTITY[i], "Cancellation identity")
    lower = shear(shear(IDENTITY, [0, 0, 1]), [0, 1, -1])
    require(degrees(lower) == (1, 3, 2), "Cancellation of highest shear coefficient")
    return {"explicit_realizations": rows, "reduced_two_block_words": words,
            "three_block_degrees": [9, 27, 18], "inverse_word_degree": 15,
            "cancellation_degrees": [[1, 1, 1], [1, 3, 2]]}


def hyperbolic(h, n, lam):
    C = power(add(one, scale(x, lam)), n + 1)
    U = scale(mul(add(one, scale(x, lam)),
                  add(*(scale(power(x, j), comb(2 * n, j + 1) * lam ** j)
                        for j in range(2 * n)))), F(1, 2 * n))
    numerator = add(mul(U, U), scale(one, -1))
    require(all(a >= 1 for a, b, e in numerator), "U^2-1 divisible by x")
    quotient = {(a - 1, b, e): coefficient for (a, b, e), coefficient in numerator.items()}
    return mul(mul(power(t, h), mul(C, power(U, h - 1))),
               add(scale(y, 2), scale(mul(power(t, 2), quotient), -1)))


def hyperbolic_controls():
    rows = []
    inputs = [(1, 1, 1, 1, 1, False), (1, 1, -2, 1, 2, False),
              (1, 2, 1, 1, 2, False), (3, 1, 1, 1, 1, False),
              (3, 1, -2, 1, 2, False), (1, 1, 1, 3, 1, False),
              (1, 1, 1, 1, 2, True)]
    for h, n, lam, r, q, flipped in inputs:
        triple = theta(r, q, flipped)
        a, c, b = degrees(triple)
        expected = (2 * h + 3) * n * a + (h + 2) * b
        actual = degree(evaluate(hyperbolic(h, n, F(lam)), triple))
        require(actual == expected, "Literal hyperbolic transformed degree")
        rows.append({"h": h, "n": n, "lambda": lam, "r": r, "q": q,
                     "swapped": flipped, "degree": actual})
    return rows


def divisor_test(D, alpha, beta):
    for r in range(1, D + 1, 2):
        if D % r:
            continue
        v = D // r
        if v >= alpha + beta and ((v - alpha) % beta == 0 or (v + alpha) % (2 * alpha + beta) == 0):
            return True
    return False


def isprime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def arithmetic_controls():
    cutoff = 300
    rows = []
    for h in (1, 3, 5, 7):
        for n in range(1, 5):
            alpha, beta = (2 * h + 3) * n, h + 2
            parameters = set()
            for r in range(1, cutoff // (alpha + beta) + 1, 2):
                for q in range(1, cutoff + 1):
                    for value in (r * (alpha + beta * q), r * ((2 * alpha + beta) * q - alpha)):
                        if value <= cutoff:
                            parameters.add(value)
            triples = set()
            for a in range(1, cutoff // alpha + 1):
                for b in range(1, cutoff // beta + 1):
                    value = alpha * a + beta * b
                    if value <= cutoff and permitted(a, 2 * b - a, b):
                        triples.add(value)
            divisors = {D for D in range(1, cutoff + 1) if divisor_test(D, alpha, beta)}
            require(parameters == triples == divisors, "Three independent finite spectrum enumerations")
            require(all(D % gcd(n, beta) == 0 for D in divisors), "Common-divisor obstruction")
            for prime in filter(isprime, range(1, cutoff + 1)):
                predicted = prime >= alpha + beta and ((prime - alpha) % beta == 0 or (prime + alpha) % (2 * alpha + beta) == 0)
                require((prime in divisors) == predicted, "Prime membership criterion")
                if gcd(n, beta) > 1:
                    require(prime not in divisors, "No prime for nontrivial gcd")
            rows.append({"h": h, "n": n, "cutoff": cutoff, "number_of_degrees": len(divisors)})
    require(not divisor_test(13, 5, 3) and not divisor_test(19, 5, 3) and divisor_test(73, 5, 3), "Smallest-orbit examples")
    return {"parameter_pairs": rows, "integer_tests": len(rows) * cutoff,
            "smallest_orbit_degrees_through_100": [D for D in range(1, 101) if divisor_test(D, 5, 3)]}


def negative_controls():
    rows = []

    def reject(false_claim, name):
        require(not false_claim, "Negative control survived: " + name)
        rows.append(name)

    good = shear(IDENTITY, [0, 1])
    malformed = [good[0], add(good[1], scale(y, -1)), good[2]]
    reject(not nf(add(mul(malformed[2], malformed[2]), scale(mul(malformed[0], malformed[1]), 2), scale(one, -1))),
           "Printed missing-v shear coordinate")
    bad_coefficient = [x, add(y, scale(mul(t, x), -1), scale(power(x, 3), -1)), add(t, power(x, 2))]
    reject(not nf(add(mul(bad_coefficient[2], bad_coefficient[2]), scale(mul(bad_coefficient[0], bad_coefficient[1]), 2), scale(one, -1))),
           "Changed quadratic shear coefficient")
    require(3 + 7 == 2 * 5 and all(d % 2 for d in (3, 7)), "Negative-control prerequisites")
    reject(permitted(3, 7, 5), "Parity and quadric degree relation incorrectly declared sufficient")
    require(degrees(theta(1, 2, True)) == (3, 1, 2), "Reflected positive control")
    reject(2 % 3 == 0, "Omitted reflected branch by requiring deg(x) to divide deg(t)")
    reject(divisor_test(16, 5, 3), "Even multiplier incorrectly allowed in the spectrum")
    cancelled = shear(shear(IDENTITY, [0, 1]), [0, -1])
    reject(max(degrees(cancelled)) == 9, "Unreduced word treated as a reduced degree product")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=HERE.parents[2])
    parser.add_argument("--output", type=Path, default=HERE / RECEIPT)
    args = parser.parse_args()
    sources = check_sources(args.repo_root.resolve())
    formal = formal_shear_controls()
    coordinates = coordinate_controls()
    hyperbolic_rows = hyperbolic_controls()
    arithmetic = arithmetic_controls()
    negatives = negative_controls()
    receipt = {
        "schema": 1, "status": "PASS", "python": platform.python_version(),
        "command": "python3 research/foundations/extended-geometry/" + CHECKER + " --repo-root .",
        "arithmetic": "Exact rational arithmetic; independent t^2 and xy quotient reductions.",
        "source_ledger_sha256": sha(HERE / LEDGER), "sources": sources,
        "formal_shear_controls": formal, "coordinate_controls": coordinates,
        "literal_hyperbolic_controls": hyperbolic_rows, "arithmetic_spectrum_controls": arithmetic,
        "quotient_equalities_checked_twice": checks, "negative_controls_rejected": negatives,
        "assurance": ASSURANCE,
        "scope": "Bounded exact controls. The classical generator theorem is imported, and the unrestricted all-word degree and field-extension assertions are proved in the written note.",
    }
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print("PASS: %d realizations; %d two-block words; three-block/inverse/cancellation controls; %d literal hyperbolic degrees; %d arithmetic spectrum tests; %d negative controls."
          % (len(coordinates["explicit_realizations"]), len(coordinates["reduced_two_block_words"]),
             len(hyperbolic_rows), arithmetic["integer_tests"], len(negatives)))
    print(args.output)


if __name__ == "__main__":
    main()
