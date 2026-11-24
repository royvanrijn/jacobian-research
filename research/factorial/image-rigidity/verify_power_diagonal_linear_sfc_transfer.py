#!/usr/bin/env python3
"""Exact finite controls for the written common-power linear SFC transfer.

Uses only the standard library. Actual multivariate powers are expanded
independently of formal series inversion by composition. Finite controls
do not prove inverse rigidity or extrapolate the unrestricted quantifiers.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys


DIRECTORY = Path("research/factorial/image-rigidity")
INSPECTED_BASE = "19eb9912fb60142f8b620dac24896931d4ff3e46"
SOURCE_PINS = {
    "POWER_DIAGONAL_LINEAR_SFC_TRANSFER.md": "40278040a10d4813c473e8256b868e4fb6d42f20912972b801136398436609d6",
    "RIGIDITY_CRITICAL_VALUES.tex": "8a2b6faa7996398d8c0474b34550a56cab6a771436a49fb2d3dff7a393e9c4dd",
    "sources/rigidity.typ": "8a6e152c797c227529e88572db140fcfb5b6c3ad61650c2a94fc6304cb334af6",
    "sources/edo-van-den-essen-1304.3956v2.pdf": "c6fd11479a33c02b5260dd54d82613530ef2448fc532336f653685a4761ea1e9",
    "sources/rocks-2015-author-preprint.pdf": "cee089adf02a4cc8add4bded8804e89bd51317108939c612ff2e3710bd43ecdc",
    "EVENTUAL_BINOMIAL_ALGEBRAIC_CONSEQUENCES.md": "1f200d3a3f060c41fae993a65ba2dc5376af81456a9b16a952ec1ca153901830",
    "UNRESTRICTED_SFC_MECHANISMS.tex": "f8b3dfc5947601fb002aee0ba942b0fc20b7ef46ca8e0d4f69b5a3e5e1361785",
}
COEFFICIENT_CASES = (
    (1,), (-2,), (0, -2), (1, -1), (2, -3),
    (0, 2, -3), (1, -2, 3), (1, 1, -2, 0),
    (0, 0, 0), (1, -1, 1, -1), (2, 0, 0, 0),
)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(repo):
    found = {}
    for name, expected in SOURCE_PINS.items():
        path = repo / DIRECTORY / name
        require(path.is_file(), "Missing written proof/source: " + name)
        actual = file_hash(path)
        require(actual == expected, "Changed written proof/source: " + name)
        found[name] = actual
    return found


def trim(poly):
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def multiply(a, b, degree):
    result = [0] * min(len(a) + len(b) - 1, degree + 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b[:max(0, degree + 1 - i)]):
                if bj:
                    result[i + j] += ai * bj
    return trim(result)


def compose(a, b, degree):
    result = [0]
    for coefficient in reversed(a):
        result = multiply(result, b, degree)
        result[0] += coefficient
    return trim(result)


def coefficient(poly, degree):
    return poly[degree] if degree < len(poly) else 0


def inverse_by_composition(a, degree):
    """Solve a(b)=z successively; no Lagrange formula is used here."""
    require(a[:2] == [0, 1], "Normalized inverse input")
    b = [0, 1] + [0] * (degree - 1)
    for k in range(2, degree + 1):
        b[k] = -coefficient(compose(a, b[:k], k), k)
    require(compose(a, b, degree) == [0, 1], "a composed with b is z")
    require(compose(b, a, degree) == [0, 1], "b composed with a is z")
    return b


def inverse_polynomial(q, coefficients):
    require(type(q) is int and q >= 1, "Positive integer common power")
    require(len(coefficients) >= 1, "At least one ambient variable")
    product = [1]
    for c in coefficients:
        product = multiply(product, [1] + [0] * (q - 1) + [-c],
                           q * len(coefficients))
    return [0] + product


def sparse_multiply(a, b):
    result = {}
    for alpha, ca in a.items():
        for beta, cb in b.items():
            gamma = tuple(x + y for x, y in zip(alpha, beta))
            result[gamma] = result.get(gamma, 0) + ca * cb
    return {gamma: c for gamma, c in result.items() if c}


def phase(common_powers, coefficients):
    require(len(common_powers) == len(coefficients) >= 1, "Phase dimensions")
    result = {}
    for j, c in enumerate(coefficients):
        if c:
            exponent = list(common_powers)
            exponent[j] += 1
            result[tuple(exponent)] = c
    return result


def actual_moments(polynomial, variables, maximum):
    power = {(0,) * variables: 1}
    moments = []
    for m in range(maximum + 1):
        moments.append(sum(c * math.prod(math.factorial(e) for e in alpha)
                           for alpha, c in power.items()))
        if m != maximum:
            power = sparse_multiply(power, polynomial)
    return moments


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for k in range(total + 1):
            for tail in compositions(total - k, length - 1):
                yield (k,) + tail


def direct_factorial_sum(q, coefficients, m):
    result = 0
    for multiindex in compositions(m, len(coefficients)):
        denominator = math.prod(math.factorial(k) for k in multiindex)
        numerator = math.factorial(m)
        require(numerator % denominator == 0, "Exact multinomial integer")
        result += (numerator // denominator
                   * math.prod(math.factorial(q * m + k)
                               for k in multiindex)
                   * math.prod(c ** k for c, k in zip(coefficients, multiindex)))
    return result


def family_checks():
    rows = []
    identities = structural = inactive = nonzero_windows = 0
    maximum = 5
    for q in (1, 2, 3, 4):
        for coefficients in COEFFICIENT_CASES:
            r = len(coefficients)
            active = tuple(c for c in coefficients if c)
            s = len(active)
            a = inverse_polynomial(q, coefficients)
            require(len(a) - 1 == q * s + 1, "Actual active polynomial degree")
            b = inverse_by_composition(a, q * maximum + 1)
            moments = actual_moments(phase([q] * r, coefficients), r, maximum)
            for ell, value in enumerate(b):
                if ell % q != 1 % q:
                    require(value == 0, "Exact structural inverse coefficient zero")
                    structural += 1
            for m, moment in enumerate(moments):
                require(moment == direct_factorial_sum(q, coefficients, m),
                        "Actual expansion agrees with factorial sum")
                factor = math.factorial(m) * math.factorial(q * m) ** r * (q * m + 1)
                require(moment == factor * b[q * m + 1],
                        "Actual moment agrees with independently inverted coefficient")
                identities += 1
            require(moments[0] == b[1] == 1, "m=0 identity including the zero phase")
            if 0 < s < r:
                reduced = actual_moments(phase([q] * s, active), s, maximum)
                for m in range(maximum + 1):
                    require(moments[m] == math.factorial(q * m) ** (r - s) * reduced[m],
                            "Inactive variables contribute exactly their own factorials")
                    inactive += 1
            if s:
                for start in range(1, maximum - s + 2):
                    require(any(moments[start:start + s]), "Finite actual nonzero window control")
                    nonzero_windows += 1
            else:
                require(moments[1:] == [0] * maximum, "All-zero phase boundary")
            rows.append({"q": q, "coefficients": list(coefficients), "active_terms": s,
                         "actual_polynomial_degree": len(a) - 1,
                         "actual_moments_through_five": moments,
                         "inverse_coefficients": b})
    return {"parameter_cases": len(rows), "moment_identities": identities,
            "structural_zero_coefficients": structural,
            "inactive_variable_identities": inactive,
            "finite_nonzero_window_controls": nonzero_windows,
            "cases": rows}


def window_checks():
    rows = []
    for q, s, n in ((1, 1, 1), (1, 4, 3), (2, 1, 1),
                    (2, 3, 1), (3, 2, 4), (5, 4, 7)):
        first, last = q * (n - 1) + 2, q * (n + s)
        selected = {q * m + 1 for m in range(n, n + s)}
        forced = {ell for ell in range(first, last + 1)
                  if ell % q != 1 % q or ell in selected}
        require(forced == set(range(first, last + 1)), "Every inverse degree is covered")
        require(first >= 2 and len(forced) == q * (s + 1) - 1 >= q * s,
                "Rigidity window starts at two and has sufficient length")
        rows.append({"q": q, "s": s, "start_moment": n,
                     "first_inverse_degree": first, "last_inverse_degree": last,
                     "forced_count": len(forced), "required_count": q * s})
    return rows


def cyclic_product(a, b):
    length = len(a)
    result = [0] * length
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[(i + j) % length] += ai * bj
    return result


def reduce_monic(poly, modulus):
    poly = trim(list(poly))
    while len(poly) >= len(modulus):
        shift, leading = len(poly) - len(modulus), poly[-1]
        for j, c in enumerate(modulus):
            poly[shift + j] -= leading * c
        trim(poly)
    return poly


def sharpness_checks():
    # Direct phase multiplication over Z[zeta]/Phi_s(zeta), s=2,3,4.
    cyclotomic = {2: [1, 1], 3: [1, 1, 1], 4: [1, 0, 1]}
    rows = []
    for q in (1, 2, 3):
        for s, modulus in cyclotomic.items():
            terms = {}
            for j in range(s):
                exponent, c = [q] * s, [0] * s
                exponent[j] += 1
                c[j] = 1
                terms[tuple(exponent)] = c
            power = {(0,) * s: [1] + [0] * (s - 1)}
            moments = []
            for m in range(1, s + 1):
                updated = {}
                for alpha, ca in power.items():
                    for beta, cb in terms.items():
                        gamma = tuple(x + y for x, y in zip(alpha, beta))
                        product = cyclic_product(ca, cb)
                        old = updated.setdefault(gamma, [0] * s)
                        for j in range(s):
                            old[j] += product[j]
                power = updated
                moment = [0] * s
                for alpha, c in power.items():
                    factor = math.prod(math.factorial(e) for e in alpha)
                    for j in range(s):
                        moment[j] += factor * c[j]
                reduced = reduce_monic(moment, modulus)
                expected = (0 if m < s else
                            math.factorial(s) * math.factorial(q * s) ** s * (q * s + 1))
                require(reduced == [expected], "Actual root-of-unity sharpness identity")
                moments.append(expected)
            rows.append({"q": q, "s": s, "cyclotomic_polynomial": modulus,
                         "positive_moments_through_s": moments})
    return rows


def adjunction_checks():
    rows = []
    for q, coefficients, lambdas, mark in (
            (2, (1, -1), (2,), (3,)),
            (3, (1, 0, -2), (2, -3), (1, 2))):
        r, t, maximum = len(coefficients), len(lambdas), 4
        f = phase([q] * r, coefficients)
        g = {alpha + (0,) * t: c for alpha, c in f.items()}
        for j, c in enumerate(lambdas):
            exponent = [0] * (r + t)
            exponent[r + j] = 1
            g[tuple(exponent)] = c
        f_moments = actual_moments(f, r, maximum)
        g_moments = actual_moments(g, r + t, maximum)
        denominator = [1]
        for c in lambdas:
            denominator = multiply(denominator, [1, -c], t)
        for m in range(maximum + 1):
            recovered = sum(Fraction(denominator[j] * g_moments[m - j],
                                     math.factorial(m - j))
                            for j in range(min(m, t) + 1))
            require(recovered == Fraction(f_moments[m], math.factorial(m)),
                    "Actual independent-linear finite filter")
        marked = {alpha + mark: c for alpha, c in g.items()}
        marked_moments = actual_moments(marked, r + t + len(mark), maximum)
        require(all(marked_moments[m] == g_moments[m] *
                    math.prod(math.factorial(exponent * m) for exponent in mark)
                    for m in range(maximum + 1)), "Actual independent monomial factor")
        rows.append({"q": q, "coefficients": coefficients, "linear_coefficients": lambdas,
                     "independent_factor_exponents": mark,
                     "checked_orders": list(range(maximum + 1))})
    return rows


def negative_controls():
    rejected = []

    def reject(name, action):
        try:
            action()
        except (ValueError, TypeError):
            rejected.append(name)
        else:
            raise ValueError("Invalid control accepted: " + name)

    reject("q-zero-in-positive-q-identity", lambda: inverse_polynomial(0, (1, 2)))
    reject("noninteger-q", lambda: inverse_polynomial(Fraction(3, 2), (1, 2)))
    reject("boolean-q", lambda: inverse_polynomial(True, (1, 2)))
    reject("empty-ambient-variable-list", lambda: inverse_polynomial(2, ()))
    q, coefficients, m = 2, (2, -3), 1
    a = inverse_polynomial(q, coefficients)
    b = inverse_by_composition(a, 7)
    actual = actual_moments(phase([q] * 2, coefficients), 2, m)[m]
    wrong = math.factorial(m) * math.factorial(q * m) ** 2 * (m + 1) * b[q * m + 1]
    reject("m-plus-one-instead-of-qm-plus-one", lambda: require(actual == wrong, "Wrong normalization"))
    altered = b.copy()
    altered[q + 1] += 1
    reject("altered-inverse-coefficient", lambda: require(compose(a, altered, 7) == [0, 1],
                                                         "Wrong inverse coefficient"))
    full = actual_moments(phase([2, 2, 2], (1, -2, 0)), 3, 1)[1]
    reduced = actual_moments(phase([2, 2], (1, -2)), 2, 1)[1]
    reject("omitted-inactive-variable-factor", lambda: require(full == reduced, "Missing factorial factor"))
    unequal = actual_moments(phase([1, 2], (1, 1)), 2, 1)[1]
    uniform = actual_moments(phase([1, 1], (1, 1)), 2, 1)[1]
    reject("unsupported-unequal-common-exponents", lambda: require(unequal == uniform,
                                                                  "Different actual support"))
    return rejected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    sources = check_sources(repo)
    family = family_checks()
    windows = window_checks()
    sharpness = sharpness_checks()
    adjunction = adjunction_checks()
    negatives = negative_controls()
    require(sources == check_sources(repo), "Written inputs unchanged during replay")
    receipt = {
        "schema": 1, "status": "PASS", "inspected_retained_base": INSPECTED_BASE,
        "theorem_scope": "Every common monomial power times a linear form, any active support size, every positive starting index, every characteristic-zero field",
        "standard_library_only": True, "python_version": sys.version.split()[0],
        "source_sha256": sources, "checker_sha256": file_hash(Path(__file__)),
        "finite_identity_controls": family, "window_index_controls": windows,
        "root_of_unity_sharpness_controls": sharpness,
        "independent_adjunction_controls": adjunction, "rejected_controls": negatives,
        "proof_boundary": {
            "universal_theorem_from_written_proof": True,
            "unrestricted_SFC_proved": False,
            "general_variable_identification_proved": False,
            "ordinary_binomial_eventual_theorem_used": False,
            "rigidity_proof_replaced_by_finite_tests": False,
            "formal_verification": False, "external_peer_review": False,
            "literature_priority_claim": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "PASS", "parameter_cases": family["parameter_cases"],
                      "moment_identities": family["moment_identities"],
                      "window_controls": len(windows), "sharpness_controls": len(sharpness),
                      "negative_controls": len(negatives), "output": str(args.output)}))


if __name__ == "__main__":
    main()
