#!/usr/bin/env python3
"""Bounded exact controls for the written full Hopf-orbit degree theorem.

The unrestricted automorphism and field assertions are written proofs, not
conclusions of these deterministic rational examples. Standard library only.
"""
import argparse
import hashlib
import importlib.util
import json
import platform
from fractions import Fraction as Q
from math import comb, factorial, prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
PREFIX = Path("research/foundations/extended-geometry")
NOTE = "HOPF_AUTOMORPHISM_ORBIT_SPECTRA.md"
CHECKER = "verify_hopf_automorphism_orbit_spectra.py"
LEDGER = "hopf_automorphism_orbit_spectra_source_ledger.json"
RECEIPT = "hopf_automorphism_orbit_spectra_verification.json"
SOURCES = {
    "written_proof": (NOTE, "Self-contained orbit spectra and geometric multiplicity classification within the displayed higher-start family; explicitly dependent Image interpretation and coarse boundary."),
    "exact_checker": (CHECKER, "Bounded exact controls only; deterministic cases are retained in this file."),
    "arithmetic_helper": ("verify_hopf_lnd_orbit.py", "Rational sparse arithmetic, one quotient normal form, formal moments, literal original lifts."),
    "least_lift_and_power_descent": ("HOPF_LND_ORBIT_DEGREE_PROPAGATION_AND_POWER_DESCENT.md", "Sections 4-5: least homogenization, radial primitivity, scalar powers in the full original ring."),
    "quartic_moments": ("TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md", "Theorem 1.1 and sphere normalization: all-order quartic purity and fixed-marker survival."),
    "sextic_moments": ("OFFCENTER_CUBIC_ALL_MULTIPLES_OF_THREE_IMAGE_FAMILY.md", "N1-N4 and N6-N8: regular seeds, all-order purity and fixed-marker survival."),
    "hyperbolic_moments": ("HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md", "H1-H9 and Section 6, H18-H22: cubic-height and every odd starting exponent; all-order purity, odd marked vanishing and even marked survival."),
    "automorphism_trace": ("RANK_ONE_TRACE_AUTOMORPHISM_INVARIANCE_ALL_IMAGE_PAIRS.md", "Sections 3-5: all polynomial automorphisms preserve the normalized trace, all characteristic-zero extensions."),
    "linear_image_safety": ("PRIME_BALANCED_IMAGE_FAILURES_ARE_ABSOLUTELY_NONPOWER.md", "LinearSampledVanishingNilpotenceAndSafety: every pure balanced linear form is safe for all original markers."),
    "cubic_image_safety": ("TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md", "Opening complete-cubic theorem and 'Complete classification and the original Image conclusion': all original markers and characteristic-zero fields; lower balanced degrees covered there."),
    "higher_pair_witnesses": ("ALL_BALANCED_DEGREES_IMAGE_WITNESSES_AT_LEAST_THREE_PAIRS.md", "AllBalancedDegreesThreeOrMorePairs: explicit primitive witnesses for every n>=3,D>=2; N1 and Sections 1-4."),
}
ASSURANCE = {
    "bounded_exact_controls_only": True,
    "finite_controls_prove_universal_quantifiers": False,
    "independent_replay": False,
    "independent_mathematical_review": False,
    "formal_verification": False,
    "external_review": False,
    "worldwide_priority": False,
}
PRIMARY = {
    "authors": "Jeremy Blanc and Adrien Dubouloz",
    "title": "Automorphisms of A1-fibered affine surfaces",
    "publication": "Transactions of the AMS 363 (2011), 5887-5924",
    "url": "https://algebra.dmi.unibas.ch/blanc/articles/autoa1fibrations.pdf",
    "locator": "Theorem 5.4.5(7) and its proof, printed pp. 5915-5918",
    "scope": "Inspected amalgamated-product description for uv=P(w), deg P>=2, including the smooth split quadric. Context only; no imported logical premise or novelty claim.",
    "local_pdf_download_or_hash_claim": False,
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(repo):
    ledger = json.loads((HERE / LEDGER).read_text(encoding="utf-8"))
    require(set(ledger) == {"schema", "sources", "primary_context", "assurance"}, "Ledger keys")
    require(type(ledger["schema"]) is int and ledger["schema"] == 1, "Ledger schema")
    require(ledger["assurance"] == ASSURANCE, "Assurance changed")
    require(ledger["primary_context"] == PRIMARY, "Primary context changed")
    rows = ledger["sources"]
    require(isinstance(rows, list) and len(rows) == len(SOURCES), "Source cardinality")
    require({row["role"] for row in rows} == set(SOURCES), "Exact source roles")
    for row in rows:
        require(set(row) == {"role", "path", "scope", "sha256"}, "Source row keys")
        name, scope = SOURCES[row["role"]]
        require(row["path"] == str(PREFIX / name) and row["scope"] == scope, "Source path/scope")
        path = repo / row["path"]
        require(path.is_file() and sha(path) == row["sha256"], "Stale source: " + row["role"])
    return rows


spec = importlib.util.spec_from_file_location("hopf_exact_arithmetic", HERE / "verify_hopf_lnd_orbit.py")
h = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h)
add, mul, scale, power, sub = h.add, h.mul, h.scale, h.power, h.sub
x, y, t, one = h.X, h.Y, h.T, h.ONE
relation = add(power(t, 2), scale(mul(x, y), 2), scale(one, -1))
checks = 0
negative_controls = []


def equal(a, b, message):
    """Check quotient equality through TWO different leading-term reductions."""
    global checks
    delta = add(a, scale(b, -1))
    require(not h.nf(delta) and not reduce_t_squared(delta), message)
    checks += 1


def reduce_t_squared(p):
    """Independent normal form t^2=1-2xy, unlike helper's xy reduction."""
    out = {}
    for (a, b, ell), c in p.items():
        k, parity = divmod(ell, 2)
        for r in range(k + 1):
            key = (a + r, b + r, parity)
            out[key] = out.get(key, Q(0)) + c * comb(k, r) * (-2) ** r
    return {key: c for key, c in out.items() if c}


def degree(p):
    value = max((sum(m) for m in reduce_t_squared(p)), default=-1)
    require(value == h.degree(p), "Independent quotient degrees disagree")
    return value


def divided_by_x(p):
    require(all(m[0] >= 1 for m in p), "Not coefficientwise divisible by x")
    return {(a - 1, b, c): v for (a, b, c), v in p.items()}


def quartic():
    return mul(add(one, x), add(y, scale(power(t, 2), -1),
                               scale(mul(x, power(t, 2)), Q(-1, 2))))


def sextic_family(j):
    v = add(t, one, scale(power(x, j), -1))
    numerator = add(mul(power(add(one, x), 4), power(v, 3)),
                    scale(mul(power(add(one, x), 2), v), -4))
    at_zero = mul(add(power(t, 2), scale(one, -1)), add(t, scale(one, 3)))
    regular = add(divided_by_x(add(numerator, scale(at_zero, -1))),
                  scale(mul(y, add(t, scale(one, 3))), -2))
    return regular, numerator


def hyperbolic_data(n, lam=Q(1)):
    base = add(one, scale(x, lam))
    c = power(base, n + 1)
    s = add(*(scale(power(x, r), comb(2 * n, r + 1) * lam ** r) for r in range(2 * n)))
    u = scale(mul(base, s), Q(1, 2 * n))
    return c, u


def hyperbolic(n, lam=Q(1)):
    c, u = hyperbolic_data(n, lam)
    regular = mul(mul(t, c), add(scale(y, 2),
                    scale(mul(power(t, 2), divided_by_x(add(power(u, 2), scale(one, -1)))), -1)))
    numerator = mul(mul(t, c), add(one, scale(mul(power(t, 2), power(u, 2)), -1)))
    return regular, numerator


def shear(q):
    return [x, add(y, scale(mul(t, q), -1), scale(mul(x, power(q, 2)), Q(-1, 2))),
            add(t, mul(x, q))]


def compose(left, right):
    """Images for left after right as algebra homomorphisms."""
    return [h.nf(sub(v, left)) for v in right]


def reject(ok, label):
    require(not ok, "Corrupted control accepted: " + label)
    negative_controls.append(label)


def eval_curve(p, at):
    return sum((c * ((1 - at * at) / 2) ** b * at ** ell
                for (a, b, ell), c in p.items()), Q(0))


def curve_derivative(p, at):
    result = Q(0)
    yy = (1 - at * at) / 2
    for (a, b, ell), c in p.items():
        if b:
            result -= c * b * yy ** (b - 1) * at ** (ell + 1)
        if ell:
            result += c * ell * yy ** b * at ** (ell - 1)
    return result


def weak_compositions(total, n):
    if n == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for tail in weak_compositions(total - first, n - 1):
                yield (first,) + tail


def run_controls():
    f = quartic()
    g, _ = sextic_family(1)
    p, _ = hyperbolic(1)
    equal(scale(mul(x, f), 2), add(one, x, scale(mul(power(add(one, x), 3), power(t, 2)), -1)),
          "Quartic numerator identity")
    for j in range(1, 5):
        gj, numerator = sextic_family(j)
        equal(mul(x, gj), numerator, "Regular sextic numerator")
        equal(sub(g, shear(add(one, scale(power(x, j - 1), -1)))), gj, "One sextic orbit")
        require(degree(gj) == 3 * j + 3, "Sextic attained degree")
    for n in range(1, 4):
        for lam in (Q(1), Q(2)):
            pn, numerator = hyperbolic(n, lam)
            equal(mul(x, pn), numerator, "Regular hyperbolic numerator")
            require(degree(pn) == 5 * n + 3, "Hyperbolic minimum attained")
        p1, _ = hyperbolic(n, Q(1))
        p2, _ = hyperbolic(n, Q(2))
        equal(sub(p1, [scale(x, 2), scale(y, Q(1, 2)), t]), scale(p2, Q(1, 2)),
              "Hyperbolic parameter scaling")
    xx, yy, tt, q, r = [h.var(i, 5) for i in range(5)]

    def formal_shear(v):
        return [xx, add(yy, scale(mul(tt, v), -1), scale(mul(xx, power(v, 2)), Q(-1, 2))),
                add(tt, mul(xx, v)), q, r]
    form = add(power(tt, 2), scale(mul(xx, yy), 2))
    require(sub(form, formal_shear(q)) == form, "Universal formal shear relation")
    require([sub(v, formal_shear(q)) for v in formal_shear(scale(q, -1))] == [xx, yy, tt, q, r],
            "Universal formal inverse")
    require([sub(v, formal_shear(q)) for v in formal_shear(r)] == formal_shear(add(q, r)),
            "Universal formal shear composition")
    identity, swap, sh = [x, y, t], [y, x, t], shear(x)
    sh_inv = shear(scale(x, -1))
    maps = {
        "identity": (identity, identity),
        "shear_q=x": (sh, sh_inv),
        "swap": (swap, swap),
        "shear_swap_shear": (compose(sh, compose(swap, sh)),
                              compose(sh_inv, compose(swap, sh_inv))),
    }
    degree_rows = []
    for name, (sigma, inverse) in maps.items():
        equal(sub(relation, sigma), {}, name + " relation")
        for v, original in zip(compose(sigma, inverse), identity):
            equal(v, original, name + " inverse")
        a, b = degree(sigma[0]), degree(sigma[2])
        hh = degree(sub(add(t, one, scale(x, -1)), sigma))
        require(min(a, b, hh) >= 1, "Positive coordinate degrees")
        for seed_name, seed, expected in (("quartic", f, 2 * (a + b)),
                                         ("sextic", g, 3 * (a + hh)),
                                         ("hyperbolic_n=1", p, 5 * a + 3 * b)):
            actual = degree(sub(seed, sigma))
            require(actual == expected, "Full-orbit degree formula " + name + seed_name)
            degree_rows.append({"map": name, "seed": seed_name, "degree": actual,
                                "a": a, "b": b, "h1": hh})
    for rvalue in range(1, 5):
        require(degree(sub(f, shear(power(x, rvalue - 1)))) == 2 * rvalue + 2,
                "Every displayed even degree attained")
    g2, _ = sextic_family(2)
    cancellation_map = shear(add(x, scale(one, -1)))
    equal(sub(g2, cancellation_map), g, "Degree-nine member descends to degree six")
    actual_h = degree(sub(add(t, one, scale(power(x, 2), -1)), cancellation_map))
    guessed_h = max(degree(cancellation_map[2]), 2 * degree(cancellation_map[0]))
    require(actual_h == 1 and guessed_h == 2, "Cancellation control was not active")
    reject(3 * (1 + guessed_h) == degree(sub(g2, cancellation_map)),
           "max(b,ja) incorrectly substituted for actual h_j")
    moment_rows = []
    for name, seed in (("quartic", f), ("sextic", g), ("hyperbolic_n=1", p)):
        for map_name in ("identity", "shear_q=x", "swap"):
            sigma = maps[map_name][0]
            transformed, marker = h.nf(sub(seed, sigma)), sigma[0]
            d, e = degree(transformed), degree(marker)
            literal, literal_marker = h.literal_lift(transformed, d), h.literal_lift(marker, e)
            for m in (1, 2):
                raised = power(transformed, m)
                pure, marked = h.lam(raised), h.lam(mul(marker, raised))
                require(pure == 0, "Pure bounded moment")
                require(marked == h.lam(mul(x, power(seed, m))), "Transported bounded marker moment")
                if name == "quartic":
                    expected = Q(factorial(m), prod(range(1, 2 * m + 2, 2)))
                    require(marked == expected, "Quartic formula")
                elif name == "sextic":
                    expected = (-8) ** m * Q(2 ** m * factorial(m), prod(m + 2 * r + 1 for r in range(m + 1)))
                    require(marked == expected, "Sextic formula")
                else:
                    require(marked == (Q(0) if m == 1 else Q(-4, 105)), "Hyperbolic parity/survival formula")
                literal_power = h.balpow(literal, m)
                require(h.contract(literal_power, m * d) == 0, "Original pure contraction")
                require(h.contract(h.balmul(literal_marker, literal_power), m * d + e)
                        == factorial(m * d + e + 1) * marked, "Original marked normalization")
                moment_rows.append({"seed": name, "map": map_name, "m": m,
                                    "marked": str(marked), "literal_lift_degree": d})
    simple_roots = [("quartic", f, Q(1, 2), Q(-4)),
                    ("sextic", g, Q(0), Q(-16))]
    simple_roots += [("hyperbolic_n=" + str(n), hyperbolic(n)[0], Q(0), Q(2 ** (n + 1))) for n in range(1, 4)]
    for name, seed, root, derivative in simple_roots:
        require(eval_curve(seed, root) == 0 and curve_derivative(seed, root) == derivative != 0,
                "Simple-root nonpower control: " + name)
    inflated = add(f, mul(relation, power(x, 7)))
    require(max(map(sum, inflated)) == 9 and degree(inflated) == 4, "Representative-degree trap")
    radial_cases = 0
    for pairs in range(1, 5):
        for d in range(4):
            for alpha in weak_compositions(d, pairs):
                for k in range(4):
                    direct = sum((Q(factorial(k), prod(factorial(c) for c in gamma))
                                  * prod(factorial(a + c) for a, c in zip(alpha, gamma))
                                  for gamma in weak_compositions(k, pairs)), Q(0))
                    expected = prod(range(pairs + d, pairs + d + k)) * prod(factorial(a) for a in alpha)
                    require(direct == expected, "Matching original radial identity")
                    radial_cases += 1
                if pairs > 1 and d > 0:
                    beta = tuple(alpha[(i + 1) % pairs] for i in range(pairs))
                    if beta != alpha:
                        for gamma in weak_compositions(2, pairs):
                            require(tuple(a + c for a, c in zip(alpha, gamma))
                                    != tuple(b + c for b, c in zip(beta, gamma)), "Nonmatching radial term")
    wrong_shear = [x, add(y, mul(t, x), scale(power(x, 3), Q(-1, 2))), add(t, power(x, 2))]
    reject(not reduce_t_squared(sub(relation, wrong_shear)), "Changed shear sign")
    reject(not reduce_t_squared(add(scale(mul(x, add(f, x)), 2),
               scale(add(one, x, scale(mul(power(add(one, x), 3), power(t, 2)), -1)), -1))),
           "Changed seed term in numerator identity")
    reject(degree(f) == 2 * (degree(x) + degree(t)) + 1, "Incorrect degree offset")
    square = power(f, 2)
    reject(eval_curve(square, Q(1, 2)) == 0 and curve_derivative(square, Q(1, 2)) != 0,
           "False simple-root certificate for a square")
    return {
        "quotient_equalities_checked_by_two_reductions": checks,
        "universal_formal_shear_identities": 3,
        "degree_formula_cases": degree_rows,
        "quotient_cancellation_case": {"before": 9, "after": 6, "actual_h2": actual_h, "incorrect_max": guessed_h},
        "moment_and_original_contraction_cases": moment_rows,
        "simple_root_cases": len(simple_roots),
        "radial_matching_monomial_cases": radial_cases,
        "negative_controls_rejected": negative_controls.copy(),
    }


def univariate(p):
    """Convert an x-only sparse polynomial to rational coefficients."""
    require(all(b == ell == 0 for a, b, ell in p), "Expected x-only polynomial")
    out = [Q(0)] * (max((a for a, b, ell in p), default=-1) + 1)
    for (a, b, ell), c in p.items():
        out[a] += c
    return trim(out)


def trim(coefficients):
    coefficients = list(coefficients)
    while coefficients and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def univariate_divmod(f, g):
    f, g = trim(f), trim(g)
    require(bool(g), "Division by zero polynomial")
    quotient = [Q(0)] * max(0, len(f) - len(g) + 1)
    while len(f) >= len(g):
        degree_difference = len(f) - len(g)
        c = f[-1] / g[-1]
        quotient[degree_difference] += c
        for i, v in enumerate(g):
            f[i + degree_difference] -= c * v
        f = trim(f)
    return trim(quotient), f


def univariate_gcd(f, g):
    f, g = trim(f), trim(g)
    while g:
        _, remainder = univariate_divmod(f, g)
        f, g = g, remainder
    return [v / f[-1] for v in f] if f else []


def univariate_valuation(f, divisor):
    require(bool(f) and len(divisor) >= 2, "Invalid polynomial valuation")
    order = 0
    while True:
        quotient, remainder = univariate_divmod(f, divisor)
        if remainder:
            return order, f
        f, order = quotient, order + 1


def higher_hyperbolic(start_h, n, lam=Q(1)):
    require(start_h >= 1 and start_h % 2 == 1 and n >= 1 and lam != 0,
            "Odd-start parameter domain")
    c, u = hyperbolic_data(n, lam)
    seed = mul(power(mul(t, u), start_h - 1), hyperbolic(n, lam)[0])
    numerator = mul(mul(power(t, start_h), mul(c, power(u, start_h - 1))),
                    add(one, scale(mul(power(t, 2), power(u, 2)), -1)))
    return seed, numerator


def clear_branch_value(p, u, sign):
    """Compute U^degree_t(p) * p(x, sign/U) without rational functions."""
    require(all(b == 0 for a, b, ell in p), "Branch chart has y terms")
    height = max(ell for a, b, ell in p)
    return add(*(scale(mul(power(x, a), power(u, height - ell)), c * sign ** ell)
                 for (a, b, ell), c in p.items()))


def recover_parameters(multiplicities):
    distinct = sorted(set(multiplicities), reverse=True)
    require(len(distinct) >= 2, "Need two distinct geometric multiplicities")
    return distinct[1], distinct[0] - distinct[1]


def family_divisor_parameters(minimum_degree):
    require(minimum_degree >= 2, "Orbit-count domain is D>=2")
    m = 2 * minimum_degree - 1
    divisors_one_mod_four = [v for v in range(1, m + 1, 4) if m % v == 0]
    pairs = {((v - 3) // 2, (m // v - 1) // 2)
             for v in divisors_one_mod_four if 5 <= v < m}
    count = len(divisors_one_mod_four) - 1 - int(minimum_degree % 2 == 1)
    require(len(pairs) == count, "Divisor exclusions disagree with orbit count")
    return pairs, count, len(divisors_one_mod_four)


def orbit_count_controls():
    # Independent complete parameter enumeration, using only the degree formula.
    # n>=1 implies D>=3h+5; for each h the degree increases strictly with n.
    limit = 512
    parameter_rows = {d: set() for d in range(2, limit + 1)}
    for start_h in range(1, (limit - 5) // 3 + 1, 2):
        n = 1
        while (d := (2 * start_h + 3) * n + start_h + 2) <= limit:
            parameter_rows[d].add((start_h, n))
            n += 1
    enumeration = []
    for d, pairs in parameter_rows.items():
        divisor_pairs, count, tau_one = family_divisor_parameters(d)
        require(divisor_pairs == pairs and count == len(pairs),
                "Divisor bijection disagrees with independent parameter enumeration")
        enumeration.append({"D": d, "tau_one": tau_one, "count": count,
                            "parameters_h_n": [list(pair) for pair in sorted(pairs)]})
    reject(family_divisor_parameters(13)[2] - 1 == len(parameter_rows[13]),
           "Omit the full divisor correction at odd minimum degree")
    samples = {2, 3, 5, 7, 8, 13, 23, 63, 313, 512}
    return {
        "minimum_degree_interval_inclusive": [2, limit],
        "degrees_checked": len(enumeration),
        "complete_parameter_bound": "Odd 1<=h<=floor((512-5)/3), n>=1 until (2h+3)n+h+2>512.",
        "parameter_pairs_checked": sum(len(pairs) for pairs in parameter_rows.values()),
        "nonempty_minimum_degrees": sum(bool(pairs) for pairs in parameter_rows.values()),
        "sample_counts": [row for row in enumeration if row["D"] in samples],
        "enumeration_sha256": hashlib.sha256(json.dumps(enumeration, sort_keys=True).encode()).hexdigest(),
        "scope": "Complete finite integer enumeration only; equation (21) is proved for all D>=2 in the note.",
    }


def run_higher_start_controls():
    before_checks, before_negatives = checks, len(negative_controls)
    degree_rows, divisor_rows, moment_rows = [], [], []
    parameters = [(start_h, n, lam) for start_h in (1, 3, 5)
                  for n in (1, 2) for lam in (Q(1), Q(2))]
    parameters.append((1, 4, Q(1)))
    maps = {"identity": [x, y, t], "shear_q=x": shear(x), "swap": [y, x, t]}
    signatures = {}
    for start_h, n, lam in parameters:
        seed, direct_numerator = higher_hyperbolic(start_h, n, lam)
        c, u = hyperbolic_data(n, lam)
        equal(mul(x, seed), direct_numerator, "Higher-start regular numerator")
        for name, sigma in maps.items():
            a, b = degree(sigma[0]), degree(sigma[2])
            actual = degree(sub(seed, sigma))
            expected = (2 * start_h + 3) * n * a + (start_h + 2) * b
            require(actual == expected, "Higher-start full-Aut degree law")
            degree_rows.append({"h": start_h, "n": n, "lambda": str(lam),
                                "map": name, "a": a, "b": b, "degree": actual})
        # Obtain the chart numerator from the REGULAR seed, not its proposed factors.
        chart = h.chart(mul(x, seed))
        require(all(a >= 0 for a, ell in chart), "Uncleared chart denominator")
        numerator = {(a, 0, ell): v for (a, ell), v in chart.items()}
        require(numerator == direct_numerator, "Independent chart factorization")
        order_t = min(ell for a, b, ell in numerator)
        in_u_coordinate = sub(numerator, [scale(add(x, scale(one, -1)), 1 / lam), y, t])
        order_at_u_zero = min(a for a, b, ell in in_u_coordinate)
        u_coefficients = univariate(u)
        derivative = [Q(i) * v for i, v in enumerate(u_coefficients) if i > 0]
        require(univariate_gcd(u_coefficients, derivative) == [Q(1)], "Repeated root of U")
        require(u_coefficients[0] == 1, "A root of U lies on the excluded chart boundary")
        v, remainder = univariate_divmod(u_coefficients, [Q(1), lam])
        require(not remainder and len(v) - 1 == 2 * n - 1, "Wrong remaining root polynomial")
        require(univariate_gcd(v, [Q(1), lam]) == [Q(1)], "Overlapping root families")
        t_coefficients = []
        for ell in sorted({ell for a, b, ell in numerator}):
            part = {(a, 0, 0): value for (a, b, e), value in numerator.items() if e == ell}
            t_coefficients.append(univariate(part))
        measured_orders = [univariate_valuation(coefficient, v)[0] for coefficient in t_coefficients]
        root_order = min(measured_orders)
        low_order, residual = univariate_valuation(t_coefficients[0], v)
        require(low_order == root_order and univariate_gcd(v, residual) == [Q(1)],
                "Geometric roots have nonuniform multiplicities")
        for sign in (-1, 1):
            require(not clear_branch_value(numerator, u, sign), "Missing branch factor")
            derivative_value = clear_branch_value(h.deriv(numerator, 2), u, sign)
            require(derivative_value == scale(mul(c, power(u, start_h + 1)), -2)
                    and bool(derivative_value), "Branch multiplicity is not one")
            boundary = sub(seed, [{}, y, scale(one, sign)])
            slope = u_coefficients[1]
            require(boundary == scale(add(y, scale(one, -slope)), 2 * sign)
                    and bool(boundary), "Missing vertical zero divisor")
        multiplicities = [order_at_u_zero, order_t, 1, 1]
        multiplicities += [root_order] * (len(v) - 1) if root_order else []
        multiplicities.sort(reverse=True)
        require(recover_parameters(multiplicities) == (start_h, n), "Multiplicity parameter recovery")
        signatures[(start_h, n, lam)] = multiplicities
        divisor_rows.append({"h": start_h, "n": n, "lambda": str(lam),
                             "geometric_multiplicities": multiplicities,
                             "remaining_simple_U_roots": len(v) - 1,
                             "branch_multiplicities": [1, 1], "vertical_components": 0})
    p1 = higher_hyperbolic(3, 2, Q(1))[0]
    p2 = higher_hyperbolic(3, 2, Q(2))[0]
    equal(sub(p1, [scale(x, 2), scale(y, Q(1, 2)), t]), scale(p2, Q(1, 2)),
          "Higher-start lambda scaling")
    # One additional higher-start moment/lift case with a nonlinear transport.
    seed = higher_hyperbolic(3, 1)[0]
    for map_name in ("identity", "shear_q=x"):
        sigma = maps[map_name]
        transformed, marker = h.nf(sub(seed, sigma)), sigma[0]
        d, e = degree(transformed), degree(marker)
        literal, literal_marker = h.literal_lift(transformed, d), h.literal_lift(marker, e)
        for m in (1, 2):
            raised = power(transformed, m)
            marked = h.lam(mul(marker, raised))
            require(h.lam(raised) == 0 and marked == (Q(0) if m == 1 else Q(-4, 693)),
                    "Higher-start exact moment formula")
            literal_power = h.balpow(literal, m)
            require(h.contract(literal_power, m * d) == 0, "Higher-start original pure contraction")
            require(h.contract(h.balmul(literal_marker, literal_power), m * d + e)
                    == factorial(m * d + e + 1) * marked, "Higher-start original marked normalization")
            moment_rows.append({"h": 3, "n": 1, "map": map_name, "m": m,
                                "marked": str(marked), "literal_lift_degree": d})
    same_degree = [degree(higher_hyperbolic(1, 4)[0]), degree(higher_hyperbolic(3, 2)[0])]
    require(same_degree == [23, 23], "Equal-minimum example")
    reject(signatures[(1, 4, Q(1))] == signatures[(3, 2, Q(1))],
           "Equal orbit minima incorrectly treated as sufficient for equivalence")
    squared_u = univariate(power(hyperbolic_data(2)[1], 2))
    squared_derivative = [Q(i) * v for i, v in enumerate(squared_u) if i > 0]
    reject(univariate_gcd(squared_u, squared_derivative) == [Q(1)],
           "Repeated U roots incorrectly treated as simple")
    reject(bool(sub(mul(x, seed), [{}, y, one])),
           "Vertical zero component omitted after corrupt multiplication by x")
    counts = orbit_count_controls()
    same_minimum_families = []
    for rvalue in range(2, 7):
        common_degree = (5 ** rvalue + 1) // 2
        pairs = [((5 ** j - 3) // 2, (5 ** (rvalue - j) - 1) // 2)
                 for j in range(1, rvalue)]
        require(len(set(pairs)) == rvalue - 1, "Repeated parameter pair in common-degree construction")
        require(all(start_h >= 1 and start_h % 2 == 1 and n >= 1
                    and (2 * start_h + 3) * n + start_h + 2 == common_degree
                    for start_h, n in pairs), "Invalid common-degree construction")
        divisor_pairs, exact_count, _ = family_divisor_parameters(common_degree)
        require(set(pairs) == divisor_pairs and exact_count == rvalue - 1,
                "Power-of-five family is not the exact orbit list")
        same_minimum_families.append({"r": rvalue, "minimum_degree": common_degree,
                                     "exact_displayed_family_orbits": exact_count,
                                     "distinct_parameters_h_n": [list(pair) for pair in pairs]})
    return {
        "quotient_equalities_checked_by_two_reductions": checks - before_checks,
        "degree_formula_cases": degree_rows,
        "geometric_divisor_controls": divisor_rows,
        "moment_and_original_contraction_cases": moment_rows,
        "equal_minimum_inequivalent_pair": {
            "parameters_h_n": [[1, 4], [3, 2]], "minimum_degrees": same_degree,
            "multiplicities": [signatures[(1, 4, Q(1))], signatures[(3, 2, Q(1))]],
        },
        "common_minimum_integer_parameter_controls": same_minimum_families,
        "exact_orbit_count_controls": counts,
        "negative_controls_rejected": negative_controls[before_negatives:],
        "scope": "Finite exact factorization/multiplicity controls. Smoothness, geometric component completeness and all-parameter orbit classification are written proofs.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=HERE.parents[2])
    parser.add_argument("--output", type=Path, default=HERE / RECEIPT)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    sources = check_sources(repo)
    results = run_controls()
    higher_results = run_higher_start_controls()
    receipt = {
        "schema": 1, "status": "PASS", "python": platform.python_version(),
        "command": "python3 research/foundations/extended-geometry/" + CHECKER + " --repo-root .",
        "arithmetic": "Exact rational arithmetic; two independent quotient reductions and literal original-variable contractions.",
        "source_ledger_sha256": sha(HERE / LEDGER), "sources": sources,
        "checks": results, "higher_odd_start_checks": higher_results, "assurance": ASSURANCE,
        "scope": "Finite identity and example controls, not universal proofs. Orbit degree and displayed-family geometric multiplicity arguments are self-contained; all-order Image interpretation and coarse boundary have the separate pinned inputs.",
    }
    args.output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print("PASS: %d degree cases; %d moment/lift cases; %d radial cases; %d corrupted controls rejected."
          % (len(results["degree_formula_cases"]), len(results["moment_and_original_contraction_cases"]),
             results["radial_matching_monomial_cases"], len(results["negative_controls_rejected"])))
    print("PASS higher odd starts: %d degree cases; %d divisor profiles; %d moment/lift cases; %d corrupted controls rejected."
          % (len(higher_results["degree_formula_cases"]), len(higher_results["geometric_divisor_controls"]),
             len(higher_results["moment_and_original_contraction_cases"]), len(higher_results["negative_controls_rejected"])))
    print(args.output)


if __name__ == "__main__":
    main()
