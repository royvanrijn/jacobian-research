#!/usr/bin/env python3
"""Exact support for the written classical Gamma/product-annulus proof.

Finite algebra controls do not replace its all-parameter argument.
Standard library only. No canonical files are written.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
PROOF = "GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md"
CHECKER = "verify_gamma_product_sharp_lower_annulus.py"
INPUTS = "gamma_product_sharp_lower_annulus_inputs.json"
LEDGER = "gamma_product_sharp_lower_annulus_source_ledger.json"
RECEIPT = "gamma_product_sharp_lower_annulus_verification.json"
ROLE_SPEC = {
    "written_proof": ("local", PROOF, "written_argument"),
    "exact_checker": ("local", CHECKER, "computational_support"),
    "exact_fixture_input": ("local", INPUTS, "computational_support"),
    "classical_bessel_primary": ("local", "sources/varga-generalized-bessel-1981-I-II.pdf", "classical_context_only"),
    "gamma_upper_disk_proof": ("repo", "research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex", "mathematical_input_upper_factor_only"),
    "closed_disc_grace_proof": ("repo", "research/factorial/image-rigidity/GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md", "mathematical_input_closed_disc_only"),
}
SCOPE = {
    "positive_real_shapes": True,
    "all_positive_integer_degrees_and_gaps": True,
    "factor_simplicity_proved": True,
    "product_simplicity_claimed": False,
    "fixed_lower_exponent_all_order_sfc_claimed": False,
    "criterion_necessity_claimed": False,
    "classical_factor_circle_priority_claimed": False,
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def repo_root(value):
    if value is not None:
        root = Path(value).resolve()
        require((root / "research/MATH_STATUS.json").is_file(), "invalid --repo-root")
        return root
    for candidate in (Path.cwd(), *Path.cwd().parents, HERE, *HERE.parents):
        if (candidate / "research/MATH_STATUS.json").is_file():
            return candidate
    raise ValueError("provide --repo-root when the package is outside the repository")


def validate_ledger(data, root):
    require(type(data) is dict, "ledger must be an object")
    require(type(data.get("schema")) is int and data["schema"] == 1, "ledger schema must be integer1")
    require(type(data.get("sources")) is list, "ledger sources must be a list")
    records = data["sources"]
    require(len(records) == len(ROLE_SPEC), "exactly six source roles required")
    require(all(type(row) is dict for row in records), "source records must be objects")
    roles = [row.get("role") for row in records]
    require(all(type(role) is str for role in roles), "source roles must be strings")
    require(len(set(roles)) == len(roles) and set(roles) == set(ROLE_SPEC), "source-role set mismatch")
    pins = {}
    for row in records:
        require(set(row) == {"role", "location", "path", "sha256", "relation"}, "source record fields mismatch")
        role = row["role"]
        location, path, relation = ROLE_SPEC[role]
        require(row["location"] == location and row["path"] == path and row["relation"] == relation,
                "source role/path/relation mismatch: " + role)
        expected = row["sha256"]
        require(type(expected) is str and re.fullmatch(r"[0-9a-f]{64}", expected), "invalid source SHA256")
        actual_path = (HERE if location == "local" else root) / path
        require(actual_path.is_file(), "missing source: " + role)
        require(digest(actual_path) == expected, "stale or corrupt source: " + role)
        pins[role] = expected
    primary = data.get("primary_inspection")
    require(type(primary) is dict and primary.get("inspected") is True, "actual primary inspection required")
    require(primary.get("theorem") == "Paper I, printed3, Theorem3.1; printed2 (2.1),(2.3)", "primary locator mismatch")
    require(primary.get("original_docev_paper_inspected") is False, "do not imply original Docev inspection")
    flags = data.get("assurance_flags")
    require(type(flags) is dict and set(flags) == {"formal_verification", "external_review", "independent_replay"}
            and all(v is False for v in flags.values()), "assurance flags must remain false")
    return pins


def trim(poly):
    out = [F(c) for c in poly]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def add(p, q):
    out = [F(0)] * max(len(p), len(q))
    for i, x in enumerate(p):
        out[i] += x
    for i, x in enumerate(q):
        out[i] += x
    return trim(out)


def scale(p, c):
    return trim([c * x for x in p])


def mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return trim(out)


def deriv(p):
    return trim([i * p[i] for i in range(1, len(p))] or [0])


def power(p, n):
    out = [F(1)]
    for _ in range(n):
        out = mul(out, p)
    return out


def rem(p, q):
    p, q = trim(p), trim(q)
    require(q != [0], "zero polynomial divisor")
    while p != [0] and len(p) >= len(q):
        shift, c = len(p) - len(q), p[-1] / q[-1]
        for j, v in enumerate(q):
            p[shift + j] -= c * v
        p = trim(p)
    return p


def gcd_poly(p, q):
    p, q = trim(p), trim(q)
    while q != [0]:
        p, q = q, rem(p, q)
        if q != [0]:
            q = scale(q, 1 / q[-1])
    return scale(p, 1 / p[-1]) if p != [0] else p


def rising(gamma, n):
    out = F(1)
    for j in range(n):
        out *= gamma + j
    return out


def g_vector(n, gamma):
    return [F(comb(n, j)) * rising(gamma, n - j) for j in range(n + 1)]


def t_vector(n, gamma, gap):
    return [F(comb(n, k)) * rising(gamma, gap * k) for k in range(n + 1)]


def compose(p, q):
    require(len(p) == len(q), "composition degrees differ")
    n = len(p) - 1
    return [p[k] * q[k] / comb(n, k) for k in range(n + 1)]


def sparse_add(p, q):
    r = dict(p)
    for k, c in q.items():
        r[k] = r.get(k, F(0)) + c
        if r[k] == 0:
            del r[k]
    return r


def sparse_scale(p, c):
    return {k: c * v for k, v in p.items() if c * v}


def sparse_mul(p, q):
    r = {}
    for a, c in p.items():
        for b, v in q.items():
            k = tuple(x + y for x, y in zip(a, b))
            r[k] = r.get(k, F(0)) + c * v
    return {k: v for k, v in r.items() if v}


def sparse_deriv(p, var):
    r = {}
    for k, v in p.items():
        if k[var]:
            j = list(k)
            j[var] -= 1
            r[tuple(j)] = v * k[var]
    return r


def check_ode(max_n):
    one, z, gamma = {(0, 0): F(1)}, {(1, 0): F(1)}, {(0, 1): F(1)}
    count = 0
    for n in range(1, max_n + 1):
        g = {}
        for k in range(n + 1):
            r = [F(1)]
            for j in range(k):
                r = mul(r, [F(j), F(1)])
            for e, c in enumerate(r):
                if c:
                    g[(n - k, e)] = c * comb(n, k)
        gp = sparse_deriv(g, 0)
        a = sparse_add(sparse_scale(one, 1 - n), sparse_scale(gamma, -1))
        a = sparse_add(a, sparse_scale(z, -1))
        identity = sparse_add(sparse_mul(z, sparse_deriv(gp, 0)), sparse_mul(a, gp))
        identity = sparse_add(identity, sparse_scale(g, n))
        require(not identity, "formal gamma ODE failed")
        # The actual next-to-leading coefficient is n*gamma, fixing Vieta.
        vieta = {e: c for (ze, e), c in g.items() if ze == n - 1}
        require(vieta == {1: F(n)}, "formal Vieta identity failed")
        count += 1
    # Cleared real-part identity: D-2(1-x)=x*x+y*y-1.
    x, y = {(1, 0): F(1)}, {(0, 1): F(1)}
    d = sparse_add(sparse_mul(sparse_add(one, sparse_scale(x, -1)),
                               sparse_add(one, sparse_scale(x, -1))), sparse_mul(y, y))
    lhs = sparse_add(d, sparse_scale(sparse_add(one, sparse_scale(x, -1)), -2))
    rhs = sparse_add(sparse_add(sparse_mul(x, x), sparse_mul(y, y)), sparse_scale(one, -1))
    require(lhs == rhs, "cleared complex inequality identity failed")
    return count


def check_composition(data):
    count = 0
    for row in data["composition_fixtures"]:
        n, gap, gamma = row["n"], row["d"], F(row["gamma"])
        c = None
        for j in range(gap):
            g = g_vector(n, (gamma + j) / gap)
            c = g if c is None else compose(c, g)
        actual = [c[k] * gap ** (gap * (n - k)) for k in range(n + 1)]
        expected = list(reversed(t_vector(n, gamma, gap)))
        require(actual == expected, "full multiplication/dilation vector failed")
        for k in range(n + 1):
            product = F(gap ** (gap * k))
            for j in range(gap):
                product *= rising((gamma + j) / gap, k)
            require(product == rising(gamma, gap * k), "rising multiplication failed")
        # Deliberately nonmonic, non-unit constant factors test reversal.
        p = [F((k + 2) * (n + 3), k + 1) for k in range(n + 1)]
        q = [F((k + 3) ** 2, k + 2) for k in range(n + 1)]
        require(list(reversed(compose(p, q))) == compose(list(reversed(p)), list(reversed(q))),
                "reversal/composition commutation failed")
        require(expected[0] and expected[-1] == 1, "exact degree or constant lost")
        count += 1
    return count


def check_shape_gcds(data):
    count, equality, newly_covered = 0, 0, 0
    specs = data["shape_gcd_domain"]
    for d in range(1, specs["d_max"] + 1):
        for r in range(1, specs["r_max"] + 1):
            for s in range(r + 1, specs["s_max"] + 1):
                for gamma in map(F, specs["gammas"]):
                    for offset in map(F, specs["threshold_offsets"]):
                        delta = d * (r - 1) + offset
                        p, q = t_vector(r, gamma, d), t_vector(s, gamma + delta, d)
                        require(gcd_poly(p, q) == [1], "shape threshold QQ gcd failed")
                        count += 1
                        equality += int(offset == 0)
                        newly_covered += int(delta < d * (r + s - 2))
    return count, equality, newly_covered


def check_moments(data):
    pair_count, normalized, annulus = 0, 0, 0
    specs = data["moment_domain"]
    for d in range(1, specs["d_max"] + 1):
        for n in range(2, specs["N_max"] + 1):
            a = d * (n - 2)
            for beta in map(F, specs["betas"]):
                polynomials = {j: t_vector(j, F(a * j) + beta, d) for j in range(1, n + 1)}
                for r, s in combinations(range(1, n + 1), 2):
                    require(a * (s - r) >= d * (r - 1), "pairwise maximum inequality failed")
                    require(gcd_poly(polynomials[r], polynomials[s]) == [1], "actual-shape QQ gcd failed")
                    pair_count += 1
                for j, p in polynomials.items():
                    # Direct expansion of the original phase X^a(1+t X^d),
                    # then apply L_beta(X^k)=(beta)_k. This preserves each
                    # actual moment's own shape a*j+beta.
                    direct = [F(comb(j, k)) * rising(beta, a * j + d * k) for k in range(j + 1)]
                    expected = scale(p, rising(beta, a * j))
                    require(direct == expected, "actual moment shape normalization failed")
                    normalized += 1
    for d in range(1, specs["d_max"] + 1):
        for m in range(1, specs["N_max"] + 1):
            for a in map(F, data["ordinary_annulus_a"]):
                lower = rising(a * m + 1, d) / m ** d
                upper = rising(a * m + 1 + d * (m - 1), d) / m ** d
                require(a ** d < lower <= upper <= (a + d) ** d, "ordinary endpoint bound failed")
                annulus += 1
    return pair_count, normalized, annulus


def check_cusp_cutoff(data):
    count = 0
    for d in range(3, data["cusp_cutoff_d_max"] + 1, 2):
        a = F((d - 1) ** 2, 4)
        for n in range(2, data["cusp_cutoff_N_max"] + 1):
            if a < d * (n - 2):
                require(d * d - (4 * n - 6) * d + 1 < 0, "cusp quadratic rearrangement failed")
                require(d < 4 * n - 6 and F(d, n - 1) < 4, "cusp hard complement failed")
                count += 1
    return count


def mathematical_negative_controls():
    # Wrong sign on zG' fails for an actual nonconstant Gamma polynomial.
    n, gamma = 3, F(5, 2)
    g = g_vector(n, gamma)
    correct = add(add([0] + deriv(deriv(g)),
                      mul([F(1) - gamma - n, F(-1)], deriv(g))), scale(g, n))
    wrong = add(add([0] + deriv(deriv(g)),
                    mul([F(1) - gamma - n, F(1)], deriv(g))), scale(g, n))
    require(correct == [0] and wrong != [0], "wrong ODE sign control did not fail")
    # n=2 roots are -gamma +/- i*sqrt(gamma): gamma+n-1 is not a lower bound.
    gamma = F(3, 2)
    require(gamma ** 2 + gamma < (gamma + 1) ** 2, "false larger lower circle not rejected")
    # Omitting the d^d dilation gives a genuinely different full vector.
    n, d, gamma = 3, 2, F(7, 3)
    wrong = compose(g_vector(n, gamma / d), g_vector(n, (gamma + 1) / d))
    require(wrong != list(reversed(t_vector(n, gamma, d))), "missing dilation not rejected")
    # The sufficient shift condition is not necessary: this actual pair
    # lies below threshold and still has no common root.
    r, s, d, gamma, delta = 2, 3, 2, F(1), F(0)
    require(delta < d * (r - 1) and gcd_poly(t_vector(r, gamma, d), t_vector(s, gamma + delta, d)) == [1],
            "criterion-necessity negative control failed")
    return 4


def source_negative_controls(ledger, root, inputs):
    count = 0
    for index in range(len(ROLE_SPEC)):
        bad = deepcopy(ledger)
        del bad["sources"][index]
        try:
            validate_ledger(bad, root)
        except (ValueError, OSError):
            count += 1
        else:
            raise ValueError("a mandatory source omission was accepted")
        bad = deepcopy(ledger)
        bad["sources"][index]["sha256"] = "0" * 64
        try:
            validate_ledger(bad, root)
        except (ValueError, OSError):
            count += 1
        else:
            raise ValueError("a bad mandatory source pin was accepted")
    for value in (True, "1", 2):
        bad = deepcopy(ledger)
        bad["schema"] = value
        try:
            validate_ledger(bad, root)
        except ValueError:
            count += 1
        else:
            raise ValueError("a wrong schema type/version was accepted")
    for key in ("product_simplicity_claimed", "fixed_lower_exponent_all_order_sfc_claimed", "criterion_necessity_claimed"):
        bad = deepcopy(inputs)
        bad["scope"][key] = True
        try:
            validate_inputs(bad)
        except ValueError:
            count += 1
        else:
            raise ValueError("a false scope extension was accepted")
    return count


def validate_inputs(data):
    require(type(data) is dict and type(data.get("schema")) is int and data["schema"] == 1, "input schema")
    require(type(data.get("scope")) is dict and set(data["scope"]) == set(SCOPE), "input scope fields")
    require(all(type(v) is bool for v in data["scope"].values()) and data["scope"] == SCOPE, "input scope mismatch")
    require(type(data.get("ode_max_n")) is int and data["ode_max_n"] == 16, "ODE fixture coverage changed")
    for domain, names in (("shape_gcd_domain", ("d_max", "r_max", "s_max")),
                          ("moment_domain", ("d_max", "N_max"))):
        require(type(data.get(domain)) is dict and all(type(data[domain].get(k)) is int for k in names),
                "domain bounds must be exact integers")
    require(data.get("shape_gcd_domain") == {"d_max": 4, "r_max": 4, "s_max": 6,
            "gammas": ["1/3", "1", "5/2"], "threshold_offsets": ["0", "1/7"]}, "shape domain changed")
    require(data.get("moment_domain") == {"d_max": 4, "N_max": 7, "betas": ["1/2", "1"]}, "moment domain changed")
    require(data.get("ordinary_annulus_a") == ["1/7", "1", "3", "100"], "annulus domain changed")
    require(type(data.get("cusp_cutoff_d_max")) is int and data["cusp_cutoff_d_max"] == 101, "cusp d coverage changed")
    require(type(data.get("cusp_cutoff_N_max")) is int and data["cusp_cutoff_N_max"] == 60, "cusp N coverage changed")
    expected = [{"n": n, "d": d, "gamma": gamma} for n in range(1, 7) for d in range(1, 6)
                for gamma in ("1/3", "2", "11/2")]
    require(type(data.get("composition_fixtures")) is list and all(type(row) is dict
            and set(row) == {"n", "d", "gamma"} and type(row["n"]) is int
            and type(row["d"]) is int and type(row["gamma"]) is str for row in data["composition_fixtures"]),
            "composition input types mismatch")
    require(data.get("composition_fixtures") == expected, "composition coverage changed")


def main():
    if not __debug__ or sys.flags.optimize:
        raise SystemExit("optimized Python is refused; run without -O/-OO")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root")
    parser.add_argument("--output", type=Path, default=HERE / RECEIPT)
    args = parser.parse_args()
    root = repo_root(args.repo_root)
    ledger = load_json(HERE / LEDGER)
    pins = validate_ledger(ledger, root)
    data = load_json(HERE / INPUTS)
    validate_inputs(data)
    ode = check_ode(data["ode_max_n"])
    composition = check_composition(data)
    shape = check_shape_gcds(data)
    moments = check_moments(data)
    cusp = check_cusp_cutoff(data)
    math_neg = mathematical_negative_controls()
    source_neg = source_negative_controls(ledger, root, data)
    receipt = {
        "schema": 1,
        "status": "PASS",
        "scope": SCOPE,
        "source_sha256": pins,
        "source_ledger_sha256": digest(HERE / LEDGER),
        "checks": {
            "formal_gamma_ode_and_vieta_identities": ode,
            "universal_cleared_complex_real_part_identity": 1,
            "full_product_dilation_and_reversal_vectors": composition,
            "exact_shape_threshold_QQ_gcds": shape[0],
            "shape_threshold_equality_cases": shape[1],
            "shape_cases_not_covered_by_old_threshold": shape[2],
            "actual_shape_pairwise_QQ_gcds": moments[0],
            "full_original_gamma_moment_vectors": moments[1],
            "ordinary_normalized_annulus_endpoints": moments[2],
            "cusp_hard_complement_controls": cusp,
            "rejected_mathematical_extensions": math_neg,
            "rejected_source_schema_and_scope_controls": source_neg,
        },
        "universal_proof": "Written minimum-modulus/strictness argument, closed-disc composition and reversal; finite fixtures are support only.",
        "classical_comparison": "1981 paper I (2.1),(2.3),Theorem3.1 and following Docev-circle attribution inspected; classical factor circle is not claimed new.",
        "assurance_flags": {"formal_verification": False, "external_review": False, "independent_replay": False},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "checks": receipt["checks"]}, sort_keys=True))


if __name__ == "__main__":
    main()
