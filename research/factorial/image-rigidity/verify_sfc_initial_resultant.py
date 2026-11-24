#!/usr/bin/env python3
"""Exact checks for initial resultants and two complete ternary window identities.

Python standard library only.  The all-support theorem is a hand induction;
the interpolation checks below certify polynomial identities, not sampled
parameter claims.  Run from any directory.  No floating-point arithmetic.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import factorial, gcd, lcm
from pathlib import Path
import sys

if not __debug__:
    raise SystemExit("Refusing optimized Python: verification must run without -O.")

HERE = Path(__file__).resolve().parent
NOTE = "SFC_INITIAL_RESULTANT_AND_SHIFTED_WINDOWS.md"
EXPECTED = {
    NOTE: "9df960f24475ad06f6447739dc30c8d25f9a00b629f10b9f32c88a00b2631751",
    "sfc_ternary_orders_2_4_witness.json": "d3c25d8280d238533abbba7eb417445e27d92ac18a458742137f2e77e1b93639",
    "sfc_ternary_orders_3_5_witness.json": "4226c5203c2ab09835cc511dca8d40d0f5e95a2c839394a3b011efbb0eb91e46",
}


def need(ok, message):
    if not ok:
        raise ArithmeticError(message)


def pc(dim, c=1):
    return {(0,) * dim: Q(c)} if c else {}


def pv(dim, i):
    e = [0] * dim
    e[i] = 1
    return {tuple(e): Q(1)}


def pa(*polys):
    ans = {}
    for poly in polys:
        for e, c in poly.items():
            ans[e] = ans.get(e, Q(0)) + c
    return {e: c for e, c in ans.items() if c}


def ps(poly, c):
    return {e: Q(c) * v for e, v in poly.items() if c * v}


def pm(a, b):
    ans = {}
    for e, c in a.items():
        for f, d in b.items():
            g = tuple(x + y for x, y in zip(e, f))
            ans[g] = ans.get(g, Q(0)) + c * d
    return {e: c for e, c in ans.items() if c}


def pp(poly, n, dim):
    out = pc(dim)
    for _ in range(n):
        out = pm(out, poly)
    return out


def pr(poly, n, dim):
    """Rising factorial of a polynomial."""
    out = pc(dim)
    for j in range(n):
        out = pm(out, pa(poly, pc(dim, j)))
    return out


def monomial_shift(poly, index, power):
    out = {}
    for e, c in poly.items():
        f = list(e)
        f[index] += power
        out[tuple(f)] = c
    return out


def degree(poly, index):
    return max((e[index] for e in poly), default=-1)


def evaluate(poly, values):
    return sum(c * product(Q(x) ** k for x, k in zip(values, e))
               for e, c in poly.items())


def product(values):
    out = 1
    for value in values:
        out *= value
    return out


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for first in range(total, -1, -1):
            for tail in compositions(total - first, length - 1):
                yield (first,) + tail


def read_terms(rows, dim):
    ans = {}
    for row in rows:
        need(len(row) == dim + 1, "wrong witness row length")
        e, c = tuple(row[:-1]), Q(row[-1])
        need(e not in ans and c and all(isinstance(x, int) and x >= 0 for x in e),
             "noncanonical witness term")
        ans[e] = c
    return ans


def bareiss(matrix):
    a = [list(row) for row in matrix]
    n = len(a)
    if not n:
        return 1
    need(all(len(row) == n for row in a), "nonsquare determinant")
    sign, old = 1, 1
    for k in range(n - 1):
        if not a[k][k]:
            swap = next((j for j in range(k + 1, n) if a[j][k]), None)
            if swap is None:
                return 0
            a[k], a[swap] = a[swap], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                value = pivot * a[i][j] - a[i][k] * a[k][j]
                need(value % old == 0, "Bareiss division was not exact")
                a[i][j] = value // old
            a[i][k] = 0
        old = pivot
    return sign * a[-1][-1]


def determinant(matrix):
    integers, factor = [], 1
    for row in matrix:
        den = lcm(*(Q(v).denominator for v in row)) if row else 1
        integers.append([int(Q(v) * den) for v in row])
        factor *= den
    return Q(bareiss(integers), factor)


def sylvester(f, g):
    """Standard resultant; inputs are ascending, with fixed degrees retained."""
    m, n = len(f) - 1, len(g) - 1
    a = [[0] * (m + n) for _ in range(m + n)]
    for j in range(n):
        a[j][j:j + m + 1] = list(reversed(f))
    for j in range(m):
        a[n + j][j:j + n + 1] = list(reversed(g))
    return determinant(a)


def at_parameter(poly, values, final_degree):
    out = [Q(0)] * (final_degree + 1)
    for e, c in poly.items():
        out[e[-1]] += c * product(Q(x) ** k for x, k in zip(values, e[:-1]))
    return out


def ternary_components(m, fixed_a=None):
    """Coefficient polynomials in w for the actual formula (16)."""
    dim = 3 if fixed_a is None else 2
    av = pv(dim, 0) if fixed_a is None else pc(dim, fixed_a)
    bv, vi = pv(dim, dim - 2), dim - 1
    alpha, beta = pa(ps(av, m), pc(dim)), pa(ps(bv, m), pc(dim))
    out = [{} for _ in range(m // 2 + 1)]
    for k in range(m + 1):
        for j in range(k // 2 + 1):
            term = pm(pr(alpha, k - j, dim), pr(beta, m - k, dim))
            term = monomial_shift(term, vi, k - 2 * j)
            term = ps(term, Q((-1) ** j,
                             factorial(k - 2 * j) * factorial(j) * factorial(m - k)))
            out[j] = pa(out[j], term)
    return out


def resultant_factor_value(data, values, dim):
    out = Q(data["resultant_scale"])
    for factor in data["factors"]:
        out *= evaluate(read_terms(factor["terms"], dim), values) ** factor["power"]
    return out


def factor_degree_bound(data, dim):
    return tuple(sum(f["power"] * degree(read_terms(f["terms"], dim), i)
                     for f in data["factors"]) for i in range(dim))


def verify_ternary_n2(data):
    dim = 3
    av, bv, vv = pv(dim, 0), pv(dim, 1), pv(dim, 2)
    aa, bb = pa(ps(av, 2), pc(dim)), pa(ps(bv, 2), pc(dim))
    numerator = pa(pm(aa, pm(pa(av, pc(dim)), pp(vv, 2, dim))),
                   pm(aa, pm(bb, vv)), pm(bb, pa(bv, pc(dim))))
    c3, c4 = ternary_components(3), ternary_components(4)
    f = ps(pa(pm(aa, c3[0]), pm(numerator, c3[1])), 6)
    g = ps(pa(pm(pp(aa, 2, dim), c4[0]), pm(pm(aa, numerator), c4[1]),
              pm(pp(numerator, 2, dim), c4[2])), 24)
    need(f == read_terms(data["F_terms"], 3), "n2 F not derived from actual moments")
    need(g == read_terms(data["G_terms"], 3), "n2 G not derived from actual moments")
    need(all(c.denominator == 1 for c in list(f.values()) + list(g.values())),
         "nonintegral elimination coefficient")
    need(tuple(degree(f, i) for i in range(3)) == (4, 3, 3), "n2 F degree")
    need(tuple(degree(g, i) for i in range(3)) == (6, 4, 4), "n2 G degree")
    bounds = (4 * degree(f, 0) + 3 * degree(g, 0),
              4 * degree(f, 1) + 3 * degree(g, 1))
    need(bounds == (34, 24), "n2 determinant degree bound")
    need(all(x <= y for x, y in zip(factor_degree_bound(data, 2), bounds)),
         "n2 RHS degree too large")
    q70 = read_terms(data["factors"][-1]["terms"], 2)
    need(len(q70) == 70 and q70[(0, 0)] == 432 and all(c > 0 for c in q70.values()),
         "n2 positive-factor certificate failed")
    need(max(map(sum, q70)) == 12, "n2 Q total degree")
    count = 0
    for a in range(bounds[0] + 1):
        for b in range(bounds[1] + 1):
            lhs = sylvester(at_parameter(f, (a, b), 3), at_parameter(g, (a, b), 4))
            need(lhs == resultant_factor_value(data, (a, b), 2),
                 f"n2 exact polynomial identity failed at {a},{b}")
            count += 1
    return {"complete_interpolation_grid": [35, 25], "exact_determinants": count,
            "degree_bounds": list(bounds), "positive_Q_terms": len(q70)}


def verify_ternary_n3(data):
    dim = 2
    c3, c4, c5 = (ternary_components(m, 1) for m in (3, 4, 5))
    numerator, denominator = ps(c3[0], 2), ps(c3[1], -2)
    f = ps(pa(*(pm(c4[j], pm(pp(numerator, j, dim),
                            pp(denominator, 2 - j, dim))) for j in range(3))), 24)
    g = ps(pa(*(pm(c5[j], pm(pp(numerator, j, dim),
                            pp(denominator, 2 - j, dim))) for j in range(3))), 120)
    for key, poly in (("N_terms", numerator), ("J_terms", denominator),
                      ("F_terms", f), ("G_terms", g)):
        need(poly == read_terms(data[key], 2), f"n3 {key} not actual moment formula")
    need(tuple(degree(f, i) for i in range(2)) == (6, 6), "n3 F degree")
    need(tuple(degree(g, i) for i in range(2)) == (7, 7), "n3 G degree")
    bound = 7 * degree(f, 0) + 6 * degree(g, 0)
    need(bound == 84 and factor_degree_bound(data, 1)[0] <= bound, "n3 degree bound")
    q22 = read_terms(data["factors"][-1]["terms"], 1)
    need(len(q22) == 23 and degree(q22, 0) == 22 and all(c > 0 for c in q22.values()),
         "n3 positive-factor certificate failed")
    for b in range(bound + 1):
        lhs = sylvester(at_parameter(f, (b,), 6), at_parameter(g, (b,), 7))
        need(lhs == resultant_factor_value(data, (b,), 1),
             f"n3 exact polynomial identity failed at {b}")
    # Degree-three identity after substituting the sole denominator root.
    for b in range(4):
        need(evaluate(numerator, (b, -Q(3 * b + 1, 5))) ==
             Q((b + 2) * (3 * b + 1) * (3 * b + 11), 25),
             "n3 excluded-denominator identity")
    need(degree(numerator, 0) <= 3 and max(map(sum, numerator)) == 3,
         "n3 excluded-denominator interpolation degree")
    return {"complete_interpolation_grid": [85], "exact_determinants": 85,
            "degree_bound": bound, "positive_Q_terms": len(q22)}


def numeric_moment_polynomial(r, m):
    out = {}
    for e in compositions(m, len(r)):
        c = Q(1)
        for ri, ki in zip(r, e):
            for j in range(ki):
                c *= Q(m * ri + 1 + j, j + 1)
        if c:
            out[e] = c
    return out


def macaulay_resultant(polys, degrees):
    """Exact Macaulay numerator/extraneous-minor formula, fixed coordinate order."""
    s = len(degrees)
    level = sum(degrees) - s + 1
    basis = list(compositions(level, s))
    indices = {e: j for j, e in enumerate(basis)}
    matrix, extra = [], []
    for row, e in enumerate(basis):
        eligible = [i for i, d in enumerate(degrees) if e[i] >= d]
        need(eligible, "Macaulay degree did not cover a monomial")
        i = eligible[0]
        shift = list(e)
        shift[i] -= degrees[i]
        vector = [Q(0)] * len(basis)
        for alpha, c in polys[i].items():
            beta = tuple(x + y for x, y in zip(alpha, shift))
            vector[indices[beta]] += c
        matrix.append(vector)
        if len(eligible) > 1:
            extra.append(row)
    den = determinant([[matrix[i][j] for j in extra] for i in extra])
    need(den, "finite Macaulay control hit a zero extraneous minor")
    return determinant(matrix) / den, len(basis), len(extra)


def subset_product(r):
    s, out = len(r), Q(1)
    for k in range(1, s + 1):
        for indices in combinations(range(s), k):
            out *= (1 + sum(r[i] for i in indices)) ** (factorial(k - 1) * factorial(s - k))
    return out


def series_product(a, b, level, dim):
    out = [{} for _ in range(level + 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j <= level:
                out[i + j] = pa(out[i + j], pm(x, y))
    return out


def symbolic_moment(s, m):
    dim, ans = 2 * s, {}
    for e in compositions(m, s):
        term = pc(dim)
        for i, ki in enumerate(e):
            alpha = pa(ps(pv(dim, i), m), pc(dim))
            factor = ps(pr(alpha, ki, dim), Q(1, factorial(ki)))
            term = pm(term, monomial_shift(factor, s + i, ki))
        ans = pa(ans, term)
    return ans


def verify_triangular(s):
    dim, level = 2 * s, s
    a = [pc(dim)] + [{} for _ in range(level)]
    for i in range(s):
        coeffs = []
        falling = pc(dim)
        for k in range(level + 1):
            if k:
                falling = pm(falling, pa(pv(dim, i), pc(dim, -(k - 1))))
            coeffs.append(monomial_shift(ps(falling, Q((-1) ** k, factorial(k))), s + i, k))
        a = series_product(a, coeffs, level, dim)
    a = [{}] + a[:-1]
    d = [pc(dim)] + [{} for _ in range(level)]
    for k in range(1, s + 1):
        for indices in combinations(range(s), k):
            term = pa(pc(dim), *(pv(dim, i) for i in indices))
            for i in indices:
                term = monomial_shift(term, s + i, 1)
            d[k] = pa(d[k], ps(term, (-1) ** k))
    power = [pc(dim)] + [{} for _ in range(level)]
    total = [pc(dim)] + [{} for _ in range(level)]
    for m in range(1, level + 1):
        power = series_product(power, a, level, dim)
        gm = symbolic_moment(s, m)
        total = [pa(old, pm(gm, new)) for old, new in zip(total, power)]
    result = series_product(d, total, level, dim)
    need(result == [pc(dim)] + [{} for _ in range(level)],
         f"formal triangular generating identity failed for support {s}")
    return {"support": s, "formal_order": level, "parameter_variables": dim}


def verify_initial():
    formal = [verify_triangular(s) for s in range(1, 5)]
    recursion = []
    for s in range(1, 10):
        sign_exponent = factorial(s - 1) * (s + 1 + s * (s - 1) // 2)
        need(sign_exponent % 2 == 0, "resultant sign")
        for k in range(1, s):
            need((s - k) * factorial(k - 1) * factorial(s - k - 1) ==
                 factorial(k - 1) * factorial(s - k), "subset exponent recursion")
        recursion.append(s)
    determinant_controls = []
    for r in ((0,), (1,), (0, 0), (1, 2), (Q(-1, 3), Q(2, 5)),
              (0, 0, 0), (1, 1, 2), (1, 2, 3), (0, 0, 0, 0), (1, 1, 1, 1)):
        polys = [numeric_moment_polynomial(r, m) for m in range(1, len(r) + 1)]
        value, rows, extra = macaulay_resultant(polys, list(range(1, len(r) + 1)))
        need(value == subset_product(r), f"Macaulay initial resultant failed at {r}")
        determinant_controls.append({"r": list(map(str, r)), "resultant": str(value),
                                     "matrix_size": rows, "extraneous_minor_size": extra})
    return {"formal_triangular_controls": formal, "recursion_controls_supports": recursion,
            "actual_Macaulay_controls": determinant_controls}


def polynomial_determinant(matrix, dim):
    n, ans = len(matrix), {}
    for perm in permutations(range(n)):
        sign = (-1) ** sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
        term = pc(dim, sign)
        for i, j in enumerate(perm):
            term = pm(term, matrix[i][j])
        ans = pa(ans, term)
    return ans


def verify_binary_shift():
    dim, rv, sv = 2, pv(2, 0), pv(2, 1)
    def coeffs(m):
        return [ps(pm(pr(pa(ps(rv, m), pc(dim)), k, dim),
                      pr(pa(ps(sv, m), pc(dim)), m - k, dim)),
                   Q(1, factorial(k) * factorial(m - k))) for k in range(m + 1)]
    f, g = coeffs(2), coeffs(3)
    matrix = [[{} for _ in range(5)] for _ in range(5)]
    for j in range(3):
        matrix[j][j:j + 3] = list(reversed(f))
    for j in range(2):
        matrix[3 + j][j:j + 4] = list(reversed(g))
    lhs = polynomial_determinant(matrix, dim)
    q = {(3, 2): Q(81), (3, 1): Q(74), (3, 0): Q(17),
         (2, 3): Q(81), (2, 2): Q(257), (2, 1): Q(181), (2, 0): Q(37),
         (1, 3): Q(74), (1, 2): Q(181), (1, 1): Q(115), (1, 0): Q(22),
         (0, 3): Q(17), (0, 2): Q(37), (0, 1): Q(22), (0, 0): Q(4)}
    rhs = ps(pm(pm(pa(rv, pc(dim)), pa(sv, pc(dim))),
                pm(pp(pa(rv, sv, pc(dim)), 2, dim), q)), Q(1, 4))
    need(lhs == rhs, "binary nonlinear shifted resultant identity")
    need([evaluate(x, (Q(-1, 3), Q(-1, 3))) for x in f] ==
         [Q(2, 9), Q(1, 9), Q(2, 9)], "shifted countercontrol g2")
    need(all(evaluate(x, (Q(-1, 3), Q(-1, 3))) == 0 for x in g),
         "shifted countercontrol g3")
    need(subset_product((Q(-1, 3), Q(-1, 3))) == Q(4, 27),
         "countercontrol unexpectedly has initial resonance")
    return {"symbolic_resultant_verified": True, "nonlinear_factor_terms": len(q),
            "negative_parameter_countercontrol": ["-1/3", "-1/3"],
            "initial_resultant_at_countercontrol": "4/27"}


def verify_boundaries():
    # These checks use polynomial identities, including parameter a.
    dim, av, vv = 2, pv(2, 0), pv(2, 1)
    w = pm(pa(av, pc(dim)), pp(vv, 2, dim))
    def two_moment(m):
        ans = {}
        for j in range(m // 2 + 1):
            factor = pr(pa(ps(av, m), pc(dim)), m - j, dim)
            factor = ps(factor, Q((-1) ** j, factorial(m - 2 * j) * factorial(j)))
            ans = pa(ans, pm(monomial_shift(factor, 1, m - 2 * j), pp(w, j, dim)))
        return ans
    need(not two_moment(2), "n2 boundary g2 substitution")
    expected = ps(pm(pm(pa(ps(av, 3), pc(dim)), pa(ps(av, 3), pc(dim, 2))),
                     pa(av, pc(dim))), Q(-1, 2))
    expected = monomial_shift(expected, 1, 3)
    need(two_moment(3) == expected, "n2 boundary g3 identity")
    # Exact a=1 forms in v,w.
    for m, wanted in ((3, {(3, 0): Q(20), (1, 1): Q(-20)}),
                      (4, {(4, 0): Q(70), (2, 1): Q(-105), (0, 2): Q(15)})):
        got = {}
        for j in range(m // 2 + 1):
            c = Q((-1) ** j, factorial(m - 2 * j) * factorial(j))
            c *= product(range(m + 1, 2 * m - j + 1))
            got[(m - 2 * j, j)] = c
        need(got == wanted, "n3 boundary formula")
    # Roots-of-unity resonance: formal coefficient for every m in a finite control.
    resonance_count = 0
    for k in range(1, 9):
        for m in range(1, 25):
            if m % k == 0:
                ell = m // k
                coeff = Q(1)
                exponent = Q(m, k) - 1
                for j in range(ell):
                    coeff *= Q(exponent - j, j + 1)
                need(coeff == 0, "minimal resonance coefficient")
            resonance_count += 1
    return {"exact_boundary_identities": 4,
            "minimal_resonance_coefficient_controls": resonance_count}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", type=Path, default=HERE / NOTE)
    parser.add_argument("--output", type=Path,
                        default=HERE / "sfc_initial_resultant_verification.json")
    args = parser.parse_args()
    pins = {}
    for name, expected in EXPECTED.items():
        path = args.note if name == NOTE else HERE / name
        need(path.is_file(), f"missing proof or witness: {path}")
        actual = sha256(path.read_bytes()).hexdigest()
        need(actual == expected, f"source pin mismatch: {name}")
        pins[name] = actual
    n2 = json.loads((HERE / "sfc_ternary_orders_2_4_witness.json").read_text())
    n3 = json.loads((HERE / "sfc_ternary_orders_3_5_witness.json").read_text())
    receipt = {"status": "PASS", "arithmetic": "exact rational/integer; Python standard library",
               "scope": "all-support initial formula by written induction; two uniform shifted "
                        "ternary polynomial identities by complete degree-bounded interpolation; "
                        "no arbitrary-start SFC closure",
               "source_sha256": pins,
               "initial_resultant": verify_initial(),
               "ternary_orders_2_4": verify_ternary_n2(n2),
               "ternary_orders_3_5": verify_ternary_n3(n3),
               "shifted_countercontrol": verify_binary_shift(),
               "boundaries_and_resonances": verify_boundaries()}
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print("PASS: initial-resultant controls; 960 complete-interpolation determinants; "
          "two uniform ternary windows; shifted countercontrol.")


if __name__ == "__main__":
    try:
        main()
    except (ArithmeticError, OSError, ValueError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        raise SystemExit(1)
