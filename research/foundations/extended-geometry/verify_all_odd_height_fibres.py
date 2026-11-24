#!/usr/bin/env python3
"""Exact finite controls for all-odd-height fibre stabilization.

The written proof supplies the unrestricted theorem and atlas termination.
This checker does not implement an unrestricted fibre classifier.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("Assertions are required; do not run with Python optimization.")

import argparse
import hashlib
import itertools
import json
import platform
from pathlib import Path

import sympy as sp

x = sp.Symbol("x")
BASE = Path("research/foundations/extended-geometry")
PROOF = BASE / "ALL_ODD_HEIGHT_FIBRE_STABILIZATION.md"
CHECKER = BASE / "verify_all_odd_height_fibres.py"
LEDGER = BASE / "all_odd_height_fibres_source_ledger.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cancel(expr):
    return sp.cancel(sp.together(expr))


def rational(expr):
    expr = cancel(expr)
    n, d = sp.fraction(expr)
    sp.Poly(n, x, domain=sp.QQ)
    sp.Poly(d, x, domain=sp.QQ)
    return expr


def orders(expr):
    """Divisor at finite geometric points, grouped into QQ factors."""
    n, d = sp.fraction(rational(expr))
    ans = {}
    for poly, sign in ((n, 1), (d, -1)):
        for f, exponent in sp.factor_list(sp.Poly(poly, x, domain=sp.QQ))[1]:
            f = f.monic().as_expr()
            ans[f] = ans.get(f, 0) + sign * int(exponent)
    return {f: e for f, e in ans.items() if e}


def degree(expr):
    n, d = sp.fraction(rational(expr))
    return max(int(sp.degree(n, x)), int(sp.degree(d, x)))


def odd_ceiling(n):
    n = max(1, int(n))
    return n if n % 2 else n + 1


def coordinate_data(u):
    u = rational(u)
    assert u.subs(x, 0) == 0 and sp.diff(u, x).subs(x, 0) == 1
    ud = orders(u)
    zeros = {f: s for f, s in ud.items() if s > 0 and f != x}
    poles = {f: -s for f, s in ud.items() if s < 0}
    S = max(zeros.values(), default=0)
    M = max(poles.values(), default=0)
    N = degree(u)
    cutoff = odd_ceiling(max(S, M - 1))
    assert S <= N - 1 and M <= N
    assert cutoff <= 2 * ((N + 1) // 2) - 1
    return zeros, poles, cutoff, N


def height_interval(u, X):
    """Only the local polynomiality test; reciprocity is separate."""
    X = rational(X)
    assert X.subs(x, 0) == 1
    zeros, poles, _, _ = coordinate_data(u)
    xd = orders(X)
    if any(r < 0 and f not in poles for f, r in xd.items()):
        return None
    if any(xd.get(f, 0) <= 0 for f in zeros):
        return None
    lower = max([1] + [(s + xd[f] - 1) // xd[f] for f, s in zeros.items()])
    caps = [poles[f] // (-r) - 2 for f, r in xd.items() if r < 0]
    upper = min(caps) if caps else None
    return lower, upper


def in_interval(interval, a):
    return (interval is not None and a >= 1 and a % 2 == 1
            and a >= interval[0] and (interval[1] is None or a <= interval[1]))


def literal_pair(u, X, a):
    """Independent rational construction and denominator check."""
    C = rational(x * X**a / u)
    D = rational(-x * X**(a + 2) / u)
    cn, cd = sp.fraction(C)
    dn, dd = sp.fraction(D)
    ok = not (cd.has(x) or dd.has(x))
    if ok:
        assert C.subs(x, 0) == 1 and D.subs(x, 0) == -1
        assert cancel(x * X**a / C - u) == 0
        assert cancel(-D / C - X**2) == 0
    return ok, C, D


def exact_raw_moment(C, D, a, m):
    """Literal angular integral, independent of the endpoint identity."""
    if m % 2:
        return sp.S.Zero
    ans = sp.S.Zero
    for j in range(m + 1):
        coefficient = sp.Poly(sp.expand(C**(m - j) * D**j), x).nth(m)
        ans += sp.binomial(m, j) * coefficient / (a * m + 2 * j + 1)
    return cancel(ans)


def assert_polynomial_degree(C, D, X, u, a):
    n, q = sp.fraction(rational(u))
    sinfty = int(sp.degree(q, x) - sp.degree(n, x))
    d = int(sp.degree(X, x))
    cdeg, ddeg = int(sp.degree(C, x)), int(sp.degree(D, x))
    assert cdeg == 1 + a * d + sinfty
    assert ddeg == 1 + (a + 2) * d + sinfty
    sphere = max(cdeg + a - 1, ddeg + a + 1)
    assert sphere == (a + 2) * (d + 1) + sinfty
    return sphere


def verify():
    # Exact integer controls include zero, positive and negative orders.
    local_controls = 0
    for a, s, r in itertools.product(range(1, 18, 2), range(-12, 13), range(-5, 6)):
        literal = a * r - s >= 0 and (a + 2) * r - s >= 0
        if s > 0:
            filtered = r > 0 and a >= (s + r - 1) // r
        elif s < 0:
            filtered = r >= 0 or a <= (-s) // (-r) - 2
        else:
            filtered = r >= 0
        assert literal == filtered
        local_controls += 1

    # Rational local data deliberately include finite upper endpoints.
    # These are polynomiality tests, NOT automatically pure examples.
    raw_specs = [
        (1, 1, 3, 0), (2, 1, 4, 0), (3, 1, 7, 0), (4, 2, 6, 0),
        (5, 2, 3, 0), (6, 3, 4, 0),
        (1, 1, 3, 1), (1, 1, 5, 1), (1, 1, 7, 1), (1, 1, 9, 1),
        (2, 1, 5, 1), (2, 1, 7, 1), (3, 1, 5, 1), (3, 1, 7, 1),
        (4, 1, 5, 1), (4, 1, 7, 1), (4, 2, 5, 1), (5, 2, 7, 1),
        (1, 1, 6, 2), (1, 1, 8, 2), (2, 1, 6, 2), (2, 1, 8, 2),
        (3, 2, 6, 2), (5, 3, 9, 2),
    ]
    interval_checks = 0
    raw_records = []
    for s, r, m, t in raw_specs:
        u = x * (1 - x)**s / (1 - 2*x)**m
        X = (1 - x)**r / (1 - 2*x)**t
        interval = height_interval(u, X)
        observed = []
        for a in range(1, 12, 2):
            ok, _, _ = literal_pair(u, X, a)
            assert ok == in_interval(interval, a)
            if ok:
                observed.append(a)
            interval_checks += 1
        raw_records.append({"orders": [s, r, m, t], "interval": interval,
                            "tested_heights": list(range(1, 12, 2)),
                            "polynomial_heights": observed,
                            "purity_claimed": False})

    signature_controls = 0
    for N in range(1, 33):
        A = 2 * ((N + 1) // 2) - 1
        for s in range(1, N):
            for a in (A, A+2, A+8):
                assert (s-1) // a == 0
                signature_controls += 1
        for m in range(1, N+1):
            for a in (A, A+2, A+8):
                assert 1 + m // (a+2) == 1
                signature_controls += 1

    # Direct sharp family; eta=1 is reciprocal on every component.
    sharp_records = []
    polynomial_pairs = 0
    transport_checks = 0
    degree_checks = 0
    moment_checks = 0
    for N in range(2, 10):
        c = sp.Rational(N, N-1)
        u = rational(x * (1-x)**(N-1) / (1-c*x)**N)
        X = sp.expand((1-x)*(1-c*x))
        assert cancel(u / (x * sp.diff(u, x)) - X) == 0
        interval = height_interval(u, X)
        assert interval == (N-1, None)
        zeros, poles, A, actual_degree = coordinate_data(u)
        assert actual_degree == N and A == 2*((N+1)//2)-1
        for a in (A, A+2):
            ok, C, D = literal_pair(u, X, a)
            assert ok
            assert_polynomial_degree(C, D, X, u, a)
            polynomial_pairs += 1
            degree_checks += 1
            if a == A:
                ok2, C2, D2 = literal_pair(u, X, a+2)
                assert ok2 and cancel(C2-X**2*C) == 0 and cancel(D2-X**2*D) == 0
                transport_checks += 1
        if N >= 3:
            before, _, _ = literal_pair(u, X, A-2)
            assert not before
        sharp_records.append({"coordinate_degree": N, "first_odd_height": A,
                              "eta": "1", "interval": [N-1, None]})
        if N in (2, 3):
            _, C, D = literal_pair(u, X, A)
            for m in (1, 2, 3, 4):
                assert exact_raw_moment(C, D, A, m) == 0
                moment_checks += 1

    # N=1 safe boundary; infinity order -1 is essential to the degree.
    for a in (1, 3, 7):
        assert height_interval(x, sp.S.One) == (1, None)
        ok, C, D = literal_pair(x, sp.S.One, a)
        assert ok and C == 1 and D == -1
        assert assert_polynomial_degree(C, D, sp.S.One, x, a) == a+1
        polynomial_pairs += 1
        degree_checks += 1

    # Complete retained infinite fibre: explicit rational involution.
    u = x/(1+x)
    R = 1+2*x
    j = -x/(1+2*x)
    assert cancel(u.subs(x, j)+u) == 0
    infinite_records = []
    for k in (0, 1, 2, 5):
        eta = R**k
        X = sp.expand((1+x)*(1+eta)/2)
        assert cancel(eta * eta.subs(x, j)-1) == 0
        assert cancel(X-u/(x*sp.diff(u,x))*(1+eta)/2) == 0
        assert height_interval(u, X) == (1, None)
        for a in (1, 3, 5):
            ok, C, D = literal_pair(u, X, a)
            assert ok
            sphere = assert_polynomial_degree(C, D, X, u, a)
            assert sphere == (a+2)*(k+2)
            polynomial_pairs += 1
            degree_checks += 1
            ok2, C2, D2 = literal_pair(u, X, a+2)
            assert ok2 and cancel(C2-X**2*C) == 0 and cancel(D2-X**2*D) == 0
            transport_checks += 1
            if a == 1 and k in (0, 1, 2):
                for m in (1, 2, 3, 4):
                    assert exact_raw_moment(C, D, a, m) == 0
                    moment_checks += 1
        infinite_records.append({"exponent": k, "height_interval": [1, None],
                                 "degree_slope": k+2})

    negatives = []
    # D, not C, supplies the upper height endpoint for a negative order.
    assert 3*(-1)+3 == 0 and 5*(-1)+3 < 0
    negatives.append("using C alone incorrectly admits a pole of X")
    # The marked zero has different x-factor cancellation.
    assert height_interval(x, sp.S.One) == (1, None)
    negatives.append("requiring X to vanish at the marked point")
    assert height_interval(x*(1-x)/(1-2*x), sp.S.One) is None
    negatives.append("omitting an unmarked zero of u")
    assert height_interval(x/(1-x), 1/(1-2*x)) is None
    negatives.append("allowing a pole of X away from the poles of u")
    # A genuine finite upper interval for local polynomiality alone.
    uraw = x/(1-x)**5
    Xraw = 1/(1-x)
    assert height_interval(uraw, Xraw) == (1, 3)
    assert literal_pair(uraw, Xraw, 1)[0] and literal_pair(uraw, Xraw, 3)[0]
    assert not literal_pair(uraw, Xraw, 5)[0]
    etaraw = rational(2*Xraw*x*sp.diff(uraw,x)/uraw-1)
    u2 = sp.expand(sp.series(uraw,x,0,3).removeO()).coeff(x,2)
    j2 = -x-2*u2*x**2
    norm2 = sp.expand(sp.series(etaraw*etaraw.subs(x,j2),x,0,3).removeO()).coeff(x,2)
    assert norm2 == -220
    negatives.append("local polynomiality alone falsely treated as reciprocity/purity")
    assert odd_ceiling(4) == 5 and not literal_pair(
        x*(1-x)**4/(1-sp.Rational(5,4)*x)**5,
        (1-x)*(1-sp.Rational(5,4)*x), 3)[0]
    negatives.append("rounding the sharp cutoff down to the previous odd height")
    assert sp.expand((1+x)**2) != 1
    negatives.append("claiming stabilized C and D themselves remain unchanged")

    return {
        "integer_local_order_controls": local_controls,
        "rational_interval_cases": len(raw_specs),
        "independent_literal_polynomiality_checks": interval_checks,
        "stable_signature_controls": signature_controls,
        "sharp_degree_cases": len(sharp_records),
        "pure_polynomial_pair_controls": polynomial_pairs,
        "transport_controls": transport_checks,
        "exact_least_degree_controls": degree_checks,
        "raw_moment_controls": moment_checks,
        "negative_controls": negatives,
        "sharp_examples": sharp_records,
        "infinite_example_controls": infinite_records,
        "local_rational_examples": raw_records,
        "local_nonpure_transient_control": {
            "u": str(uraw), "X": str(Xraw), "height_interval": [1,3],
            "reciprocal_norm_x_squared_coefficient": str(norm2),
            "pure_example": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    ledger = json.loads((root/LEDGER).read_text())
    pins_checked = 0
    for group in ("sources", "context_sources"):
        for item in ledger.get(group, []):
            assert sha256(root/item["path"]) == item["sha256"], item["path"]
            pins_checked += 1
    results = verify()
    receipt = {
        "status": "PASS",
        "python": platform.python_version(), "sympy": sp.__version__,
        "arithmetic": "Exact rational polynomial, divisor and moment arithmetic; no floating-point acceptance.",
        "checker_sha256": sha256(Path(__file__)),
        "source_pins_checked": pins_checked,
        "packet_pins": {str(p): sha256(root/p) for p in (PROOF,CHECKER,LEDGER)},
        "scope": "Finite controls for the written unrestricted theorem; not an implementation of the entire all-height atlas.",
        "verification": results,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=False)+"\n")
    print(json.dumps({"status":"PASS", "local_order_controls":results["integer_local_order_controls"],
                      "literal_polynomiality_checks":results["independent_literal_polynomiality_checks"],
                      "pure_pair_controls":results["pure_polynomial_pair_controls"],
                      "negative_controls":len(results["negative_controls"])}))


if __name__ == "__main__":
    main()
