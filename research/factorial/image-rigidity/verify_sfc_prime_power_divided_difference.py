#!/usr/bin/env python3
"""Exact finite controls for the uniform divided-difference SFC proof.

Only the Python standard library is used. The unrestricted theorem is
proved in the pinned note; these finite controls are not an extrapolation.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, factorial, isqrt
from pathlib import Path


if not __debug__:
    raise SystemExit("Refusing optimized Python: assertions must remain enabled.")


PROOF_DIR = Path("research/factorial/image-rigidity")
SOURCE_PINS = {
    "SFC_PRIME_POWER_DIVIDED_DIFFERENCE_CLOSURE.md":
        "bd7d60241cc455581bd90c8eb82be4ee8bcb2d6638be6bb5117ec5de6fed3dc6",
    "SFC_PRIME_POWER_SCALAR_TRANSFER.md":
        "9fa457f230397c670c7b4152c0510df93da24df03fc08112f1c5bf5c73bc0349",
    "SFC_PRIME_POWER_UPPER_RESPONSE.md":
        "c32d3c76f705aa5fbb5e046b23b1ca4c18fdc91568b5f3e5e36b3fd4774dbdf2",
}
INF = 10**30


def rational(x):
    x = F(x)
    return {"numerator": x.numerator, "denominator": x.denominator}


def vp(x, p):
    x = F(x)
    if not x:
        return INF
    a, b = abs(x.numerator), x.denominator
    out = 0
    while a % p == 0:
        a //= p
        out += 1
    while b % p == 0:
        b //= p
        out -= 1
    return out


def vfac(n, p):
    out = 0
    while n:
        n //= p
        out += n
    return out


def vbinom(n, k, p):
    assert 0 <= k <= n
    return vfac(n, p) - vfac(k, p) - vfac(n-k, p)


def rise(x, length):
    out = F(1)
    for j in range(length):
        out *= x+j
    return out


def trim(a):
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a


def pmul(a, b):
    c = [F(0)] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return trim(c)


def peval(a, x):
    out = F(0)
    for c in reversed(a):
        out = out*x+c
    return out


def root_product(roots):
    out = [F(1)]
    for x in roots:
        out = pmul(out, [-x, F(1)])
    return out


def rising_poly(h):
    out = [F(1)]
    for i in range(h):
        out = pmul(out, [F(i), F(1)])
    return out


def divided_difference(coeff, nodes):
    """Polynomial divided difference, without dividing node differences."""
    m = len(nodes)-1
    if len(coeff)-1 < m:
        return F(0)
    degree = len(coeff)-1-m
    complete = [F(1)] + [F(0)]*degree
    for node in nodes:
        for i in range(1, degree+1):
            complete[i] += node*complete[i-1]
    return sum((coeff[m+i]*complete[i] for i in range(degree+1)), F(0))


def direct_functionals(poly, q, n, gamma):
    """Use the original rising-jet triangular basis, not formula (2)."""
    c = gamma+q-1
    C = gamma+q*(n+1)
    basis = []
    falling = [F(1)]
    for k in range(q-1):
        basis.append(pmul([F(q*n), F(1)], falling))
        falling = pmul(falling, [c-k, F(-1)])
    rem = list(poly)+[F(0)]*(q-len(poly))
    rem[0] -= peval(poly, -q*n)
    d = [F(0)]*(q-1)
    for k in range(q-2, -1, -1):
        d[k] = rem[k+1]/basis[k][-1]
        for i, value in enumerate(basis[k]):
            rem[i] -= d[k]*value
    assert not any(rem)
    return [F(0)]+[
        sum((d[k]*rise(C+q-k, h-q+k)
             for k in range(max(0, q-h), q-1)), F(0))
        for h in range(1, q+1)
    ]


def bound(p, h, m=3):
    out = 0
    h //= p
    while h:
        out += max(h-m, 0)
        h //= p
    return out


def R_value(p, h):
    out = 0
    h //= p
    while h:
        out += min(h, 3)
        h //= p
    return out


def check_functionals():
    identities = 0
    inequalities = 0
    digest = sha256()
    cases = []
    for p, q in [(2, 2), (3, 3), (2, 4), (5, 5), (7, 7),
                 (2, 8), (3, 9), (11, 11), (2, 16), (5, 25)]:
        shapes = [F(1), F(7), F(1, 3 if p == 2 else 2)]
        for n, gamma in [(1, shapes[0]), (2, shapes[1]), (5, shapes[2])]:
            assert vp(gamma, p) >= 0
            c, C = gamma+q-1, gamma+q*(n+1)
            r = (c.numerator*pow(c.denominator, -1, q)) % q
            nodes = [c-j for j in range(q)]
            Qr = root_product([nodes[j] for j in range(q) if j != r])
            qvals = direct_functionals(Qr, q, n, gamma)
            D = [INF]*(q+1)
            for h in range(1, q+1):
                expected = (-1)**q*divided_difference(
                    rising_poly(h), [C, C-q, F(r-q+1)])
                assert qvals[h] == expected
                identities += 1
                digest.update(f"Q|{q}|{n}|{gamma}|{h}|{expected}\n".encode())
            for j in range(q):
                if j == r:
                    continue
                Gj = root_product([nodes[i] for i in range(q)
                                   if i != r and i != j])
                gvals = direct_functionals(Gj, q, n, gamma)
                for h in range(1, q+1):
                    expected = (-1)**(q+1)*divided_difference(
                        rising_poly(h), [C, C-q, F(r-q+1), F(j-q+1)])
                    assert gvals[h] == expected
                    assert vp(expected, p) >= bound(p, h)
                    identities += 1
                    inequalities += 1
                    D[h] = min(D[h], vp(expected, p))
                    digest.update(f"G|{q}|{n}|{gamma}|{j}|{h}|{expected}\n".encode())
            assert qvals[1] == 0
            assert qvals[2] == (-1)**q
            if q >= 3:
                assert qvals[3] == (-1)**q*(2*C-2*q+r+4)
                assert D[3] == 0
            cases.append({"p": p, "q": q, "n": n, "gamma": rational(gamma),
                          "omitted_node": r,
                          "functional_minimum_valuations":
                          [None if x == INF else x for x in D[1:]],
                          "proved_lower_bounds": [bound(p, h) for h in range(1, q+1)]})
    return {"identity_count": identities, "valuation_inequality_count": inequalities,
            "ordered_identity_stream_sha256": digest.hexdigest(), "cases": cases}


def check_blocks():
    controls = 0
    minima = []
    for p in [2, 3, 5, 7]:
        for h in range(3, 25):
            minimum = INF
            for l0 in range(h-2):
                for l1 in range(h-2-l0):
                    for l2 in range(h-2-l0-l1):
                        l3 = h-3-l0-l1-l2
                        lengths = [l0, l1, l2, l3]
                        assert sum(lengths) == h-3
                        value = sum(vfac(l, p) for l in lengths)
                        assert value >= bound(p, h)
                        P = p
                        while P <= h:
                            assert sum(l//P for l in lengths) >= max(h//P-3, 0)
                            P *= p
                        controls += 1
                        minimum = min(minimum, value)
            minima.append({"p": p, "h": h, "exact_composition_minimum": minimum,
                           "bound": bound(p, h)})
    return {"composition_count": controls, "minima": minima}


def check_power_cancellation():
    controls = 0
    rows = []
    for p, max_s in [(2, 7), (3, 5), (5, 3), (7, 2), (11, 2)]:
        for s in range(1, max_s+1):
            q = p**s
            for w in range(1, s+1):
                a = p**w
                delta = int(p == 2 and w == s)
                expected = s-vfac(a, p)-delta
                margin = expected+2*s+bound(p, a)
                assert margin >= 3*(s-w)+2-delta > 0
                assert R_value(p, a) <= 3*w-2
                for beta in sorted({1, a, max(1, a//2)}):
                    assert vbinom(a-1, beta-1, p) == 0
                    for N in [2, 3, p**(s+1), 10**20+1]:
                        C = (a+q)*N+beta-a
                        actual = s+vp(N, p)-sum(vp(C+j, p) for j in range(a))
                        assert actual == expected
                        controls += 1
                rows.append({"p": p, "s": s, "w": w, "a": a, "q": q,
                             "kappa_valuation": expected,
                             "functional_bound": bound(p, a),
                             "strict_criterion_margin": margin,
                             "coarse_margin_without_functional_bound": expected+2*s})
    assert any(r["coarse_margin_without_functional_bound"] <= 0 for r in rows)
    return {"exact_cancellation_count": controls, "families": rows}


def check_general_threshold():
    controls = 0
    rows = []
    for p, s in [(2, 5), (3, 4), (5, 3), (7, 2), (11, 2)]:
        q = p**s
        for a in range(1, q):
            w = vp(a, p)
            A = a//p**w
            K = 3*s-R_value(p, a)
            assert K >= 3
            assert vfac(a, p)-R_value(p, a) == bound(p, a)
            for N in [2, 3, p**2+1, 2*p**3+7]:
                C = (a+q)*N-a+1
                actual = s+vp(N, p)-sum(vp(C+j, p) for j in range(a))
                V = vbinom((A+p**(s-w))*N-1, A-1, p)
                assert actual == s-vfac(a, p)-V
                assert actual+2*s+bound(p, a) == K-V
                controls += 1
            if a in [1, q-1, q//p, max(1, q//2)]:
                rows.append({"p": p, "s": s, "a": a, "K": K,
                             "old_factorial_threshold": 3*s-vfac(a, p)})
    residues = []
    for p, s, a in [(5, 1, 3), (5, 2, 15), (7, 2, 21)]:
        q = p**s
        w = vp(a, p)
        A = a//p**w
        K = 3*s-R_value(p, a)
        M = p**K
        multiplier = A+p**(s-w)
        exceptional = {j*pow(multiplier, -1, M) % M for j in range(1, A)}
        # Full residue system, not a sample, for these bounded controls.
        found = set()
        for n_res in range(M):
            N = n_res if n_res >= 2 else n_res+M
            V = vbinom(multiplier*N-1, A-1, p)
            if V >= K:
                found.add(n_res)
        assert found == exceptional
        residues.append({"p": p, "s": s, "a": a, "K": K,
                         "modulus": M, "exceptional_terminal_classes": sorted(found),
                         "scope": "failure of this sufficient criterion only"})
    return {"ordinary_identity_count": controls, "threshold_examples": rows,
            "complete_residue_controls": residues}


def mod_trim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def mod_divrem(a, b, prime):
    a = list(a)
    assert any(b)
    q = [0]*max(1, len(a)-len(b)+1)
    inv = pow(b[-1], -1, prime)
    while any(a) and len(a) >= len(b):
        i = len(a)-len(b)
        c = a[-1]*inv % prime
        q[i] = c
        for j, v in enumerate(b):
            a[i+j] = (a[i+j]-c*v) % prime
        mod_trim(a)
    return mod_trim(q), mod_trim(a)


def mod_add(a, b, prime):
    return mod_trim([((a[i] if i < len(a) else 0)+
                      (b[i] if i < len(b) else 0)) % prime
                     for i in range(max(len(a), len(b)))])


def mod_mul(a, b, prime):
    c = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] = (c[i+j]+x*y) % prime
    return mod_trim(c)


def mod_neg(a, prime):
    return [(-x) % prime for x in a]


def bezout(a, b, prime):
    r0, r1 = list(a), list(b)
    s0, s1, t0, t1 = [1], [0], [0], [1]
    while any(r1):
        quotient, rem = mod_divrem(r0, r1, prime)
        r0, r1 = r1, rem
        s0, s1 = s1, mod_add(s0, mod_neg(mod_mul(quotient, s1, prime), prime), prime)
        t0, t1 = t1, mod_add(t0, mod_neg(mod_mul(quotient, t1, prime), prime), prime)
    assert len(r0) == 1 and r0[0]
    inv = pow(r0[0], -1, prime)
    s0 = [(x*inv) % prime for x in s0]
    t0 = [(x*inv) % prime for x in t0]
    assert mod_add(mod_mul(a, s0, prime), mod_mul(b, t0, prime), prime) == [1]
    return s0, t0


def moment_mod(n, gamma, q, prime):
    out = []
    product = 1
    old = 0
    for k in range(n+1):
        for j in range(old, q*k):
            product = product*(gamma+j) % prime
        old = q*k
        out.append(comb(n, k)*product % prime)
    assert out[0] == 1 and out[-1] != 0
    return out


def check_actual_pairs():
    prime = 1000003
    assert all(prime % j for j in range(2, isqrt(prime)+1))
    rows = []
    for p, s, w in [(2, 1, 1), (2, 3, 2), (2, 6, 6), (2, 7, 6),
                    (3, 1, 1), (3, 4, 4), (3, 5, 4), (5, 2, 2),
                    (5, 3, 2), (7, 2, 2), (11, 2, 2), (13, 1, 1)]:
        a, q = p**w, p**s
        for beta in sorted({1, a, max(1, a//2)}):
            for n in [1, 2, 3, 5]:
                assert a*(n+1)+beta+q*(n+1)-1 < prime
                lower = moment_mod(n, a*n+beta, q, prime)
                upper = moment_mod(n+1, a*(n+1)+beta, q, prime)
                U, V = bezout(lower, upper, prime)
                rows.append({"p": p, "s": s, "w": w, "n": n, "beta": beta,
                             "lower": lower, "upper": upper,
                             "bezout_lower": U, "bezout_upper": V})
    return {"certificate_prime": prime, "pair_count": len(rows),
            "coefficient_orientation": "ascending",
            "bezout_identity": "bezout_lower*lower + bezout_upper*upper = 1",
            "witnesses": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repo_root is None:
        candidates = [p for p in Path(__file__).resolve().parents
                      if (p/PROOF_DIR).is_dir() and (p/"AGENTS.md").is_file()]
        assert candidates, "Repository root not found"
        repo = candidates[0]
    else:
        repo = args.repo_root.resolve()
    pins = {}
    for name, expected in SOURCE_PINS.items():
        path = repo/PROOF_DIR/name
        actual = sha256(path.read_bytes()).hexdigest()
        assert actual == expected, f"Reviewed source pin mismatch: {name}"
        pins[str(PROOF_DIR/name)] = actual
    receipt = {
        "status": "PASS",
        "scope": "Exact finite controls; unrestricted quantifiers are supplied by the pinned proofs.",
        "source_sha256": pins,
        "functional_identities": check_functionals(),
        "contiguous_block_valuation": check_blocks(),
        "uniform_power_cancellation": check_power_cancellation(),
        "general_ordinary_threshold": check_general_threshold(),
        "actual_moment_pair_controls": check_actual_pairs(),
        "limits": ["No assertion of unrestricted SFC.",
                   "A surviving valuation class is not a counterexample.",
                   "Gamma application here requires integer 1 <= beta <= a."],
    }
    receipt["checker_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    output = args.output or repo/PROOF_DIR/"sfc_prime_power_divided_difference_verification.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True)+"\n")
    print("PASS: "
          f"{receipt['functional_identities']['identity_count']} functional identities; "
          f"{receipt['contiguous_block_valuation']['composition_count']} block controls; "
          f"{receipt['uniform_power_cancellation']['exact_cancellation_count']} cancellations; "
          f"{receipt['actual_moment_pair_controls']['pair_count']} actual-pair Bezout certificates.")
    print(output)


if __name__ == "__main__":
    main()
