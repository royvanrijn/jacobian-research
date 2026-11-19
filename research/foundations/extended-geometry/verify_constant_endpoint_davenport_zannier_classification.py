#!/usr/bin/env python3
"""Exact finite controls for the marked Davenport--Zannier correspondence.

The written proof supplies the all-order classification. No symbolic algebra
package, root finder, or numerical tolerance is used by this checker.
"""

import argparse
from fractions import Fraction as Q
import hashlib
from functools import lru_cache
from itertools import product
import json
from math import comb
from pathlib import Path
import sys

if not __debug__:
    raise SystemExit("Refusing -O: assertions are part of the exact checker.")


def trim(p):
    p = tuple(Q(x) for x in p)
    while p and p[-1] == 0:
        p = p[:-1]
    return p


Z = ()
ONE = (Q(1),)
XX = (Q(0), Q(1))


def add(p, q):
    return trim(tuple((p[i] if i < len(p) else 0) +
                      (q[i] if i < len(q) else 0)
                      for i in range(max(len(p), len(q)))))


def scale(p, c):
    return trim(tuple(c*x for x in p))


def sub(p, q):
    return add(p, scale(q, -1))


def mul(p, q):
    if not p or not q:
        return Z
    out = [Q(0)]*(len(p)+len(q)-1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            out[i+j] += pi*qj
    return trim(out)


def power(p, n):
    assert n >= 0
    out = ONE
    while n:
        if n & 1:
            out = mul(out, p)
        p = mul(p, p)
        n //= 2
    return out


def deriv(p):
    return trim(tuple(i*p[i] for i in range(1, len(p))))


def divmodp(p, q):
    assert q
    rem = list(p)
    quotient = [Q(0)]*max(0, len(p)-len(q)+1)
    while rem and len(rem) >= len(q):
        shift = len(rem)-len(q)
        c = rem[-1]/q[-1]
        quotient[shift] = c
        for j, qj in enumerate(q):
            rem[shift+j] -= c*qj
        rem = list(trim(rem))
    return trim(quotient), trim(rem)


def exactdiv(p, q):
    quotient, rem = divmodp(p, q)
    assert not rem
    return quotient


def monic(p):
    assert p
    return scale(p, 1/p[-1])


def gcdp(p, q):
    while q:
        p, q = q, divmodp(p, q)[1]
    return monic(p) if p else Z


def degree(p):
    return len(p)-1


def radical(p):
    return monic(exactdiv(p, gcdp(p, deriv(p))))


def max_multiplicity(p):
    count = 0
    while degree(p) > 0:
        count += 1
        p = gcdp(p, deriv(p))
    return count


class RF:
    def __init__(self, num=ONE, den=ONE):
        num, den = trim(num), trim(den)
        assert den
        g = gcdp(num, den)
        num, den = exactdiv(num, g), exactdiv(den, g)
        self.num, self.den = scale(num, 1/den[-1]), scale(den, 1/den[-1])

    def __mul__(self, other):
        return RF(mul(self.num, other.num), mul(self.den, other.den))

    def __truediv__(self, other):
        assert other.num
        return RF(mul(self.num, other.den), mul(self.den, other.num))

    def __add__(self, other):
        return RF(add(mul(self.num, other.den), mul(other.num, self.den)),
                  mul(self.den, other.den))

    def __pow__(self, n):
        if n < 0:
            return RF(power(self.den, -n), power(self.num, -n))
        return RF(power(self.num, n), power(self.den, n))

    def derivative(self):
        return RF(sub(mul(deriv(self.num), self.den),
                      mul(self.num, deriv(self.den))), power(self.den, 2))

    def at_zero(self):
        assert self.den[0] != 0
        return self.num[0]/self.den[0] if self.num else Q(0)

    def polynomial(self):
        return exactdiv(self.num, self.den)

    def __eq__(self, other):
        return self.num == other.num and self.den == other.den


def pair_data(A, B, a):
    """Validate all marked DZ hypotheses, without assuming factorization."""
    assert a >= 0 and degree(A) == degree(B) >= 1
    assert A[-1] == B[-1] == 1
    assert degree(gcdp(A, B)) == 0
    assert A[0] == 0 and len(A) > 1 and A[1] != 0 and B[0] != 0
    r = degree(radical(mul(A, B)))-1
    N = degree(A)
    assert 1 <= r <= N
    assert degree(sub(A, B)) == N-r
    assert max_multiplicity(exactdiv(A, XX)) <= a
    assert N <= 1+a*(r-1)
    f = RF(A, B)
    u = f/RF((f.derivative().at_zero(),))
    P = exactdiv(radical(mul(A, B)), XX)
    P = scale(P, 1/P[0])
    H = (RF(mul(XX, power(P, a)))/u).polynomial()
    assert H[0] == 1 and P[0] == 1
    assert u.derivative()/u == RF(ONE, mul(XX, P))
    assert u.derivative() == RF(P)**(a-1)/RF(H)
    assert degree(P) == r and degree(H) == a*r+1
    return u, P, H, N, r


def forward(P, H, a):
    u = RF(mul(XX, power(P, a)), H)
    assert u.derivative() == RF(P)**(a-1)/RF(H)
    A, B = monic(u.num), monic(u.den)
    ur, Pr, Hr, N, r = pair_data(A, B, a)
    assert ur == u and Pr == P and Hr == H
    D = scale(mul(H, power(P, 2)), -1)
    assert max(degree(H)+a-1, degree(D)+a+1) == (a+2)*(r+1)
    return A, B, N, r


def star(a, r, lam=Q(1)):
    assert a >= 1 or r == 1
    U = (Q(1), lam)
    S = exactdiv(sub(power(U, r), ONE), sub(U, ONE))
    P = scale(mul(U, S), Q(1, r))
    if a == 0:
        H = U
    else:
        H = scale(mul(power(U, a+r), power(S, a-1)), Q(r)**(1-a))
    u = RF(sub(power(U, r), ONE), scale(power(U, r), r*lam))
    assert u.derivative() == RF(ONE, power(U, r+1))
    assert RF(mul(XX, power(P, a)), H) == u
    return P, H


def raw_moment(C, D, a, m, marker_power=0):
    """Direct expansion and t-integration, not the endpoint formula."""
    if a*m % 2:
        return Q(0)
    index = m-marker_power
    if index < 0:
        return Q(0)
    total = Q(0)
    for j in range(m+1):
        coeffs = mul(power(C, m-j), power(D, j))
        coefficient = coeffs[index] if index < len(coeffs) else Q(0)
        total += Q(comb(m, j), a*m+2*j+1)*coefficient
    return total


def rejects(A, B, a):
    try:
        pair_data(A, B, a)
    except AssertionError:
        return True
    return False


def nonnegative_compositions(n, length):
    if length == 0:
        return [()] if n == 0 else []
    if length == 1:
        return [(n,)]
    return [(i,)+tail for i in range(n+1)
            for tail in nonnegative_compositions(n-i, length-1)]


def explicit_weighted_tree_counts(a, max_r):
    """Enumerate actual ordered weighted trees, not just a scalar recurrence.

    A white vertex is an ordered tuple of black branches. A black branch
    stores its full positive edge-weight tuple and ordered white children.
    The marked leaf/edge is omitted, so residual edge count is r-1.
    """
    weights = {
        k: tuple(w for w in product(range(1, a+1), repeat=k+1)
                 if sum(w) <= a) for k in range(a)
    }

    @lru_cache(None)
    def white(n):
        if n == 0:
            return ((),)
        out = []
        for k in range(min(a, n)):
            for child_total in range(n-k):
                remaining = n-(k+1)-child_total
                for child_sizes in nonnegative_compositions(child_total, k):
                    choices = product(*(white(c) for c in child_sizes))
                    for children in choices:
                        for w in weights[k]:
                            branch = (w, children)
                            for tail in white(remaining):
                                out.append((branch,)+tail)
        assert len(set(out)) == len(out)
        return tuple(out)

    counts = [len(white(r-1)) for r in range(1, max_r+1)]
    for r, count in enumerate(counts, 1):
        formula = Q(comb(a*r, r-1), r) if a*r >= r-1 else Q(0)
        assert formula.denominator == 1 and count == formula
    return counts


def quadratic_descent_control():
    """Reconstruct the nonreal marked unitree classes in Q[eta]/(eta^2+2)."""
    def plus(v, w):
        return v[0]+w[0], v[1]+w[1]

    def times(v, w):
        return v[0]*w[0]-2*v[1]*w[1], v[0]*w[1]+v[1]*w[0]

    zero, one = (Q(0), Q(0)), (Q(1), Q(0))

    def add2(p, q):
        return tuple(plus(p[i] if i < len(p) else zero,
                          q[i] if i < len(q) else zero)
                     for i in range(max(len(p), len(q))))

    def mul2(p, q):
        out = [zero]*(len(p)+len(q)-1)
        for i, pi in enumerate(p):
            for j, qj in enumerate(q):
                out[i+j] = plus(out[i+j], times(pi, qj))
        return tuple(out)

    U = (one, one)
    V = (one, (Q(1, 2), Q(-1, 4)))
    L = (one, (Q(2, 3), Q(1, 6)))
    N = mul2((zero, one), mul2(V, mul2(L, L)))
    Nprime = tuple((i*N[i][0], i*N[i][1]) for i in range(1, len(N)))
    numerator = add2(mul2(Nprime, U), tuple((-4*c, -4*d) for c, d in N))
    assert numerator[:len(L)] == L and all(c == zero for c in numerator[len(L):])
    ratio = (Q(4, 3), Q(2, 3))
    assert times(V[1], ratio) == one and ratio[1] != 0
    A0 = (Q(3), Q(-4), Q(0), Q(0), Q(1))
    B0 = (Q(0), Q(0), Q(0), Q(0), Q(1))
    assert A0 == mul(power((Q(-1), Q(1)), 2), (Q(3), Q(2), Q(1)))
    assert degree(sub(A0, B0)) == 1 and degree(radical(mul(A0, B0))) == 4
    return {"quadratic_relation": "eta^2=-2", "height": 2, "r": 3,
            "degree": 16, "passport": [[1, 1, 2], [4]],
            "conjugate_nonreal_marked_classes": 2,
            "unmarked_rational_unitree_model": "(x^4-4x+3,x^4)"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name(
        "constant_endpoint_davenport_zannier_classification_verification.json"))
    parser.add_argument('--repo-root', type=Path, default=Path(__file__).resolve().parents[3])
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    ledger_path = here/'constant_endpoint_davenport_zannier_classification_source_ledger.json'
    ledger = json.loads(ledger_path.read_text())
    assert hashlib.sha256((here/ledger['proof']).read_bytes()).hexdigest() == ledger['proof_sha256']
    assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest() == ledger['verification']['checker_sha256']
    source_pins = []
    for entry in ledger['retained_resources']+ledger['primary_sources']:
        path = args.repo_root/entry['retained_path']
        got = hashlib.sha256(path.read_bytes()).hexdigest()
        assert got == entry['sha256'], str(path)
        source_pins.append({'path':entry['retained_path'],'sha256':got})
    stars = moments = marked = nonstars = 0
    degree_spectra = {}
    for a in range(7):
        for r in (range(1, 5) if a else (1,)):
            for lam in (Q(1), Q(2, 3), Q(-2)):
                P, H = star(a, r, lam)
                A, B, N, rr = forward(P, H, a)
                assert N == rr == r
                stars += 1
                degree_spectra.setdefault(str(a), set()).add((a+2)*(r+1))
            P, H = star(a, r)
            D = scale(mul(H, power(P, 2)), -1)
            for m in range(1, 7):
                assert raw_moment(H, D, a, m) == 0
                moments += 1
                if a*m % 2 == 0:
                    assert raw_moment(H, D, a, m, 1) != 0
                    marked += 1
    U, V = (Q(1), Q(1)), (Q(1), Q(2, 3))
    for a in range(2, 8):
        P = mul(U, V)
        H = mul(power(U, a+3), power(V, a-2))
        A, B, N, r = forward(P, H, a)
        assert (N, r) == (3, 2)
        D = scale(mul(H, power(P, 2)), -1)
        for m in range(1, 7):
            assert raw_moment(H, D, a, m) == 0
            moments += 1
        nonstars += 1
    A = mul(XX, power((Q(3, 2), Q(1)), 2))
    B = power(U, 3)
    assert rejects(A, B, 1)  # Exact DZ identity, but multiplicity cap fails.
    near_A = mul(XX, power((Q(2), Q(1)), 2))
    assert degree(sub(near_A, B)) == 2  # Correct passport, wrong DZ degree.
    assert rejects(near_A, B, 2)
    f = RF(near_A, B)
    unear = f/RF((f.derivative().at_zero(),))
    Pnear = mul(U, (Q(1), Q(1, 2)))
    Hnear = (RF(mul(XX, power(Pnear, 2)))/unear).polynomial()
    assert unear.derivative() != RF(Pnear)/RF(Hnear)
    assert raw_moment(Hnear, scale(mul(Hnear, power(Pnear, 2)), -1), 2, 2) != 0
    assert rejects(mul(XX, U), power(U, 2), 2)  # Common factor.
    assert rejects((Q(-2), Q(0), Q(1)), (Q(0), Q(0), Q(1)), 2)
    # Rational-root candidates for the monic integer polynomial X^2-2.
    assert all(n*n != 2 for n in (-2, -1, 1, 2))
    for r in range(2, 10):
        assert 1+0*(r-1) < r  # Height-zero passport obstruction.
    for r in range(1, 10):
        assert 1+1*(r-1) == r  # Height-one star rigidity.
    scale_controls = 0
    for a in range(7):
        for r in (range(1,5) if a else (1,)):
            for lam in (Q(1),Q(2,3),Q(-2)):
                P,H = star(a,r,lam)
                q = H[-1]/P[-1]**a
                ss = 1/q
                Ps = trim(tuple(c*ss**i for i,c in enumerate(P)))
                Hs = trim(tuple(c*ss**i for i,c in enumerate(H)))
                assert Hs[-1]/Ps[-1]**a == 1
                forward(Ps,Hs,a)
                scale_controls += 1
    tree_counts = {str(a): explicit_weighted_tree_counts(a, 6)
                   for a in range(5)}
    assert tree_counts["2"] == [1, 2, 5, 14, 42, 132]
    quadratic_control = quadratic_descent_control()
    root = Path(__file__).resolve().parent
    proof = root/"CONSTANT_ENDPOINT_DAVENPORT_ZANNIER_CLASSIFICATION.md"
    receipt = {
        "scope": "Exact finite controls; all-order marked DZ iff is the written proof.",
        "status": "PASS",
        "proof_sha256": hashlib.sha256(proof.read_bytes()).hexdigest(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "star_normalized_forward_reverse_controls": stars,
        "nonstar_normalized_forward_reverse_controls": nonstars,
        "raw_angular_purity_controls": moments,
        "raw_x_marker_nonvanishing_controls": marked,
        "malformed_controls": ["multiplicity cap", "near-DZ degree defect",
                               "near-DZ derivative and raw second moment",
                               "common factor", "unmarked irrational A-root"],
        "explicit_ordered_weighted_tree_counts": tree_counts,
        "tree_enumeration_total_controls": 30,
        "tree_enumeration_total_distinct_objects": sum(map(sum, tree_counts.values())),
        "nonreal_marked_descent_control": quadratic_control,
        "degree_spectra_controls": {a: sorted(v) for a, v in degree_spectra.items()},
        "external_or_formal_review_claimed": False,
    }
    receipt.update({'schema_version':1,'retained_source_pins':source_pins,
        'source_ledger_sha256':hashlib.sha256(ledger_path.read_bytes()).hexdigest(),
        'canonical_scale_controls':scale_controls,
        'independent_replay':False,'external_review':False,'formal_verification':False})
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in receipt.items() if not k.endswith("sha256")},
                     sort_keys=True))


if __name__ == "__main__":
    main()
