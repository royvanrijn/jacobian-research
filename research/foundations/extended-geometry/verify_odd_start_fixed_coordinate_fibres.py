#!/usr/bin/env python3
"""Exact bounded controls for the written fixed-coordinate fibre theorem.

Standard library only.  This does not enumerate all endpoint coordinates or
prove the universal finite-translate theorem by finite testing.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as Q
from hashlib import sha256
from math import comb
from pathlib import Path
import json
import platform


def trim(p):
    p = list(map(Q, p))
    while p and not p[-1]:
        p.pop()
    return tuple(p)


ZERO, ONE, X = (), (Q(1),), (Q(0), Q(1))


def add(p, q):
    return trim((p[i] if i < len(p) else 0) +
                (q[i] if i < len(q) else 0)
                for i in range(max(len(p), len(q))))


def scale(p, c):
    return trim(v * c for v in p)


def sub(p, q):
    return add(p, scale(q, -1))


def mul(p, q, bound=None):
    if not p or not q:
        return ZERO
    n = len(p) + len(q) - 1
    if bound is not None:
        n = min(n, bound + 1)
    out = [Q(0)] * n
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            if i + j < n:
                out[i + j] += a * b
    return trim(out)


def power(p, n, bound=None):
    assert n >= 0
    out = ONE
    while n:
        if n & 1:
            out = mul(out, p, bound)
        p = mul(p, p, bound)
        n //= 2
    return out


def divmod_poly(p, q):
    assert q
    rem = list(p)
    out = [Q(0)] * max(0, len(p) - len(q) + 1)
    while len(rem) >= len(q):
        k = len(rem) - len(q)
        c = rem[-1] / q[-1]
        out[k] += c
        for j, a in enumerate(q):
            rem[j + k] -= c * a
        while rem and not rem[-1]:
            rem.pop()
    return trim(out), tuple(rem)


def gcd_poly(p, q):
    while q:
        p, q = q, divmod_poly(p, q)[1]
    return scale(p, 1 / p[-1]) if p else ZERO


def derivative(p):
    return trim(i * p[i] for i in range(1, len(p)))


def evaluate(p, value):
    out = Q(0)
    for a in reversed(p):
        out = out * value + a
    return out


def compose_series(p, q, bound):
    out = ZERO
    for a in reversed(p):
        out = add(mul(out, q, bound), (a,))
    return trim(out[:bound + 1])


@dataclass(frozen=True)
class Rat:
    n: tuple
    d: tuple = ONE

    def __post_init__(self):
        n, d = trim(self.n), trim(self.d)
        assert d
        if not n:
            object.__setattr__(self, "n", ZERO)
            object.__setattr__(self, "d", ONE)
            return
        g = gcd_poly(n, d)
        n, rn = divmod_poly(n, g)
        d, rd = divmod_poly(d, g)
        assert not rn and not rd
        c = 1 / d[-1]
        object.__setattr__(self, "n", scale(n, c))
        object.__setattr__(self, "d", scale(d, c))

    @staticmethod
    def cast(v):
        return v if isinstance(v, Rat) else Rat((Q(v),))

    def __add__(self, other):
        other = self.cast(other)
        return Rat(add(mul(self.n, other.d), mul(other.n, self.d)),
                   mul(self.d, other.d))

    __radd__ = __add__

    def __neg__(self):
        return Rat(scale(self.n, -1), self.d)

    def __sub__(self, other):
        return self + -self.cast(other)

    def __rsub__(self, other):
        return self.cast(other) - self

    def __mul__(self, other):
        other = self.cast(other)
        return Rat(mul(self.n, other.n), mul(self.d, other.d))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.cast(other)
        assert other.n
        return Rat(mul(self.n, other.d), mul(self.d, other.n))

    def __rtruediv__(self, other):
        return self.cast(other) / self

    def __pow__(self, n):
        if n < 0:
            assert self.n
            return Rat(power(self.d, -n), power(self.n, -n))
        return Rat(power(self.n, n), power(self.d, n))

    def diff(self):
        return Rat(sub(mul(derivative(self.n), self.d),
                       mul(self.n, derivative(self.d))), power(self.d, 2))

    def compose(self, other):
        other = self.cast(other)

        def ev(p):
            out = Rat(ZERO)
            for c in reversed(p):
                out = out * other + c
            return out

        return ev(self.n) / ev(self.d)

    def series(self, bound):
        assert self.d[0]
        out = [Q(0)] * (bound + 1)
        for i in range(bound + 1):
            v = self.n[i] if i < len(self.n) else Q(0)
            for j in range(1, min(i, len(self.d) - 1) + 1):
                v -= self.d[j] * out[i - j]
            out[i] = v / self.d[0]
        return trim(out)

    def polynomial(self):
        return len(self.d) == 1

    def as_poly(self):
        assert self.polynomial()
        return scale(self.n, 1 / self.d[0])

    def at(self, value):
        return evaluate(self.n, value) / evaluate(self.d, value)

    def degree_at_infinity(self):
        assert self.n
        return len(self.n) - len(self.d)


RX, RONE = Rat(X), Rat(ONE)


def inverse_series(u, bound):
    coeffs = u.series(bound)
    assert len(coeffs) >= 2 and coeffs[0] == 0 and coeffs[1] == 1
    v = X
    for n in range(2, bound + 1):
        trial = compose_series(coeffs, v, n)
        c = trial[n] if len(trial) > n else 0
        v = add(v, (Q(0),) * n + (-c,))
    assert compose_series(coeffs, v, bound) == X
    return v


def marked_return(u, bound):
    return compose_series(inverse_series(u, bound),
                          scale(u.series(bound), -1), bound)


def raw_moment(c, d, a, m, marker=0):
    if a * m % 2 or m < marker:
        return Q(0)
    n = m - marker
    out = Q(0)
    for j in range(m + 1):
        term = mul(power(c[:n + 1], m - j, n),
                   power(d[:n + 1], j, n), n)
        coeff = term[n] if len(term) > n else 0
        out += Q(comb(m, j), a * m + 2 * j + 1) * coeff
    return out


def kappa(a, m):
    return sum((Q((-1) ** j * comb(m, j), a * m + 2 * j + 1)
                for j in range(m + 1)), Q(0))


def jet_holds(eta, modulus):
    # eta == -1 in the local quotient; its denominator must be a unit.
    if modulus == ONE:
        return True
    if len(gcd_poly(eta.d, modulus)) != 1:
        return False
    return not divmod_poly(add(eta.n, eta.d), modulus)[1]


def construct(u, eta, a):
    b = u / (RX * u.diff())
    xx = b * (1 + eta) / 2
    c = RX * xx ** a / u
    d = -RX * xx ** (a + 2) / u
    return b, xx, c, d


def verify_pair(u, eta, a, expected_degree, counts, return_rational=None):
    b, xx, c, d = construct(u, eta, a)
    assert c.polynomial() and d.polynomial()
    cp, dp = c.as_poly(), d.as_poly()
    assert cp[0] == 1 and dp[0] == -1 and xx.at(0) == 1
    assert RX * xx ** a / c == u
    assert -(d / c) == xx ** 2
    assert max(len(cp) + a - 2, len(dp) + a) == expected_degree
    assert (2 * xx / b - 1) == eta
    tt = b / xx
    if return_rational is not None:
        j = return_rational
        assert u.compose(j) == -u
        assert eta * eta.compose(j) == RONE
        assert tt + tt.compose(j) == Rat((2,))
        counts["rational_return_identities"] += 1
    else:
        bound = 10
        j = marked_return(u, bound)
        assert mul(eta.series(bound),
                   compose_series(eta.series(bound), j, bound), bound) == ONE
        assert add(tt.series(bound),
                   compose_series(tt.series(bound), j, bound)) == (Q(2),)
        counts["nonrational_return_series_controls"] += 1
    for m in (2, 4, 6):
        assert raw_moment(cp, dp, a, m) == 0
        counts["raw_even_moments"] += 1
    counts["polynomial_pairs"] += 1
    return cp, dp


def run_controls():
    counts = {k: 0 for k in (
        "polynomial_pairs", "raw_even_moments", "rational_return_identities",
        "nonrational_return_series_controls", "marker_controls",
        "finite_empty_fibre_candidates", "jet_exponent_cases",
        "higher_height_rejections", "negative_controls")}
    samples = []

    # Entire degree-one-coordinate fibres, including k=0 and k=1.
    for lam in (Q(1), Q(-2)):
        u = RX / (1 + lam * RX)
        rr = 1 + 2 * lam * RX
        jj = -RX / (1 + 2 * lam * RX)
        for a in (1, 3, 5):
            for k in range(6):
                cp, dp = verify_pair(u, rr ** k, a, (a + 2) * (k + 2),
                                     counts, jj)
                if k in (0, 1):
                    for m in (2, 4, 6):
                        if k == 0:
                            assert raw_moment(cp, dp, a, m, 1) == \
                                kappa(a, m) * lam ** (m - 1)
                        else:
                            assert raw_moment(cp, dp, a, m, 1) == 0
                            assert raw_moment(cp, dp, a, m, 2) == \
                                kappa(a, m) * lam ** (m - 2)
                        counts["marker_controls"] += 1

    # Safe polynomial coordinate and a nonlinear polynomial-coordinate failure.
    verify_pair(RX, RONE, 3, 4, counts, -RX)
    u_bad = RX + RX ** 2
    _, _, c_bad, d_bad = construct(u_bad, RONE, 3)
    assert not c_bad.polynomial() and not d_bad.polynomial()
    counts["negative_controls"] += 1

    # Hyperbolic reciprocal quotient: exact congruence k == 2 (mod 4).
    rr = 1 + RX
    u = (rr - 1 / rr) / 2
    jj = -RX / (1 + RX)
    modulus = (rr ** 2 + 1).as_poly()
    assert (rr ** 4 - 1).n and not divmod_poly((rr ** 4 - 1).n, modulus)[1]
    assert divmod_poly((rr ** 2 - 1).n, modulus)[1]
    for a in (1, 3, 5):
        for k in range(-1, 15):
            eta = rr ** k
            accept = k % 4 == 2
            assert jet_holds(eta, modulus) == accept
            _, _, c, d = construct(u, eta, a)
            assert (c.polynomial() and d.polynomial()) == accept
            counts["jet_exponent_cases"] += 1
            if accept:
                verify_pair(u, eta, a, (a + 2) * k + a + 1, counts, jj)

    # A genuinely nonrational original return, with a polynomial quotient R.
    rr = 1 + RX + RX ** 2 / 8
    u = 2 * (rr - 1) / (rr + 1)
    assert u == 2 * RX * (RX + 8) / (RX + 4) ** 2
    assert u / (RX * u.diff()) == (RX + 4) * (RX + 8) / 32
    for a in (1, 3, 5):
        for k in range(5):
            verify_pair(u, rr ** k, a, (a + 2) * (2 * k + 3), counts)
    assert rr * rr.compose(-RX / (1 + RX)) != RONE
    counts["negative_controls"] += 1

    # Higher-height local cancellation: the same eta fails at a=1.
    u = (rr - 1 / rr) / 2
    critical_modulus = (rr ** 2 + 1).as_poly()
    height_one_modulus = mul(critical_modulus, (Q(4), Q(1)))
    for k in (2, 6, 10):
        eta = rr ** k
        assert jet_holds(eta, critical_modulus)
        assert not jet_holds(eta, height_one_modulus)
        _, _, c1, d1 = construct(u, eta, 1)
        assert not c1.polynomial() or not d1.polynomial()
        counts["higher_height_rejections"] += 1
        for a in (3, 5):
            verify_pair(u, eta, a, 2 * (a + 2) * k + a, counts)

    # A finite empty fibre: marked return is y=-x, not the other xy=1 factor.
    u = RX / (1 - RX ** 2)
    modulus = (1 + RX ** 2).as_poly()
    for k in (-1, 0, 1):
        eta = ((1 + RX) / (1 - RX)) ** k
        assert eta * eta.compose(-RX) == RONE
        assert not jet_holds(eta, modulus)
        for a in (1, 3, 5):
            _, _, c, d = construct(u, eta, a)
            assert not c.polynomial() or not d.polynomial()
        counts["finite_empty_fibre_candidates"] += 1

    # A root-of-unity value alone misses the multiple-root jet obstruction.
    modulus = power((Q(-1), Q(1)), 2)
    rr_jet = 1 - 2 * RX
    assert rr_jet.at(1) == -1
    for k in (2, 4, 6, 8, 10, 12):
        assert divmod_poly((rr_jet ** k - 1).n, modulus)[1]
    counts["negative_controls"] += 1

    # An off-base component can satisfy reciprocity while the marked one fails.
    u = RX * (RX - 1) * (RX - 2) / 2
    h = u * (RX - 1) ** 2
    eta = (1 + h) / (1 - h)
    assert u.at(0) == 0 and u.diff().at(0) == 1 and eta.at(0) == 1
    assert u.compose(2 - RX) == -u
    assert eta * eta.compose(2 - RX) == RONE
    j = marked_return(u, 10)
    assert mul(eta.series(10), compose_series(eta.series(10), j, 10), 10) != ONE
    counts["negative_controls"] += 1

    # A polynomial, normalized corrupted phase must fail actual moments too.
    u = RX / (1 + RX)
    eta = (1 + 2 * RX) ** 2 + RX
    _, _, c, d = construct(u, eta, 1)
    assert c.polynomial() and d.polynomial() and eta.at(0) == 1
    assert eta * eta.compose(-RX / (1 + 2 * RX)) != RONE
    corrupt_moments = [raw_moment(c.as_poly(), d.as_poly(), 1, m)
                       for m in (2, 4, 6)]
    assert any(corrupt_moments)
    counts["negative_controls"] += 1
    samples.append({"control": "corrupted_coefficient_raw_moments",
                    "orders": [2, 4, 6],
                    "values": list(map(str, corrupt_moments))})

    # Wrong reciprocal constant and wrong normalization are separate failures.
    assert (-RONE).at(0) != 1
    assert (1 + 2 * RX) * (-(1 + 2 * RX)).compose(-RX / (1 + 2 * RX)) != RONE
    counts["negative_controls"] += 1
    return counts, samples


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path,
                        default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    directory = Path("research/foundations/extended-geometry")
    ledger_path = directory / "odd_start_fixed_coordinate_fibres_source_ledger.json"
    ledger = json.loads((repo / ledger_path).read_text())
    for source in ledger["sources"]:
        actual = digest(repo / source["path"])
        assert actual == source["sha256"], (source["path"], actual)
    counts, samples = run_controls()
    proof = directory / "ODD_START_FIXED_COORDINATE_FIBRE_CLASSIFICATION.md"
    checker = directory / Path(__file__).name
    result = {
        "schema": "odd_start_fixed_coordinate_fibres_verification_v1",
        "status": "PASS",
        "arithmetic": "fractions.Fraction; Python standard library only",
        "python": platform.python_version(),
        "counts": counts,
        "samples": samples,
        "source_pins_checked": len(ledger["sources"]),
        "packet_pins": {str(p): digest(repo / p)
                        for p in (proof, checker, ledger_path)},
        "scope": "Bounded exact controls for the displayed fibres and failure cases. "
                 "The written proof supplies the universal classification, "
                 "Dickson finiteness and arithmetic-progression conclusions. "
                 "No enumeration of all endpoint coordinates or all fibre bases.",
        "assurance": "Internal author replay; not external or proof-assistant verification."
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "counts": counts,
                      "source_pins_checked": result["source_pins_checked"],
                      "receipt": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
