#!/usr/bin/env python3
"""Exact bounded controls for the global infinite-fibre classification.

Standard library only. The written proof supplies the unrestricted theorem.
These controls do not enumerate all rational coordinates or prove exhaustion.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("Run this exact checker without -O or -OO.")

import argparse
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import gcd, lcm
from pathlib import Path
import json
import platform
import sys


PREFIX = Path("research/foundations/extended-geometry")
PROOF = PREFIX / "GLOBAL_CYCLOTOMIC_INFINITE_FIBRE_CLASSIFICATION.md"
LEDGER = PREFIX / "global_cyclotomic_infinite_fibres_source_ledger.json"
HELPER = PREFIX / "verify_odd_start_fixed_coordinate_fibres.py"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def quadratic_field(radicand):
    """Exact Q[s]/(s^2-radicand), used only for radicands -1 and 2."""
    assert radicand in (-1, 2)

    @dataclass(frozen=True, eq=False)
    class Element:
        a: Fraction
        b: Fraction = Fraction(0)

        @staticmethod
        def cast(value):
            if isinstance(value, Element):
                return value
            if isinstance(value, (int, Fraction)):
                return Element(Fraction(value))
            return NotImplemented

        def __bool__(self):
            return bool(self.a or self.b)

        def __eq__(self, other):
            other = self.cast(other)
            return NotImplemented if other is NotImplemented else (
                self.a == other.a and self.b == other.b)

        def __hash__(self):
            return hash(self.a) if not self.b else hash((radicand, self.a, self.b))

        def __add__(self, other):
            other = self.cast(other)
            return NotImplemented if other is NotImplemented else Element(
                self.a + other.a, self.b + other.b)

        __radd__ = __add__

        def __neg__(self):
            return Element(-self.a, -self.b)

        def __sub__(self, other):
            other = self.cast(other)
            return NotImplemented if other is NotImplemented else self + -other

        def __rsub__(self, other):
            other = self.cast(other)
            return NotImplemented if other is NotImplemented else other + -self

        def __mul__(self, other):
            other = self.cast(other)
            return NotImplemented if other is NotImplemented else Element(
                self.a * other.a + radicand * self.b * other.b,
                self.a * other.b + self.b * other.a)

        __rmul__ = __mul__

        def inverse(self):
            norm = self.a * self.a - radicand * self.b * self.b
            assert norm
            return Element(self.a / norm, -self.b / norm)

        def __truediv__(self, other):
            other = self.cast(other)
            return NotImplemented if other is NotImplemented else self * other.inverse()

        def __rtruediv__(self, other):
            other = self.cast(other)
            return NotImplemented if other is NotImplemented else other * self.inverse()

        def __pow__(self, exponent):
            assert isinstance(exponent, int)
            if exponent < 0:
                return self.inverse() ** (-exponent)
            out, base = Element(Fraction(1)), self
            while exponent:
                if exponent & 1:
                    out *= base
                base *= base
                exponent //= 2
            return out

    def coefficient(n=0, d=None):
        if d is None:
            value = Element.cast(n)
            assert value is not NotImplemented
            return value
        return Element(Fraction(n, d))

    radical = Element(Fraction(0), Fraction(1))
    assert radical ** 2 == radicand
    for a, b in product(range(-2, 3), repeat=2):
        v = Element(Fraction(a), Fraction(b))
        if v:
            assert v * v.inverse() == 1
    return coefficient, radical


def load_arithmetic(repo, radicand=None):
    """Reuse unchanged polynomial/Rat formulae with an exact coefficient field."""
    name = "global_fibre_arithmetic_" + str(radicand)
    spec = spec_from_file_location(name, repo / HELPER)
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    radical = None
    if radicand is not None:
        module.Q, radical = quadratic_field(radicand)
        module.ZERO = ()
        module.ONE = (module.Q(1),)
        module.X = (module.Q(0), module.Q(1))
        module.RX = module.Rat(module.X)
        module.RONE = module.Rat(module.ONE)
    return module, radical


def equal(left, right):
    assert left == right


def rejected(name, action):
    try:
        action()
    except (AssertionError, ValueError, ZeroDivisionError):
        return name
    raise AssertionError("Negative control accepted: " + name)


def root_multiplicity_gcd(t, p):
    """Characteristic-zero squarefree decomposition, without a factor oracle."""
    assert len(p) > 1
    g = t.gcd_poly(p, t.derivative(p))
    w, rem = t.divmod_poly(p, g)
    assert not rem
    i, result = 1, 0
    while len(w) > 1:
        y = t.gcd_poly(w, g)
        z, rem = t.divmod_poly(w, y)
        assert not rem
        if len(z) > 1:
            result = gcd(result, i)
        w = y
        g, rem = t.divmod_poly(g, y)
        assert not rem
        i += 1
    assert len(g) <= 1
    return result


def remainder(t, rational, modulus):
    """Divisibility of a rational numerator; denominator must be a local unit."""
    assert len(t.gcd_poly(rational.d, modulus)) == 1
    return t.divmod_poly(rational.n, modulus)[1]


def order_at(t, rational, value):
    linear = (-t.Q(value), t.Q(1))

    def multiplicity(p):
        assert p
        n = 0
        while len(p) > 1:
            q, r = t.divmod_poly(p, linear)
            if r:
                break
            n += 1
            p = q
        return n

    assert rational.n
    return multiplicity(rational.n) - multiplicity(rational.d)


def pair_control(t, u, eta, a, counts):
    b, xhat, c, d = t.construct(u, eta, a)
    assert u.at(0) == 0 and u.diff().at(0) == 1
    assert c.polynomial() and d.polynomial()
    cp, dp = c.as_poly(), d.as_poly()
    assert cp[0] == 1 and dp[0] == -1 and xhat.at(0) == 1
    assert t.RX * xhat ** a / c == u and -d / c == xhat ** 2
    for m in range(1, 7):
        assert t.raw_moment(cp, dp, a, m) == 0
        counts["raw_moments"] += 1
    counts["polynomial_pairs"] += 1
    return b, xhat, max(len(cp) + a - 2, len(dp) + a)


def local_and_cyclotomic_controls():
    zero_cases = regular_cases = pole_cases = 0
    for a in (1, 3, 5, 7, 9):
        for e, s in product(range(1, 9), repeat=2):
            q = (e * s + a - 1) // a - 1
            assert (q > 0) == (e * s > a)
            assert (q <= e) == (e * (s - a) <= a)
            zero_cases += 1
        for e, f in product(range(1, 9), repeat=2):
            q = e * f - 1
            assert (q <= e) == (e * (f - 1) <= 1)
            if f > 1 and q <= e:
                assert (e, f) == (1, 2)
            regular_cases += 1
        for e, m, h in product(range(1, 9), range(1, 13), range(6)):
            bound = Fraction(m, a + 2) + Fraction(1, e)
            assert (e * h <= 1 + e * m // (a + 2)) == (h <= bound)
            pole_cases += 1

    def v2(n):
        result = 0
        while n % 2 == 0:
            n //= 2
            result += 1
        return result

    crt_cases = 0
    for length, upper in ((1, 16), (2, 12), (3, 8)):
        for orders in product(range(1, upper + 1), repeat=length):
            period = lcm(*orders)
            actual = [k for k in range(period)
                      if all(n % 2 == 0 and k % n == n // 2 for n in orders)]
            criterion = all(n % 2 == 0 for n in orders) and len(
                {v2(n) for n in orders}) == 1
            assert bool(actual) == criterion
            assert actual == ([period // 2] if criterion else [])
            crt_cases += 1
    return {
        "zero_ramification_cases": zero_cases,
        "regular_ramification_cases": regular_cases,
        "pole_bound_cases": pole_cases,
        "simultaneous_untwisted_order_cases": crt_cases,
    }


def algebra_controls(repo):
    counts = {"polynomial_pairs": 0, "raw_moments": 0,
              "period_residue_cases": 0, "thick_jet_cases": 0}
    negatives, samples = [], []
    t, imaginary = load_arithmetic(repo, -1)
    x, z = t.RX, t.RX
    for e in (2, 3, 4):
        r = imaginary + (1 - imaginary) * (1 + x) ** e
        U = (2 / (e * (1 - imaginary))) * (z - 1) / (z + 1)
        u = U.compose(r)
        assert U.compose(1 / z) == -U
        assert root_multiplicity_gcd(t, r.as_poly()) == 1
        assert len(t.gcd_poly((r - 1).n, t.derivative((r - 1).n))) == 1
        J = ((1 + x) ** (e - 1)).as_poly()
        assert order_at(t, u.diff(), -1) == e - 1
        assert not remainder(t, r - imaginary, J)
        for k in range(8):
            assert (not remainder(t, r ** k + 1, J)) == (k % 4 == 2)
            counts["period_residue_cases"] += 1
        assert not remainder(t, r ** 4 - 1, J)
        assert remainder(t, r ** 4 - 1, ((1 + x) ** (e + 1)).as_poly())
        counts["thick_jet_cases"] += 1
        for a in (1, 3):
            b, _, degree = pair_control(t, u, r ** 2, a, counts)
            assert not b.polynomial() and order_at(t, b, -1) == -(e - 1)
            samples.append({"family": "Qi_thick_jet", "e": e, "a": a,
                            "exponent": 2, "least_sphere_degree": degree})
    negatives.append(rejected("cyclotomic residues with a jet made too thick",
        lambda: equal(remainder(t, r ** 4 - 1, ((1 + x) ** (e + 1)).as_poly()), ())))

    t, root2 = load_arithmetic(repo, 2)
    x, z = t.RX, t.RX
    alpha = root2 - 1
    H = (z - alpha) / (1 - alpha * z)
    critical = (z ** 2 - root2 * z + 1).as_poly()
    assert H.compose(1 / z) * H == t.RONE and H.at(1) == 1 and H.at(-1) == -1
    assert 1 + H * z ** 2 == (z + 1) * (z ** 2 - root2 * z + 1) / (1 - alpha * z)
    assert not remainder(t, z ** 8 - 1, critical)
    assert all(remainder(t, z ** j - 1, critical) for j in range(1, 8))
    for e in (2, 3, 4):
        r = -1 + 2 * (1 + x) ** e
        U = ((1 - root2) / (2 * e)) * (z ** 2 - 1) / (z ** 2 - 2 * root2 * z + 1)
        u = U.compose(r)
        assert U.compose(1 / z) == -U
        assert root_multiplicity_gcd(t, r.as_poly()) == 1
        assert order_at(t, r + 1, -1) == e
        assert len(t.gcd_poly((r - alpha).n, t.derivative((r - alpha).n))) == 1
        assert len(t.gcd_poly((r - 1 / alpha).n, t.derivative((r - 1 / alpha).n))) == 1
        assert U.diff() == ((1 - root2) / (2 * e)) * (-2 * root2) * (
            z ** 2 - root2 * z + 1) / (z ** 2 - 2 * root2 * z + 1) ** 2
        for a in (1, 3, 5):
            for n in (-1, 0, 1):
                survivors = []
                for k in range(8):
                    candidate = H ** n * z ** k
                    critical_ok = not remainder(t, candidate + 1, critical)
                    zero_ok = a >= e or candidate.at(-1) == -1
                    if critical_ok and zero_ok:
                        survivors.append(k)
                    counts["period_residue_cases"] += 1
                expected = [6] if n == -1 else [2] if n == 1 else [4] if a >= e else []
                assert survivors == expected
        # Nontrivial twists at every displayed e; one conjugate inverse-twist
        # branch and the newly admitted pure branch are checked separately.
        cases = [(1, 1, 2), (3, 1, 2)]
        if e == 2:
            cases += [(1, -1, 6), (3, 0, 4)]
        for a, n, k in cases:
            eta = (H ** n * z ** k).compose(r)
            _, xx, degree = pair_control(t, u, eta, a, counts)
            if n == 1 and k == 2:
                explicit = ((r - 1) / x) * (r - alpha) * (
                    1 + x) ** (e + 1) / (2 * e * root2 * alpha)
                assert xx == explicit and explicit.polynomial()
            samples.append({"family": "Qsqrt2_pole_twist", "e": e, "a": a,
                            "twist": n, "exponent": k, "least_sphere_degree": degree})
        if e == 2:
            # Pure k=4 satisfies the ordinary critical values but fails the
            # active -1 value at height1. A square pole twist passes all value
            # tests at height3 but violates the simple-pole allowance.
            bad_pure = t.construct(u, r ** 4, 1)
            assert not (bad_pure[2].polynomial() and bad_pure[3].polynomial())
            bad_pole = t.construct(u, (H ** 2).compose(r), 3)
            assert not (bad_pole[2].polynomial() and bad_pole[3].polynomial())
            negatives += [
                rejected("drop the incompatible active value -1",
                         lambda: equal(bad_pure[2].polynomial() and bad_pure[3].polynomial(), True)),
                rejected("permit a pole twist beyond the exact bound",
                         lambda: equal(bad_pole[2].polynomial() and bad_pole[3].polynomial(), True)),
            ]

    t, _ = load_arithmetic(repo)
    x, z = t.RX, t.RX
    w = (z - 1) / (z + 1)
    U = 2 * (w + Fraction(2, 3) * w ** 3 + Fraction(1, 5) * w ** 5)
    r = 1 + x
    u = U.compose(r)
    assert U.compose(1 / z) == -U and U.diff().at(1) == 1
    assert len(t.gcd_poly(U.diff().n, t.derivative(U.diff().n))) > 1
    modulus = ((r ** 2 + 1) ** 2).as_poly()
    eta = r ** 2
    assert not remainder(t, eta + 1, (r ** 2 + 1).as_poly())
    assert remainder(t, eta + 1, modulus)
    bad = t.construct(u, eta, 1)
    assert not (bad[2].polynomial() and bad[3].polynomial())
    negatives.append(rejected("ignore higher quotient ramification",
        lambda: equal(bad[2].polynomial() and bad[3].polynomial(), True)))

    # Marked z=1 and nonunit z=0 controls use canonical nonpower quotients.
    marked = 1 + x * (1 + x) ** 2
    assert root_multiplicity_gcd(t, marked.as_poly()) == 1
    assert order_at(t, marked - 1, -1) == 2
    assert all((marked ** k).at(-1) == 1 for k in range(8))
    negatives.append(rejected("ignore the active marked value 1",
                              lambda: equal((marked ** 3).at(-1), -1)))
    zero = (1 + x) ** 2 * (1 + 2 * x)
    assert root_multiplicity_gcd(t, zero.as_poly()) == 1
    assert order_at(t, zero, -1) == 2
    assert all((zero ** k - 1).at(-1) == -1 for k in range(1, 9))
    negatives.append(rejected("treat zero as a torsion unit",
                              lambda: equal((zero ** 4 - 1).at(-1), 0)))
    for m in (2, 3, 4):
        r = (1 + x) ** m
        assert root_multiplicity_gcd(t, r.as_poly()) == m
        u = Fraction(2, m) * (r - 1) / (r + 1)
        j = -x / (1 + x)
        assert u.compose(j) == -u and (1 + x).compose(j) * (1 + x) == t.RONE
    negatives.append(rejected("a proper-power quotient claimed canonical",
                              lambda: equal(root_multiplicity_gcd(t, r.as_poly()), 1)))
    negatives.append(rejected("use an unramified preimage for the pole bound",
                              lambda: equal(2 <= 1 + 2 // 3, True)))
    negatives.append(rejected("distinct 2-adic orders claimed compatible",
                              lambda: equal(any(k % 2 == 1 and k % 4 == 2 for k in range(4)), True)))
    return counts, samples, negatives


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    ledger = json.loads((repo / LEDGER).read_text())
    for source in ledger["sources"]:
        assert digest(repo / source["path"]) == source["sha256"], source["path"]
    local = local_and_cyclotomic_controls()
    counts, samples, negatives = algebra_controls(repo)
    result = {
        "status": "PASS",
        "scope": "Bounded exact ramification, cyclic-jet, pole-twist, construction and moment controls. The written proof supplies the global classification and exhaustion.",
        "arithmetic": "Exact integers and rational arithmetic in Q, Q(i), and Q(sqrt(2)); standard library only; active assertions.",
        "python": platform.python_version(),
        "proof_sha256": digest(repo / PROOF),
        "checker_sha256": digest(Path(__file__)),
        "source_ledger_sha256": digest(repo / LEDGER),
        "source_pins_checked": len(ledger["sources"]),
        "local_and_cyclotomic_controls": local,
        "algebraic_controls": counts,
        "samples": samples,
        "negative_controls": negatives,
        "not_claimed": ["Finite tests proving the universal theorem", "A complete degree atlas of coordinates",
                        "A classification of all finite fibres", "Formal verification or external review"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS", json.dumps({"local": local, "algebra": counts,
                              "negative_controls": len(negatives)}, sort_keys=True))


if __name__ == "__main__":
    main()
