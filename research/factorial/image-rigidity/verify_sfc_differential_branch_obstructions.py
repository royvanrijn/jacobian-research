#!/usr/bin/env python3
"""Exact two-parameter obstruction controls; no universal SFC claim.

Standard library only. The pinned retained connection checker certifies
minimal differential operators over each exact quadratic field. Separate
inverse-branch expansions certify that the residue sequences are actual
algebraic conjugates. The all-order implications are proved in the note.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from math import comb
from pathlib import Path

if not __debug__:
    raise SystemExit("REJECT: optimized execution is not an approved replay")

HERE = Path(__file__).resolve().parent
SOURCE = HERE
NOTE = HERE / "SFC_DIFFERENTIAL_BRANCH_OBSTRUCTIONS.md"
NOTE_PIN = "f1a3ebedcc15a90443ed196b4b8e798341720a17f6b581b0da97a609a949ea2b"
PINS = {
    "verify_sfc_monomial_linear_contact.py":
        "c80ea0aa4ade828a89d18e028072c65fc0525c6700afc661183b5cb14c2e1055",
    "verify_sfc_uniform_rational_window.py":
        "384f05608a1047cdaa87bdad5647fde67b3e265877e1d1dcbc45501533682829",
    "SFC_MONOMIAL_LINEAR_CONTACT_BOUND.md":
        "409ca1ed8a1ef3218fb07673d14830d04cf10562f5714c88bd8de8e9a681b417",
}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def field(q0, q1, q2):
    """Exact Q[u]/(q2*u^2+q1*u+q0); both inputs are irreducible."""
    aa, bb = Fraction(-q0, q2), Fraction(-q1, q2)
    need(q1*q1-4*q0*q2 < 0, "quadratic must have nonreal roots")

    class E:
        def __init__(self, a=0, b=None):
            if isinstance(a, E):
                self.a, self.b = a.a, a.b
            else:
                self.a, self.b = Fraction(a), Fraction(0)
            if b is not None:
                value = self / E(b)
                self.a, self.b = value.a, value.b

        @classmethod
        def pair(cls, a, b):
            value = cls()
            value.a, value.b = Fraction(a), Fraction(b)
            return value

        def __add__(self, other):
            other = E(other)
            return E.pair(self.a+other.a, self.b+other.b)
        __radd__ = __add__

        def __neg__(self):
            return E.pair(-self.a, -self.b)

        def __sub__(self, other):
            return self + -E(other)

        def __rsub__(self, other):
            return E(other) + -self

        def __mul__(self, other):
            other = E(other)
            return E.pair(self.a*other.a+aa*self.b*other.b,
                          self.a*other.b+self.b*other.a+bb*self.b*other.b)
        __rmul__ = __mul__

        def __truediv__(self, other):
            other = E(other)
            norm = other.a*(other.a+bb*other.b)-aa*other.b*other.b
            need(norm != 0, "division by zero field element")
            return self * E.pair((other.a+bb*other.b)/norm, -other.b/norm)

        def __rtruediv__(self, other):
            return E(other)/self

        def __pow__(self, exponent):
            if exponent < 0:
                return (1/self)**(-exponent)
            value, base = E(1), self
            while exponent:
                if exponent & 1:
                    value = value*base
                base = base*base
                exponent //= 2
            return value

        def __eq__(self, other):
            try:
                other = E(other)
            except (ValueError, TypeError):
                return False
            return self.a == other.a and self.b == other.b

        def __bool__(self):
            return bool(self.a or self.b)

        def __repr__(self):
            return f"({self.a})+({self.b})u"

    return E


def serialize(value):
    if hasattr(value, "a") and hasattr(value, "b"):
        return [str(value.a), str(value.b)]
    if isinstance(value, (tuple, list)):
        return [serialize(v) for v in value]
    if isinstance(value, dict):
        return {str(k): serialize(v) for k, v in value.items()}
    if isinstance(value, Fraction):
        return str(value)
    return value


def run_case(modulus, gap_start, h, q):
    E = field(*modulus)
    h.F = q.F = E
    u = E.pair(0, 1)
    need(modulus[2]*u*u+modulus[1]*u+modulus[0] == 0, "field relation")
    need(u != 0 and u != 1 and u != 2, "active distinct nonsingular parameters")
    data = q.check_case((1, 2), (E(1), u), h)
    a, D = data["a"], data["D_phase"]

    # Full binomial residues through order 4, authenticated independently
    # by locally inverting the actual polynomial at both simple fibers.
    x, y = [], []
    for m in range(5):
        x.append(sum(E(comb(2*m-k, m-k)*comb(2*m+k, k))*u**k
                     for k in range(m+1)))
        y.append(sum(E((-1)**(k+1)*comb(2*m-k, m-k)*comb(2*m+k, k))*u**k
                     /(1-u)**(2*m+1+k) for k in range(m+1)))
    for alpha, expected in ((E(0), x), (E(1), y)):
        shifted = q.compose(a, [alpha, E(1)])
        need(shifted[0] == 0 and shifted[1] != 0, "simple zero fiber")
        normalized = q.scale(shifted, 1/shifted[1])
        inverse = h.inverse(normalized, 4)
        inverse = [E(c)/shifted[1]**j for j, c in enumerate(inverse)]
        need(h.compose(shifted, inverse, 4) == [0, 1], "actual local inverse")
        rho = h.reciprocal(h.compose(q.compose(D, [alpha, E(1)]), inverse, 4), 4)
        need(all(E(h.coeff(rho, m)) == expected[m] for m in range(5)),
             "residue equals actual conjugate coefficient")

    ell = [y[gap_start]*xx-x[gap_start]*yy for xx, yy in zip(x, y)]
    need(ell[0] != 0, "normalizable residue period")
    normalized_ell = [v/ell[0] for v in ell]
    need(ell[gap_start] == ell[gap_start+1] == 0, "exact forbidden period window")
    need(ell[gap_start+2] != 0, "nonzero post-gap return")
    if gap_start == 1:
        need(ell[0] == E(Fraction(158, 625))-Fraction(14, 125)*u,
             "first unnormalized constant")
        need(normalized_ell[3] == E(Fraction(-3984, 625))+Fraction(9360, 625)*u,
             "first exact return")
    else:
        need(ell[1] == -5*u**3*(u-2)/(4*(u-1)**6), "second linear term")
        need(normalized_ell[1] == E(Fraction(25, 101))+Fraction(300, 707)*u,
             "second normalized linear term")
        need(normalized_ell[4] == E(Fraction(-5000125, 349056))
             + Fraction(2253875, 203616)*u, "second exact return")

    # Exact trace-H computation: Newton sums of z and z^2 do not depend
    # on the constant term -t of a(z)-t. This is a finite control for the
    # all-degree trace argument in the note, not a replacement for it.
    monic = q.scale(a, 1/a[-1])
    s1 = -monic[-2]
    s2 = -monic[-2]*s1-2*monic[-3]
    trace_h = s1-u*s2/2
    need(trace_h == 1+1/u-u/2, "constant trace of actual H")

    def origin(poly):
        return next(i for i, c in enumerate(poly) if c)

    def primitive(operator):
        common = next(poly for poly in operator if poly != [0])
        for poly in operator:
            common = q.gcd(common, poly)
        reduced = [q.exact(poly, common) for poly in operator]
        scale = 1/reduced[-1][-1]
        return [q.scale(poly, scale) for poly in reduced], common

    def euler(operator):
        lower = min(origin(poly)-i for i, poly in enumerate(operator) if poly != [0])
        upper = max(q.deg(poly)-i for i, poly in enumerate(operator) if poly != [0])
        falling = [[E(1)]]
        for i in range(1, len(operator)):
            falling.append(q.mul(falling[-1], [E(1-i), E(1)]))
        polys = []
        for j in range(upper-lower+1):
            value = [E(0)]
            for i, poly in enumerate(operator):
                index = j+i+lower
                if 0 <= index < len(poly):
                    value = q.add(value, q.scale(falling[i], poly[index]))
            polys.append(value)
        return polys, lower, upper

    reduced, common = primitive(data["operator"])
    ep, lo, hi = euler(reduced)
    need(data["minimal_order_q"] == 3 and hi-lo == 3, "actual minimal order and width")
    need(q.gcd(reduced[0], reduced[1]) == [1], "primitive low coefficients")
    # Add precisely the period truncation. Its image under the actual
    # minimal operator is a pure power, so the displayed left multiple
    # is an exact annihilator without recurrence guessing.
    truncation = normalized_ell[:gap_start]
    image = [E(0)]
    derivative = truncation
    for coefficient in reduced:
        image = q.add(image, q.mul(coefficient, derivative))
        derivative = q.der(derivative)
    need(q.deg(image) == gap_start and origin(image) == gap_start,
         "minimal operator image of truncation is a pure power")
    # t*(partial-gap_start/t)*L, followed by exact scalar-gcd removal.
    enlarged = []
    for i in range(len(reduced)+1):
        current = reduced[i] if i < len(reduced) else [E(0)]
        previous = reduced[i-1] if i else [E(0)]
        enlarged.append(q.add(q.sub(q.mul([0, 1], q.der(current)),
                                    q.scale(current, gap_start)),
                              q.mul([0, 1], previous)))
    enlarged, common2 = primitive(enlarged)
    ep2, lo2, hi2 = euler(enlarged)
    need(len(enlarged)-1 == 4 and hi2-lo2 == 2, "exact narrow left multiple")
    need(q.evaluate(ep2[0], E(gap_start+2)) == 0, "forced head resonance")
    head = [E(1)]
    for root in (0, 1, Fraction(1, 2), gap_start+2):
        head = q.mul(head, [E(-root), E(1)])
    scale = ((3-8*u)/128 if gap_start == 1 else (7-19*u)/320)
    need(ep2[0] == q.scale(head, scale), "complete narrow head factorization")
    top = [E(1)]
    for root in (gap_start-1, Fraction(-1, 2), Fraction(-3, 4), Fraction(-5, 4)):
        top = q.mul(top, [E(-root), E(1)])
    need(ep2[-1] == top, "complete polynomial-degree obstruction")
    # Full polynomial identity: the operator annihilates the displayed
    # truncation, including the nonconstant truncation in the second case.
    zero = [E(0)]
    for j, poly in enumerate(ep2):
        term = [E(c)*q.evaluate(poly, E(i)) for i, c in enumerate(truncation)]
        zero = q.add(zero, [E(0)]*j+term)
    need(zero == [0], "exact polynomial-kernel witness")

    return serialize({
        "modulus_ascending": modulus,
        "profile": [1, 2], "coefficients": [E(1), u],
        "residue_x_orders_0_to_4": x,
        "residue_y_orders_0_to_4": y,
        "normalized_period_orders_0_to_4": normalized_ell,
        "gap_start": gap_start, "forced_head_root": gap_start+2,
        "polynomial_kernel_witness": truncation,
        "trace_H": trace_h,
        "minimal_order": 3, "minimal_euler_width": hi-lo,
        "raw_common_factor_degree": q.deg(common),
        "raw_common_factor_origin_order": origin(common),
        "minimal_operator": reduced, "minimal_euler": ep,
        "truncation_image": image,
        "left_multiple_common_factor": common2,
        "left_multiple_common_factor_origin_order": origin(common2),
        "narrow_order": 4, "narrow_euler_width": hi2-lo2,
        "narrow_operator": enlarged, "narrow_euler": ep2,
        "independent_derivative_rows": data["independent_rows"],
        "all_coordinate_cofactor_relation": "PASS",
        "actual_factorial_orders_verified": list(range(9)),
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", type=Path, default=NOTE)
    parser.add_argument("--source-dir", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path,
                        default=HERE/"sfc_differential_branch_obstructions_verification.json")
    args = parser.parse_args()
    need(args.note.is_file(), "missing note")
    note_hash = sha256(args.note.read_bytes()).hexdigest()
    need(note_hash == NOTE_PIN, "note SHA-256 mismatch")
    for name, digest in PINS.items():
        path = args.source_dir/name
        need(path.is_file(), f"missing pinned source: {name}")
        need(sha256(path.read_bytes()).hexdigest() == digest,
             f"source SHA-256 mismatch: {name}")
    h = load("euler_contact_helper", args.source_dir/"verify_sfc_monomial_linear_contact.py")
    q = load("euler_uniform_helper", args.source_dir/"verify_sfc_uniform_rational_window.py")
    results = [run_case((2, -2, 5), 1, h, q), run_case((5, -5, 12), 2, h, q)]
    receipt = {
        "status": "PASS", "note_sha256": note_hash, "source_sha256": PINS,
        "exact_cases": 2, "residue_branch_coefficients_compared": 20,
        "actual_factorial_coefficients_compared": 18,
        "scope": "Two exact algebraic-parameter controls, complete cofactor identities, and polynomial witnesses; universal implications are proved in the note.",
        "claims_not_made": ["SFC counterexample", "full SFC closure", "universal recurrence-width bound"],
        "cases": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True)+"\n")
    print("PASS: 2 exact residue countercontrols; 20 branch coefficients; 18 factorial coefficients; 2 minimal operators and 2 narrow left multiples.")


if __name__ == "__main__":
    main()
