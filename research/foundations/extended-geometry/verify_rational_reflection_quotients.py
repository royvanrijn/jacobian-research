#!/usr/bin/env python3
"""Exact rational reflection extractor and bounded fibre controls.

The Extractor class implements the finite linear-space and common-fibre-gcd
construction over Q or a supplied number field.  The receipt covers only the
listed exact inputs.  Universal extraction, completeness and the cardinality
bound are established by the separate written proof, not by these controls.
"""

from __future__ import annotations

if not __debug__:
    raise RuntimeError("Assertions must be enabled.")

import argparse
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import platform
import time

import sympy as sp
from sympy.polys.matrices import DomainMatrix


x, y, T, z = sp.symbols("x y T z")
PREFIX = Path("research/foundations/extended-geometry")
PROOF = PREFIX / "RATIONAL_REFLECTION_QUOTIENT_AND_FINITE_FIBRES.md"
LEDGER = PREFIX / "rational_reflection_quotients_source_ledger.json"


class Extractor:
    """Exact extraction over QQ or QQ(radical), with no preselected quotient.

    ``radical`` may be a presented primitive algebraic number.  SymPy supplies
    exact polynomial factorization and field linear algebra; no series cutoff,
    floating point decision or numerical root selection is used.
    """

    def __init__(self, radical=None):
        self.radical = radical
        self.field = sp.QQ if radical is None else sp.QQ.algebraic_field(radical)
        self.rational = self.field.frac_field(x)

    def cancel(self, value):
        value = sp.sympify(value)
        variables = sorted(value.free_symbols, key=str)
        if not variables:
            return self.field.to_sympy(self.field.from_sympy(value))
        field = self.field.frac_field(*variables)
        return field.from_sympy(value).as_expr()

    def same(self, left, right):
        return self.cancel(left - right) == 0

    def fraction(self, value):
        return self.cancel(value).as_numer_denom()

    def degree(self, value):
        numerator, denominator = self.fraction(value)
        return max(sp.degree(numerator, x), sp.degree(denominator, x))

    def nullspace(self, columns, variables):
        polynomials = [sp.Poly(value, *variables, domain=self.field)
                       for value in columns]
        monomials = sorted(set().union(*(set(p.monoms()) for p in polynomials)))
        if not monomials:
            return sp.eye(len(columns))
        data = [[p.coeff_monomial(monomial) for p in polynomials]
                for monomial in monomials]
        matrix = DomainMatrix.from_list_sympy(
            len(data), len(columns), data).convert_to(self.field)
        return matrix.nullspace().to_Matrix()

    def marked(self, numerator, denominator):
        correspondence = sp.Poly(
            numerator*denominator.subs(x, y)
            + numerator.subs(x, y)*denominator,
            x, y, domain=self.field)
        factors = correspondence.factor_list()[1]
        choices = [factor for factor, _ in factors
                   if factor.eval({x: 0, y: 0}) == 0]
        assert len(choices) == 1
        marked = choices[0].monic()
        assert marked.diff(x).eval({x: 0, y: 0}) != 0
        assert marked.diff(y).eval({x: 0, y: 0}) != 0
        return marked

    def remainder(self, numerator, marked):
        return sp.reduced(
            sp.expand(numerator), [marked.as_expr()],
            x, y, domain=self.field)[1]

    def anti(self, value, marked):
        numerator, denominator = self.fraction(value + value.subs(x, y))
        assert self.remainder(denominator, marked) != 0
        return self.remainder(numerator, marked) == 0

    def normalize_generator(self, generator, coordinate):
        numerator, denominator = self.fraction(generator)
        if denominator.subs(x, 0) == 0:
            w = self.cancel(1/generator)
        else:
            w = self.cancel(generator - numerator.subs(x, 0)/denominator.subs(x, 0))
        derivative = self.cancel(sp.diff(w, x)).subs(x, 0)
        assert derivative != 0
        w = self.cancel(w/derivative)
        assert self.same(w.subs(x, 0), 0)
        assert self.same(sp.diff(w, x).subs(x, 0), 1)
        b = self.cancel((sp.diff(w, x, 2).subs(x, 0)
                         - sp.diff(coordinate, x, 2).subs(x, 0))/2)
        h = self.cancel(w/(1 + b*w))
        return h, w, b

    def outer(self, coordinate, h):
        numerator, denominator = self.fraction(coordinate)
        hn, hd = self.fraction(h)
        n, g = int(self.degree(coordinate)), int(self.degree(h))
        assert n % g == 0
        bound = n//g
        powers = [sp.expand(hn**i*hd**(bound-i)) for i in range(bound+1)]
        columns = [sp.expand(denominator*p) for p in powers]
        columns += [sp.expand(-numerator*p) for p in powers]
        kernel = self.nullspace(columns, (x,))
        assert kernel.rows == 1
        coefficients = list(kernel.row(0))
        a = sum(coefficients[i]*z**i for i in range(bound+1))
        b = sum(coefficients[bound+1+i]*z**i for i in range(bound+1))
        assert b != 0
        result = self.cancel(a/b)
        assert self.same(result.subs(z, h), coordinate)
        assert self.same(result.subs(z, -z), -result)
        assert self.same(result.subs(z, 0), 0)
        assert self.same(sp.diff(result, z).subs(z, 0), 1)
        return result

    def extract(self, coordinate):
        coordinate = self.cancel(coordinate)
        assert self.same(coordinate.subs(x, 0), 0)
        assert self.same(sp.diff(coordinate, x).subs(x, 0), 1)
        numerator, denominator = self.fraction(coordinate)
        n = int(self.degree(coordinate))
        marked = self.marked(numerator, denominator)
        common = sp.expand(numerator*denominator)
        common_y = common.subs(x, y)
        columns = [
            self.remainder(x**i*common_y + y**i*common, marked)
            for i in range(2*n+1)
        ]
        kernel = self.nullspace(columns, (x, y))
        assert kernel.rows > 0
        basis = []
        for row in range(kernel.rows):
            polynomial = sum(kernel[row, i]*x**i for i in range(2*n+1))
            value = self.cancel(polynomial/common)
            assert self.anti(value, marked)
            basis.append(value)

        minimal = None
        for value in basis:
            a, b = self.fraction(value)
            relation = sp.Poly(
                a.subs(x, T) - value*b.subs(x, T),
                T, domain=self.rational)
            if relation.is_zero:
                continue
            minimal = relation.monic() if minimal is None else minimal.gcd(relation).monic()
        assert minimal is not None and minimal.degree() > 0
        assert self.same(minimal.as_expr().subs(T, x), 0)
        generators = [coefficient for coefficient in minimal.all_coeffs()[1:]
                      if not self.same(sp.diff(coefficient, x), 0)]
        assert generators
        h, w, b = self.normalize_generator(generators[0], coordinate)
        assert self.anti(h, marked)
        assert self.same(h.subs(x, 0), 0)
        assert self.same(sp.diff(h, x).subs(x, 0), 1)
        assert self.degree(h) == minimal.degree()
        outer = self.outer(coordinate, h)
        # Every nonconstant minimal-polynomial coefficient must give the same
        # canonical normalized generator, not merely the same degree.
        for generator in generators[1:]:
            other, _, _ = self.normalize_generator(generator, coordinate)
            assert self.same(other, h)
        return {
            "coordinate": coordinate, "N": n, "marked_component": marked.as_expr(),
            "linear_space_dimension": kernel.rows, "basis": basis,
            "minimal_polynomial": minimal.as_expr(), "h": h,
            "h_degree": int(self.degree(h)), "w": w, "mobius_b": b,
            "outer": outer,
        }


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def polynomial(extractor, value):
    _, denominator = extractor.fraction(value)
    return not denominator.has(x)


def polynomial_direction(extractor, h):
    """Apply Q21 directly, including the entire-fibre condition."""
    a, b = extractor.fraction(h)
    ap, bp = (sp.Poly(value, x, domain=extractor.field) for value in (a, b))
    if ap.degree() != bp.degree():
        return None
    c = extractor.cancel(ap.LC()/bp.LC())
    difference = extractor.cancel(a-c*b)
    if c == 0 or difference == 0 or difference.has(x):
        return None
    result = extractor.cancel((c+h)/(c-h))
    assert polynomial(extractor, result)
    assert extractor.same(result.subs(x, 0), 1)
    return c, result


def jet_holds(extractor, eta, modulus):
    numerator, denominator = extractor.fraction(eta)
    modulus = sp.Poly(modulus, x, domain=extractor.field)
    if modulus.degree() == 0:
        return True
    denominator = sp.Poly(denominator, x, domain=extractor.field)
    if denominator.gcd(modulus).degree() != 0:
        return False
    return sp.Poly(numerator+denominator.as_expr(), x,
                   domain=extractor.field).rem(modulus).is_zero


def construct(extractor, u, eta, height):
    assert height > 0 and height % 2
    xx = extractor.cancel(u*(1+eta)/(2*x*sp.diff(u, x)))
    c = extractor.cancel(x*xx**height/u)
    d = extractor.cancel(-x*xx**(height+2)/u)
    return xx, c, d


def multiply_truncated(a, b, bound):
    result = [sp.S.Zero]*(bound+1)
    for i, ai in enumerate(a[:bound+1]):
        for j, bj in enumerate(b[:bound+1-i]):
            result[i+j] += ai*bj
    return result


def power_truncated(a, power, bound):
    result = [sp.S.One]+[sp.S.Zero]*bound
    for _ in range(power):
        result = multiply_truncated(result, a, bound)
    return result


def raw_moment(extractor, c, d, height, moment):
    """Direct coefficient/integral formula for this finite positive moment."""
    assert polynomial(extractor, c) and polynomial(extractor, d)
    cp, dp = (sp.Poly(value, x, domain=extractor.field) for value in (c, d))
    ca = [cp.nth(i) for i in range(moment+1)]
    da = [dp.nth(i) for i in range(moment+1)]
    value = 0
    for j in range(moment+1):
        term = multiply_truncated(power_truncated(ca, moment-j, moment),
                                  power_truncated(da, j, moment), moment)
        value += sp.Rational(comb(moment, j), height*moment+2*j+1)*term[moment]
    return extractor.cancel(value)


def order_at(extractor, value, point):
    numerator, denominator = extractor.fraction(value)
    orders = []
    for expression in (numerator, denominator):
        p = sp.Poly(expression, x, domain=extractor.field)
        linear = sp.Poly(x-point, x, domain=extractor.field)
        multiplicity = 0
        while not p.is_zero:
            q, r = p.div(linear)
            if not r.is_zero:
                break
            multiplicity += 1
            p = q
        assert not p.is_zero
        orders.append(multiplicity)
    return orders[0]-orders[1]


def rejected(name, action):
    try:
        action()
    except (AssertionError, ValueError, ZeroDivisionError):
        return name
    raise AssertionError("Negative control accepted: " + name)


def require(value):
    assert value


def examples():
    quadratic = 1+x+x**2/8
    thick_r = sp.I+(1-sp.I)*(1+x)**3
    thick_u = 2/(3*(1-sp.I))*(thick_r-1)/(thick_r+1)
    dz = x*(1-x)**2/(1-sp.Rational(3, 2)*x)**3
    return [
        ("linear", None, x, x, None),
        ("rational_odd", None, x/(1-x**2), x, None),
        ("polynomial_nonrational_return", None, x+x**2, x+x**2, None),
        ("weighted_DZ_singleton", None, dz, dz, None),
        ("mobius_polynomial_generator", None, x/(1+x), x/(1+x), 1+2*x),
        ("quadratic_nonrational_return", None,
         2*(quadratic-1)/(quadratic+1), 2*(quadratic-1)/(quadratic+1), quadratic),
        ("Qi_thick_jet", sp.I, thick_u, thick_u, thick_r),
    ]


def extraction_controls():
    records, outputs = [], {}
    for name, radical, coordinate, expected, expected_r in examples():
        extractor = Extractor(radical)
        result = extractor.extract(coordinate)
        assert extractor.same(result["h"], expected), name
        direction = polynomial_direction(extractor, result["h"])
        if expected_r is None:
            assert direction is None, name
        else:
            assert direction is not None, name
            assert extractor.same(direction[1], expected_r), name
            assert extractor.anti((direction[1]-1)/(direction[1]+1),
                                  sp.Poly(result["marked_component"], x, y,
                                          domain=extractor.field))
        record = {key: value if isinstance(value, int) else str(value)
                  for key, value in result.items() if key != "basis"}
        record["basis"] = [str(value) for value in result["basis"]]
        record["name"] = name
        record["coefficient_field"] = str(extractor.field)
        record["polynomial_direction"] = None if direction is None else {
            "c": str(direction[0]), "R": str(direction[1])}
        records.append(record)
        outputs[name] = (extractor, result)
        print("extracted", name, "N="+str(result["N"]),
              "degree_h="+str(result["h_degree"]), flush=True)
    return records, outputs


def fibre_controls(outputs):
    counts = {"finite_pole_box_candidates": 0, "polynomial_pairs": 0,
              "raw_even_moments": 0, "nonpolynomial_reconstructions": 0,
              "thick_jet_tests": 0}
    negative, samples = [], []

    def valid_pair(extractor, u, eta, height):
        xx, c, d = construct(extractor, u, eta, height)
        assert polynomial(extractor, c) and polynomial(extractor, d)
        assert extractor.same(c.subs(x, 0), 1)
        assert extractor.same(d.subs(x, 0), -1)
        assert extractor.same(xx.subs(x, 0), 1)
        assert extractor.same(x*xx**height/c, u)
        assert extractor.same(-d/c, xx**2)
        for moment in (2, 4, 6):
            assert raw_moment(extractor, c, d, height, moment) == 0
            counts["raw_even_moments"] += 1
        counts["polynomial_pairs"] += 1
        return {"height": height, "eta": str(eta), "X": str(xx),
                "C": str(c), "D": str(d)}

    ex, out = outputs["linear"]
    for height in (1, 3, 5):
        sample = valid_pair(ex, x, sp.S.One, height)
        assert ex.same(sp.sympify(sample["C"]), 1)
        assert ex.same(sp.sympify(sample["D"]), -1)
    samples.append({"name": "linear_complete_singleton", "heights": [1, 3, 5]})

    ex, out = outputs["rational_odd"]
    u, h = out["coordinate"], out["h"]
    marked = sp.Poly(out["marked_component"], x, y, domain=ex.field)
    assert ex.same(marked.as_expr(), x+y)
    numerator, denominator = ex.fraction(u)
    full = sp.Poly(numerator*denominator.subs(x, y)
                   + numerator.subs(x, y)*denominator, x, y, domain=ex.field)
    assert ex.remainder(x+y, marked) == 0
    assert ex.remainder(x+y, full) != 0
    negative.append(rejected("replace the marked component by the full correspondence",
                             lambda: require(ex.remainder(x+y, full) == 0)))
    negative.append(rejected("normalize by an arbitrary Mobius transform",
                             lambda: require(ex.anti(h/(1+h), marked))))
    assert ex.anti(h+h**3, marked)
    assert sp.diff(h+h**3, x).subs(x, 0) == 1
    negative.append(rejected("omit the field-generator requirement",
                             lambda: require(ex.degree(h+h**3) == out["h_degree"])))
    outer_x = out["outer"].subs(z, x)
    caps = {}
    for value in (-1, 1):
        pole_order = -order_at(ex, outer_x, value)
        local_degree = order_at(ex, h-value, value)
        assert pole_order == local_degree == 1
        caps[value] = int(sp.floor(sp.Rational(pole_order, 3)
                                   +sp.Rational(1, local_degree)))
        assert caps[value] == 1
    box = []
    for exponent in range(-caps[1], caps[-1]+1):
        eta = ex.cancel(((1-x)/(1+x))**exponent)
        assert ex.same(eta*eta.subs(x, -x), 1)
        assert order_at(ex, eta, 1) == exponent
        assert order_at(ex, eta, -1) == -exponent
        assert not jet_holds(ex, eta, 1+x*x)
        _, c, d = construct(ex, u, eta, 1)
        assert not (polynomial(ex, c) and polynomial(ex, d))
        box.append({"exponent": exponent, "eta": str(eta), "jet_passed": False})
        counts["finite_pole_box_candidates"] += 1
        counts["nonpolynomial_reconstructions"] += 1
    samples.append({"name": "empty_three_candidate_fibre", "candidates": box})
    negative.append(rejected("omit a quotient-denominator unit condition",
                             lambda: require(jet_holds(ex, -1+1/(1+x*x), 1+x*x))))

    ex, out = outputs["polynomial_nonrational_return"]
    u = out["coordinate"]
    assert not jet_holds(ex, sp.S.One, 1+2*x)
    _, c, d = construct(ex, u, sp.S.One, 1)
    assert not (polynomial(ex, c) and polynomial(ex, d))
    counts["nonpolynomial_reconstructions"] += 1
    negative.append(rejected("ignore an ordinary critical jet",
                             lambda: require(polynomial(ex, c) and polynomial(ex, d))))
    # The marked quadratic is not linear in y; the original return is not a
    # rational function of x even though its canonical reflection quotient is.
    marked_y = sp.Poly(out["marked_component"], y, domain=sp.QQ.frac_field(x))
    assert marked_y.degree() == 2
    discriminant = sp.Poly(sp.discriminant(marked_y.as_expr(), y), x)
    assert discriminant.gcd(discriminant.diff()).degree() == 0
    assert discriminant.degree() == 2

    ex, out = outputs["weighted_DZ_singleton"]
    u = out["coordinate"]
    for height in (3, 5, 7):
        samples.append({"name": "weighted_DZ_singleton",
                        **valid_pair(ex, u, sp.S.One, height)})
    _, c, d = construct(ex, u, sp.S.One, 1)
    assert not (polynomial(ex, c) and polynomial(ex, d))
    counts["nonpolynomial_reconstructions"] += 1
    negative.append(rejected("accept the weighted singleton at height one",
                             lambda: require(polynomial(ex, c) and polynomial(ex, d))))

    ex, out = outputs["mobius_polynomial_generator"]
    u = out["coordinate"]
    c, r = polynomial_direction(ex, out["h"])
    assert c == 1 and ex.same(r, 1+2*x)
    hn, hd = ex.fraction(out["h"])
    assert not ex.cancel(hn-c*hd).has(x)
    for exponent in range(5):
        samples.append({"name": "degree_one_polynomial_direction",
                        **valid_pair(ex, u, r**exponent, 1)})

    ex, out = outputs["quadratic_nonrational_return"]
    r = polynomial_direction(ex, out["h"])[1]
    samples.append({"name": "quadratic_polynomial_direction",
                    **valid_pair(ex, out["coordinate"], r**2, 1)})

    ex, out = outputs["Qi_thick_jet"]
    r = polynomial_direction(ex, out["h"])[1]
    assert order_at(ex, r-sp.I, -1) == 3
    assert order_at(ex, sp.diff(out["coordinate"], x), -1) == 2
    for exponent in range(8):
        assert jet_holds(ex, r**exponent, (x+1)**2) == (exponent % 4 == 2)
        counts["thick_jet_tests"] += 1
    assert ex.same((r**2).subs(x, -1), -1)

    # Maximum finite ramification, as opposed to a generic degree-one fibre,
    # controls the allowed quotient pole order.
    ex = Extractor()
    h = x+x*x
    u = h/(1-16*h*h)
    eta = (1-4*h)/(1+4*h)
    assert order_at(ex, h+sp.Rational(1, 4), -sp.Rational(1, 2)) == 2
    assert order_at(ex, u, -sp.Rational(1, 2)) == -2
    assert order_at(ex, eta, -sp.Rational(1, 2)) == -2
    allowed_source_pole = 1+2//3
    assert allowed_source_pole == 1
    assert sp.floor(sp.Rational(1, 3)+sp.Rational(1, 2)) == 0
    negative.append(rejected("use unramified allowance at a ramified preimage",
                             lambda: require(2 <= allowed_source_pole)))
    samples.append({"name": "maximum_finite_ramification",
                    "h": str(h), "pole_value": "-1/4", "local_degree": 2,
                    "height": 1, "outer_pole_order": 1, "quotient_pole_cap": 0})
    return counts, samples, negative


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path,
                        default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    ledger = json.loads((repo / LEDGER).read_text())
    for source in ledger["sources"]:
        assert digest(repo / source["path"]) == source["sha256"], source["path"]
    assert ledger["proof_sha256"] == digest(repo / PROOF)
    assert ledger["checker_sha256"] == digest(Path(__file__))
    started = time.monotonic()
    records, outputs = extraction_controls()
    counts, samples, negative = fibre_controls(outputs)
    receipt = {
        "status": "PASS",
        "scope": "Implemented exact quotient extraction on seven specified inputs, and bounded exact pole-box, polynomial reconstruction, full-jet and moment controls. Universal extraction, exhaustion and the finite-fibre bound rely on the separate written proof.",
        "arithmetic": "Exact rational and Q(i) polynomial/rational-function arithmetic; SymPy factorization, field linear algebra and univariate gcd; active assertions.",
        "python": platform.python_version(), "sympy": sp.__version__,
        "proof_sha256": digest(repo / PROOF),
        "checker_sha256": digest(Path(__file__)),
        "source_ledger_sha256": digest(repo / LEDGER),
        "source_pins_checked": len(ledger["sources"]),
        "extractions": records, "fibre_control_counts": counts,
        "samples": samples, "negative_controls": negative,
        "omitted_control": "The optional degree-four Q(sqrt(2)) twist-family extraction exceeded a 90-second exploratory cap and is not included in this receipt. No output or completion is claimed for that run.",
        "not_claimed": ["A full arbitrary-number-field fibre classifier implementation",
                        "A runtime bound for number-field factorization",
                        "Finite controls proving the unrestricted classification",
                        "Formal verification or external human review"],
        "elapsed_seconds": round(time.monotonic()-started, 6),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2)+"\n")
    print("PASS", json.dumps({"extractions": len(records), **counts,
                              "negative_controls": len(negative)}, sort_keys=True))


if __name__ == "__main__":
    main()
