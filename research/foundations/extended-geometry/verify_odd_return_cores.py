#!/usr/bin/env python3
"""Exact controls for odd-return descent and the bounded Cayley diagnostic.

The unrestricted reduction is a written theorem. The finite computation
covers exactly the stated six values of n and all reciprocal Q-divisors;
for each divisor, exact polynomial factorization covers every integer m.
It does not exclude arbitrary odd-return cores or nonrational divisors.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("Assertions are required; Python optimization is unsupported.")

import argparse
import hashlib
import itertools
import json
import platform
from pathlib import Path

import sympy as sp

x, z, r, m = sp.symbols("x z r m")
BASE = Path("research/foundations/extended-geometry")
PROOF = BASE / "ODD_RETURN_CORE_AND_CAYLEY_OBSTRUCTION.md"
CHECKER = BASE / "verify_odd_return_cores.py"
LEDGER = BASE / "odd_return_cores_source_ledger.json"
N_VALUES = (12, 15, 18, 20, 24, 30)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cancel(expr):
    return sp.cancel(sp.together(expr))


def polynomial(expr, variable=x):
    numerator, denominator = sp.fraction(cancel(expr))
    sp.Poly(numerator, variable, domain=sp.QQ)
    sp.Poly(denominator, variable, domain=sp.QQ)
    return not denominator.has(variable)


def order_at(expr, point, variable=x):
    numerator, denominator = sp.fraction(cancel(expr))

    def multiplicity(poly):
        poly = sp.Poly(poly, variable, domain=sp.QQ)
        assert not poly.is_zero
        divisor = sp.Poly(variable-point, variable, domain=sp.QQ)
        count = 0
        while poly.eval(point) == 0:
            poly, remainder = sp.div(poly, divisor)
            assert remainder.is_zero
            count += 1
        return count

    return multiplicity(numerator)-multiplicity(denominator)


def raw_moment(C, phase_D, height, exponent):
    if exponent % 2:
        return sp.S.Zero
    return cancel(sum(
        sp.binomial(exponent, j)
        * sp.Poly(sp.expand(C**(exponent-j)*phase_D**j), x).nth(exponent)
        / (height*exponent+2*j+1)
        for j in range(exponent+1)))


def pullback(Cbar, Dbar, h, height):
    derivative = sp.diff(h, x)
    C = cancel(Cbar.subs(z, h)*(h/x)**(height-1)/derivative**height)
    phase_D = cancel(Dbar.subs(z, h)*(h/x)**(height+1)/derivative**(height+2))
    return C, phase_D


def cleared_middle(Afactor):
    poly = sp.Poly(Afactor, r, domain=sp.QQ)
    degree = poly.degree()
    assert degree % 2 == 0 and poly.LC() == poly.TC() == 1
    assert sp.expand(r**degree*Afactor.subs(r, 1/r)-Afactor) == 0
    half = degree//2
    answer = sp.S.Zero
    for j in range(degree+1):
        distance = abs(j-half)
        falling = sp.prod(m-half-i for i in range(distance))
        remainder = sp.prod(m-half+i for i in range(distance+1, half+1))
        answer += poly.nth(j)*falling*remainder
    answer = sp.Poly(sp.expand(answer), m, domain=sp.QQ)
    assert not answer.is_zero and answer.degree() <= half
    return answer, half


def literal_middle(Afactor, pole_order):
    half = int(sp.degree(Afactor, r))//2
    power = 2*pole_order-2*half
    if pole_order < 1 or power < 0:
        raise ValueError("Positive pole order and nonnegative Cayley power are required.")
    # Literal coefficient extraction, independent of the factorial-ratio code.
    return sp.Poly(sp.expand((r+1)**power*Afactor), r).nth(pole_order)


def binomial_middle(Afactor, pole_order):
    half = int(sp.degree(Afactor, r))//2
    power = 2*pole_order-2*half
    if pole_order < 1 or power < 0:
        raise ValueError("Positive pole order and nonnegative Cayley power are required.")
    poly = sp.Poly(Afactor, r)
    return sum(poly.nth(j)*sp.binomial(power, pole_order-j)
               for j in range(2*half+1))


def exact_roots(poly):
    """Full QQ factorization, with its identity and irreducibility checked."""
    content, factors = sp.factor_list(poly)
    rebuilt = sp.Poly(content, m, domain=sp.QQ)
    roots = {}
    for factor, multiplicity in factors:
        assert factor.is_irreducible
        rebuilt *= factor**multiplicity
        if factor.degree() == 1:
            root = -factor.nth(0)/factor.nth(1)
            assert poly.eval(root) == 0
            roots[root] = multiplicity
    assert rebuilt == poly
    # Cross-check the rational-root interface against the linear factors.
    assert roots == sp.polys.polytools.ground_roots(poly)
    return roots, {
        "content": str(content),
        "irreducible_factors": [
            {"polynomial": str(f.as_expr()), "multiplicity": exponent}
            for f, exponent in factors],
    }


def cayley_fixture(n, pole_order, Afactor):
    degree = int(sp.degree(Afactor, r))
    power = 2*pole_order-degree
    assert pole_order >= 1 and power >= 0
    assert sp.rem(r**n+1, Afactor, r) == 0 and Afactor.subs(r, -1) != 0
    assert sp.expand(r**degree*Afactor.subs(r, 1/r)-Afactor) == 0
    constant = sp.Rational(1, 2**(power+1))/Afactor.subs(r, 1)
    numerator = sp.Poly(sp.expand((r+1)**power*Afactor), r)
    residue = constant*numerator.nth(pole_order)
    derivative = constant*r**(-pole_order-1)*numerator.as_expr()
    assert cancel(derivative.subs(r, 1/r)-r**2*derivative) == 0
    assert derivative.subs(r, 1) == sp.Rational(1, 2)
    if residue != 0:
        return {"n": n, "m": pole_order, "Afactor": str(Afactor),
                "logarithmic_residue": str(residue), "integrable": False}
    primitive = sum(constant*coefficient*r**(j-pole_order)/(j-pole_order)
                    for (j,), coefficient in numerator.terms() if j != pole_order)
    primitive = cancel(primitive-primitive.subs(r, 1))
    assert cancel(sp.diff(primitive, r)-derivative) == 0
    assert cancel(primitive.subs(r, 1/r)+primitive) == 0
    cayley = (1+x)/(1-x)
    u = cancel(primitive.subs(r, cayley))
    eta = cayley**n
    X = cancel(u*(1+eta)/(2*x*sp.diff(u, x)))
    C = cancel(x*X/u)
    phase_D = cancel(-x*X**3/u)
    assert u.subs(x, 0) == 0 and sp.diff(u, x).subs(x, 0) == 1
    assert cancel(u.subs(x, -x)+u) == 0
    assert cancel(eta*eta.subs(x, -x)-1) == 0
    assert C.subs(x, 0) == 1 and phase_D.subs(x, 0) == -1
    assert cancel(sp.diff(u, x)-(1/C+1/C.subs(x, -x))/2) == 0
    observed = [order_at(C, 1), order_at(phase_D, 1),
                order_at(C, -1), order_at(phase_D, -1)]
    assert observed == [pole_order+1-n, pole_order-3*(n-1),
                        pole_order+1, pole_order+3]
    return {"n": n, "m": pole_order, "Afactor": str(Afactor),
            "logarithmic_residue": "0", "integrable": True,
            "u": str(u), "X": str(X), "C": str(C), "phase_D": str(phase_D),
            "C_polynomial": polynomial(C), "phase_D_polynomial": polynomial(phase_D),
            "orders_C_D_at_plus_and_minus_one": observed,
            "admissible_pole_bound": pole_order >= 3*(n-1),
            "pure_polynomial_core_claimed": False}


def repeated_pole_controls():
    records = []
    jet_controls = 0
    for pole_order, scale in itertools.product((1, 2, 7), (sp.S.One, sp.Integer(4))):
        alpha = sp.S.One if scale == 1 else sp.Rational(1, 2)
        u = x/(1-scale*x*x)**pole_order
        critical = sp.Poly(1+(2*pole_order-1)*scale*x*x, x, domain=sp.QQ)
        assert cancel(sp.diff(u, x)-critical.as_expr()/(1-scale*x*x)**(pole_order+1)) == 0
        phi = (alpha-x)/(alpha+x)
        assert cancel(phi*phi.subs(x, -x)-1) == 0

        def remainder(expr):
            numerator, denominator = sp.fraction(cancel(expr))
            assert sp.gcd(sp.Poly(denominator, x, domain=sp.QQ), critical).degree() == 0
            return sp.rem(sp.Poly(numerator, x, domain=sp.QQ), critical).as_expr()

        trace = sp.Rational(2*(pole_order-1), pole_order)
        assert remainder(phi+1/phi-trace) == 0
        if pole_order in (1, 2):
            root_order = 4 if pole_order == 1 else 6
            assert remainder(phi**root_order-1) == 0
            assert all(remainder(phi**j-1) != 0 for j in range(1, root_order))
        else:
            root_order = None
            assert trace.q > 1
        for height in (1, 3):
            cap = 1+pole_order//(height+2)
            remainders = []
            for exponent in range(-cap, cap+1):
                rem = remainder(phi**exponent+1)
                assert rem != 0
                remainders.append({"exponent": exponent, "nonzero_jet_numerator": str(rem)})
                jet_controls += 1
            records.append({"m": pole_order, "lambda": str(scale), "a": height,
                            "critical_polynomial": str(critical.as_expr()),
                            "critical_unit_trace": str(trace), "torsion_order_when_present": root_order,
                            "exact_exponent_cap": cap, "rejected_jet_values": remainders})
    assert len(records) == 12 and jet_controls == 48
    return records, jet_controls


def verify():
    negatives = []
    divisor_records = []
    binomial_controls = 0
    integer_root_count = 0
    clearing_exception_count = 0
    admissible = []

    for n in N_VALUES:
        factors = [f.as_expr() for f, exponent in
                   sp.factor_list(sp.Poly(r**n+1, r, domain=sp.QQ))[1]
                   if f.as_expr() != r+1]
        excluded = r+1 if n % 2 else sp.S.One
        assert sp.expand(excluded*sp.prod(factors)-(r**n+1)) == 0
        assert all(sp.degree(f, r) % 2 == 0 for f in factors)
        for selection in itertools.product((0, 1), repeat=len(factors)):
            Afactor = sp.expand(sp.prod(f for f, use in zip(factors, selection) if use))
            assert sp.rem(r**n+1, Afactor, r) == 0 and Afactor.subs(r, -1) != 0
            poly, half = cleared_middle(Afactor)
            assert poly.LC() == Afactor.subs(r, 1) > 0
            for value in (max(1, half), 3*(n-1), 3*(n-1)+7):
                denominator = sp.prod(value-half+i for i in range(1, half+1))
                central = sp.binomial(2*value-2*half, value-half)
                assert denominator > 0
                assert binomial_middle(Afactor, value)*denominator == central*poly.eval(value)
                binomial_controls += 1
            roots, factorization = exact_roots(poly)
            integer_records = []
            for root in sorted(roots):
                if root.q != 1:
                    continue
                value = int(root)
                integer_root_count += 1
                denominator = sp.prod(value-half+i for i in range(1, half+1))
                clearing_exception = denominator == 0
                clearing_exception_count += int(clearing_exception)
                in_domain = value >= max(1, half)
                direct = binomial_middle(Afactor, value) if in_domain else None
                if in_domain:
                    assert denominator > 0 and direct == 0
                passes = in_domain and value >= 3*(n-1)
                if passes:
                    admissible.append({"n": n, "m": value, "Afactor": str(Afactor)})
                integer_records.append({
                    "m": value, "multiplicity": roots[root],
                    "cayley_power": 2*value-2*half,
                    "positive_pole_and_factorial_domain": in_domain,
                    "clearing_denominator": str(denominator),
                    "clearing_exception": clearing_exception,
                    "original_middle_coefficient": str(direct) if direct is not None else None,
                    "admissible_pole_bound": value >= 3*(n-1),
                    "admissible_candidate": passes})
            divisor_records.append({
                "n": n,
                "selected_irreducible_factors": [str(f) for f, use in zip(factors, selection) if use],
                "Afactor": str(Afactor), "degree_Afactor": 2*half,
                "cleared_polynomial": str(poly.as_expr()), "factorization": factorization,
                "all_rational_roots": [{"root": str(q), "multiplicity": roots[q]} for q in sorted(roots)],
                "all_integer_roots": integer_records})
    assert len(divisor_records) == 44 and binomial_controls == 132
    assert integer_root_count == 91 and clearing_exception_count == 0
    assert admissible == []

    small = [(r**2+1, 2*(m-1), [1]),
             (r**2-r+1, m-2, [2]),
             (r**4-r**3+r**2-r+1, (m-3)*(m-4), [3, 4]),
             (r**6-r**3+1, m**3-21*m**2+92*m-120, [])]
    for Afactor, expected, integers in small:
        poly, half = cleared_middle(Afactor)
        assert sp.expand(poly.as_expr()-expected) == 0
        roots, _ = exact_roots(poly)
        assert sorted(int(q) for q in roots if q.q == 1) == integers
        for value in (max(1, half), max(12, half)):
            assert literal_middle(Afactor, value) == binomial_middle(Afactor, value)
    for Afactor, _, _ in small[:3]:
        poly, half = cleared_middle(Afactor)
        value = max(12, half)
        denominator = sp.prod(value-half+i for i in range(1, half+1))
        central = sp.binomial(2*value-2*half, value-half)
        assert literal_middle(Afactor, value)*denominator != central*(poly.eval(value)+1)
        negatives.append("Corrupted cleared polynomial rejected for "+str(Afactor))

    # An artificial reciprocal polynomial exercises the forbidden-domain
    # branch; it is deliberately NOT a divisor from the bounded computation.
    synthetic = r**4+2*r**3+2*r+1
    poly, half = cleared_middle(synthetic)
    assert half == 2 and poly.eval(1) == 0
    assert sp.prod(1-half+i for i in range(1, half+1)) == 0
    try:
        literal_middle(synthetic, 1)
    except ValueError:
        negatives.append("A cleared root outside m>=degree(Afactor)/2 was rejected.")
    else:
        raise AssertionError("A negative Cayley power was incorrectly accepted.")

    cayley_records = [cayley_fixture(2, 1, r**2+1),
                      cayley_fixture(5, 3, r**4-r**3+r**2-r+1),
                      cayley_fixture(2, 3, r**2+1)]
    assert cayley_records[0]["C_polynomial"] and not cayley_records[0]["phase_D_polynomial"]
    assert not cayley_records[1]["C_polynomial"] and not cayley_records[1]["phase_D_polynomial"]
    assert not cayley_records[2]["integrable"]
    negatives.extend(["Rational integration did not suffice for polynomial phase D.",
                      "A rationally integrable fixture with a pole in C was rejected.",
                      "A nonzero logarithmic residue was rejected."])
    malformed = r**2+2
    assert sp.rem(r**2+1, malformed, r) != 0
    negatives.append("A proposed Cayley factor not dividing r^n+1 was rejected.")
    repeated_records, repeated_jet_controls = repeated_pole_controls()

    local_records = []
    for height, e in itertools.product((1, 3), (1, 2, 3)):
        h = (-1)**e*x*(x-1)**e
        C, phase_D = pullback(sp.S.One, -sp.S.One, h, height)
        observed = [order_at(C, 1), order_at(phase_D, 1)]
        assert observed == [height-e, height+2-e]
        local_records.append({"case": "unmarked zero", "a": height, "e": e,
                              "orders": observed})
    for height, e, cdegree, ddegree in (
            (1, 2, 0, 0), (1, 2, 1, 2), (3, 2, 2, 3),
            (1, 3, 1, 1), (3, 3, 2, 4), (3, 3, 3, 4)):
        h = ((x-1)**e-(-1)**e)/(e*(-1)**(e-1))
        C, phase_D = pullback((1-e*z)**cdegree, -(1-e*z)**ddegree, h, height)
        observed = [order_at(C, 1), order_at(phase_D, 1)]
        assert observed == [e*cdegree-height*(e-1), e*ddegree-(height+2)*(e-1)]
        local_records.append({"case": "nonzero ramified value", "a": height,
                              "e": e, "outer_orders": [cdegree, ddegree], "orders": observed})
    for height, e, trivial in itertools.product((1, 3), (1, 2), (False, True)):
        cdegree, ddegree = (0, 0) if trivial else (height+1, height+3)
        h = x/(1-x)**e
        C, phase_D = pullback((1+z)**cdegree, -(1+z)**ddegree, h, height)
        observed = [order_at(C, 1), order_at(phase_D, 1)]
        assert observed == [height-e*(cdegree-1), height+2-e*(ddegree-1)]
        local_records.append({"case": "finite pole", "a": height, "e": e,
                              "outer_degrees": [cdegree, ddegree], "orders": observed})
    assert len(local_records) == 20
    assert any(min(case["orders"]) < 0 for case in local_records)
    negatives.append("Ramified pullbacks with an actual finite coefficient pole were rejected.")

    pure_records = []
    moment_controls = 0
    for degree, height in ((2, 1), (2, 3), (3, 3), (4, 3)):
        scalar = sp.Rational(degree, degree-1)
        u = x*(1-x)**(degree-1)/(1-scalar*x)**degree
        X = (1-x)*(1-scalar*x)
        assert cancel(u/(x*sp.diff(u, x))-X) == 0
        C, phase_D = pullback(sp.S.One, -sp.S.One, u, height)
        assert cancel(C-x*X**height/u) == 0
        assert cancel(phase_D+x*X**(height+2)/u) == 0
        assert polynomial(C) and polynomial(phase_D)
        for exponent in range(1, 5):
            assert raw_moment(C, phase_D, height, exponent) == 0
            moment_controls += 1
        pure_records.append({"family": "DZ trivial outer", "degree_u": degree,
                             "a": height, "u": str(cancel(u)), "X": str(X),
                             "degree_C": int(sp.degree(C, x)),
                             "degree_phase_D": int(sp.degree(phase_D, x)),
                             "control_scope": "Trivial-outer pullback and polynomial purity; no claim that this source lacks a polynomial direction."})
    for height in (1, 3):
        source_degree = 2
        h = ((1+x)**source_degree-1)/source_degree
        Cbar = (1+source_degree*z)**(height+1)
        Dbar = -(1+source_degree*z)**(height+3)
        C, phase_D = pullback(Cbar, Dbar, h, height)
        assert polynomial(C) and polynomial(phase_D)
        assert sp.degree(C, x) == source_degree*(sp.degree(Cbar, z)-1)+1
        assert sp.degree(phase_D, x) == source_degree*(sp.degree(Dbar, z)-1)+1
        u = cancel(h/(1+source_degree*h))
        X = cancel(u/(x*sp.diff(u, x)))
        assert cancel(C-x*X**height/u) == 0
        assert cancel(phase_D+x*X**(height+2)/u) == 0
        for exponent in range(1, 5):
            assert raw_moment(C, phase_D, height, exponent) == 0
            moment_controls += 1
        pure_records.append({"family": "polynomial pullback degree identity", "a": height,
                             "degree_h": source_degree, "degree_C": int(sp.degree(C, x)),
                             "degree_phase_D": int(sp.degree(phase_D, x)),
                             "outer_endpoint_is_odd": False,
                             "control_scope": "Pullback degree and polynomial purity; not a nontrivial odd-return core."})
    assert len(pure_records) == 6 and moment_controls == 24

    # Essential-hypothesis control: polynomial upstairs, rational outer pair.
    height = 1
    h = x/(1+x)
    Cbar, Dbar = 1/(1-z), -1/(1-z)**3
    C, phase_D = pullback(Cbar, Dbar, h, height)
    assert cancel(C-(1+x)**3) == 0 and cancel(phase_D+(1+x)**7) == 0
    assert polynomial(C) and polynomial(phase_D)
    assert not polynomial(Cbar, z) and not polynomial(Dbar, z)
    assert cancel((1+h)/(1-h)-(1+2*x)) == 0
    negatives.append("Dropping the no-polynomial-direction premise invalidates outer polynomiality.")

    return {
        "bounded_n_values": list(N_VALUES), "Q_divisor_cases": len(divisor_records),
        "all_integer_roots_recorded": integer_root_count,
        "actual_clearing_exceptions": clearing_exception_count,
        "binomial_coefficient_controls": binomial_controls,
        "small_cleared_polynomial_controls": len(small),
        "literal_small_coefficient_controls": 2*len(small),
        "admissible_Cayley_candidates": admissible,
        "Cayley_construction_controls": cayley_records,
        "repeated_pole_family_controls": repeated_records,
        "repeated_pole_nonzero_jet_controls": repeated_jet_controls,
        "source_pullback_local_controls": local_records,
        "pure_polynomial_pair_controls": pure_records,
        "moment_controls_total": moment_controls,
        "direct_even_raw_moment_computations": moment_controls//2,
        "odd_moment_parity_controls": moment_controls//2,
        "negative_controls": negatives,
        "divisor_records": divisor_records,
        "limits": ["Only the displayed n-list and rational reciprocal factors are enumerated.",
                   "No nontrivial pure odd-return polynomial core is exhibited or excluded universally.",
                   "The unrestricted descent and converse rely on the written proof."]}


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
        "status": "PASS", "python": platform.python_version(), "sympy": sp.__version__,
        "arithmetic": "Exact rational functions, polynomial factorization and integer roots; no floating-point acceptance.",
        "checker_sha256": sha256(Path(__file__)), "source_pins_checked": pins_checked,
        "packet_pins": {str(path): sha256(root/path) for path in (PROOF, CHECKER, LEDGER)},
        "scope": "Finite exact controls for the written reduction, and all-integer-m coverage for exactly 44 rational-divisor cases.",
        "verification": results}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=False)+"\n")
    print(json.dumps({"status": "PASS", "Q_divisor_cases": results["Q_divisor_cases"],
                      "all_integer_roots": results["all_integer_roots_recorded"],
                      "admissible_candidates": len(results["admissible_Cayley_candidates"]),
                      "pure_pair_controls": len(results["pure_polynomial_pair_controls"]),
                      "negative_controls": len(results["negative_controls"])}))


if __name__ == "__main__":
    main()
