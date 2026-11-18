#!/usr/bin/env python3
"""Arithmetic-certificate proof of the rational degree12 scalar affine exclusion.

Default execution rederives the necessary systems and replays five retained
arithmetic circuits. It makes no Groebner-basis call and requires no Singular.
The uniform reduction and theorem scope are supplied by the accompanying note.
"""
from __future__ import annotations
import argparse
import hashlib
import math
import json
import platform
from pathlib import Path
import sympy as sp
from sympy.polys.rings import ring

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
w, c, u, q0, q1, q2, lam, z = sp.symbols("w c u q0 q1 q2 lam z")
q = w**3 + q2*w**2 + q1*w + q0
EXCEPTIONS = {sp.Integer(0), -sp.Rational(1, 4), -sp.Rational(1, 3),
              -sp.Rational(1, 2), -sp.Integer(1)}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def primitive(expression, variables, allowed_denominators=()):
    numerator, denominator = sp.fraction(sp.cancel(expression))
    assert denominator.free_symbols <= set(allowed_denominators)
    return sp.Poly(numerator, *variables, domain=sp.QQ).primitive()[1].as_expr()


def denominator_certificates(bundle, retained):
    """Check every factor product and projective root-free residue witness."""
    evidence = json.loads((HERE/"degree12_scalar_affine_denominator_certificates.json").read_text())
    factors = []
    factor_roots = []
    witnesses = []
    for item in evidence["factors"]:
        coefficients = item["coefficients_descending"]
        assert all(isinstance(i, int) for i in coefficients)
        assert coefficients[0] > 0
        assert math.gcd(*coefficients) == 1
        polynomial = sp.Poly.from_list(coefficients, c, domain=sp.QQ)
        assert polynomial.degree() == item["degree"] >= 1
        factors.append(polynomial.monic())
        if polynomial.degree() == 1:
            root = -sp.Rational(coefficients[1], coefficients[0])
            assert str(root) == item["rational_root"] and root in EXCEPTIONS
            factor_roots.append({root: 1})
        else:
            prime = item["prime"]
            assert prime >= 2 and all(prime % d for d in range(2, math.isqrt(prime)+1))
            assert coefficients[0] % prime == item["leading_coefficient_mod_prime"] != 0
            values = []
            for x in range(prime):
                value = 0
                for coefficient in coefficients:
                    value = (value*x+coefficient) % prime
                values.append(value)
            assert values == item["values"] and all(values)
            factor_roots.append({})
            witnesses.append({"degree": polynomial.degree(), "prime": prime,
                              "leading_coefficient_mod_prime": coefficients[0] % prime,
                              "all_finite_values_nonzero": True})

    polynomials = []
    lookup = {}
    for item in evidence["polynomials"]:
        target = sp.Poly.from_list(item["monic_coefficients_descending"], c, domain=sp.QQ)
        assert target.LC() == 1
        product = sp.Poly(1, c, domain=sp.QQ)
        roots = {}
        for index, exponent in item["factorization"]:
            assert 0 <= index < len(factors) and isinstance(exponent, int) and exponent > 0
            product *= factors[index]**exponent
            for root, multiplicity in factor_roots[index].items():
                roots[root] = roots.get(root, 0)+exponent*multiplicity
        assert product == target
        polynomials.append(target)
        key = tuple(map(str, target.all_coeffs()))
        assert key not in lookup
        lookup[key] = roots

    maps = evidence["circuit_denominator_polynomial_indices"]
    cases = {"generic": bundle["generic"], **bundle["specializations"]}
    assert set(maps) == set(cases)
    for label, certificate in cases.items():
        assert len(maps[label]) == len(certificate["coefficients"])
        for row, index in zip(certificate["coefficients"], maps[label]):
            denominator = sp.Poly(sum(sp.Rational(v)*c**int(e)
                                      for e, v in row["denominator"]), c, domain=sp.QQ)
            assert denominator and denominator.monic() == polynomials[index]
    assert len(evidence["triangular_delta_polynomial_indices"]) == 7
    for text, index in zip(retained["deltas_j_0_to_6"],
                           evidence["triangular_delta_polynomial_indices"]):
        assert sp.Poly(sp.sympify(text), c, domain=sp.QQ).monic() == polynomials[index]
    return lookup, {"exact_factor_product_count": len(polynomials),
                    "linear_factor_roots": sorted(str(next(iter(r)))
                                                  for r in factor_roots if r),
                    "nonlinear_root_free_witnesses": witnesses,
                    "no_factorization_or_irreducibility_oracle_used": True}


def generic_system(retained, denominator_lookup):
    """Independently reconstruct the triangular trace/square reduction."""
    rho = 1 + c*w
    ss, pp = 2*rho + w, rho**2
    C = [sp.Integer(0), sp.Integer(1)]
    T = [sp.Integer(2), ss]
    for j in range(2, 7):
        C.append(sp.expand(ss*C[-1] - pp*C[-2]))
        T.append(sp.expand(ss*T[-1] - pp*T[-2]))
    D = [sp.expand(T[j] - w*C[j]) for j in range(7)]
    deltas = [sp.factor(sp.Poly(D[j], w).nth(j)) for j in range(7)]
    assert [str(t) for t in deltas] == retained["deltas_j_0_to_6"]
    rational_delta_roots = set()
    for delta in deltas:
        key = tuple(map(str, sp.Poly(delta, c, domain=sp.QQ).monic().all_coeffs()))
        rational_delta_roots.update(denominator_lookup[key])
    assert rational_delta_roots == {sp.Integer(0), -sp.Integer(1)}
    assert D[0] == 2
    assert all(sp.cancel(D[j].subs(w, -1/c)) == 0 for j in range(1, 7))

    domain = sp.QQ.frac_field(c).poly_ring(q0, q1, q2)
    residual = sp.Poly(q*q, w, domain=domain)
    coefficients = [sp.Integer(0)]*7
    for j in range(6, -1, -1):
        basis = sp.Poly(D[j], w, domain=domain)
        coefficients[j] = sp.cancel(residual.nth(j)/basis.nth(j))
        residual -= sp.Poly(coefficients[j]*D[j], w, domain=domain)
        assert residual.is_zero or residual.degree() < j
    assert residual.is_zero
    assert all(sp.cancel(actual - sp.sympify(saved)) == 0
               for actual, saved in zip(
                   coefficients, retained["N_coefficients_ascending"]))
    baseC = sp.Poly(sum(coefficients[j]*C[j] for j in range(7)),
                    w, domain=domain)
    divisor = sp.Poly(q, w, domain=domain)
    crem = baseC.rem(divisor, auto=False)
    erem = (baseC + 2*sp.Poly(w, w, domain=domain)*baseC.diff()).rem(
        divisor, auto=False)
    equations = [primitive(e, (q0, q1, q2, c), (c,)) for e in [
        crem.nth(1), crem.nth(2), erem.nth(0), erem.nth(1), erem.nth(2)]]
    assert [str(e) for e in equations] == retained["equations"]
    assert sp.cancel(crem.nth(0)-sp.sympify(
        retained["C_remainder_constant"])) == 0
    return equations + [z*q0*sp.discriminant(q, w)-1], deltas


def exceptional_system(value):
    """Direct quadratic trace chart; no triangular denominator inversion."""
    rho = 1 + value*w
    quadratic = u*u - (2*rho+w)*u + rho*rho
    critical = sp.Poly(sp.resultant(q, quadratic, w), u)
    assert critical.degree() == 6 and critical.LC() == 1
    numerator = sp.Poly(sum(11*critical.nth(j)*u**j/sp.Integer(2*j-1)
                            for j in range(7)), u)
    assert sp.expand(2*u*sp.diff(numerator.as_expr(), u) -
                     numerator.as_expr() - 11*critical.as_expr()) == 0
    rem = sp.rem(lam*numerator.as_expr() - (u-rho), quadratic, u)
    equations = []
    for i in range(2):
        coeffrem = sp.Poly(sp.rem(sp.Poly(rem, u).nth(i), q, w), w)
        equations.extend(primitive(coeffrem.nth(j), (q0, q1, q2, lam))
                         for j in range(3))
    normrho = sp.resultant(q, rho, w)
    normdelta = sp.resultant(q, w*(w+4*rho), w)
    open_product = primitive(
        lam*q0*sp.discriminant(q, w)*normrho*normdelta,
        (q0, q1, q2, lam))
    return equations + [z*open_product-1]


def replay(certificate, expected_inputs, allowed_roots, denominator_lookup):
    """Verify ideal membership by arithmetic only; no polynomial division."""
    domain = sp.QQ.frac_field(c)
    R, *_ = ring(",".join(certificate["variables"]), domain, order="grevlex")
    parameter_ring = domain.field.ring
    coefficients = []
    denominator_controls = []
    seen = set()
    for item in certificate["coefficients"]:
        def parse(rows):
            return parameter_ring.from_dict({
                (int(power),): sp.QQ(coefficient)
                for power, coefficient in rows})
        numerator, denominator = parse(item["numerator"]), parse(item["denominator"])
        assert denominator
        coefficients.append(domain.field.new(numerator, denominator))
        key = str(denominator.monic())
        if key not in seen:
            seen.add(key)
            polynomial = sp.Poly(denominator.as_expr(), c, domain=sp.QQ)
            root_key = tuple(map(str, polynomial.monic().all_coeffs()))
            roots = denominator_lookup[root_key]
            assert set(roots) <= allowed_roots
            denominator_controls.append({
                "degree": polynomial.degree(),
                "rational_roots": {str(r): int(m) for r, m in roots.items()}})

    def polynomial(terms):
        exponents = [tuple(map(int, ex)) for ex, _ in terms]
        assert len(set(exponents)) == len(exponents)
        assert all(len(ex) == R.ngens and all(i >= 0 for i in ex)
                   for ex in exponents)
        return R.from_dict({ex: coefficients[int(index)]
                            for ex, (_, index) in zip(exponents, terms)})
    inputs = [polynomial(p) for p in certificate["inputs"]]
    assert inputs == [R.from_expr(e) for e in expected_inputs]
    nodes = []
    for index, node in enumerate(certificate["nodes"]):
        def parent(i):
            assert 0 <= i < index
            return nodes[i]
        op = node["op"]
        if op == "input":
            value = inputs[node["index"]]
        elif op == "zero":
            value = R.zero
        elif op == "scale":
            value = parent(node["parent"]).mul_ground(
                coefficients[node["coefficient"]])
        elif op == "monomial_difference":
            left, right = tuple(node["left_monomial"]), tuple(node["right_monomial"])
            assert len(left) == len(right) == R.ngens
            assert all(isinstance(i, int) and i >= 0 for i in left+right)
            value = (parent(node["left"]).mul_monom(left) -
                     parent(node["right"]).mul_monom(right))
        elif op == "remainder_identity":
            value = parent(node["parent"])
            for term in node["subtract"]:
                value = value-polynomial(term["multiplier"])*parent(term["parent"])
        else:
            raise AssertionError(op)
        nodes.append(value)
    assert nodes[certificate["final_node"]] == R.one
    return {"exact_final_identity": "1", "nodes": len(nodes),
            "field_coefficients": len(coefficients),
            "distinct_coefficient_denominators": len(seen),
            "denominator_controls": denominator_controls}


def run():
    system_path = HERE/"degree12_scalar_affine_reduced_system.json"
    cert_path = HERE/"degree12_scalar_affine_arithmetic_certificates.json"
    retained = json.loads(system_path.read_text())
    certificates = json.loads(cert_path.read_text())
    denominator_lookup, denominator_report = denominator_certificates(certificates, retained)
    generic_inputs, deltas = generic_system(retained, denominator_lookup)
    generic = certificates["generic"]
    assert generic["input_system_sha256"] == digest(system_path)
    result = replay(generic, generic_inputs, EXCEPTIONS, denominator_lookup)
    special = {}
    assert set(certificates["specializations"]) == {"-1", "-1/2", "-1/3", "-1/4"}
    for key, certificate in certificates["specializations"].items():
        assert certificate["fixed_parameter"] == key
        special[key] = replay(certificate, exceptional_system(sp.Rational(key)), set(), denominator_lookup)

    # Critical-pair derivative simplification for arbitrary affine rho.
    rr = sp.symbols("rr")
    rho = rr + c*w
    delta = w*(w+4*rho)
    assert sp.expand((w+4*rho)*sp.diff(delta, w) -
                     2*(2*c+1)*delta - 4*rho*(w+4*rho)) == 0

    paths = [
        HERE/"DEGREE_TWELVE_SCALAR_AFFINE_OBSTRUCTION.md",
        Path(__file__).resolve(), system_path, cert_path,
        HERE/"degree12_scalar_affine_denominator_certificates.json",
        HERE/"PAIRED_CRITICAL_TRACE_CHART.md",
        HERE/"SCALAR_SECANT_PROPORTIONAL_SCALE_OBSTRUCTION.md"]
    return {
        "schema_version": 1, "result": "PASS",
        "python": platform.python_version(), "sympy": sp.__version__,
        "claim": "No good rational degree12 even-numerator paired pencil has scalar secant and affine reciprocal scale.",
        "uniform_written_reduction_required": True,
        "exact_normalization_requires_no_radical_extension": True,
        "generic_certificate": result,
        "denominator_arithmetic_certificates": denominator_report,
        "generic_rational_exceptional_parameters": ["-1", "-1/2", "-1/3", "-1/4", "0"],
        "triangular_denominators": list(map(str, deltas)),
        "triangular_rational_exceptions": ["-1", "0"],
        "exceptional_direct_certificates": special,
        "c_zero": "Excluded by frozen both-scalar Proposition2.",
        "no_groebner_call_in_replay": True,
        "no_rational_denominator_restriction": True,
        "scope_limit": "Scalar secant with affine scale in degree12 only; unrestricted nonreciprocal models remain open.",
        "sha256": {str(path.relative_to(ROOT)): digest(path) for path in paths},
        "assurance": {"formal_verification": False, "independent_external_replay": False,
                      "external_human_review": False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    result = run()
    text = json.dumps(result, indent=2, sort_keys=True)+"\n"
    if args.write_receipt:
        args.write_receipt.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()

