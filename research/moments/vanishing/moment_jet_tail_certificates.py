"""Exact local moment-tail certificates at the exceptional cubic shear.

The general finite-jet/Nakayama theorem is in the companion proof. This
backend derives its inputs from the actual Gaussian or normalized S^2
functional at P*=W+ZY-Z^3/2. It proves local/completed ideal equality,
not global membership and not purity from an unproved base hypothesis.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import json
from math import factorial
from pathlib import Path

import sympy as sp

Z, W, Y, T = sp.symbols("Z W Y t")
SPATIAL = (Z, W, Y)
PSTAR = W + Z*Y - Z**3/2


def monomials(count, degree, exact=False):
    """Parameter exponent tuples in increasing degree/lexicographic order."""
    if count < 1 or degree < 0:
        return []
    def compositions(total, length):
        if length == 1:
            yield (total,)
        else:
            for first in range(total+1):
                for tail in compositions(total-first, length-1):
                    yield (first,)+tail
    degrees = [degree] if exact else range(degree+1)
    return [powers for total in degrees for powers in compositions(total, count)]


def monomial(variables, powers):
    return sp.prod(variable**power for variable, power in zip(variables, powers))


def truncate(expression, parameters, order):
    if order < 0:
        return sp.Integer(0)
    return sp.expand(sum(coefficient*monomial(parameters, powers)
                         for powers, coefficient in sp.Poly(expression, *parameters).terms()
                         if sum(powers) <= order))


def encode_polynomial(expression, parameters):
    return [{"powers": list(powers), "coefficient": str(coefficient)}
            for powers, coefficient in sp.Poly(sp.expand(expression), *parameters,
                                               domain=sp.QQ).terms() if coefficient]


def decode_polynomial(terms, parameters):
    answer = sp.Integer(0)
    seen = set()
    for term in terms:
        powers = tuple(term["powers"])
        if len(powers) != len(parameters) or any(type(v) is not int or v < 0 for v in powers):
            raise ValueError("Invalid nonnegative parameter exponent tuple")
        if powers in seen:
            raise ValueError("Duplicate polynomial monomial")
        seen.add(powers)
        coefficient = sp.Rational(term["coefficient"])
        answer += coefficient*monomial(parameters, powers)
    return sp.expand(answer)


def decode_instance(instance):
    if instance.get("functional") not in ("gaussian", "sphere_3"):
        raise ValueError("Only Gaussian and normalized three-coordinate sphere functionals are built in")
    names = instance["parameters"]
    if not names or len(set(names)) != len(names) or any(not isinstance(n, str) or not n for n in names):
        raise ValueError("Supply distinct nonempty parameter names")
    # Internal names avoid collisions with spatial coordinates and t.
    parameters = tuple(sp.Symbol("parameter_"+str(i)) for i in range(len(names)))
    perturbation = sp.Integer(0)
    degree = 0
    for term in instance["perturbation"]:
        pp, xx = tuple(term["parameter_powers"]), tuple(term["spatial_powers"])
        if (len(pp) != len(parameters) or len(xx) != 3
                or any(type(v) is not int or v < 0 for v in pp+xx) or not sum(pp)):
            raise ValueError("Perturbations must vanish at the parameter origin")
        coefficient = sp.Rational(term["coefficient"])
        if coefficient:
            perturbation += coefficient*monomial(parameters, pp)*monomial(SPATIAL, xx)
            degree = max(degree, sum(xx))
    if degree > 3:
        raise ValueError("This cubic backend requires spatial perturbation degree at most three")
    prefix, nilpotence = instance["prefix_N"], instance["nilpotence_nu"]
    if type(prefix) is not int or prefix < 0 or type(nilpotence) is not int or nilpotence < 1:
        raise ValueError("N must be nonnegative and nu positive")
    relations = [decode_polynomial(terms, parameters)
                 for terms in instance.get("slice_relations", [])]
    if any(value.subs(dict.fromkeys(parameters, 0)) != 0 for value in relations):
        raise ValueError("Slice relations must pass through the parameter origin")
    return parameters, sp.expand(perturbation), degree, relations


def encode_instance(perturbation, parameters, functional, prefix, nilpotence, relations=()):
    """Build a sparse, executable-code-free instance from exact expressions."""
    terms = []
    for powers, coefficient in sp.Poly(sp.expand(perturbation), *parameters, Z, W, Y,
                                       domain=sp.QQ).terms():
        if coefficient:
            terms.append({"parameter_powers": list(powers[:len(parameters)]),
                          "spatial_powers": list(powers[len(parameters):]),
                          "coefficient": str(coefficient)})
    return {"functional": functional, "parameters": [str(p) for p in parameters],
            "perturbation": terms, "prefix_N": prefix, "nilpotence_nu": nilpotence,
            "slice_relations": [encode_polynomial(value, parameters) for value in relations]}


@lru_cache(maxsize=None)
def odd_factorial(last):
    return sp.prod(range(1, last+1, 2)) if last >= 1 else sp.Integer(1)


def linear_null_period(marker, functional):
    """L(marker exp(tW)), retaining all original spatial coefficients."""
    answer = sp.Integer(0)
    for (a, b, j), coefficient in sp.Poly(sp.expand(marker), Z, W, Y).terms():
        if a < b or j % 2:
            continue
        value = sp.Rational(factorial(a), factorial(a-b))*odd_factorial(j-1)
        if functional == "sphere_3":
            value /= odd_factorial(2*a+j+1)
        answer += coefficient*value*T**(a-b)
    return sp.expand(answer)


@lru_cache(maxsize=16)
def _shear_moment_jet_items(perturbation, parameters, order, functional):
    """All nonzero moment jets through parameter order `order`.

    The omitted moments are zero at this jet order by the written sharp
    cutoff. This uses the original inverse shear and exact monomial L.
    """
    transformed = sp.expand(perturbation.subs(
        {Z: Z, W: W-Z*Y-Z**3/2, Y: Y+Z**2}, simultaneous=True))
    exponential_jet = sp.Integer(1)
    power = sp.Integer(1)
    for j in range(1, order+1):
        power = truncate(sp.expand(power*transformed), parameters, order)
        if power == 0:
            break
        exponential_jet += T**j*linear_null_period(power, functional)/factorial(j)
    polynomial = sp.Poly(sp.expand(exponential_jet), T)
    return tuple((n, sp.expand(coefficient*factorial(n)))
                 for (n,), coefficient in polynomial.terms() if n > 0 and coefficient)


def shear_moment_jets(perturbation, parameters, order, functional):
    return dict(_shear_moment_jet_items(perturbation, tuple(parameters), order, functional))


def _terms(expression, parameters):
    return {powers: sp.Rational(coefficient)
            for powers, coefficient in sp.Poly(expression, *parameters,
                                               domain=sp.QQ).terms() if coefficient}


def _add_scaled(target, source, scalar):
    for key, value in source.items():
        updated = target.get(key, 0)+scalar*value
        if updated:
            target[key] = updated
        else:
            target.pop(key, None)


class JetIdealBasis:
    """Sparse exact elimination with provenance in ORIGINAL generators."""
    def __init__(self, generators, parameters, order):
        self.generators = generators
        self.parameters = parameters
        self.order = order
        self.rows = {}
        for index, generator in enumerate(generators):
            generator = truncate(generator, parameters, order)
            if generator == 0:
                continue
            initial_order = min(sum(p) for p in _terms(generator, parameters))
            for powers in monomials(len(parameters), order-initial_order):
                value = truncate(generator*monomial(parameters, powers), parameters, order)
                vector = _terms(value, parameters)
                provenance = {(index, powers): sp.Integer(1)}
                while vector:
                    pivot = min(vector, key=lambda p: (sum(p), p))
                    if pivot not in self.rows:
                        scale = vector[pivot]
                        vector = {p: c/scale for p, c in vector.items()}
                        provenance = {p: c/scale for p, c in provenance.items()}
                        self.rows[pivot] = (vector, provenance)
                        break
                    row, identity = self.rows[pivot]
                    scale = vector[pivot]
                    _add_scaled(vector, row, -scale)
                    _add_scaled(provenance, identity, -scale)

    def witness(self, target):
        vector = _terms(truncate(target, self.parameters, self.order), self.parameters)
        answer = {}
        while vector:
            pivot = min(vector, key=lambda p: (sum(p), p))
            if pivot not in self.rows:
                return None
            row, identity = self.rows[pivot]
            scale = vector[pivot]
            _add_scaled(vector, row, -scale)
            _add_scaled(answer, identity, scale)
        multipliers = [sp.Integer(0) for _ in self.generators]
        for (index, powers), coefficient in answer.items():
            multipliers[index] += coefficient*monomial(self.parameters, powers)
        return [sp.expand(value) for value in multipliers]

    def normal_form(self, target):
        vector = _terms(truncate(target, self.parameters, self.order), self.parameters)
        remainder = {}
        while vector:
            pivot = min(vector, key=lambda p: (sum(p), p))
            if pivot not in self.rows:
                remainder[pivot] = vector.pop(pivot)
            else:
                row, _ = self.rows[pivot]
                _add_scaled(vector, row, -vector[pivot])
        return remainder

    def nonmembership_certificate(self, target_powers):
        target = monomial(self.parameters, target_powers)
        remainder = self.normal_form(target)
        if not remainder:
            return None
        selected = min(remainder, key=lambda p: (sum(p), p))
        scale = remainder[selected]
        dual = []
        for powers in monomials(len(self.parameters), self.order):
            value = self.normal_form(monomial(self.parameters, powers)).get(selected, 0)/scale
            if value:
                dual.append({"powers": list(powers), "coefficient": str(value)})
        return {"target": list(target_powers), "dual_coefficients": dual}


def _jet_input(instance):
    parameters, perturbation, degree, relations = decode_instance(instance)
    nu, prefix = instance["nilpotence_nu"], instance["prefix_N"]
    jets = shear_moment_jets(perturbation, parameters, nu, instance["functional"])
    if any(value.subs(dict.fromkeys(parameters, 0)) != 0 for value in jets.values()):
        raise ValueError("Moment jets must vanish at the known pure base")
    labels = [n for n in range(1, min(prefix, (3*degree+1)*nu)+1) if jets.get(n, 0) != 0]
    generators = [jets[n] for n in labels]+relations
    return parameters, degree, jets, labels, generators


def produce_certificate(instance):
    """Construct explicit finite witnesses, or return an exact failed test."""
    parameters, degree, jets, labels, generators = _jet_input(instance)
    nu, prefix = instance["nilpotence_nu"], instance["prefix_N"]
    basis = JetIdealBasis(generators, parameters, nu)
    nilpotence = []
    for powers in monomials(len(parameters), nu, exact=True):
        identity = basis.witness(monomial(parameters, powers))
        if identity is None:
            return {"status": "rejected_nilpotence", "target": list(powers)}
        nilpotence.append({"target": list(powers),
                           "multipliers": [encode_polynomial(v, parameters) for v in identity]})
    bound = (3*degree+1)*(nu-1)
    lower_basis = JetIdealBasis(generators, parameters, nu-1)
    intervening = []
    for n in range(prefix+1, bound+1):
        identity = lower_basis.witness(jets.get(n, sp.Integer(0)))
        if identity is None:
            return {"status": "rejected_intervening_moment", "moment": n,
                    "jet": encode_polynomial(truncate(jets.get(n, 0), parameters, nu-1), parameters)}
        intervening.append({"moment": n,
                            "multipliers": [encode_polynomial(v, parameters) for v in identity]})
    packet = {"schema_version": 1, "instance": instance,
            "prefix_generator_moments": labels,
            "nilpotence_witnesses": nilpotence,
            "intervening_witnesses": intervening}
    if instance.get("certify_exact_nilpotence", False):
        for powers in monomials(len(parameters), nu-1, exact=True):
            dual = lower_basis.nonmembership_certificate(powers)
            if dual is not None:
                packet["nilpotence_minimality"] = dual
                break
        else:
            return {"status": "rejected_exact_nilpotence", "reason": "m^(nu-1) already vanishes in the jet quotient"}
    return packet


def verify_certificate(certificate):
    """Replay polynomial identities; no producer status or Groebner flag is trusted."""
    if certificate.get("schema_version") != 1:
        raise ValueError("Unknown certificate schema")
    instance = certificate["instance"]
    parameters, degree, jets, labels, generators = _jet_input(instance)
    if certificate["prefix_generator_moments"] != labels:
        raise ValueError("Prefix generator labels differ from actual moment jets")
    nu, prefix = instance["nilpotence_nu"], instance["prefix_N"]

    def check_identity(target, encoded_multipliers, order):
        if len(encoded_multipliers) != len(generators):
            raise ValueError("Every original generator must have a multiplier")
        multipliers = [decode_polynomial(terms, parameters) for terms in encoded_multipliers]
        residual = target-sum(a*g for a, g in zip(multipliers, generators))
        if truncate(sp.expand(residual), parameters, order) != 0:
            raise ValueError("False truncated polynomial membership identity")

    expected = set(monomials(len(parameters), nu, exact=True))
    seen = set()
    for identity in certificate["nilpotence_witnesses"]:
        powers = tuple(identity["target"])
        if powers not in expected or powers in seen:
            raise ValueError("Invalid or duplicate nilpotence target")
        check_identity(monomial(parameters, powers), identity["multipliers"], nu)
        seen.add(powers)
    if seen != expected:
        raise ValueError("Missing degree-nu monomial witness")

    bound = (3*degree+1)*(nu-1)
    required = set(range(prefix+1, bound+1))
    seen = set()
    for identity in certificate["intervening_witnesses"]:
        n = identity["moment"]
        if type(n) is not int or n not in required or n in seen:
            raise ValueError("Invalid or duplicate intervening moment")
        check_identity(jets.get(n, sp.Integer(0)), identity["multipliers"], nu-1)
        seen.add(n)
    if seen != required:
        raise ValueError("Missing intervening moment witness")
    exact_nilpotence = False
    if "nilpotence_minimality" in certificate:
        dual_packet = certificate["nilpotence_minimality"]
        powers = tuple(dual_packet["target"])
        if powers not in set(monomials(len(parameters), nu-1, exact=True)):
            raise ValueError("Minimality target must have parameter degree nu-1")
        dual = _terms(decode_polynomial(dual_packet["dual_coefficients"], parameters), parameters)
        if any(sum(p) > nu-1 for p in dual):
            raise ValueError("Dual certificate lies outside the prescribed jet algebra")
        def pairing(expression):
            return sum(coefficient*dual.get(p, 0)
                       for p, coefficient in _terms(truncate(expression, parameters, nu-1), parameters).items())
        if pairing(monomial(parameters, powers)) != 1:
            raise ValueError("Dual certificate does not detect its target")
        for generator in generators:
            for multiplier in monomials(len(parameters), nu-1):
                if pairing(generator*monomial(parameters, multiplier)) != 0:
                    raise ValueError("Dual certificate fails to annihilate the full truncated ideal")
        exact_nilpotence = True
    if instance.get("certify_exact_nilpotence", False) and not exact_nilpotence:
        raise ValueError("Requested exact nilpotence lacks its nonmembership certificate")
    return {"status": "passed", "functional": instance["functional"],
            "base_phase": "W+Z*Y-Z^3/2", "known_pure_base": True,
            "prefix_N": prefix, "nilpotence_nu": nu,
            "nilpotence_exponent_is_exact": exact_nilpotence,
            "perturbation_spatial_degree": degree, "mixed_cutoff_C": 3,
            "last_required_moment": bound, "intervening_moments_checked": len(required),
            "conclusion": "I_infinity=I_N in the local ring at the parameter origin and its completion",
            "global_ideal_equality_claimed": False,
            "cutoff_source": "MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md, Sections 1-3"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--produce", action="store_true", help="Read an instance and export explicit witnesses")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    value = json.loads(args.certificate.read_text(encoding="utf-8"))
    result = produce_certificate(value) if args.produce else verify_certificate(value)
    if args.produce and result.get("schema_version"):
        verify_certificate(result)
    report = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(report, encoding="utf-8", newline="\n")
    print(report, end="")
    if result.get("status", "passed") != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
