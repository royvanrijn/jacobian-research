"""Replay exact QQ local moment-tail certificates at the cubic shear.

Default execution replays retained polynomial identities and dual
annihilators. --generate rebuilds them from the literal exact inputs below.
The written mixed-cutoff/Nakayama theorem supplies the infinite reduction;
these finite certificates assert local/completed, rather than global,
ideal equality. The historical finite.tex input is pinned as context only.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
from math import factorial, prod
from pathlib import Path
import platform

import sympy as sp

import moment_jet_tail_certificates as backend


if not __debug__:
    raise RuntimeError("Run this exact checker without Python optimization (-O)")


DIRECTORY = Path(__file__).resolve().parent
ROOT = DIRECTORY.parents[2]
FIXTURES = DIRECTORY / "moment_jet_tail_certificate_examples.json"
PROOF = DIRECTORY / "MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md"
CONTEXT = ROOT / "unpublished/cubic-sphere-image/finite.tex"


def literal_instances():
    """The source phases are expressions over QQ, never executable input text."""
    e, f, u, v = sp.symbols("eps eta u v")
    profile_e = sp.Symbol("e")
    z, w, y = backend.SPATIAL

    def instance(phase, parameters, functional, prefix, nu, relations=()):
        result = backend.encode_instance(phase, parameters, functional,
                                         prefix, nu, relations)
        result["certify_exact_nilpotence"] = True
        return result

    profile = (
        u*(y**2-sp.Rational(1, 3)) + v*z*(y**2-sp.Rational(1, 5)) + profile_e*z**2*y
        - sp.Rational(2, 15)*u**2*z
        - (sp.Rational(124, 945)*u**3+sp.Rational(8, 35)*u*v)*z**2
        - (sp.Rational(2312, 14175)*u**4+sp.Rational(188, 525)*u**2*v
           + sp.Rational(8, 75)*v**2)*z**3
    )
    positives = {
        "gaussian_immersed_nu2": instance(e*(y-z**2), (e,), "gaussian", 2, 2),
        "sphere_immersed_nu2": instance(e*(y-z**2), (e,), "sphere_3", 2, 2),
        "sphere_three_parameter_profile_nu7": instance(profile, (u, v, profile_e), "sphere_3", 7, 7),
        "gaussian_two_parameter_nu2": instance(e*(y-z**2)+f, (e, f), "gaussian", 2, 2),
        "gaussian_local_not_global": instance(
            e*(1-e)*(y-z**2)+e**3*(y**2-1), (e,), "gaussian", 2, 2),
    }
    negatives = {}
    for name, functional, center, coefficient in [
        ("gaussian_intervening_failure", "gaussian", sp.Integer(1), sp.Integer(120)),
        ("sphere_intervening_failure", "sphere_3", sp.Rational(1, 3), sp.Rational(8, 63)),
    ]:
        value = backend.encode_instance(e*(y**2-center), (e,), functional, 2, 2)
        negatives[name] = (value, coefficient)
    return positives, negatives


def negative_packet(instance):
    """A valid nilpotence witness followed by an intentionally false M5 identity."""
    parameters, degree, jets, labels, generators = backend._jet_input(instance)
    nu = instance["nilpotence_nu"]
    basis = backend.JetIdealBasis(generators, parameters, nu)
    nilpotence = []
    for powers in backend.monomials(len(parameters), nu, exact=True):
        identity = basis.witness(backend.monomial(parameters, powers))
        assert identity is not None
        nilpotence.append({"target": list(powers), "multipliers": [
            backend.encode_polynomial(value, parameters) for value in identity]})
    bound = (3*degree+1)*(nu-1)
    return {
        "schema_version": 1, "instance": instance,
        "prefix_generator_moments": labels,
        "nilpotence_witnesses": nilpotence,
        "intervening_witnesses": [
            {"moment": n, "multipliers": [[] for _ in generators]}
            for n in range(instance["prefix_N"]+1, bound+1)],
    }


def least_prefix_packet(instance):
    """Produce a coefficient dual separating actual F7 from the six-moment ideal."""
    parameters, perturbation, _, relations = backend.decode_instance(instance)
    order, prefix, moment = 7, 6, 7
    jets = backend.shear_moment_jets(perturbation, parameters, order, instance["functional"])
    generators = [jets.get(n, sp.Integer(0)) for n in range(1, prefix+1)]+relations
    basis = backend.JetIdealBasis(generators, parameters, order)
    remainder = basis.normal_form(jets.get(moment, sp.Integer(0)))
    assert remainder, "F7 has no nonmembership dual at the requested jet order"
    selected = min(remainder, key=lambda powers: (sum(powers), powers))
    scale = remainder[selected]
    dual = []
    for powers in backend.monomials(len(parameters), order):
        value = basis.normal_form(backend.monomial(parameters, powers)).get(selected, 0)/scale
        if value:
            dual.append({"powers": list(powers), "coefficient": str(value)})
    return {"source_instance": instance, "jet_order": order, "prefix_N": prefix,
            "nonmember_moment": moment, "dual_coefficients": dual}


def generate(path):
    positives, negatives = literal_instances()
    packet = {
        "schema": "moment-jet-tail-examples-v1", "coefficient_field": "QQ",
        "fixed_base_phase": "W+Z*Y-Z^3/2",
        "scope": "Explicit local/completed membership and exact nilpotence certificates; no general global ideal claim.",
        "positive": [], "negative": [],
    }
    for name, instance in positives.items():
        certificate = backend.produce_certificate(instance)
        assert certificate.get("schema_version") == 1, (name, certificate)
        result = backend.verify_certificate(certificate)
        assert result["nilpotence_exponent_is_exact"] is True
        packet["positive"].append({"id": name, "certificate": certificate})
    for name, (instance, coefficient) in negatives.items():
        packet["negative"].append({
            "id": name, "certificate": negative_packet(instance),
            "first_nonmember_moment": 5,
            "linear_coefficient": str(coefficient),
        })
    packet["least_prefix_witness"] = least_prefix_packet(
        positives["sphere_three_parameter_profile_nu7"])
    path.write_text(json.dumps(packet, indent=2)+"\n", encoding="utf-8", newline="\n")


def require_rejection(certificate, label, expected_text):
    try:
        backend.verify_certificate(certificate)
    except ValueError as error:
        if expected_text not in str(error):
            raise AssertionError((label, str(error))) from error
        return {"case": label, "rejected": True, "reason": str(error)}
    raise AssertionError("Corrupted certificate was accepted: "+label)


def verify_negative(entry, expected_instance, coefficient):
    certificate = entry["certificate"]
    assert certificate["instance"] == expected_instance
    assert entry["first_nonmember_moment"] == 5
    assert sp.Rational(entry["linear_coefficient"]) == coefficient
    parameters, degree, jets, labels, generators = backend._jet_input(expected_instance)
    assert certificate["prefix_generator_moments"] == labels
    nu = expected_instance["nilpotence_nu"]
    expected_targets = set(backend.monomials(len(parameters), nu, exact=True))
    seen = set()
    for identity in certificate["nilpotence_witnesses"]:
        powers = tuple(identity["target"])
        assert powers in expected_targets and powers not in seen
        seen.add(powers)
        multipliers = [backend.decode_polynomial(value, parameters)
                       for value in identity["multipliers"]]
        assert len(multipliers) == len(generators)
        residual = backend.monomial(parameters, powers)-sum(
            a*g for a, g in zip(multipliers, generators))
        assert backend.truncate(sp.expand(residual), parameters, nu) == 0
    assert seen == expected_targets
    required = list(range(expected_instance["prefix_N"]+1, (3*degree+1)*(nu-1)+1))
    assert [item["moment"] for item in certificate["intervening_witnesses"]] == required
    # The WHOLE truncated prefix ideal is zero in degree<=1. A nonzero M5
    # there proves nonmembership, independently of producer rejection flags.
    assert all(backend.truncate(g, parameters, nu-1) == 0 for g in generators)
    first_nonzero = next(n for n in required
                         if backend.truncate(jets.get(n, 0), parameters, nu-1) != 0)
    assert first_nonzero == 5
    assert backend.truncate(jets[5], parameters, nu-1) == coefficient*parameters[0]
    rejection = require_rejection(certificate, entry["id"],
                                  "False truncated polynomial membership identity")
    return {"case": entry["id"], "nilpotence_witnesses_replayed": len(seen),
            "first_nonmember_moment": first_nonzero,
            "nonzero_linear_coefficient": str(coefficient), **rejection}


def corruptions(certificates):
    answer = []
    original = certificates["gaussian_immersed_nu2"]
    value = deepcopy(original)
    value["nilpotence_witnesses"].pop()
    answer.append(require_rejection(value, "missing_nilpotence_monomial",
                                    "Missing degree-nu monomial witness"))
    value = deepcopy(original)
    value["intervening_witnesses"].pop()
    answer.append(require_rejection(value, "missing_intervening_moment",
                                    "Missing intervening moment witness"))
    value = deepcopy(original)
    term = next(term for multiplier in value["nilpotence_witnesses"][0]["multipliers"]
                for term in multiplier)
    term["coefficient"] = str(sp.Rational(term["coefficient"])+1)
    answer.append(require_rejection(value, "corrupted_identity_coefficient",
                                    "False truncated polynomial membership identity"))
    value = deepcopy(original)
    value.pop("nilpotence_minimality")
    answer.append(require_rejection(value, "missing_requested_dual",
                                    "Requested exact nilpotence lacks"))
    value = deepcopy(original)
    value["nilpotence_minimality"]["dual_coefficients"][0]["coefficient"] = "0"
    answer.append(require_rejection(value, "corrupted_dual",
                                    "Dual certificate does not detect its target"))
    value = deepcopy(original)
    value["prefix_generator_moments"].pop()
    answer.append(require_rejection(value, "wrong_prefix_generator_labels",
                                    "Prefix generator labels differ"))
    value = deepcopy(original)
    value["instance"]["perturbation"].append({
        "parameter_powers": [0], "spatial_powers": [0, 0, 0], "coefficient": "1"})
    answer.append(require_rejection(value, "base_origin_violation",
                                    "Perturbations must vanish at the parameter origin"))
    value = deepcopy(certificates["gaussian_two_parameter_nu2"])
    value["nilpotence_witnesses"] = [item for item in value["nilpotence_witnesses"]
                                     if item["target"] != [1, 1]]
    answer.append(require_rejection(value, "missing_mixed_degree_nu_monomial",
                                    "Missing degree-nu monomial witness"))
    return answer


def direct_moment(phase, parameters, n, order, functional):
    """Independent original-coordinate contraction; no inverse-shear formula."""
    result = sp.Integer(0)
    for powers, coefficient in sp.Poly(sp.expand(phase**n),
                                       *backend.SPATIAL, *parameters).terms():
        a, b, c = powers[:3]
        pp = powers[3:]
        if sum(pp) > order or a != b or c % 2:
            continue
        value = sp.Rational(factorial(a)*factorial(c), 2**(c//2)*factorial(c//2))
        if functional == "sphere_3":
            value /= prod(range(1, 2*a+c+2, 2))
        result += coefficient*value*backend.monomial(parameters, pp)
    return sp.expand(result)


def direct_controls():
    e, f = sp.symbols("eps eta")
    z, w, y = backend.SPATIAL
    perturbation = e*z*w+f*y**2+e**2*z**3
    reports = []
    for functional in ("gaussian", "sphere_3"):
        jets = backend.shear_moment_jets(perturbation, (e, f), 2, functional)
        for n in range(1, 7):
            actual = direct_moment(backend.PSTAR+perturbation, (e, f), n, 2, functional)
            assert sp.expand(actual-jets.get(n, 0)) == 0
        reports.append({"functional": functional, "moments": [1, 6],
                        "parameters": 2, "jet_order": 2,
                        "method": "Independent original-coordinate Wick or sphere contraction"})
    local_phase = backend.PSTAR+e*(1-e)*(y-z**2)+e**3*(y**2-1)
    m2 = direct_moment(local_phase, (e,), 2, 9, "gaussian")
    m3 = direct_moment(local_phase, (e,), 3, 9, "gaussian")
    assert sp.expand(m2-e**2*((1-e)**2+2*e**4)) == 0
    assert sp.expand(m3-2*e**5*(3*(1-e)**2+4*e**4)) == 0
    remainder = sp.div(m3, m2, e)[1]
    assert remainder != 0
    return reports, {"case": "gaussian_local_not_global",
                     "M2": str(sp.factor(m2)), "M3": str(sp.factor(m3)),
                     "M3_remainder_mod_global_M2": str(remainder),
                     "global_prefix_equality": False,
                     "local_prefix_equality": True}


def verify_least_prefix(packet, source_instance):
    assert packet["source_instance"] == source_instance
    assert (packet["jet_order"], packet["prefix_N"], packet["nonmember_moment"]) == (7, 6, 7)
    parameters, perturbation, _, relations = backend.decode_instance(source_instance)
    order = packet["jet_order"]
    jets = backend.shear_moment_jets(perturbation, parameters, order, source_instance["functional"])
    dual = backend.decode_polynomial(packet["dual_coefficients"], parameters)
    dual_terms = {powers: coefficient for powers, coefficient
                  in sp.Poly(dual, *parameters, domain=sp.QQ).terms() if coefficient}
    assert all(sum(powers) <= order for powers in dual_terms)

    def pairing(expression):
        expression = backend.truncate(sp.expand(expression), parameters, order)
        return sum(coefficient*dual_terms.get(powers, 0)
                   for powers, coefficient in sp.Poly(expression, *parameters,
                                                      domain=sp.QQ).terms())

    assert pairing(jets.get(7, sp.Integer(0))) == 1
    generators = [jets.get(n, sp.Integer(0)) for n in range(1, 7)]+relations
    multipliers = backend.monomials(len(parameters), order)
    for generator in generators:
        for powers in multipliers:
            assert pairing(generator*backend.monomial(parameters, powers)) == 0
    return {"case": "sphere_three_parameter_profile_least_local_prefix",
            "jet_order": order, "separated_moment": 7, "smaller_prefix": 6,
            "target_pairing": "1", "ideal_generator_count": len(generators),
            "monomial_multipliers_per_generator": len(multipliers),
            "all_ideal_pairings": "0", "least_local_prefix": 7,
            "global_least_prefix_claimed": False}


def replay(path):
    packet = json.loads(path.read_text(encoding="utf-8"))
    assert packet["schema"] == "moment-jet-tail-examples-v1"
    assert packet["coefficient_field"] == "QQ"
    assert packet["fixed_base_phase"] == "W+Z*Y-Z^3/2"
    positives, negatives = literal_instances()
    assert {item["id"] for item in packet["positive"]} == set(positives)
    assert len(packet["positive"]) == len(positives)
    assert {item["id"] for item in packet["negative"]} == set(negatives)
    assert len(packet["negative"]) == len(negatives)
    reports, certificates = [], {}
    for entry in packet["positive"]:
        name, certificate = entry["id"], entry["certificate"]
        assert certificate["instance"] == positives[name]
        result = backend.verify_certificate(certificate)
        assert result["nilpotence_exponent_is_exact"] is True
        assert result["global_ideal_equality_claimed"] is False
        reports.append({"case": name, **result})
        certificates[name] = certificate
    negative_reports = [verify_negative(entry, *negatives[entry["id"]])
                        for entry in packet["negative"]]
    least_prefix = verify_least_prefix(packet["least_prefix_witness"],
                                       positives["sphere_three_parameter_profile_nu7"])
    corruption_reports = corruptions(certificates)
    corrupted_dual = deepcopy(packet["least_prefix_witness"])
    for term in corrupted_dual["dual_coefficients"]:
        term["coefficient"] = str(2*sp.Rational(term["coefficient"]))
    try:
        verify_least_prefix(corrupted_dual, positives["sphere_three_parameter_profile_nu7"])
    except AssertionError:
        corruption_reports.append({"case": "corrupted_least_prefix_dual",
                                   "rejected": True, "reason": "Target pairing changed from one"})
    else:
        raise AssertionError("Corrupted least-prefix dual was accepted")
    direct_reports, local_report = direct_controls()
    return {"status": "passed", "coefficient_field": "QQ",
            "default_mode": "Replay retained identities and dual annihilators; no certificate producer is used",
            "positive_certificates": reports, "negative_cases": negative_reports,
            "corrupted_packets_rejected": corruption_reports,
            "least_local_prefix_certificate": least_prefix,
            "independent_direct_controls": direct_reports,
            "local_versus_global_control": local_report,
            "scope": "Finite exact certificates plus the written mixed-cutoff/Nakayama proof; the generic conclusion is local/completed ideal equality.",
            "independent_replay": False, "formal_verification": False,
            "external_review": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate", action="store_true", help="Regenerate explicit retained witnesses from literal QQ inputs")
    parser.add_argument("--fixtures", type=Path, default=FIXTURES)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.generate:
        generate(args.fixtures)
    result = replay(args.fixtures)
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    result.update(python=platform.python_version(), sympy=sp.__version__,
                  proof_sha256=digest(PROOF), module_sha256=digest(Path(backend.__file__)),
                  checker_sha256=digest(Path(__file__)), fixture_sha256=digest(args.fixtures),
                  historical_input_context_sha256=digest(CONTEXT),
                  historical_input_context_role="Context only; no historical ideal-membership or nilpotence assertion is trusted")
    report = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(report, encoding="utf-8", newline="\n")
    print(report, end="")


if __name__ == "__main__":
    main()
