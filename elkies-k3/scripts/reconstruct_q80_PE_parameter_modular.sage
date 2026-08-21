#!/usr/bin/env sage
"""Recover exact centered P(u),E(u) on the rational q80 surface branch.

The exact quartic ideal has a five-element lexicographic Groebner basis with
one equation monic linear in P and two equations linear in E.  Substitution of
the certified D(u),Q(u) map creates large extraneous factors.  We cancel those
modulo deterministic 61-bit primes, CRT/rational-reconstruct the reduced
degree-10 P and degree-15 E functions, and then verify all five Groebner-basis
relations as exact polynomial identities over QQ.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


parser = argparse.ArgumentParser()
parser.add_argument(
    "--ideal",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-quartic-ideal.json",
)
parser.add_argument(
    "--dq-parameter",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-DQ-parameter.json",
)
parser.add_argument(
    "--output",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-PDQE-parameter.json",
)
parser.add_argument("--maximum-primes", type=int, default=1000)
parser.add_argument("--reconstruct-every", type=int, default=50)
parser.add_argument("--minimum-primes", type=int, default=650)
arguments = parser.parse_args()


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


ideal_path = Path(arguments.ideal)
dq_path = Path(arguments.dq_parameter)
ideal_payload = json.loads(ideal_path.read_text())
dq_payload = json.loads(dq_path.read_text())
if ideal_payload.get("schema") != "q80-cm24-qq-quartic-ideal-v1":
    raise ValueError("unexpected quartic-ideal schema")
if dq_payload.get("schema") != "q80-cm24-qq-DQ-parameter-v1":
    raise ValueError("unexpected D,Q parameter schema")

source_ring = PolynomialRing(QQ, names=("D", "P", "Q", "E"))
D0, P0, Q0, E0 = source_ring.gens()
lex_ring = PolynomialRing(QQ, names=("P", "E", "D", "Q"), order="lex")
P, E, D, Q = lex_ring.gens()
relations = tuple(
    lex_ring(source_ring(value)(D=D, P=P, Q=Q, E=E))
    for value in ideal_payload["relations"]
)
groebner_basis = tuple(lex_ring.ideal(relations).groebner_basis())
if len(groebner_basis) != 5:
    raise ArithmeticError("expected a five-element lex Groebner basis")
if groebner_basis[0].lm() != P:
    raise ArithmeticError("first lex relation is not monic linear in P")

dq_ring = PolynomialRing(QQ, names=("D", "Q"))
d, q = dq_ring.gens()


def coefficient_in_PE(polynomial, p_exponent, e_exponent):
    return dq_ring(
        sum(
            coefficient*d**exponents[2]*q**exponents[3]
            for exponents, coefficient in polynomial.dict().items()
            if exponents[0] == p_exponent and exponents[1] == e_exponent
        )
    )


P_relation = groebner_basis[0]
E_relation = groebner_basis[2]
E_coefficient = coefficient_in_PE(E_relation, 0, 1)
E_constant = coefficient_in_PE(E_relation, 0, 0)
P_E_coefficient = coefficient_in_PE(P_relation, 0, 1)
P_constant = coefficient_in_PE(P_relation, 0, 0)
if any(
    exponents[0] or exponents[1] > 1
    for exponents in E_relation.dict()
):
    raise ArithmeticError("selected E relation is not linear in E")
if coefficient_in_PE(P_relation, 1, 0) != 1:
    raise ArithmeticError("selected P relation is not monic in P")

parameter_ring = PolynomialRing(QQ, "u")
u = parameter_ring.gen()
D_pair = tuple(
    parameter_ring(dq_payload["functions"]["D"][key])
    for key in ("numerator", "denominator")
)
Q_pair = tuple(
    parameter_ring(dq_payload["functions"]["Q"][key])
    for key in ("numerator", "denominator")
)

graph_polynomials = (
    E_coefficient, E_constant, P_E_coefficient, P_constant
)
all_denominators = tuple(
    {
        coefficient.denominator()
        for polynomial in graph_polynomials+D_pair+Q_pair
        for coefficient in polynomial.coefficients()
    }
)
candidate_primes = []
prime = next_prime(ZZ(2)**60 + 7654321)
while len(candidate_primes) < arguments.maximum_primes:
    if all(denominator % prime for denominator in all_denominators):
        candidate_primes.append(ZZ(prime))
    prime = next_prime(prime+2)
full_modulus = prod(candidate_primes)


def rational_residue(coefficient):
    return (
        coefficient.numerator()
        * inverse_mod(coefficient.denominator(), full_modulus)
    ) % full_modulus


def pre_reduce_univariate(polynomial):
    return tuple(
        rational_residue(polynomial[index])
        for index in range(polynomial.degree()+1)
    )


def pre_reduce_bivariate(polynomial):
    return {
        "degrees": (
            polynomial.degree(polynomial.parent().gen(0)),
            polynomial.degree(polynomial.parent().gen(1)),
        ),
        "coefficients": {
            exponents: rational_residue(coefficient)
            for exponents, coefficient in polynomial.dict().items()
        },
    }


D_pair_mod_full = tuple(pre_reduce_univariate(value) for value in D_pair)
Q_pair_mod_full = tuple(pre_reduce_univariate(value) for value in Q_pair)
graph_mod_full = tuple(
    pre_reduce_bivariate(value) for value in graph_polynomials
)


def evaluate_bivariate_mod_prime(
    record, prime, finite, finite_parameter_ring, finite_D_pair, finite_Q_pair,
):
    D_degree, Q_degree = record["degrees"]
    numerator = finite_parameter_ring.zero()
    for exponents, coefficient in record["coefficients"].items():
        numerator += (
            finite(coefficient % prime)
            * finite_D_pair[0]**exponents[0]
            * finite_D_pair[1]**(D_degree-exponents[0])
            * finite_Q_pair[0]**exponents[1]
            * finite_Q_pair[1]**(Q_degree-exponents[1])
        )
    denominator = finite_D_pair[1]**D_degree*finite_Q_pair[1]**Q_degree
    return numerator, denominator


expected_degrees = {"P": 10, "E": 15}
crt_modulus = ZZ.one()
residues = {
    "P_numerator": [ZZ.zero()]*(expected_degrees["P"]+1),
    "P_denominator": [ZZ.zero()]*(expected_degrees["P"]+1),
    "E_numerator": [ZZ.zero()]*(expected_degrees["E"]+1),
    "E_denominator": [ZZ.zero()]*(expected_degrees["E"]+1),
}
accepted_primes = []
bad_primes = []
reconstructed = None


def reconstruct_residues():
    result = {}
    for name, values in residues.items():
        result[name] = parameter_ring(
            [value.rational_reconstruction(crt_modulus) for value in values]
        )
    return result


def normalize_mod_prime(numerator, denominator):
    common = numerator.gcd(denominator)
    numerator //= common
    denominator //= common
    scale = ~denominator.leading_coefficient()
    return numerator*scale, denominator*scale, common.degree()


for prime in candidate_primes:
    finite = GF(prime)
    finite_parameter_ring = PolynomialRing(finite, "u")
    try:
        finite_D_pair = tuple(
            finite_parameter_ring(
                [finite(value % prime) for value in polynomial]
            )
            for polynomial in D_pair_mod_full
        )
        finite_Q_pair = tuple(
            finite_parameter_ring(
                [finite(value % prime) for value in polynomial]
            )
            for polynomial in Q_pair_mod_full
        )
        E_A_num, E_A_den = evaluate_bivariate_mod_prime(
            graph_mod_full[0], prime, finite, finite_parameter_ring,
            finite_D_pair, finite_Q_pair,
        )
        E_B_num, E_B_den = evaluate_bivariate_mod_prime(
            graph_mod_full[1], prime, finite, finite_parameter_ring,
            finite_D_pair, finite_Q_pair,
        )
        E_num, E_den, E_cancel_degree = normalize_mod_prime(
            -E_B_num*E_A_den, E_A_num*E_B_den
        )
        P_C_num, P_C_den = evaluate_bivariate_mod_prime(
            graph_mod_full[2], prime, finite, finite_parameter_ring,
            finite_D_pair, finite_Q_pair,
        )
        P_B_num, P_B_den = evaluate_bivariate_mod_prime(
            graph_mod_full[3], prime, finite, finite_parameter_ring,
            finite_D_pair, finite_Q_pair,
        )
        P_num, P_den, P_cancel_degree = normalize_mod_prime(
            -(P_C_num*E_num*P_B_den + P_B_num*P_C_den*E_den),
            P_C_den*E_den*P_B_den,
        )
    except (ZeroDivisionError, ArithmeticError):
        bad_primes.append({"prime": int(prime), "reason": "denominator_or_gcd_failure"})
        continue
    if (
        (P_num.degree(), P_den.degree()) != (10, 10)
        or (E_num.degree(), E_den.degree()) != (15, 15)
    ):
        bad_primes.append(
            {
                "prime": int(prime),
                "reason": "degree_drop",
                "P_degrees": [int(P_num.degree()), int(P_den.degree())],
                "E_degrees": [int(E_num.degree()), int(E_den.degree())],
            }
        )
        continue

    inverse_old_modulus = inverse_mod(crt_modulus % prime, prime)
    modular_pairs = {
        "P_numerator": P_num,
        "P_denominator": P_den,
        "E_numerator": E_num,
        "E_denominator": E_den,
    }
    for name, polynomial in modular_pairs.items():
        for index in range(len(residues[name])):
            value = ZZ(polynomial[index])
            correction = (
                (value-residues[name][index]) % prime*inverse_old_modulus
            ) % prime
            residues[name][index] += crt_modulus*correction
    crt_modulus *= prime
    accepted_primes.append(prime)
    accepted_count = len(accepted_primes)
    if (
        accepted_count < arguments.minimum_primes
        or accepted_count % arguments.reconstruct_every
    ):
        continue
    try:
        candidate = reconstruct_residues()
    except ArithmeticError:
        print(
            "Q80PEMOD|stage=reconstruct|"
            f"accepted_primes={accepted_count}|modulus_bits={crt_modulus.nbits()}|"
            "status=insufficient_modulus",
            flush=True,
        )
        continue
    reconstructed = candidate
    print(
        "Q80PEMOD|stage=reconstruct|"
        f"accepted_primes={accepted_count}|modulus_bits={crt_modulus.nbits()}|"
        "status=candidate_available",
        flush=True,
    )
    break

if reconstructed is None:
    raise ArithmeticError("P,E functions did not rationally reconstruct")


def primitive_pair(numerator, denominator):
    coefficient_denominator = lcm(
        coefficient.denominator()
        for polynomial in (numerator, denominator)
        for coefficient in polynomial.coefficients()
    )
    numerator = parameter_ring(coefficient_denominator*numerator)
    denominator = parameter_ring(coefficient_denominator*denominator)
    content = gcd(
        tuple(
            ZZ(coefficient)
            for polynomial in (numerator, denominator)
            for coefficient in polynomial.coefficients()
        )
    )
    numerator //= content
    denominator //= content
    if denominator.leading_coefficient() < 0:
        numerator = -numerator
        denominator = -denominator
    return numerator, denominator


P_pair = primitive_pair(
    reconstructed["P_numerator"], reconstructed["P_denominator"]
)
E_pair = primitive_pair(
    reconstructed["E_numerator"], reconstructed["E_denominator"]
)

# Exact global verification against the full lex Groebner basis.  Each
# rational substitution is cross-multiplied over one common denominator.
function_pairs = (P_pair, E_pair, D_pair, Q_pair)


def exact_relation_residual(relation):
    variable_degrees = tuple(
        relation.degree(variable) for variable in lex_ring.gens()
    )
    powers = []
    for (numerator, denominator), maximum_degree in zip(
        function_pairs, variable_degrees
    ):
        powers.append(
            tuple(
                numerator**exponent*denominator**(maximum_degree-exponent)
                for exponent in range(maximum_degree+1)
            )
        )
    residual = parameter_ring.zero()
    for exponents, coefficient in relation.dict().items():
        term = parameter_ring(coefficient)
        for coordinate, exponent in enumerate(exponents):
            term *= powers[coordinate][exponent]
        residual += term
    return residual


exact_residuals = tuple(
    exact_relation_residual(relation) for relation in groebner_basis
)
if any(exact_residuals):
    first = next(index for index, value in enumerate(exact_residuals) if value)
    raise ArithmeticError(f"reconstructed parameter misses Groebner relation {first}")

original_P_pair = primitive_pair(
    P_pair[0]+QQ(9)/4*P_pair[1], P_pair[1]
)
original_E_pair = primitive_pair(
    E_pair[0]-QQ(27)/32*E_pair[1], E_pair[1]
)


def function_record(pair):
    return {
        "numerator": str(pair[0]),
        "denominator": str(pair[1]),
        "value": f"({pair[0]})/({pair[1]})",
        "degrees": [int(value.degree()) for value in pair],
    }


output_payload = {
    "schema": "q80-cm24-qq-PDQE-parameter-v1",
    "scope": "exact_global_parameterization_of_the_reconstructed_unmarked_surface_branch",
    "claim_boundary": {
        "proved": "exact substitution into the reconstructed coefficient ideal",
        "not_proved": [
            "globalization of the three marked Mordell-Weil sections",
            "identification with X(6,79)/<w_474>",
            "generic Mordell-Weil rank three",
        ],
    },
    "status": "exact_lex_groebner_basis_substitution",
    "slope": "8/87",
    "parameter": "u",
    "centered_functions": {
        "D": function_record(D_pair),
        "P": function_record(P_pair),
        "Q": function_record(Q_pair),
        "E": function_record(E_pair),
    },
    "original_functions": {
        "d": function_record((D_pair[0]-QQ(1)/2*D_pair[1], D_pair[1])),
        "p": function_record(original_P_pair),
        "q": function_record((Q_pair[0]-QQ(9)/4*Q_pair[1], Q_pair[1])),
        "e": function_record(original_E_pair),
    },
    "degrees": {"D": [5, 5], "P": [10, 10], "Q": [10, 10], "E": [15, 15]},
    "modular_reconstruction": {
        "accepted_primes": [int(value) for value in accepted_primes],
        "bad_primes": bad_primes,
        "crt_modulus_bits": int(crt_modulus.nbits()),
        "expected_generic_cancellations": {"P": 230, "E": 215},
    },
    "groebner_basis": [str(value) for value in groebner_basis],
    "checks": {
        "groebner_residuals": [str(value) for value in exact_residuals],
        "relations_checked": len(groebner_basis),
    },
    "inputs": [
        {"path": str(ideal_path), "sha256": sha256(ideal_path)},
        {"path": str(dq_path), "sha256": sha256(dq_path)},
    ],
}
output_path = Path(arguments.output)
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text(json.dumps(output_payload, indent=2, sort_keys=True)+"\n")

print(
    "Q80PEMOD|D_degrees=5,5|P_degrees=10,10|Q_degrees=10,10|"
    "E_degrees=15,15|"
    f"accepted_primes={len(accepted_primes)}|modulus_bits={crt_modulus.nbits()}|"
    f"groebner_residuals=0,0,0,0,0|output={output_path}|"
    "status=PASS_EXACT_PDQE_PARAMETERIZATION",
    flush=True,
)
