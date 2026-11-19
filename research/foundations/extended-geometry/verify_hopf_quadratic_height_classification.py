"""Exact finite controls for HOPF_LIFT_CLASSIFICATION.md Section8.

Formal uniqueness, the universal ordinary-ideal identities and the unrestricted
polynomial obstruction are proved in the written note, not by these fixtures.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import comb, factorial, prod
from pathlib import Path
import platform

if not __debug__:
    raise RuntimeError("This exact certificate requires assertions enabled")


def trim(a):
    answer = list(map(Q, a))
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer or [Q(0)]


def coefficient(a, n):
    return a[n] if 0 <= n < len(a) else Q(0)


def add(a, b, scale=1):
    return trim([coefficient(a, j)+scale*coefficient(b, j)
                 for j in range(max(len(a), len(b)))])


def multiply(a, b, order=None):
    length = len(a)+len(b)-1
    if order is not None:
        length = min(length, order+1)
    answer = [Q(0)]*length
    for i, av in enumerate(a):
        for j, bv in enumerate(b):
            if i+j < length:
                answer[i+j] += av*bv
    return trim(answer)


def power(a, n, order=None):
    answer = [Q(1)]
    for _ in range(n):
        answer = multiply(answer, a, order)
    return answer


def inverse_series(a, order):
    assert a[0]
    answer = [Q(1)/a[0]]
    for n in range(1, order+1):
        answer.append(-sum((coefficient(a, j)*answer[n-j]
                            for j in range(1, n+1)), Q(0))/a[0])
    return answer


def formal_solution(C, order):
    B = [Q(1-i)*c for i, c in enumerate(C)]
    numerator = power(C, 3, order)
    denominator_inverse = power(inverse_series(B, order), 2, order)
    return [-c for c in multiply(numerator, denominator_inverse, order)]


def moment(C, D, m, marker_power=0):
    index = m-marker_power
    if index < 0:
        return Q(0)
    return sum((Q(comb(m, j), 2*j+1)*coefficient(
        multiply(power(C, m-j, index), power(D, j, index), index), index)
                for j in range(m+1)), Q(0))


def tau(c0, m):
    return c0**(m-1)*Q(2**(m-1)*factorial(m), prod(range(1, 2*m+2, 2)))


def triangular_solution(C, order):
    D = [-C[0]]+[Q(0)]*order
    for m in range(1, order+1):
        residual = moment(C, D, m)
        D[m] = -residual/tau(C[0], m)
        assert moment(C, D, m) == 0
    return trim(D)


def polynomial_division(numerator, denominator):
    numerator, denominator = trim(numerator), trim(denominator)
    assert denominator != [0]
    quotient = [Q(0)]*max(1, len(numerator)-len(denominator)+1)
    while numerator != [0] and len(numerator) >= len(denominator):
        shift = len(numerator)-len(denominator)
        value = numerator[-1]/denominator[-1]
        quotient[shift] = value
        subtract = [Q(0)]*shift+[value*c for c in denominator]
        numerator = add(numerator, subtract, -1)
    return trim(quotient), numerator


def require_zero(value):
    assert value == 0


def rejects(check):
    try:
        check()
    except AssertionError:
        return True
    return False


def delta_moment_polynomial(C, star, m):
    """Actual f_m in Q[delta1,...,delta_m], by (x,t,delta) multiplication."""
    zero = (0,)*m
    terms = []
    for j in range(m+1):
        if coefficient(C,j):
            terms.append((j,0,zero,coefficient(C,j)))
        if coefficient(star,j):
            terms.append((j,2,zero,coefficient(star,j)))
        if j:
            variable = tuple(int(i==j-1) for i in range(m))
            terms.append((j,2,variable,Q(1)))
    current = {(0,0,zero): Q(1)}
    for _ in range(m):
        updated = {}
        for (xdegree,tdegree,powers), value in current.items():
            for j,t,p,c in terms:
                if xdegree+j > m:
                    continue
                new_powers = tuple(a+b for a,b in zip(powers,p))
                key = xdegree+j,tdegree+t,new_powers
                updated[key] = updated.get(key,Q(0))+value*c
        current = {key:value for key,value in updated.items() if value}
    answer = {}
    for (xdegree,tdegree,powers),value in current.items():
        if xdegree==m:
            answer[powers] = answer.get(powers,Q(0))+value/Q(tdegree+1)
    return {key:value for key,value in answer.items() if value}


def delta_triangular_certificate(polynomial, c0, m):
    unit = tuple(int(i==m-1) for i in range(m))
    remainder = dict(polynomial)
    remainder[unit] = remainder.get(unit,Q(0))-tau(c0,m)
    remainder = {p:c for p,c in remainder.items() if c}
    multipliers = [{} for _ in range(m-1)]
    for powers,value in remainder.items():
        index = next((i for i in range(m-1) if powers[i]), None)
        assert index is not None, "Moment residual is not in the earlier delta ideal"
        quotient = list(powers)
        quotient[index] -= 1
        quotient = tuple(quotient)
        multipliers[index][quotient] = multipliers[index].get(quotient,Q(0))+value
    reconstructed = {unit:tau(c0,m)}
    for index,multiplier in enumerate(multipliers):
        for powers,value in multiplier.items():
            product = list(powers)
            product[index] += 1
            product = tuple(product)
            reconstructed[product] = reconstructed.get(product,Q(0))+value
    assert {p:c for p,c in reconstructed.items() if c}==polynomial
    encode = lambda p:[{"delta_powers":list(powers),"coefficient":str(c)}
                        for powers,c in sorted(p.items()) if c]
    return {"moment_index":m,"unit_coefficient_tau":str(tau(c0,m)),
            "actual_moment_polynomial":encode(polynomial),
            "earlier_delta_multipliers":[encode(p) for p in multipliers]}


def run():
    order = 10
    fixtures = [[1], [1,1], [2,-3], [1,2,3], [2,-3,0,4],
                [1,0,0,0,5], [3,1,-2,1]]
    records = []
    triangular_coefficient_controls = 0
    for raw in fixtures:
        C = trim(raw)
        star = formal_solution(C, order)
        reconstructed = triangular_solution(C, order)
        assert trim(star) == reconstructed
        assert all(moment(C, star, m) == 0 for m in range(1, order+1))
        # Coefficient-level check after all earlier deltas have vanished.
        for m in range(1, order+1):
            past = [coefficient(star, j) for j in range(m+1)]
            for u in (-2, 1, 3):
                perturbed = past[:]
                perturbed[m] += u
                assert moment(C, perturbed, m) == tau(C[0], m)*u
                triangular_coefficient_controls += 1
        B = [Q(1-i)*c for i,c in enumerate(C)]
        quotient, remainder = polynomial_division(power(C, 3), power(B, 2))
        affine = len(C) <= 2
        assert (remainder == [0]) == affine
        if affine:
            D = [-c for c in quotient]
            assert trim(D) == trim(star)
            c0, c1 = C[0], coefficient(C, 1)
            for m in range(1, order+1):
                assert moment(C, D, m) == 0
                kappa = Q(2**m*factorial(m), prod(range(1, 2*m+2, 2)))
                expected = c0**m*(c1/c0)**(m-1)*kappa
                assert moment(C, D, m, 1) == expected
        altered = [coefficient(star, j) for j in range(order+1)]
        altered[order] += 1
        assert rejects(lambda: require_zero(moment(C, altered, order)))
        records.append({"C_coefficients": [str(c) for c in C],
                        "D_star_coefficients_through_order": [str(c) for c in star],
                        "triangular_solution_equal": True,
                        "checked_moment_orders": list(range(1, order+1)),
                        "polynomial_division_remainder": [str(c) for c in remainder],
                        "affine_polynomial_solution": affine,
                        "altered_tenth_coefficient_rejected": True})

    # The singular branch is handled by exact first/second moment identities.
    singular_controls = 0
    for d in range(-4,5):
        c = -Q(d, 3)
        C, D = [Q(0), c], [Q(0), Q(d)]
        assert moment(C,D,1) == 0
        assert moment(C,D,2) == Q(4,45)*d*d
        assert (moment(C,D,2) == 0) == (d == 0)
        singular_controls += 1
    for C,D in [([0,0,1,2],[0,0,-3,0,4]), ([0],[0,0,2]),
                ([0,0,0,3],[0,0,1,-2])]:
        C,D=trim(C),trim(D)
        assert all(moment(C,D,m) == 0 for m in range(1,order+1))
        singular_controls += 1

    # Finite coefficient controls of E=B^2(D-D*) with unrestricted delta values.
    ideal_controls = 0
    for C_raw in fixtures:
        C=trim(C_raw)
        star=formal_solution(C, order)
        B=[Q(1-i)*c for i,c in enumerate(C)]
        D=[-C[0]]+[Q((-1)**j*(j+1),j+2) for j in range(1,order+1)]
        delta=add(D,star,-1)
        E=add(multiply(power(B,2),D,order),power(C,3,order))
        rhs=multiply(power(B,2),delta,order)
        assert E==rhs
        for m in range(1,order+1):
            triangular = C[0]**2*coefficient(delta,m)
            triangular += sum((coefficient(power(B,2),j)*coefficient(delta,m-j)
                               for j in range(1,m)), Q(0))
            assert coefficient(E,m)==triangular
            ideal_controls += 1
    prefix_certificates = []
    for raw in ([1,2,3],[2,-3,0,4]):
        C = trim(raw)
        star = formal_solution(C,6)
        witnesses = []
        for m in range(1,7):
            actual = delta_moment_polynomial(C,star,m)
            witnesses.append(delta_triangular_certificate(actual,C[0],m))
        corrupted = dict(actual)
        corrupted[(0,)*6] = Q(1)
        assert rejects(lambda:delta_triangular_certificate(corrupted,C[0],6))
        prefix_certificates.append({"C_coefficients":[str(c) for c in C],
                                    "D_star_coefficients":[str(c) for c in star],
                                    "ordinary_prefix_ideal_witnesses":witnesses,
                                    "injected_constant_residual_rejected":True})
    return {"status":"passed", "python":platform.python_version(),
            "arithmetic":"standard-library Fraction, exact univariate polynomial/series arithmetic",
            "formal_order":order, "fixtures":records,
            "triangular_new_coefficient_controls":triangular_coefficient_controls,
            "singular_branch_controls":singular_controls,
            "E_delta_triangular_coefficient_controls":ideal_controls,
            "actual_ordinary_prefix_ideal_certificates":prefix_certificates,
            "scope":"The listed finite fixtures and orders, including twelve explicit Q[delta] triangular ordinary-ideal certificates. Universal formal uniqueness, full ordinary-ideal equality, explicit finite prefix and unrestricted polynomial classification use the written proofs in Section8.",
            "not_claimed":["radical ordinary moment ideal", "nilpotence index",
                           "full primitive degree-five classification", "formal verification",
                           "external human review", "historical priority"]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    here=Path(__file__).resolve().parent
    root=here.parents[2]
    result=run()
    paths=[Path(__file__).resolve(),here/"HOPF_LIFT_CLASSIFICATION.md"]
    result["source_sha256"]={p.relative_to(root).as_posix():sha256(p.read_bytes()).hexdigest()
                             for p in paths}
    output=args.output or here/"hopf_quadratic_height_classification_verification.json"
    output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"fixtures":len(result["fixtures"]),
                      "order":result["formal_order"],"output":str(output)}))


if __name__=="__main__":
    main()
