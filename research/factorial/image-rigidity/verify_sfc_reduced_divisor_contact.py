#!/usr/bin/env python3
"""Exact finite controls for the reduced-divisor contact bound.

The written proof supplies all-parameter quantifiers. This replay compares
derivatives computed from inverse series with independently formed rational
numerators and checks the exact factor removed from the previous formula.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from math import comb, factorial, prod
from pathlib import Path

if not __debug__:
    raise SystemExit("REJECT: optimized execution is not an approved replay")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / "research/factorial/image-rigidity"
NOTE = HERE / "SFC_REDUCED_DIVISOR_CONTACT_BOUND.md"
NOTE_PIN = "a35a390ae1b3bb84fa52d3dede3b488ef469536388cb620cd320574ed676594b"
PINS = {
    "SFC_MONOMIAL_LINEAR_CONTACT_BOUND.md": "409ca1ed8a1ef3218fb07673d14830d04cf10562f5714c88bd8de8e9a681b417",
    "verify_sfc_monomial_linear_contact.py": "c80ea0aa4ade828a89d18e028072c65fc0525c6700afc661183b5cb14c2e1055",
    "SFC_RAMIFIED_PULLBACK_CLOSURE.md": "eb444a8f02feb92834336589024e0d920b4eba23970582291af406c11bd7f519",
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def qjson(value):
    if isinstance(value, F):
        return {"numerator": value.numerator, "denominator": value.denominator}
    raise TypeError(type(value).__name__)


def polydiv(a, b):
    """Exact polynomial division with ascending Fraction coefficients."""
    remainder = list(map(F, a))
    divisor = list(map(F, b))
    quotient = [F(0)] * max(1, len(a) - len(b) + 1)
    require(divisor != [0], "Division by zero polynomial")
    while remainder != [0] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        value = remainder[-1] / divisor[-1]
        quotient[shift] += value
        for j, coefficient in enumerate(divisor):
            remainder[j+shift] -= value * coefficient
        while len(remainder) > 1 and remainder[-1] == 0:
            remainder.pop()
    while len(quotient) > 1 and quotient[-1] == 0:
        quotient.pop()
    return quotient, remainder


def polygcd(a, b):
    left, right = list(map(F, a)), list(map(F, b))
    while right != [0]:
        _, remainder = polydiv(left, right)
        left, right = right, remainder
    require(left[0] != 0, "GCD must be a unit at zero in these controls")
    return [value/left[0] for value in left]


def exactquotient(a, b):
    quotient, remainder = polydiv(a, b)
    require(remainder == [0], "Expected polynomial divisibility")
    return quotient


def check_case(p, A, a, Eold, D, J, E, metadata, helper):
    add, sub, mul = helper.add, helper.sub, helper.mul
    scale, power, der = helper.scale, helper.power, helper.derivative
    S, T = len(p) - 1, len(J) - 1
    require(len(D) - 1 == S and D[0] == 1, "D degree/normalization")
    require(E == [0] if T == 0 else len(E) - 1 <= T - 1, "Reduced logarithmic numerator degree")
    require(mul(J, sub(mul(der(A), p), mul(A, der(p)))) == mul(E, mul(A, p)), "Reduced logarithmic identity")
    G = polygcd(D, der(D))
    K = exactquotient(D, G)
    U = exactquotient(der(D), G)
    L = polygcd(J, K)
    J0 = exactquotient(J, L)
    KoverL = exactquotient(K, L)
    d, T0 = len(K) - 1, len(J0) - 1
    require(1 <= d <= S and len(U)-1 == d-1, "Radical quotient degree")
    require(polygcd(K, der(K)) == [1], "K is squarefree")
    maximum = 16
    W = mul(A, helper.reciprocal(p, maximum), maximum)
    H = [0] + [F(value, j + 1) for j, value in enumerate(W)]
    B = helper.compose(H, helper.inverse(a, maximum + 1), maximum + 1)
    if metadata.get("exponents") is not None:
        for m in range(0, 9):
            actual = helper.actual_moment(metadata["exponents"], metadata["coefficients"], m)
            factor = factorial(m) * prod(factorial(m*r) for r in metadata["exponents"])
            require(actual == factor * (m + 1) * helper.coeff(B, m + 1), "Direct actual factorial expansion")
    if "cyclic_order" in metadata:
        q = metadata["cyclic_order"]
        for m in range(0, maximum):
            expected = 0 if m % q else comb((q + 1) * (m // q), m // q)
            require((m + 1) * helper.coeff(B, m + 1) == expected, "Cyclic inverse control")
    quotient = mul(p, helper.reciprocal(J, S), S)
    require(mul(quotient, J) == p, "P/J is polynomial")
    N = scale(der(D), -1)
    M = N[:]
    C = scale(U, -1)
    Chat = C[:]
    Bder = der(B)
    rows = []
    for k in range(2, 7):
        Bder = der(Bder)
        truncation = 7
        denominator = mul(power(W, k-1, truncation),
                          mul(power(D, 2*k-1, truncation), power(J, k-2, truncation), truncation),
                          truncation)
        evaluated = mul(helper.compose(Bder, a, truncation), denominator, truncation)
        require(evaluated == helper.trim(N[:truncation+1]), "Derivative identity from inverse series")
        upper = S - 1 + (k-2)*(S+T-1)
        require(N != [0] and len(N)-1 <= upper, "Reduced numerator degree")
        require(M == mul(power(quotient, k-2), N), "Exact removable factor in old numerator")
        radical_denominator = mul(power(W, k-1, truncation),
                                  mul(power(D, k, truncation),
                                      mul(power(K, k-1, truncation), power(J, k-2, truncation), truncation),
                                      truncation), truncation)
        evaluated_radical = mul(helper.compose(Bder, a, truncation), radical_denominator, truncation)
        require(evaluated_radical == helper.trim(C[:truncation+1]), "Radical derivative identity from inverse series")
        radical_upper = d-1+(k-2)*(d+T-1)
        require(C != [0] and len(C)-1 <= radical_upper, "Radical numerator degree")
        require(N == mul(power(G, k-1), C), "Exact repeated-critical-root factor")
        final_denominator = mul(power(W, k-1, truncation),
                                mul(power(D, k, truncation),
                                    mul(power(K, k-1, truncation), power(J0, k-2, truncation), truncation),
                                    truncation), truncation)
        evaluated_final = mul(helper.compose(Bder, a, truncation), final_denominator, truncation)
        require(evaluated_final == helper.trim(Chat[:truncation+1]), "Overlap-reduced derivative identity from inverse series")
        final_upper = d-1+(k-2)*(d+T0-1)
        require(Chat != [0] and len(Chat)-1 == final_upper, "Exact final numerator degree")
        require(C == mul(power(L, k-2), Chat), "Exact shared-divisor factor")
        rows.append({"k": k, "actual_degree": len(N)-1, "upper_degree": upper, "N": N,
                     "radical_actual_degree":len(C)-1, "radical_upper_degree":radical_upper, "C":C,
                     "final_actual_degree":len(Chat)-1, "final_upper_degree":final_upper, "C_hat":Chat})
        bracket = add(add(scale(mul(E, D), k-1), scale(mul(J, der(D)), 2*k-1)),
                      scale(mul(der(J), D), k-2))
        Nnext = sub(mul(mul(J, D), der(N)), mul(bracket, N))
        oldbracket = sub(sub(mul(der(p), D), scale(mul(Eold, D), k-1)),
                         scale(mul(p, der(D)), 2*k-1))
        M = add(mul(mul(p, D), der(M)), mul(oldbracket, M))
        radical_bracket = add(add(scale(mul(E,K), k-1), scale(mul(J,U), k)),
                              add(scale(mul(J,der(K)), k-1), scale(mul(der(J),K), k-2)))
        C = sub(mul(mul(J,K),der(C)),mul(radical_bracket,C))
        final_bracket = add(add(scale(mul(E,KoverL), k-1), scale(mul(J0,U), k)),
                            add(scale(mul(J0,der(K)), k-1), scale(mul(der(J0),K), k-2)))
        Chat = sub(mul(mul(J0,K),der(Chat)),mul(final_bracket,Chat))
        N = Nnext
    return {**metadata, "S": S, "T": T, "P": p, "A": A, "D": D, "J": J, "E": E,
            "d":d, "K":K, "G":G, "U":U, "L":L, "J0":J0, "T0":T0,
            "numerators": rows, "derivative_truncation_degree": 7,
            "contact_bound_coefficients": [d+T0, 1-T0],
            "window_bound_coefficients": [d+T0-1, 1-T0]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", type=Path, default=NOTE)
    parser.add_argument("--output", type=Path, default=HERE / "sfc_reduced_divisor_contact_verification.json")
    args = parser.parse_args()
    require(sha256(args.note.read_bytes()).hexdigest() == NOTE_PIN, "Note pin mismatch")
    for name, expected in PINS.items():
        require(sha256((SOURCE/name).read_bytes()).hexdigest() == expected, "Dependency pin mismatch: "+name)
    spec = importlib.util.spec_from_file_location("retained_contact_arithmetic", SOURCE / "verify_sfc_monomial_linear_contact.py")
    h = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(h)
    cases = [((1,), (1,)), ((4,), (-2,)), ((1,1,2), (1,2,-2)),
             ((1,1,7), (1,-2,3)), ((1,1,1,2), (1,2,-1,-3)),
             ((1,1,1,1), (1,-1,2,-2)), ((0,2,3), (1,1,-2)),
             ((0,0,2), (1,1,-2)), ((0,1,2), (1,-1,2)), ((2,3), (2,-3)),
             ((1,1,1), (1,1,1)), ((1,2), (1,1))]
    records = []
    for exponents, coefficients in cases:
        p,A,a,_,Eold,D = h.build(exponents, coefficients)
        groups = defaultdict(lambda: [0,0])
        for r,c in zip(exponents,coefficients):
            groups[c][0] += 1
            groups[c][1] += r
        defects = [(c,R-pcount) for c,(pcount,R) in groups.items() if R != pcount]
        J = [1]
        for c,_ in defects:
            J = h.mul(J,[1,-c])
        E = [0]
        for i,(c,delta) in enumerate(defects):
            term = [1]
            for j,(other,_) in enumerate(defects):
                if j != i:
                    term = h.mul(term,[1,-other])
            E = h.add(E,h.scale(term,-delta*c))
        records.append(check_case(p,A,a,Eold,D,J,E,
                                 {"exponents":exponents,"coefficients":coefficients},h))
    for q in (2,3,5):
        p = [1]+[0]*(q-1)+[-1]
        A,a = p[:],[0]+p
        Eold = h.derivative(p)
        D = h.add(p,[0]+Eold)
        records.append(check_case(p,A,a,Eold,D,[1],[0],{"cyclic_order":q},h))
    primitive = records[2]
    require(primitive["S"] == 3 and primitive["T"] == 1, "Primitive profile bookkeeping")
    require(primitive["numerators"][1]["actual_degree"] > 2*(3-1),
            "False zero-defect degree estimate must be rejected")
    require(records[10]["S"] == 3 and records[10]["d"] == 2 and records[10]["T"] == 0,
            "Actual repeated-D neutral control")
    require(records[11]["T"] == 1 and records[11]["T0"] == 0,
            "Actual J/K-overlap control")
    receipt = {
        "status":"PASS", "note_sha256":NOTE_PIN, "dependencies":PINS,
        "scope":"Finite exact identity controls only; the written induction proves the unrestricted contact bound.",
        "cases":records, "parameter_cases":len(records),
        "inverse_derivative_identities":5*len(records),
        "exact_old_numerator_factor_identities":5*len(records),
        "radical_inverse_derivative_identities":5*len(records),
        "exact_repeated_critical_root_factor_identities":5*len(records),
        "overlap_reduced_inverse_derivative_identities":5*len(records),
        "exact_shared_divisor_factor_identities":5*len(records),
        "direct_actual_moment_checks":9*len(cases),
        "cyclic_coefficient_checks":16*3,
        "negative_control":"The (1,1,2) third-derivative numerator exceeds the false T=0 degree bound.",
        "sharp_SFC_claimed":False,
    }
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True,default=qjson)+"\n")
    print(f"PASS: {len(records)} parameter cases; {15*len(records)} inverse-derivative identities; "
          f"{15*len(records)} exact removed factors; {9*len(cases)} direct factorial moments; 48 cyclic coefficients.")
    print("Receipt:",args.output)


if __name__ == "__main__":
    try:
        main()
    except (ValueError,OSError) as error:
        raise SystemExit("REJECT: "+str(error)) from error
