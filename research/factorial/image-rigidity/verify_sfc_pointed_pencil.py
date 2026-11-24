#!/usr/bin/env python3
"""Exact finite controls for the pointed-pencil note; not a gap theorem."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import sympy as sp


NOTE_SHA256 = '081d8d21cd67ded8188ffd0f6246727895da3715ef07b17f29537df7be57ea24'
DEPENDENCY_SHA256 = {'SFC_RAMIFIED_PULLBACK_CLOSURE.md': 'eb444a8f02feb92834336589024e0d920b4eba23970582291af406c11bd7f519', 'SFC_MONOMIAL_LINEAR_CONTACT_BOUND.md': '409ca1ed8a1ef3218fb07673d14830d04cf10562f5714c88bd8de8e9a681b417', 'SFC_ALGEBRAIC_UNIFORM_WINDOWS.md': '5066e9fef20d387d80bf865cb07f19563bb5907ace2c8e5c0032b29028f3b2f4'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def zero(expression):
    return sp.cancel(sp.expand(expression)) == 0


def order(expression, variable):
    polynomial = sp.Poly(expression, variable)
    return min(monomial[0] for monomial, value in polynomial.terms() if value)


def main():
    if not __debug__:
        raise SystemExit("Refusing to run with assertions disabled (-O).")
    parser = argparse.ArgumentParser()
    parser.add_argument("--note", type=Path, default=Path(__file__).with_name(
        "SFC_POINTED_PENCIL_AND_CAUSTIC_DESCENT.md"))
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name(
        "sfc_pointed_pencil_verification.json"))
    args = parser.parse_args()
    note = args.note
    require(note.is_file(), "Missing proof note: " + str(note))
    require(hashlib.sha256(note.read_bytes()).hexdigest() == NOTE_SHA256,
            "Proof note SHA-256 mismatch")
    for name, expected in DEPENDENCY_SHA256.items():
        dependency = Path(__file__).with_name(name)
        require(dependency.is_file(), "Missing pinned source: " + name)
        require(hashlib.sha256(dependency.read_bytes()).hexdigest() == expected,
                "Pinned source SHA-256 mismatch: " + name)
    z, t, eps = sp.symbols("z t eps")
    r = sp.Rational
    cases = [
        ("constant_derivative_boundary", z-z*z, z, True, 1),
        ("primitive_r112", z*(1-z)*(1-2*z)*(1+2*z)**2,
         z+z*z, True, 1),
        ("quadratic_pullback", z*(1-z)*(1-2*z)**2,
         z-z*z, True, 2),
        ("misaligned_negative_control", z-z*z/5-3*z**3/25-3*z**4/50,
         z-z*z/2, False, 1),
    ]
    receipt = {"status": "PASS", "domain": "exact rational polynomials",
               "scope": "four finite pencil controls, not a sharp-gap proof",
               "sympy_version": sp.__version__, "cases": []}
    for name, a, H, aligned, exponent in cases:
        a, H = sp.expand(a), sp.expand(H)
        W = sp.diff(H, z)
        D, rem = sp.div(sp.diff(a, z), W, z)
        require(zero(rem), name + ": polynomial critical quotient")
        require(a.subs(z, 0) == H.subs(z, 0) == 0, name + ": origin")
        require(sp.diff(a, z).subs(z, 0) == W.subs(z, 0) == 1,
                name + ": normalized derivatives")
        d, S = int(sp.degree(a, z)), int(sp.degree(D, z))
        require(d - sp.degree(H, z) == S, name + ": degree difference")
        require(zero(sp.rem(a, W, z)) == aligned, name + ": alignment control")
        T = sp.expand(a-D*H)
        require(zero(sp.diff(T, z)+H*sp.diff(D, z)), name + ": dual identity")
        f = a+eps*H-t
        fixed = sp.factor(sp.resultant(f, W, z))
        moving = sp.factor(sp.resultant(f, D+eps, z))
        full = sp.resultant(f, sp.diff(f, z), z)
        require(zero(full-fixed*moving), name + ": resultant product")
        scalar = (-1)**(S*(d+1))*sp.LC(sp.Poly(D,z))**d
        require(sp.degree(moving,t) == S, name + ": moving degree")
        require(sp.LC(sp.Poly(moving,t)) == scalar, name + ": norm scalar")
        factors = sp.factor_list(moving, t, eps)[1]
        require(len(factors) == 1 and factors[0][1] == exponent,
                name + ": exact rational factor exponent")
        reduced = sp.Poly(factors[0][0],t).monic().as_expr()
        require(zero(moving-scalar*reduced**exponent), name + ": reduced power")
        require(zero(reduced.subs({t:T, eps:-D})), name + ": caustic equation")
        require(zero((sp.diff(reduced, eps)+H*sp.diff(reduced,t)).subs(
            {t:T, eps:-D})), name + ": derivative reconstruction")
        record = {"name":name, "a":str(a), "H":str(H), "D":str(D), "T":str(T),
                  "aligned":aligned, "degree_a":d, "degree_D":S,
                  "fixed_resultant":str(fixed), "moving_resultant":str(moving),
                  "monic_reduced_caustic":str(reduced), "scalar":str(scalar),
                  "power_exponent":exponent}

        if name == "primitive_r112":
            require(math.gcd(d,int(sp.degree(H,z))) == 1,
                    "primitive degree obstruction")
            require(zero(fixed-8*(eps+4*t)), "primitive fixed line")
            require(order(1-D,z) == 2 and order(T,z) == 3,
                    "primitive marked cusp orders")
            record["geometric_irreducibility_reason"] = (
                "proved caustic-field proposition and gcd(deg a,deg H)=1")
        elif name == "quadratic_pullback":
            require(zero(a-H*(1-4*H)), "pullback a identity")
            require(zero(D-(1-8*H)), "pullback quotient identity")
            require(zero(T-4*H*H), "pullback dual identity")
            require(zero(moving-16*(16*t-(eps+1)**2)**2), "pullback exact square")
        elif name == "misaligned_negative_control":
            b_a = a-r(3,10)*a*a
            contact = sp.expand(H-b_a)
            require(order(contact,z) == 5, "negative contact order")
            require(contact.coeff(z,5) == -r(27,1250), "negative leading contact")
            require(a.subs(z,1) == r(31,50), "negative misalignment")
            require(zero(fixed-(r(31,50)+eps/2-t)), "negative fixed line")
            record["contact_with_quadratic"] = {"order":5,"n_plus_S":4,
                "exact_difference":str(sp.factor(contact))}
        else:
            require(W == 1 and fixed == 1 and S == 1,
                    "constant derivative boundary")

        # Independently compare coefficients of the inverse pencil with the
        # formal Burgers/Lagrange power expansion. Degree six and epsilon two
        # are deliberately finite controls, not all-order certificates.
        nmax, kmax = 6, 2
        B = 0
        for m in range(1,nmax+1):
            value = sp.series(W*(z/a)**m,z,0,m).removeO().expand().coeff(z,m-1)/m
            B += value*t**m
        tree = []
        for k in range(kmax+1):
            term = (-1)**k*sp.diff(sp.expand(B**(k+1)),t,k)/sp.factorial(k+1)
            tree.append(sp.Poly(term,t))
        checks = 0
        for m in range(1,nmax+1):
            direct = sp.series(W*(z/(a+eps*H))**m,z,0,m).removeO().expand()
            direct = sp.cancel(direct.coeff(z,m-1)/m)
            eps_series = sp.series(direct,eps,0,kmax+1).removeO().expand()
            for k in range(kmax+1):
                require(zero(eps_series.coeff(eps,k)-tree[k].nth(m)),
                        name+f": inverse/Burgers coefficient t^{m} eps^{k}")
                checks += 1
        record["B_through_degree_6"] = str(B)
        record["inverse_burgers_coefficient_checks"] = checks
        receipt["cases"].append(record)

    receipt["note_sha256"] = hashlib.sha256(note.read_bytes()).hexdigest()
    receipt["source_sha256"] = DEPENDENCY_SHA256
    receipt["checker_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    receipt["total_inverse_burgers_coefficient_checks"] = sum(
        c["inverse_burgers_coefficient_checks"] for c in receipt["cases"])
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("PASS: four exact pencils; 72 inverse/Burgers coefficients; "
          "misalignment counterexample retained")


if __name__ == "__main__":
    main()
