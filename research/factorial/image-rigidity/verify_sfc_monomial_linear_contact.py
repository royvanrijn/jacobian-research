#!/usr/bin/env python3
"""Exact stdlib controls for SFC_MONOMIAL_LINEAR_CONTACT_BOUND.md.

Polynomial determinant identities are checked in full. Finite derivative and
module controls support, but do not replace, the written uniform proofs.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math


NOTE_SHA = "409ca1ed8a1ef3218fb07673d14830d04cf10562f5714c88bd8de8e9a681b417"


def require(ok, message):
    if not ok:
        raise ValueError(message)


def trim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def coeff(a, k):
    return a[k] if k < len(a) else 0


def add(a, b):
    return trim([coeff(a, k) + coeff(b, k) for k in range(max(len(a), len(b)))])


def scale(a, c):
    return trim([c * x for x in a])


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b, limit=None):
    size = len(a) + len(b) - 1
    if limit is not None:
        size = min(size, limit + 1)
    result = [0] * size
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j < size:
                result[i + j] += x * y
    return trim(result)


def power(a, n, limit=None):
    result = [1]
    for _ in range(n):
        result = mul(result, a, limit)
    return result


def derivative(a):
    return trim([k * a[k] for k in range(1, len(a))]) or [0]


def compose(a, b, limit=None):
    result = [0]
    for c in reversed(a):
        result = mul(result, b, limit)
        result[0] += c
    return trim(result)


def reciprocal(a, limit):
    require(a[0] != 0, "Unit required")
    result = [F(1, a[0])]
    for k in range(1, limit + 1):
        result.append(-sum(coeff(a, j) * result[k-j] for j in range(1, k+1)) / a[0])
    return trim(result)


def inverse(a, limit):
    require(a[:2] == [0, 1], "Normalized inverse")
    b = [0, 1] + [0] * (limit - 1)
    for k in range(2, limit + 1):
        b[k] = -coeff(compose(a, b[:k], k), k)
    require(compose(a, b, limit) == [0, 1], "Composition inverse check")
    return b


def build(exponents, coefficients):
    require(len(exponents) == len(coefficients) > 0, "Nonempty dimensions")
    require(all(type(r) is int and r >= 0 for r in exponents), "Nonnegative exponents")
    require(all(c != 0 for c in coefficients), "Active coefficients")
    p, a0 = [1], [1]
    for r, c in zip(exponents, coefficients):
        p = mul(p, [1, -c])
        a0 = mul(a0, power([1, -c], r))
    e = [0]
    for i, (r, c) in enumerate(zip(exponents, coefficients)):
        term = [1]
        for j, other in enumerate(coefficients):
            if j != i:
                term = mul(term, [1, -other])
        e = add(e, scale(term, -r * c))
    dpoly = add(p, [0] + e)
    a, q = [0] + a0, [0] + p
    require(mul(p, derivative(a0)) == mul(e, a0), "Polynomial logarithmic derivative")
    require(mul(p, derivative(a)) == mul(a0, dpoly), "Derivative quotient identity")
    return p, a0, a, q, e, dpoly


def compositions(n, s):
    if s == 1:
        yield (n,)
    else:
        for k in range(n + 1):
            for rest in compositions(n-k, s-1):
                yield (k,) + rest


def actual_moment(exponents, coefficients, m):
    result = 0
    for multi in compositions(m, len(exponents)):
        multinomial = math.factorial(m) // math.prod(math.factorial(k) for k in multi)
        result += multinomial * math.prod(math.factorial(r*m+k) * c**k
                                         for r, c, k in zip(exponents, coefficients, multi))
    return result


def derivative_checks():
    cases = (((0,), (1,)), ((1,), (-2,)), ((3,), (2,)),
             ((0, 2), (1, -2)), ((1, 2), (1, 2)),
             ((2, 3), (2, -3)), ((0, 0, 2), (1, 1, -2)),
             ((1, 2, 1), (1, -1, 2)))
    rows, moments_count, differential_count, windows = [], 0, 0, 0
    maximum = 15
    for exponents, coefficients in cases:
        s = len(exponents)
        p, a0, a, q, e, dpoly = build(exponents, coefficients)
        require(len(dpoly)-1 == s and dpoly[-1] == (-1)**s * math.prod(coefficients) * (1+sum(exponents)),
                "Exact nonzero residual degree")
        hprime = mul(a0, reciprocal(p, maximum), maximum)
        h = [0] + [F(value, j+1) for j, value in enumerate(hprime)]
        bseries = compose(h, inverse(a, maximum+1), maximum+1)
        values = []
        for m in range(maximum+1):
            value = actual_moment(exponents, coefficients, m)
            factor = math.factorial(m) * math.prod(math.factorial(r*m) for r in exponents) * (m+1)
            require(value == factor * coeff(bseries, m+1), "Actual factorial moment identity")
            moments_count += 1
            values.append(value)
        require(any(values[1:s+1]), "Initial support-size window control")
        for start in (1, 2, 3):
            length = (2*s-1)*start-s+1
            require(start+length-1 <= maximum and any(values[start:start+length]), "Uniform wider window control")
            require(start+length+1 == (2*start-1)*s+2, "Window endpoint arithmetic")
            windows += 1
        mk = scale(derivative(dpoly), -1)
        bderivative = derivative(bseries)
        mk_rows = []
        for k in range(2, 6):
            bderivative = derivative(bderivative)
            denominator = mul(power(hprime, k-1, 5),
                              mul(power(dpoly, 2*k-1, 5), power(p, k-2, 5), 5), 5)
            left = mul(compose(bderivative, a, 5), denominator, 5)
            require(left == trim(mk[:6]), "Full normalized derivative recurrence through degree five")
            require(mk != [0] and len(mk)-1 <= (s-1)+(k-2)*(2*s-1), "Derivative numerator degree bound")
            differential_count += 1
            mk_rows.append({"k": k, "actual_degree": len(mk)-1,
                            "upper_degree": (s-1)+(k-2)*(2*s-1), "numerator": mk})
            bracket = sub(sub(mul(derivative(p), dpoly), scale(mul(e, dpoly), k-1)),
                          scale(mul(p, derivative(dpoly)), 2*k-1))
            mk = add(mul(mul(p, dpoly), derivative(mk)), mul(bracket, mk))
        rows.append({"exponents": exponents, "coefficients": coefficients,
                     "derivative_numerators": mk_rows, "first_five_moments": values[:6]})
    return {"cases": rows, "parameter_cases": len(rows), "moment_identities": moments_count,
            "derivative_identities": differential_count, "window_controls": windows}


def determinant(matrix):
    a = [[F(x) for x in row] for row in matrix]
    value = F(1)
    for k in range(len(a)):
        pivot = next((j for j in range(k, len(a)) if a[j][k]), None)
        if pivot is None:
            return F(0)
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            value = -value
        p = a[k][k]
        value *= p
        for j in range(k+1, len(a)):
            ratio = a[j][k] / p
            for i in range(k+1, len(a)):
                a[j][i] -= ratio * a[k][i]
    return value


def operator(poly, q, dpoly, m):
    return sub(mul(q, derivative(poly)), scale(mul(dpoly, poly), m))


def reduce_class(poly, q, dpoly, m, degree_a):
    result, s = list(poly), len(dpoly)-1
    require(len(poly)-1 <= m*degree_a+s-1, "Numerator within residue module")
    while len(result)-1 >= s:
        j = len(result)-1-s
        require(j <= m*degree_a-1, "Reducer within degree range")
        monomial = [0]*j+[1]
        exact = operator(monomial, q, dpoly, m)
        require(exact[-1] == q[-1]*(j-m*degree_a) != 0, "Nonresonant triangular leading coefficient")
        result = sub(result, scale(exact, F(result[-1], exact[-1])))
    return result


def module_checks():
    inputs = (((1,), (1,)), ((1, 2), (1, 2)), ((2, 3), (2, -3)), ((1, 2, 1), (1, -1, 2)))
    rows, exact_reductions = [], 0
    for exponents, coefficients in inputs:
        require(min(exponents) >= 1 and len(set(coefficients)) == len(coefficients), "Module hypotheses")
        p, a0, a, q, e, dpoly = build(exponents, coefficients)
        s, degree_a = len(exponents), len(a)-1
        roots = [F(1, c) for c in coefficients]
        betas = []
        cleared = [0]
        for alpha, r in zip(roots, exponents):
            beta = 1 / sum(F(k)*q[k]*alpha**(k-1) for k in range(1, len(q)))
            betas.append(beta)
            # Exact synthetic division of Q by z-alpha.
            quotient = [0]*(len(q)-1)
            quotient[-1] = q[-1]
            for k in range(len(quotient)-2, -1, -1):
                quotient[k] = q[k+1]+alpha*quotient[k+1]
            require(q[0]+alpha*quotient[0] == 0, "Exact pole-factor division")
            cleared = add(cleared, scale(quotient, beta-r))
        require(sum(betas) == -1 and cleared == sub([1], dpoly), "Complete logarithmic residue identity")
        for m in (1, 2, 3):
            for j in range(m*degree_a):
                exact = operator([0]*j+[1], q, dpoly, m)
                require(reduce_class(exact, q, dpoly, m, degree_a) == [0], "Exact differential reduces to zero")
                require(operator(mul(a, [0]*j+[1]), q, dpoly, m+1) == mul(a, exact), "Shift respects exact differentials")
                exact_reductions += 1
            columns = [reduce_class([0]*j+a, q, dpoly, m+1, degree_a) for j in range(s)]
            matrix = [[coeff(columns[j], i) for j in range(s)] for i in range(s)]
            det = determinant(matrix)
            require(det != 0, "Finite residue shift is invertible")
            rows.append({"exponents": exponents, "coefficients": coefficients, "m": m,
                         "shift_matrix": matrix, "shift_determinant": det})
    p, a0, a, q, e, dpoly = build((1,), (1,))
    counter = reduce_class([0]+operator([1], q, dpoly, 1), q, dpoly, 1, 2)
    require(counter == [F(-1, 2)], "The exact-differential subspace is not an ideal")
    return {"cases": rows, "shift_cases": len(rows), "exact_reductions": exact_reductions,
            "not_an_algebra_countercontrol": counter}


# Exact arithmetic in Q[t]/(5t^2-2t+2); no floating-point complex arithmetic.
def fadd(a, b):
    return (a[0]+b[0], a[1]+b[1])


def fneg(a):
    return (-a[0], -a[1])


def fmul(a, b):
    return (a[0]*b[0]-F(2,5)*a[1]*b[1], a[0]*b[1]+a[1]*b[0]+F(2,5)*a[1]*b[1])


def finverse(a):
    norm = a[0]**2+F(2,5)*a[0]*a[1]+F(2,5)*a[1]**2
    require(norm != 0, "Nonzero quadratic-field denominator")
    return ((a[0]+F(2,5)*a[1])/norm, -a[1]/norm)


def fpow(a, n):
    result = (F(1), F(0))
    for _ in range(n):
        result = fmul(result, a)
    return result


def feval(poly):
    result = (F(0), F(0))
    for c in reversed(poly):
        result = fadd(fmul(result, (F(0), F(1))), (F(c), F(0)))
    return result


def residue_countercontrol():
    xs, yn = [], []
    for m in range(4):
        x, y = [], [0]
        for k in range(m+1):
            c = math.comb(2*m-k, m-k)*math.comb(2*m+k, k)
            x.append(c)
            y = add(y, scale([0]*k+power([1,-1], m-k), (-1)**(k+1)*c))
        y = scale(y, (-1)**(3*m+1))
        xs.append(x)
        yn.append(y)
    require(xs == [[1], [2,3], [6,15,15], [20,70,112,84]], "Exact origin residues")
    require(yn == [[1], [-2,5], [6,-27,36], [-20,130,-312,286]], "Exact residues at one from local expansion")
    d12 = sub(mul(xs[1], yn[2]), mul(mul(xs[2], yn[1]), power([-1,1], 3)))
    d23 = sub(mul(xs[2], yn[3]), mul(mul(xs[3], yn[2]), power([-1,1], 3)))
    require(d12 == scale([0]*3+mul([-2,1], [2,-2,5]), -15), "Complete first determinant polynomial identity")
    require(d23 == scale([0]*5+mul([-2,1], [5,-5,12]), -252), "Complete second determinant polynomial identity")
    require(feval([2,-2,5]) == (0,0), "Exact singular complex parameter")
    require(feval([5,-5,12]) == (F(1,5),F(-1,5)), "Next quadratic is nonzero at the same parameter")
    vectors = [(feval(x), fmul(feval(y), finverse(fpow((-1,1),3*m+1))))
               for m,(x,y) in enumerate(zip(xs,yn))]
    det = lambda v,w: fadd(fmul(v[0],w[1]), fneg(fmul(v[1],w[0])))
    require(det(vectors[1], vectors[2]) == (0,0), "Consecutive pair degenerates")
    require(det(vectors[2], vectors[3]) != (0,0), "The next pair recovers rank")
    ell = lambda v: fadd(fmul(vectors[1][1],v[0]), fneg(fmul(vectors[1][0],v[1])))
    sequence = [ell(v) for v in vectors]
    require(sequence[0] != (0,0) and sequence[1:3] == [(0,0),(0,0)] and sequence[3] != (0,0),
            "Fixed cycle combination has a two-zero window")
    require(sequence[0] == (F(158,625), F(-14,125)), "Exact zeroth period equals (144-42i)/625")
    require(vectors[1][0] != (0,0), "This is not a factorial SFC counterexample")
    return {"origin_residue_polynomials": xs, "residue_at_one_numerators": yn,
            "first_determinant_numerator": d12, "second_determinant_numerator": d23,
            "quadratic_field_vectors": vectors, "cycle_sequence": sequence,
            "quadratic_relation": [2,-2,5]}


def correspondence_checks():
    a, h = [0,1,-2,1], [0,1,F(-1,2)]
    require(compose(a,[1,1]) == [0,0,1,1] and compose(h,[1,1]) == [F(1,2),0,F(-1,2)],
            "Shared ramification retains a cubic branch term")
    require(compose(a,[2]) == [2] and compose(h,[2]) == [0], "Zero has an additional forward image")
    return {"a_at_one_plus_w": compose(a,[1,1]), "H_at_one_plus_w": compose(h,[1,1]),
            "zero_fiber_parameter_points": [0,2], "forward_values": [0,2],
            "attracting_correspondence_family": "Uniform implicit-function proof in note; no floating-point root count is used."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    note = Path(__file__).with_name("SFC_MONOMIAL_LINEAR_CONTACT_BOUND.md")
    actual = hashlib.sha256(note.read_bytes()).hexdigest()
    require(actual == NOTE_SHA, "Pinned written proof unchanged")
    derivatives, modules = derivative_checks(), module_checks()
    countercontrol, correspondence = residue_countercontrol(), correspondence_checks()
    result = {"status": "PASS", "scope": "Uniform wider-window/contact proof support and exact obstructions; no sharp all-order SFC closure.",
              "note_sha256": actual, "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "derivatives_and_moments": derivatives, "residue_modules": modules,
              "observability_countercontrol": countercontrol, "correspondence_countercontrols": correspondence}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,default=str)+"\n")
    print(json.dumps({"status":"PASS", "parameter_cases":derivatives["parameter_cases"],
                      "moment_identities":derivatives["moment_identities"],
                      "derivative_identities":derivatives["derivative_identities"],
                      "shift_cases":modules["shift_cases"], "exact_reductions":modules["exact_reductions"],
                      "output":str(args.output)}))


if __name__ == "__main__":
    main()
