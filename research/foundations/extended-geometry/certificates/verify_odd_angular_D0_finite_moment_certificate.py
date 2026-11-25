#!/usr/bin/env python3
"""Verify the exact finite-14 obstruction for the normalized odd-angular D0 slice.

Runtime dependencies: Python standard library and python-flint. No Singular,
SymPy, network access, original experiment directory, or stored raw moments.
All arithmetic identities are over QQ; the sole modular step is a certified
univariate coprimeness argument with preserved leading degrees.
"""
from pathlib import Path
from fractions import Fraction
from collections import defaultdict
from functools import reduce
from math import factorial, gcd, isqrt, lcm
import argparse
import gzip
import hashlib
import json
import sys
import time
from flint import fmpq, fmpq_mpoly, fmpq_mpoly_ctx

CERTIFICATE_SHA256 = "f3613dd5da4a853dc5a42e8664394cd5de8185a1ca27a08152b8281d96c74dff"
NATIVE_NAMES = ("a1", "a3", "b1", "q0", "lam")
EVEN = (6, 8, 10, 12, 14)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def poly(text, context):
    return fmpq_mpoly(text.replace("**", "^"), context)


def constant(context, value):
    return context.constant(fmpq(str(value)))


def compositions(n, k):
    if k == 1:
        yield (n,)
    else:
        for j in range(n + 1):
            for tail in compositions(n - j, k - 1):
                yield (j,) + tail


def native_monomial_trace(i, j, k):
    if i != j or k % 2:
        return Fraction(0)
    a = k // 2
    return Fraction(
        2**i * factorial(i) * factorial(2 * a) * factorial(i + a),
        factorial(a) * factorial(2 * i + 2 * a + 1),
    )


def native_moments(context):
    # Each tuple is x exponent, y exponent, t exponent, parameter index.
    # None denotes one of the two fixed coefficients, both equal to one.
    phase = (
        (0, 1, 0, 0), (0, 3, 0, 1), (0, 1, 2, 2),
        (0, 3, 2, None), (0, 1, 4, 4),
        (1, 0, 2, 3), (1, 0, 4, None),
    )
    moments = {}
    for m in range(1, 15):
        answer = defaultdict(Fraction)
        for counts in compositions(m, len(phase)):
            xx = yy = tt = 0
            exponents = [0] * 5
            coefficient = Fraction(factorial(m))
            for count, (i, j, k, parameter) in zip(counts, phase):
                xx += count * i
                yy += count * j
                tt += count * k
                coefficient /= factorial(count)
                if parameter is not None:
                    exponents[parameter] += count
            value = native_monomial_trace(xx, yy, tt)
            if value:
                answer[tuple(exponents)] += coefficient * value
        moments[m] = context.from_dict({
            ex: fmpq(value.numerator, value.denominator)
            for ex, value in answer.items() if value
        })
        if m % 2:
            require(not moments[m], "An odd native moment did not vanish.")
    return moments


def pure_linear_coefficient(raw, index):
    answer = {}
    for ex, value in raw.to_dict().items():
        require(ex[index] <= 1, "The claimed pivot is not affine.")
        if ex[index] == 1:
            ee = list(map(int, ex))
            ee[index] = 0
            answer[tuple(ee)] = value
    return raw.context().from_dict(answer)


def substituted_fraction(raw, context, n1, d1, n3, d3):
    terms = raw.to_dict()
    aa = max((int(ex[0]) for ex in terms), default=0)
    bb = max((int(ex[1]) for ex in terms), default=0)
    n1p = [n1**i for i in range(aa + 1)]
    d1p = [d1**i for i in range(aa + 1)]
    n3p = [n3**i for i in range(bb + 1)]
    d3p = [d3**i for i in range(bb + 1)]
    numerator = context.constant(0)
    for ex, value in terms.items():
        i, j, b, q, ll = map(int, ex)
        term = context.from_dict({(b, q, ll): value})
        numerator += term * n1p[i] * d1p[aa-i] * n3p[j] * d3p[bb-j]
    return numerator, d1**aa * d3**bb


def only_pivot_denominator(denominator, p1, p3):
    remainder = denominator
    for pivot in (p1, p3):
        while remainder and not remainder % pivot:
            remainder = remainder // pivot
    require(bool(remainder) and all(all(int(i) == 0 for i in ex)
                                   for ex in remainder.to_dict()),
            "A generic reduced denominator has an unaccounted zero.")


def qq_dense_in_q(p, q_index):
    degree = max((int(ex[q_index]) for ex in p.to_dict()), default=0)
    answer = [Fraction(0)] * (degree + 1)
    for ex, value in p.to_dict().items():
        require(all(int(e) == 0 for i, e in enumerate(ex) if i != q_index),
                "A scalar contains a non-q variable.")
        answer[int(ex[q_index])] = Fraction(str(value))
    return answer


def convolution(a, b, modulus=None):
    answer = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            answer[i+j] += x * y
            if modulus is not None:
                answer[i+j] %= modulus
    return answer


def sum_sequences(a, b):
    return [(a[i] if i < len(a) else 0) +
            (b[i] if i < len(b) else 0)
            for i in range(max(len(a), len(b)))]


def same_sequence(a, b):
    count = max(len(a), len(b))
    return all((a[i] if i < len(a) else 0) ==
               (b[i] if i < len(b) else 0) for i in range(count))


def primitive_integer_coefficients(p, q_index):
    coefficients = qq_dense_in_q(p, q_index)
    denominator = lcm(*(x.denominator for x in coefficients))
    integers = [int(x * denominator) for x in coefficients]
    content = reduce(gcd, map(abs, integers))
    require(content != 0, "A necessary univariate polynomial is zero.")
    integers = [x // content for x in integers]
    if integers[-1] < 0:
        integers = [-x for x in integers]
    require(reduce(gcd, map(abs, integers)) == 1, "Nonprimitive polynomial.")
    return integers


def is_prime(n):
    return n >= 2 and all(n % j for j in range(2, isqrt(n) + 1))


def verify(certificate):
    started = time.monotonic()
    native = fmpq_mpoly_ctx.get(NATIVE_NAMES, ordering="lex")
    a1, a3, b1, q0, lam = native.gens()
    moments = native_moments(native)
    controls = {}
    require(moments[2] == 2*(231*a1*q0+99*a1+99*b1*q0+55*b1+55*lam*q0+35*lam)/3465,
            "The full native M2 pivot/exceptional compatibility changed.")
    p1n = 7*q0 + 3
    p3n = 1615*q0**3 + 2261*q0**2 + 1197*q0 + 231
    require(pure_linear_coefficient(moments[2], 0) == 2*p1n/105,
            "M2 pivot identity failed.")
    require(pure_linear_coefficient(moments[4], 1) == 8*p3n/4849845,
            "M4 pivot identity failed.")
    # Linear-resultant formula, evaluated by Fraction rather than a CAS
    # argument-order convention: Res(7q+3,P3)=7^3 P3(-3/7)=+2112.
    q_special = Fraction(-3, 7)
    res = 7**3 * (1615*q_special**3 + 2261*q_special**2 +
                  1197*q_special + 231)
    require(res == 2112, "The two pivot polynomials were not separated.")

    generic = fmpq_mpoly_ctx.get(("b1", "q0", "lam"), ordering="lex")
    b, q, ll = generic.gens()
    p1 = 7*q + 3
    p3 = 1615*q**3 + 2261*q**2 + 1197*q + 231
    subs = certificate["generic_substitutions"]
    n1, d1, n3, d3 = [poly(subs[key], generic)
                      for key in ("n1", "d1", "n3", "d3")]
    require(d1 == 33*p1 and d3 == 1089*p1**2*p3,
            "The two generic substitution denominators changed.")
    for m in (2, 4):
        numerator, denominator = substituted_fraction(
            moments[m], generic, n1, d1, n3, d3)
        require(not numerator and bool(denominator),
                "A solved pivot moment did not vanish identically.")
    corrupted, _ = substituted_fraction(moments[2] + 1, generic,
                                        n1, d1, n3, d3)
    controls["altered_native_M2_rejected"] = bool(corrupted)
    numerators = {}
    for m in EVEN:
        numerator, denominator = substituted_fraction(
            moments[m], generic, n1, d1, n3, d3)
        common = numerator.gcd(denominator)
        numerator, denominator = numerator//common, denominator//common
        ratio = fmpq(certificate["generic_numerator_normalization_ratios"][str(m)])
        require(bool(ratio), "A zero numerator normalization was supplied.")
        numerator, denominator = numerator*ratio, denominator*ratio
        only_pivot_denominator(denominator, p1, p3)
        numerators[m] = numerator

    scaled_context = fmpq_mpoly_ctx.get(("h", "q0", "z"), ordering="lex")
    h, qs, z = scaled_context.gens()
    ps1 = 7*qs + 3
    ps3 = 1615*qs**3 + 2261*qs**2 + 1197*qs + 231
    scaled = {}
    for m, numerator in numerators.items():
        parity = (m//2) % 2
        out = {}
        for ex, value in numerator.to_dict().items():
            ib, iq, il = map(int, ex)
            require(ib+il >= parity and (ib+il-parity) % 2 == 0,
                    "The b1=h*lambda, z=lambda^2 transformation failed.")
            out[(ib, iq, (ib+il-parity)//2)] = value
        scaled[m] = scaled_context.from_dict(out)
    scalars = []
    expected_orders = ((6,8,10,12,14), (6,8,10,14))
    lifts = certificate["generic_cleared_lifts"]
    require(len(lifts) == 2, "Exactly two cleared identities are required.")
    for index, (lift, order) in enumerate(zip(lifts, expected_orders), 1):
        require(tuple(lift["moments"]) == order, "A lift generator changed.")
        cofs = [poly(text, scaled_context) for text in lift["cofactors"]]
        require(len(cofs) == len(order), "A lift cofactor is missing.")
        answer = sum((c*scaled[m] for c, m in zip(cofs, order)),
                     scaled_context.constant(0))
        scalar = poly(lift["scalar"], scaled_context)
        require(answer == scalar and bool(scalar),
                "A complete cleared polynomial lift identity failed.")
        qq_dense_in_q(scalar, 1)
        scalars.append(scalar)
        controls["altered_generic_cofactor_%s_rejected" % index] = (
            answer + scaled[order[0]] != scalar)
    primitive = []
    for scalar, exponent in zip(scalars, (3,2)):
        pivot = ps1**7 * ps3**exponent
        require(not scalar % pivot, "A claimed pivot factor does not divide.")
        primitive.append(primitive_integer_coefficients(scalar//pivot, 1))
    cop = certificate["generic_univariate_coprimeness"]
    prime = int(cop["prime"])
    require(is_prime(prime), "The modular coprimeness modulus is not prime.")
    require([len(v)-1 for v in primitive] == list(cop["degrees"]),
            "The univariate degrees changed.")
    require(all(v[-1] % prime for v in primitive),
            "A univariate leading degree is lost modulo the prime.")
    u = list(map(int, cop["modular_Bezout_u_coefficients_ascending"]))
    v = list(map(int, cop["modular_Bezout_v_coefficients_ascending"]))
    value = sum_sequences(convolution(primitive[0], u, prime),
                          convolution(primitive[1], v, prime))
    require(value[0] % prime == 1 and
            all(x % prime == 0 for x in value[1:]),
            "The modular univariate Bezout identity failed.")
    bad_u = u[:]
    bad_u[0] += 1
    bad = sum_sequences(convolution(primitive[0], bad_u, prime),
                        convolution(primitive[1], v, prime))
    controls["altered_modular_Bezout_rejected"] = not (
        bad[0] % prime == 1 and all(x % prime == 0 for x in bad[1:]))
    controls["composite_modulus_rejected"] = not is_prime(32003*3)
    # By Gauss lemma and preserved leading degrees this modular identity
    # proves QQ-coprimeness, not a multivariate modular unit assertion.

    zero_context = fmpq_mpoly_ctx.get(("q0", "z"), ordering="lex")
    qz, zz = zero_context.gens()
    pz1 = 7*qz+3
    pz3 = 1615*qz**3+2261*qz**2+1197*qz+231
    zpolys = {}
    for m, numerator in numerators.items():
        original = numerator.subs({"lam":0})
        parity = (m//2) % 2
        out = {}
        for ex, value in original.to_dict().items():
            ib, iq, il = map(int, ex)
            require(il == 0 and ib >= parity and (ib-parity) % 2 == 0,
                    "The lambda=0, z=b1^2 transformation failed.")
            out[(iq,(ib-parity)//2)] = value
        zpolys[m] = zero_context.from_dict(out)

    def zcoefficient(p, j):
        return zero_context.from_dict({
            (int(ex[0]),0): value for ex, value in p.to_dict().items()
            if int(ex[1]) == j
        })

    require(int(zpolys[6].degrees()[1]) == 1,
            "The lambda=0 first equation is not linear in z.")
    aa, bb = zcoefficient(zpolys[6],1), zcoefficient(zpolys[6],0)
    rs = []
    for m in (8,10):
        require(int(zpolys[m].degrees()[1]) <= 2,
                "A lambda=0 resultant input is not quadratic.")
        cc, dd, ee = [zcoefficient(zpolys[m], j) for j in (2,1,0)]
        uu = -(aa*cc*zz+aa*dd-bb*cc)
        vv = aa*aa
        scalar = aa*aa*ee-aa*bb*dd+bb*bb*cc
        require(uu*zpolys[6]+vv*zpolys[m] == scalar and bool(scalar),
                "A lambda=0 division-free resultant identity failed.")
        controls["altered_lambda0_resultant_%s_rejected" % m] = (
            (uu+1)*zpolys[6]+vv*zpolys[m] != scalar)
        rs.append(scalar)

    pairs = (
        ("nonzero_b1", rs[0], rs[1], pz1**4*pz3**2),
        ("zero_b1", zpolys[8].subs({"z":0}),
         zpolys[12].subs({"z":0}), pz1**4),
    )
    for label, f, g, pivot in pairs:
        entry = certificate["lambda0_univariate_Bezout"][label]
        ua, va = ([Fraction(s) for s in entry[key]] for key in ("u","v"))
        expected = pivot*fmpq(entry["gcd_pivot_scale"])
        require(bool(expected), "A lambda=0 gcd scale is zero.")
        value = sum_sequences(convolution(qq_dense_in_q(f,0),ua),
                              convolution(qq_dense_in_q(g,0),va))
        require(same_sequence(value,qq_dense_in_q(expected,0)),
                "An exact QQ univariate extended-gcd identity failed.")
        bad = sum_sequences(value,qq_dense_in_q(f,0))
        controls["altered_lambda0_%s_Bezout_rejected" % label] = (
            not same_sequence(bad,qq_dense_in_q(expected,0)))

    p1_context = fmpq_mpoly_ctx.get(("a1","lam"), ordering="lex")
    ag, lg = p1_context.gens()
    csol = poly(certificate["exceptional_substitutions"]["P1_a3"],p1_context)
    args = (ag,csol,-10*lg/11,constant(p1_context,Fraction(-3,7)),lg)
    p1_residuals = {}
    for m in range(2,15,2):
        p1_residuals[m] = moments[m].compose(*args)
    require(not p1_residuals[2] and not p1_residuals[4],
            "The P1 exceptional substitutions are not native solutions.")
    cofs = [poly(text,p1_context)
            for text in certificate["exceptional_cofactors"]["P1"]]
    require(len(cofs) == 5, "A P1 exceptional cofactor is missing.")
    value = sum((c*p1_residuals[m] for c,m in zip(cofs,EVEN)),
                p1_context.constant(0))
    require(value == 1, "The P1 exceptional constant-one identity failed.")
    controls["altered_P1_cofactor_rejected"] = value+p1_residuals[6] != 1

    p3_context = fmpq_mpoly_ctx.get(("a3","b1","lam","q0"),ordering="lex")
    cg,bg,lg,qg = p3_context.gens()
    pp = 1615*qg**3+2261*qg**2+1197*qg+231
    asol = poly(certificate["exceptional_substitutions"]["P3_a1"],p3_context)
    args = (asol,cg,bg,qg,lg)
    p3_residuals = {
        m:moments[m].compose(*args)%pp for m in range(2,15,2)
    }
    require(not p3_residuals[2], "The P3 substitution is not a native solution.")
    cofs = [poly(text,p3_context)
            for text in certificate["exceptional_cofactors"]["P3_quotient"]]
    require(len(cofs) == 6, "A P3 exceptional quotient cofactor is missing.")
    order = (4,6,8,10,12,14)
    value = sum((c*p3_residuals[m] for c,m in zip(cofs,order)),
                p3_context.constant(0))
    require((value-1)%pp == 0,
            "The P3 exceptional quotient constant-one identity failed.")
    seventh = (1-value)//pp
    require(value+seventh*pp == 1,
            "Recovery of the seventh exact P3 cofactor failed.")
    controls["altered_P3_cofactor_rejected"] = (
        (value+p3_residuals[4]-1)%pp != 0)
    require(all(controls.values()), "A corruption-rejection control failed.")
    return {
        "status":"PASS",
        "claim":"No parameter point in the specified normalized odd-angular D0 slice has native moments M1..M14 all zero, over any characteristic-zero field.",
        "native_moment_term_counts":{str(m):len(moments[m]) for m in range(1,15)},
        "pivot_resultant_standard_order":2112,
        "generic_generator_orders":[list(order) for order in expected_orders],
        "generic_saturations_or_extra_inverses":[],
        "nonpivot_univariate_degrees":[len(a)-1 for a in primitive],
        "coprimeness_prime":prime,
        "coprimeness_argument":"Primitive integer common factor plus preserved leading degrees would give a positive-degree modular common factor, contradicting the exact modular Bezout identity.",
        "lambda_zero_chart_verified":True,
        "both_exceptional_pivot_branches_verified":True,
        "corruption_controls":controls,
        "runtime_Groebner_or_Singular":False,
        "seconds":round(time.monotonic()-started,3),
        "scope_limits":["Exactly the displayed normalized slice.","No eventual-purity, full-D0, or global degree-five assertion."],
        "assurance":{"formal_verification":False,"external_human_review":False},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate",type=Path)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    if hasattr(sys,"set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    if args.certificate is None:
        sibling = Path(__file__).with_name("odd-angular-d0-finite-moment-certificate.json.gz")
        repository_input = Path(__file__).resolve().parents[1] / "foundations" / "extended-geometry" / "certificates" / sibling.name
        args.certificate = sibling if sibling.is_file() else repository_input
    data = args.certificate.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    require(digest == CERTIFICATE_SHA256,"Compressed certificate hash mismatch.")
    certificate = json.loads(gzip.decompress(data))
    require(certificate.get("schema") == 1,"Unsupported certificate schema.")
    answer = verify(certificate)
    answer["certificate_sha256"] = digest
    answer["checker_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    answer["compressed_input_corruption_rejected"] = (
        hashlib.sha256(data+b"\0").hexdigest() != CERTIFICATE_SHA256)
    text = json.dumps(answer,indent=2,sort_keys=True)+"\n"
    if args.output:
        require(not args.output.exists(),"Refusing to overwrite an existing receipt.")
        args.output.write_text(text)
    print(text,end="")


if __name__ == "__main__":
    main()
