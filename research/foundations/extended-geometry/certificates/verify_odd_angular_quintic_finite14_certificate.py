#!/usr/bin/env python3
"""Verify the exact finite-14 obstruction for the normalized eight-monomial odd-angular quintic.

Runtime dependencies: Python standard library and python-flint. No Singular,
SymPy, network access, original experiment directory, or stored raw moments.
All arithmetic identities are over QQ. Modular arithmetic certifies large
integer resultants by rigorous CRT bounds and univariate coprimeness with
preserved leading degrees.
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
from flint import fmpq, fmpq_mpoly, fmpq_mpoly_ctx, fmpq_poly, fmpz_poly, nmod_mpoly_ctx, nmod_poly

CERTIFICATE_SHA256 = "9bb2c3ec17173f7ba4461bd5cefb6b6981b3af771921d4c38a700a0568dc5321"
NATIVE_NAMES = ("a1", "a3", "a5", "b1", "q0", "lam")
EVEN = (8, 10, 12, 14)


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
    # Eight literal monomials; the two None entries have coefficient one.
    # This reconstruction does not read stored moments or generator output.
    phase = (
        (0,1,0,0), (0,3,0,1), (0,5,0,2), (0,1,2,3),
        (0,3,2,None), (0,1,4,5), (1,0,2,4), (1,0,4,None),
    )
    moments = {}
    for m in range(1,15):
        answer = defaultdict(Fraction)
        for counts in compositions(m,8):
            xx = sum(c*term[0] for c,term in zip(counts,phase))
            yy = sum(c*term[1] for c,term in zip(counts,phase))
            tt = sum(c*term[2] for c,term in zip(counts,phase))
            if xx != yy or tt % 2:
                continue
            exponents = [0]*6
            denominator = 1
            for count,term in zip(counts,phase):
                denominator *= factorial(count)
                if term[3] is not None:
                    exponents[term[3]] += count
            coefficient = Fraction(factorial(m),denominator)
            answer[tuple(exponents)] += coefficient*native_monomial_trace(xx,yy,tt)
        moments[m] = context.from_dict({
            ex:fmpq(value.numerator,value.denominator)
            for ex,value in answer.items() if value
        })
        if m % 2:
            require(not moments[m],"An odd native moment did not vanish.")
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


def substituted_fraction(raw,context,solutions):
    data = raw.to_dict()
    maximum = {index:max((int(e[index]) for e in data),default=0)
               for index in range(3)}
    powers = {
        index:([n**j for j in range(maximum[index]+1)],
               [d**j for j in range(maximum[index]+1)])
        for index,(n,d) in enumerate(solutions)
    }
    numerator = context.constant(0)
    denominator = context.constant(1)
    for index,(n,d) in enumerate(solutions):
        denominator *= d**maximum[index]
    for e,value in data.items():
        ee = list(map(int,e))
        term = context.from_dict({tuple(ee[3:]):value})
        for index in range(3):
            nn,dd = powers[index]
            term *= nn[ee[index]]*dd[maximum[index]-ee[index]]
        numerator += term
    return numerator,denominator


def only_pivot_denominator(denominator, *pivots):
    remainder = denominator
    for pivot in pivots:
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



def primitive_polynomial(f):
    require(bool(f),"A claimed primitive polynomial is zero.")
    values = [Fraction(str(c)) for c in f.to_dict().values()]
    denominator = lcm(*(c.denominator for c in values))
    content = reduce(gcd,(abs(int(c*denominator)) for c in values))
    scale = Fraction(content,denominator)
    if f.leading_coefficient() < 0:
        scale = -scale
    prim = f/ fmpq(scale.numerator,scale.denominator)
    require(all(c.q == 1 for c in prim.to_dict().values()),
            "Primitive normalization did not give integer coefficients.")
    return prim,scale


def reduce_rational_numerator(n,d):
    common = n.gcd(d)
    n,d = n//common,d//common
    prim,scale = primitive_polynomial(n)
    return prim,d/fmpq(scale.numerator,scale.denominator)


def field_inverse(f,modulus,q_index):
    context = f.context()
    def dense(g):
        return fmpq_poly([fmpq(c.numerator,c.denominator)
                          for c in qq_dense_in_q(g,q_index)])
    gg,u,v = dense(f).xgcd(dense(modulus))
    require(gg == 1,"A claimed field coefficient is not a unit.")
    out = {}
    for degree,value in enumerate(u):
        if value:
            e = [0]*context.nvars()
            e[q_index] = degree
            out[tuple(e)] = value
    inverse = context.from_dict(out)
    require((inverse*f)%modulus == 1,"Polynomial inverse identity failed.")
    return inverse


def derivative(f,index):
    out = {}
    for e,value in f.to_dict().items():
        if e[index]:
            ee = list(map(int,e))
            power = ee[index]
            ee[index] -= 1
            out[tuple(ee)] = value*power
    return f.context().from_dict(out)


def coefficient(f,index,degree):
    return f.context().from_dict({
        tuple(0 if j == index else int(v) for j,v in enumerate(e)):c
        for e,c in f.to_dict().items() if int(e[index]) == degree
    })


def constant_linear_solution(f,index):
    c = pure_linear_coefficient(f,index)
    require(bool(c) and all(all(int(j) == 0 for j in e)
                            for e in c.to_dict()),
            "A constant exceptional pivot was not a nonzero constant.")
    return -f.subs({index:0})/c.leading_coefficient()


def cofactor_group(group,context):
    cs = [poly(s,context) for s in group["cofactors"]]
    require(all(c.q == 1 for f in cs for c in f.to_dict().values()),
            "A purported integer-scaled cofactor is not integral.")
    scalar = fmpq(group["nonzero_rational_scale"])
    require(bool(scalar),"A cofactor normalization scalar is zero.")
    return cs,context.constant(scalar)


def verify_modular_coprimeness(first,second,prime,u,v,controls,label):
    require(is_prime(prime),"A coprimeness modulus is not prime.")
    require(all(a and a[-1]%prime for a in (first,second)),
            "A full leading degree was lost in the coprimeness check.")
    require(all(reduce(gcd,map(abs,a)) == 1 for a in (first,second)),
            "A univariate integer polynomial was not primitive.")
    a,b = [[x%prime for x in arr] for arr in (first,second)]
    u,v = list(map(int,u)),list(map(int,v))
    value = sum_sequences(convolution(a,u,prime),convolution(b,v,prime))
    require(value[0]%prime == 1 and
            all(c%prime == 0 for c in value[1:]),
            "A complete modular univariate Bezout identity failed.")
    bad = u[:]
    require(bool(bad),"A modular Bezout coefficient vector is empty.")
    bad[0] += 1
    altered = sum_sequences(convolution(a,bad,prime),convolution(b,v,prime))
    controls[label+"_Bezout_corruption_rejected"] = not (
        altered[0]%prime == 1 and all(c%prime == 0 for c in altered[1:]))


def integer_univariate_evaluation(coefficients,value):
    answer = 0
    for c in reversed(coefficients):
        answer = answer*value+c
    return answer


def integer_h_specialization(f,value):
    answer = {}
    for (ih,iq,iz),c in f.to_dict().items():
        require(iz == 0 and c.q == 1,"A CRT input is not an integer h,q polynomial.")
        answer[int(ih)] = answer.get(int(ih),0)+int(c)*value**int(iq)
    return fmpz_poly([answer.get(i,0) for i in range(max(answer)+1)])


def prime_table(limit):
    sieve = bytearray(b"\x01")*(limit+1)
    sieve[:2] = b"\x00\x00"
    for p in range(2,isqrt(limit)+1):
        if sieve[p]:
            sieve[p*p:limit+1:p] = b"\x00"*(((limit-p*p)//p)+1)
    return [i for i in range(2,limit+1) if sieve[i]]


def verify_large_resultant(a,b,entry,controls,prime_cache):
    started = time.monotonic()
    da,db = int(a.degrees()[0]),int(b.degrees()[0])
    qa,qb = int(a.degrees()[1]),int(b.degrees()[1])
    degree_bound = db*qa+da*qb
    cs = [int(s) for s in entry["primitive_coefficients_ascending_q"]]
    require(cs and cs[-1] and len(cs)-1 <= degree_bound,
            "A candidate resultant violates its degree bound.")
    require(reduce(gcd,map(abs,cs)) == 1,
            "A candidate integer resultant is not primitive.")
    value = int(entry["scalar_evaluation_q"])
    av,bv = integer_h_specialization(a,value),integer_h_specialization(b,value)
    cv = integer_univariate_evaluation(cs,value)
    require(av.degree() == da and bv.degree() == db and cv,
            "The scalar evaluation loses degree or has zero denominator.")
    exact = int(av.resultant(bv))
    require(exact%cv == 0,"The resultant scalar is not an exact integer.")
    scalar = exact//cv
    require(scalar and scalar == int(entry["primitive_scalar"]),
            "The claimed primitive resultant scalar changed.")
    claimed = [scalar*c for c in cs]
    norm_a = sum(abs(int(c)) for c in a.to_dict().values())
    norm_b = sum(abs(int(c)) for c in b.to_dict().values())
    # The l1 norm is submultiplicative. Expanding the Sylvester determinant,
    # its permutation terms are a subset of all row-wise entry products.
    # Hence EVERY coefficient is bounded by this exact integer, no factorial.
    bound = norm_a**db*norm_b**da
    maximum = max(map(abs,claimed))
    require(maximum <= bound,"A claimed coefficient exceeds the determinant bound.")
    target = bound+maximum
    primes = list(map(int,entry["primes"]))
    require(primes and len(set(primes)) == len(primes),
            "CRT moduli are absent or repeated.")
    require(all(2 <= p < 2**31 for p in primes),"CRT modulus outside trial-division range.")
    product = 1
    scalar_control = False
    trial = prime_cache["trial"]
    for p in primes:
        if p not in prime_cache["verified"]:
            require(all(p%r for r in trial if r*r <= p),
                    "A CRT modulus is composite.")
            prime_cache["verified"].add(p)
        context = nmod_mpoly_ctx.get(("h","q0"),p,ordering="lex")
        def reduced(f):
            return context.from_dict({
                (int(ih),int(iq)):int(c)%p
                for (ih,iq,iz),c in f.to_dict().items() if int(c)%p
            })
        aa,bb = reduced(a),reduced(b)
        require(int(aa.degrees()[0]) == da and int(bb.degrees()[0]) == db,
                "A CRT input lost its full h degree.")
        actual = aa.resultant(bb,"h")
        dd = actual.to_dict()
        require(all(int(ih) == 0 for ih,iq in dd),"A CRT resultant is not univariate.")
        require(int(actual.degrees()[1]) <= degree_bound,
                "A modular resultant violates the Sylvester degree bound.")
        observed = [int(dd.get((0,i),0)) for i in range(degree_bound+1)]
        expected = [(claimed[i] if i < len(claimed) else 0)%p
                    for i in range(degree_bound+1)]
        require(observed == expected,"A complete CRT resultant coefficient check failed.")
        if scalar%p:
            altered = expected[:]
            altered[0] = (altered[0]+scalar)%p
            require(altered != observed,"A changed primitive constant escaped detection.")
            scalar_control = True
        product *= p
    require(product > target,"The exact CRT product does not exceed the error bound.")
    controls["CRT_%s_candidate_corruption_rejected"%entry["other"]] = scalar_control
    controls["CRT_%s_insufficient_product_rejected"%entry["other"]] = primes[0] <= target
    controls["CRT_%s_wrong_scalar_rejected"%entry["other"]] = exact != (scalar+1)*cv
    quotient,remainder = divmod(fmpz_poly(cs),
                               fmpz_poly([3,7])**int(entry["P1_exponent"]))
    require(not remainder,"The claimed P1 resultant factor division failed.")
    residual = [int(c) for c in quotient]
    content = reduce(gcd,map(abs,residual))
    require(content != 0,"A resultant residual is zero.")
    residual = [c//content for c in residual]
    return residual,{
        "other_moment":int(entry["other"]),
        "h_degrees":[da,db],"q_degrees":[qa,qb],
        "determinant_q_degree_bound":degree_bound,
        "candidate_degree":len(cs)-1,
        "true_coefficient_bound_bits":bound.bit_length(),
        "claimed_max_coefficient_bits":maximum.bit_length(),
        "prime_count":len(primes),"CRT_product_bits":product.bit_length(),
        "strict_CRT_error_bound_verified":True,
        "residual_degree":len(residual)-1,
        "seconds":round(time.monotonic()-started,3),
    }


def evaluate_mod(f,values,p):
    answer = 0
    for ex,c in f.to_dict().items():
        denominator = int(c.q)%p
        require(denominator != 0,"A sharpness moment coefficient is not p-integral.")
        term = int(c.p)%p * pow(denominator,-1,p)%p
        for value,power in zip(values,ex):
            term = term*pow(int(value),int(power),p)%p
        answer = (answer+term)%p
    return answer


def determinant_mod(matrix,p):
    a = [row[:] for row in matrix]
    answer = 1
    for i in range(len(a)):
        pivot = next((j for j in range(i,len(a)) if a[j][i]%p),None)
        if pivot is None:
            return 0
        if pivot != i:
            a[i],a[pivot] = a[pivot],a[i]
            answer = -answer
        value = a[i][i]%p
        answer = answer*value%p
        inverse = pow(value,-1,p)
        for j in range(i+1,len(a)):
            factor = a[j][i]*inverse%p
            for k in range(i,len(a)):
                a[j][k] = (a[j][k]-factor*a[i][k])%p
    return answer%p


def verify(certificate):
    started = time.monotonic()
    controls = {}
    native = fmpq_mpoly_ctx.get(NATIVE_NAMES,ordering="lex")
    a1,a3,a5,b1,q0,lam = native.gens()
    moments = native_moments(native)
    require(moments[2] == 2*(231*a1*q0+99*a1+99*b1*q0+55*b1+
                             55*lam*q0+35*lam)/3465,
            "The full native M2 changed.")
    p1n = 7*q0+3
    p3n = 1615*q0**3+2261*q0**2+1197*q0+231
    p5n = (930465*q0**5+2225025*q0**4+2314026*q0**3+
           1285570*q0**2+376805*q0+46189)
    require(pure_linear_coefficient(moments[2],0) == 2*p1n/105,
            "The a1 pivot changed.")
    require(pure_linear_coefficient(moments[4],1) == 8*p3n/4849845,
            "The a3 pivot changed.")
    require(pure_linear_coefficient(moments[6],2) == 16*p5n/300840735195,
            "The a5 pivot changed.")
    qq = Fraction(-3,7)
    require(7**3*(1615*qq**3+2261*qq**2+1197*qq+231) == 2112,
            "The P1/P3 ordered resultant changed.")
    require(7**5*(930465*qq**5+2225025*qq**4+2314026*qq**3+
                  1285570*qq**2+376805*qq+46189) == 4756480,
            "The P1/P5 ordered resultant changed.")

    generic = fmpq_mpoly_ctx.get(("b1","q0","lam"),ordering="lex")
    b,q,ll = generic.gens()
    p1 = 7*q+3
    p3 = 1615*q**3+2261*q**2+1197*q+231
    p5 = 930465*q**5+2225025*q**4+2314026*q**3+1285570*q*q+376805*q+46189
    for f,g in ((p1,p3),(p1,p5),(p3,p5)):
        field_inverse(f,g,1)
    for f in (p1,p3,p5):
        field_inverse(derivative(f,1),f,1)
    expected_denominators = (p1,p1**2*p3,p1**3*p3*p5)
    solutions = []
    for name,expected in zip(("a1","a3","a5"),expected_denominators):
        entry = certificate["solutions"][name]
        n,d = poly(entry["n"],generic),poly(entry["d"],generic)
        require(d == expected/expected.leading_coefficient(),
                "A declared monic triangular denominator changed.")
        solutions.append((n,d))
    for m in (2,4,6):
        n,d = substituted_fraction(moments[m],generic,solutions)
        require(not n and bool(d),"A native triangular solution identity failed.")
    bad,_ = substituted_fraction(moments[2]+1,generic,solutions)
    controls["altered_native_M2_rejected"] = bool(bad)
    numerators = {}
    for m in EVEN:
        n,d = substituted_fraction(moments[m],generic,solutions)
        n,d = reduce_rational_numerator(n,d)
        only_pivot_denominator(d,p1,p3,p5)
        numerators[m] = n

    # Generic q and lambda nonzero.
    gc = fmpq_mpoly_ctx.get(("h","q0","z"),ordering="lex")
    hh,qg,zz = gc.gens()
    gp1 = 7*qg+3
    gp3 = 1615*qg**3+2261*qg**2+1197*qg+231
    gp5 = (930465*qg**5+2225025*qg**4+2314026*qg**3+
           1285570*qg*qg+376805*qg+46189)
    scaled = {}
    for m,n in numerators.items():
        par = (m//2)%2
        dd = {}
        for ex,c in n.to_dict().items():
            ib,iq,il = map(int,ex)
            require(ib+il >= par and (ib+il-par)%2 == 0,
                    "A generic h,z scaling is invalid.")
            dd[(ib,iq,(ib+il-par)//2)] = c
        scaled[m] = gc.from_dict(dd)
    resultants = {}
    first_stage_report = []
    for row in certificate["generic_resultants"]["first_stage"]:
        m = int(row["m"])
        require(m in (10,12,14),"Unexpected first-stage resultant.")
        actual = scaled[8].resultant(scaled[m],"z")
        exponents = [int(v) for v in row["q_content_exponents"]]
        require(exponents == ([8,4,1] if m == 10 else [12,6,3]),
                "A q-content exponent changed.")
        factor = gp1**exponents[0]*gp3**exponents[1]*gp5**exponents[2]
        require(not actual%factor,"An exact first-stage q content does not divide.")
        primitive,scale = primitive_polynomial(actual//factor)
        require(scale == Fraction(row["constant_scale"]),
                "A first-stage integer normalization scale changed.")
        require(actual == factor*primitive*fmpq(scale.numerator,scale.denominator),
                "The complete first-stage identity failed.")
        resultants[m] = primitive
        first_stage_report.append({"other_moment":m,
            "primitive_h_q_degrees":list(map(int,primitive.degrees()[:2])),
            "q_content_exponents":exponents})
        controls["first_stage_%s_content_corruption_rejected"%m] = (
            actual+1 != factor*primitive*fmpq(scale.numerator,scale.denominator))
    require(set(resultants) == {10,12,14},"A generic resultant is missing.")
    trial = prime_table(isqrt(2**31-1))
    prime_cache = {"trial":trial,"verified":set()}
    residuals = []
    crt_report = []
    for entry in certificate["generic_resultants"]["eliminants"]:
        m = int(entry["other"])
        require(m in (12,14),"Unexpected large resultant.")
        require(int(entry["P1_exponent"]) == (120 if m == 12 else 140),
                "A large resultant pivot exponent changed.")
        residual,report = verify_large_resultant(
            resultants[10],resultants[m],entry,controls,prime_cache)
        residuals.append(residual)
        crt_report.append(report)
    require([r["other_moment"] for r in crt_report] == [12,14],
            "Large resultant order changed.")
    gauss = certificate["generic_resultants"]["gauss"]
    require([len(v)-1 for v in residuals] ==
            [int(gauss["first_degree"]),int(gauss["second_degree"])],
            "Generic residual degrees changed.")
    verify_modular_coprimeness(residuals[0],residuals[1],
        int(gauss["modulus"]),gauss["bezout_first"],gauss["bezout_second"],
        controls,"generic_Gauss")

    # Lambda zero on the generic q chart.
    zc = fmpq_mpoly_ctx.get(("q0","z"),ordering="lex")
    qz,z = zc.gens()
    zp1 = 7*qz+3
    zp3 = 1615*qz**3+2261*qz*qz+1197*qz+231
    zp5 = (930465*qz**5+2225025*qz**4+2314026*qz**3+
           1285570*qz*qz+376805*qz+46189)
    zpolys = {}
    for m,n in numerators.items():
        f = n.subs({"lam":0})
        par = (m//2)%2
        dd = {}
        for ex,c in f.to_dict().items():
            ib,iq,il = map(int,ex)
            require(il == 0 and ib >= par and (ib-par)%2 == 0,
                    "A lambda-zero scaling failed.")
            dd[(iq,(ib-par)//2)] = c
        zpolys[m] = zc.from_dict(dd)
    ldata = certificate["lambda0_generic"]
    scalar = {}
    for row in ldata["resultant_Bezout_identities"]:
        m = int(row["other_moment"])
        require(m in (10,12),"Unexpected lambda-zero resultant.")
        u = sum((poly(c,zc)*z**i for i,c in
                 enumerate(row["U_coefficients_ascending_z"])),zc.constant(0))
        v = sum((poly(c,zc)*z**i for i,c in
                 enumerate(row["V_coefficients_ascending_z"])),zc.constant(0))
        r = poly(row["scalar"],zc)
        require(u*zpolys[8]+v*zpolys[m] == r and bool(r),
                "A lambda-zero literal Sylvester identity failed.")
        qq_dense_in_q(r,0)
        actual = zpolys[8].resultant(zpolys[m],"z")
        require(r == actual or r == -actual,
                "A lambda-zero scalar is not the signed resultant.")
        scalar[m] = r
        controls["lambda0_%s_cofactor_corruption_rejected"%m] = (
            (u+1)*zpolys[8]+v*zpolys[m] != r)
    pairs = (
        ("b1_nonzero_coprimeness",scalar[10],scalar[12],zp1**8*zp3**4*zp5),
        ("b1_zero_coprimeness",zpolys[8].subs({"z":0}),
         zpolys[12].subs({"z":0}),zp1**4*zp5),
    )
    for name,f,g,pivot in pairs:
        entry = ldata[name]
        require(poly(entry["removed_common_pivot"],zc) == pivot,
                "A lambda-zero removed pivot changed.")
        require(not f%pivot and not g%pivot,
                "A lambda-zero scalar pivot does not divide exactly.")
        aa,bb = [primitive_integer_coefficients(v//pivot,0) for v in (f,g)]
        verify_modular_coprimeness(aa,bb,int(entry["prime"]),
            entry["modular_u_coefficients_ascending"],
            entry["modular_v_coefficients_ascending"],controls,"lambda0_"+name)

    # The P1 exceptional fiber, with a5 freshly retained and solved.
    c1 = fmpq_mpoly_ctx.get(("a1","a3","a5","lam"),ordering="lex")
    aa,cc,ee,lg = c1.gens()
    bspecial = -10*lg/11
    qspecial = constant(c1,Fraction(-3,7))
    f4 = moments[4].compose(aa,cc,ee,bspecial,qspecial,lg)
    csol = constant_linear_solution(f4,1)
    f6 = moments[6].compose(aa,csol,ee,bspecial,qspecial,lg)
    esol = constant_linear_solution(f6,2)
    i1 = {m:moments[m].compose(aa,csol,esol,bspecial,qspecial,lg)
          for m in range(2,15,2)}
    require(all(not i1[m] for m in (2,4,6)),
            "P1 native pivot substitutions failed.")
    cofs,unit = cofactor_group(certificate["exceptional_cofactors"]["P1"],c1)
    order = (8,10,12,14)
    require(len(cofs) == len(order),"P1 cofactor count changed.")
    answer = sum((u*i1[m] for u,m in zip(cofs,order)),c1.constant(0))
    require(answer == unit,"P1 complete native cofactor identity failed.")
    controls["P1_cofactor_corruption_rejected"] = answer+i1[8] != unit

    # The P3 exceptional fibers: retain M4 compatibility and solve a5.
    c3 = fmpq_mpoly_ctx.get(("a3","a5","b1","lam","q0"),ordering="lex")
    cg,eg,bg,lg,qg = c3.gens()
    mod3 = 1615*qg**3+2261*qg*qg+1197*qg+231
    n,d = solutions[0]
    asol = n.compose(bg,qg,lg)*field_inverse(d.compose(bg,qg,lg),mod3,4)%mod3
    f6 = moments[6].compose(asol,cg,eg,bg,qg,lg)%mod3
    coef = pure_linear_coefficient(f6,1)
    esol = -f6.subs({"a5":0})*field_inverse(coef,mod3,4)%mod3
    i3 = {m:moments[m].compose(asol,cg,esol,bg,qg,lg)%mod3
          for m in range(2,15,2)}
    require(not i3[2] and not i3[6],"P3 native solved moments failed.")
    cofs,unit = cofactor_group(certificate["exceptional_cofactors"]["P3"],c3)
    order = (4,8,10,12,14)
    require(len(cofs) == len(order),"P3 cofactor count changed.")
    answer = sum((u*i3[m] for u,m in zip(cofs,order)),c3.constant(0))
    require((answer-unit)%mod3 == 0,"P3 native quotient cofactor identity failed.")
    controls["P3_cofactor_corruption_rejected"] = (answer+i3[4]-unit)%mod3 != 0

    # The P5 fibers: native a1,a3 solutions; M6 compatibility is retained.
    c5 = fmpq_mpoly_ctx.get(("a5","b1","lam","q0"),ordering="lex")
    eg,bg,lg,qg = c5.gens()
    mod5 = (930465*qg**5+2225025*qg**4+2314026*qg**3+
            1285570*qg*qg+376805*qg+46189)
    branch = []
    for n,d in solutions[:2]:
        branch.append(n.compose(bg,qg,lg)*
                      field_inverse(d.compose(bg,qg,lg),mod5,3)%mod5)
    i5 = {m:moments[m].compose(branch[0],branch[1],eg,bg,qg,lg)%mod5
          for m in range(2,15,2)}
    require(not i5[2] and not i5[4],"P5 native solved moments failed.")
    cofs,unit = cofactor_group(certificate["exceptional_cofactors"]["P5_lambda0"],c5)
    order = (6,8,10,12,14)
    require(len(cofs) == len(order),"P5 lambda-zero cofactor count changed.")
    answer = sum((u*i5[m].subs({"lam":0}) for u,m in zip(cofs,order)),c5.constant(0))
    require((answer-unit)%mod5 == 0,"P5 lambda-zero native identity failed.")
    controls["P5_lambda0_cofactor_corruption_rejected"] = (
        answer+i5[6].subs({"lam":0})-unit)%mod5 != 0

    # P5, lambda nonzero: split EVERY h-pivot before solving z and w.
    wc = fmpq_mpoly_ctx.get(("w","h","z","q0"),ordering="lex")
    w,h,z,qw = wc.gens()
    pmod = (930465*qw**5+2225025*qw**4+2314026*qw**3+
            1285570*qw*qw+376805*qw+46189)
    gs = {}
    for m in (6,8,10,12,14):
        par = (m//2)%2
        dd = {}
        for ex,c in i5[m].to_dict().items():
            ia,ib,il,iq = map(int,ex)
            power = ia+ib+il
            require(power >= par and (power-par)%2 == 0,
                    "P5 h,z,w parity transformation failed.")
            dd[(ia,ib,(power-par)//2,iq)] = c
        gs[m] = wc.from_dict(dd)%pmod
    require(all(e[0] == 0 and e[2] <= 1 for e in gs[6].to_dict()),
            "P5 G6 is not independent of w and affine in z.")
    A,B = coefficient(gs[6],2,1),coefficient(gs[6],2,0)
    require(all(e[0] <= 1 for e in gs[8].to_dict()),"P5 G8 is not affine in w.")
    wco = coefficient(gs[8],0,1)
    require(bool(wco) and all(e[2] == 1 for e in wco.to_dict()),
            "The complete G8 w coefficient is not z times C(h).")
    C = wco//z
    E = gs[8].subs({"w":0})
    require(int(E.degrees()[2]) <= 2,"The G8 remainder is not quadratic in z.")
    E2,E1,E0 = [coefficient(E,2,j) for j in (2,1,0)]
    R8 = (E2*B*B-E1*A*B+E0*A*A)%pmod
    pivot = A*B*C%pmod
    require(gs[6] == A*z+B and gs[8] == z*C*w+E,
            "P5 structural G6/G8 identities failed.")
    pd = certificate["P5_hzw"]
    for label,f,g in (("AB",A,B),("CR",C,R8)):
        cofs,unit = cofactor_group(pd[label],wc)
        require(len(cofs) == 2,"A P5 h-pivot cofactor is missing.")
        ans = cofs[0]*f+cofs[1]*g
        require((ans-unit)%pmod == 0,"A complete P5 h-pivot unit identity failed.")
        controls["P5_%s_cofactor_corruption_rejected"%label] = (ans+f-unit)%pmod != 0
    hs = {}
    for m in (10,12,14):
        f = gs[m]
        aw,az = int(f.degrees()[0]),int(f.degrees()[2])
        ans = wc.constant(0)
        for ex,c in f.to_dict().items():
            iw,ih,iz,iq = map(int,ex)
            term = (c*R8**iw*(A*B*C)**(aw-iw)*(-B)**iz*
                    A**(az-iz)*h**ih*qw**iq)
            ans = (ans+term)%pmod
        require(all(e[0] == e[2] == 0 for e in ans.to_dict()),
                "A cleared P5 equation still contains w or z.")
        hs[m] = ans
    cofs,multiplier = cofactor_group(pd["H"],wc)
    require(len(cofs) == 3,"A P5 univariate cofactor is missing.")
    J = poly(pd["J"],wc)
    quotient = poly(pd["pivot_power_quotient"],wc)
    power = int(pd["pivot_power"])
    require(power == 3 and bool(J),"The P5 pivot-power certificate changed.")
    ans = sum((u*hs[m] for u,m in zip(cofs,(10,12,14))),wc.constant(0))
    require((ans-multiplier*J)%pmod == 0,"The P5 H combination failed.")
    require((quotient*J-pivot**power)%pmod == 0,"P5 J does not divide the pivot power.")
    controls["P5_H_cofactor_corruption_rejected"] = (
        ans+hs[10]-multiplier*J)%pmod != 0
    controls["P5_pivot_quotient_corruption_rejected"] = (
        (quotient+1)*J-pivot**power)%pmod != 0

    # Sharpness: a smooth native prefix-12 point, not a bare modular proper ideal.
    sh = certificate["sharpness"]
    require(tuple(sh["parameter_order"]) == NATIVE_NAMES,"Sharpness coordinate order changed.")
    p = int(sh["prime"])
    require(is_prime(p),"Sharpness modulus is not prime.")
    point = list(map(int,sh["point"]))
    require(len(point) == 6 and point[2]%p != 0,"The sharpness point loses a5.")
    orders = (2,4,6,8,10,12)
    require(tuple(sh["zero_even_orders"]) == orders,"Sharpness moment orders changed.")
    residues = [evaluate_mod(moments[m],point,p) for m in orders]
    require(residues == [0]*6,"The native prefix does not vanish at the sharpness point.")
    for m in (1,3,5,7,9,11,13):
        require(evaluate_mod(moments[m],point,p) == 0,"An odd sharpness moment is nonzero.")
    m14 = evaluate_mod(moments[14],point,p)
    require(m14 == int(sh["M14_residue"]) and m14 != 0,
            "The sharpness M14 residue is not a nonzero unit.")
    jac = [[evaluate_mod(derivative(moments[m],j),point,p)
            for j in range(6)] for m in orders]
    jacdet = determinant_mod(jac,p)
    require(jacdet == int(sh["native_Jacobian_determinant_residue"]) and jacdet != 0,
            "The full native sharpness Jacobian is singular.")
    bad = [row[:] for row in jac]
    for row in bad:
        row[0] = 0
    controls["sharpness_singular_Jacobian_rejected"] = determinant_mod(bad,p) == 0
    controls["sharpness_wrong_M14_residue_rejected"] = m14 != (m14+1)%p
    require(all(controls.values()),"A corruption-rejection control failed.")
    return {
        "status":"PASS",
        "claim":"For the normalized eight-monomial odd-angular quintic, no characteristic-zero parameter point has native moments M1..M14 all zero.",
        "native_moment_term_counts":{str(m):len(moments[m]) for m in range(1,15)},
        "generic_resultant_first_stage":first_stage_report,
        "generic_resultant_CRT":crt_report,
        "generic_Gauss_degrees":[len(a)-1 for a in residuals],
        "all_three_q_pivot_fibers_verified":True,
        "lambda_zero_chart_verified":True,
        "P5_lambda_nonzero_all_h_pivots_verified":True,
        "sharpness":{"prime":p,"point":point,"first13_moments_zero_mod_p":True,
                     "M14_residue":m14,"full_native_Jacobian_determinant":jacdet,
                     "a5_is_unit":True,
                     "interpretation":"Written Newton/Hensel proof gives a characteristic-zero first13-zero, M14-nonzero phase; it is not a bare modular unit inference."},
        "corruption_controls":controls,
        "seconds":round(time.monotonic()-started,3),
        "runtime_Singular_SymPy_network_or_private_paths":False,
        "scope_limits":["Exactly the displayed normalized eight-monomial family.",
                        "The written coefficient-normalization corollary requires B*C nonzero.",
                        "No full even14/S24 or global degree-five classification."],
        "assurance":{"independent_replay":False,"formal_verification":False,
                     "external_human_review":False},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate",type=Path)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    if hasattr(sys,"set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    if args.certificate is None:
        args.certificate = Path(__file__).with_name("odd-angular-quintic-finite14-certificate.json.gz")
    raw = args.certificate.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    require(digest == CERTIFICATE_SHA256,"Compressed certificate hash mismatch.")
    data = json.loads(gzip.decompress(raw))
    require(data.get("schema") == 1,"Unsupported certificate schema.")
    result = verify(data)
    result["certificate_sha256"] = digest
    result["checker_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result["compressed_input_corruption_rejected"] = (
        hashlib.sha256(raw+b"\0").hexdigest() != CERTIFICATE_SHA256)
    text = json.dumps(result,indent=2,sort_keys=True)+"\n"
    if args.output:
        require(not args.output.exists(),"Refusing to overwrite a receipt.")
        args.output.write_text(text)
    print(text,end="")


if __name__ == "__main__":
    main()
