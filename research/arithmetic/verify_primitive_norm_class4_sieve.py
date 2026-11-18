#!/usr/bin/env python3
"""Finite identities supporting the C4 class-error sieve for x^2+14*y^2.

The infinite class-error theorem requires its separate arithmetic and analytic
proof. This checker certifies bounded coordinate/local comparisons, square
envelopes, exact CRT roots, and finite rational Selberg/Euler identities only.
Standard library only.

python3 research/arithmetic/verify_primitive_norm_class4_sieve.py \
  --output research/arithmetic/primitive_norm_class4_sieve_checks.json
"""
from __future__ import annotations

import argparse
from array import array
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import gcd, isqrt, prod
from pathlib import Path
import json

PREFIXES = (1, 2, 7, 14)
D_SAMPLES = (1, 3, 5, 9, 13, 15, 65)
GOOD_PRIMES = (3, 5, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
TOY_PRIMES = (3, 5, 11, 17, 23)


def require(condition, detail="check failed"):
    if not condition:
        raise AssertionError(detail)


def rejected(action):
    try:
        action()
    except AssertionError:
        return True
    raise AssertionError("deliberate corruption was accepted")


def psi8(n):
    require(n % 2 == 1)
    return 1 if n % 8 in (1, 7) else -1


def prime_kind(p):
    """Periodic characters only: H even split, O odd split, I inert."""
    require(p not in (2, 7))
    chi = psi8(p) * (1 if p % 7 in (1, 2, 4) else -1)
    if chi == -1:
        return "I"
    return "H" if psi8(p) == 1 else "O"


def prime_tables(bound):
    flags = bytearray(b"\1") * (bound + 1)
    flags[:2] = b"\0\0"
    for p in range(2, isqrt(bound) + 1):
        if flags[p]:
            flags[p*p:bound+1:p] = b"\0" * ((bound-p*p)//p+1)
    primes = [p for p in range(2, bound + 1) if flags[p]]
    spf = array("I", range(bound + 1))
    for p in primes:
        if p*p > bound:
            break
        for k in range(p*p, bound + 1, p):
            if spf[k] == k:
                spf[k] = p
    return primes, spf


def factor(n, spf):
    out = []
    while n > 1:
        p, e = spf[n], 0
        while n % p == 0:
            n //= p
            e += 1
        out.append((p, e))
    return out


def coordinate_flags(bound):
    flags = bytearray(bound + 1)
    for y in range(isqrt(bound//14) + 1):
        for x in range(isqrt(bound-14*y*y) + 1):
            if gcd(x, y) == 1:
                n = x*x + 14*y*y
                if n:
                    flags[n] = 1
    return flags


def local_from_factors(n, factors):
    """No principal ideal membership formula is used."""
    core = n
    for p, e in factors:
        if p in (2, 7):
            if e > 1:
                return False
            core //= p**e
        elif prime_kind(p) == "I":
            return False
    return psi8(core) == 1


def square_envelope(n, factors):
    """Assert the necessary envelope for a hidden failure and return r,d,c."""
    r, d, c = 1, 1, 1
    for p, e in factors:
        if p in (2, 7):
            require(e == 1)
            r *= p
        elif prime_kind(p) == "O":
            require(e % 2 == 0, ("active odd-class prime", n, p, e))
            d *= p**(e//2)
        else:
            require(prime_kind(p) == "H")
            c *= p**e
    require(r in PREFIXES and r*d*d*c == n)
    require(gcd(d, 14*c) == 1)
    # Exact exponents in the three disjoint prime supports force uniqueness.
    for p, e in factors:
        if p in (2, 7):
            require(r % p == 0 and d % p and c % p)
        elif prime_kind(p) == "O":
            require(c % p and d % (p**(e//2)) == 0)
        else:
            require(d % p and c % (p**e) == 0)
    return r, d, c


def arithmetic_checks(bound):
    primes, spf = prime_tables(bound)
    flags = coordinate_flags(bound)
    local = bytearray(bound + 1)
    counts, examples, d_values = Counter(), [], set()
    for n in range(1, bound + 1):
        factors = factor(n, spf)
        local[n] = local_from_factors(n, factors)
        require(not flags[n] or local[n], ("global not local", n))
        if flags[n]:
            counts["principal_primitive_values"] += 1
        if local[n]:
            counts["individual_local_candidates"] += 1
        if local[n] and not flags[n]:
            r, d, c = square_envelope(n, factors)
            counts["hidden_class_failures"] += 1
            counts["envelope_prefix_" + str(r)] += 1
            d_values.add(d)
            if len(examples) < 12:
                examples.append([n, r, d, c])
    require(local[9] and not flags[9])
    require(local[450] and not flags[450])
    require(flags[15] and local[15])
    controls = []
    rejected(lambda: square_envelope(15, factor(15, spf)))
    controls.append("principal norm15 cannot be inserted as an inactive-odd-prime envelope")
    rejected(lambda: require(local[3]))
    controls.append("wrong genus admits3")
    rejected(lambda: require(not local[9]))
    controls.append("hidden class failure9 cannot be erased from the local candidates")
    rejected(lambda: require(square_envelope(9, factor(9, spf)) == (1, 1, 9)))
    controls.append("moving an O-prime square into the H-core is rejected")
    counts["distinct_square_factors_among_hidden_failures"] = len(d_values)
    return {
        "range": [1, bound],
        "inputs": "direct gcd-restricted coordinates and independent periodic chi_-56/psi_8 local conditions; no CN4 membership formula",
        "counts": dict(sorted(counts.items())),
        "first_envelopes_n_r_d_c": examples,
        "negative_controls": controls,
    }


def progression(r, s, d, target):
    g = gcd(r, s)
    if target % g:
        return None
    a, b, n = (r//g)*d*d, s//g, target//g
    require(gcd(a, b) == 1)
    c0 = 0 if b == 1 else (n * pow(a, -1, b)) % b
    v0 = (n-a*c0)//b
    require(a*c0+b*v0 == n)
    return g, a, b, n, c0, v0


def expected_nu(p, target):
    k = prime_kind(p)
    return 2-int(target % p == 0) if k == "I" else 1 if k == "O" else 0


def roots(p, target, data):
    _, a, b, _, c0, v0 = data
    kind = prime_kind(p)
    if kind == "H":
        return set()
    return {
        k for k in range(p)
        if (c0+b*k) % p == 0
        or (kind == "I" and (v0-a*k) % p == 0)
    }


def ceil_div(a, b):
    return -((-a)//b)


def progression_and_root_checks():
    for d in D_SAMPLES:
        k = d
        for q in (3, 5, 13):
            while k % q == 0:
                require(prime_kind(q) == "O")
                k //= q
        require(k == 1)
    root_cases, deep_cases, bad_gcd_cases, positive_cases = 0, 0, 0, 0
    slope_zero_cases, target_collision_cases = 0, 0
    for r, s, d in product(PREFIXES, PREFIXES, D_SAMPLES):
        g, b = gcd(r, s), s//gcd(r, s)
        for rem in range(g):
            if rem:
                require(progression(r, s, d, rem) is None)
                bad_gcd_cases += 1
        for p in GOOD_PRIMES:
            # n modulo b*p covers all combined root-position residues:
            # c0 depends on n modulo b; collisions depend on n modulo p.
            for n in range(b*p):
                target = g*n
                data = progression(r, s, d, target)
                actual = roots(p, target, data)
                require(len(actual) == expected_nu(p, target), (r,s,d,p,target,actual))
                if prime_kind(p) == "I":
                    require(data[1] % p and data[2] % p)
                    if target % p == 0:
                        target_collision_cases += 1
                if d % p == 0:
                    require(prime_kind(p) == "O")
                    require(data[1] % p == 0 and data[2] % p)
                    slope_zero_cases += 1
                root_cases += 1
            for n in (p, p*p, 2*p*p, p*p*p):
                target = g*n
                data = progression(r, s, d, target)
                require(len(roots(p,target,data)) == expected_nu(p,target))
                deep_cases += 1
        for target in (1,2,7,14,15,18,49,98,105,450,999,1000,1729,10000):
            direct = {
                (c, (target-r*d*d*c)//s)
                for c in range(1, target//(r*d*d)+1)
                if target-r*d*d*c >= s and (target-r*d*d*c) % s == 0
            }
            data = progression(r,s,d,target)
            if data is None:
                require(not direct)
            else:
                _,a,b,_,c0,v0 = data
                lo, hi = ceil_div(1-c0,b), (v0-1)//a
                parameterized = {(c0+b*k, v0-a*k) for k in range(lo,hi+1)}
                require(parameterized == direct, (r,s,d,target))
            positive_cases += 1
    # Purposeful target-collision and coefficient-divisor mutants.
    datum = progression(1,1,1,11)
    require(prime_kind(11) == "I")
    rejected(lambda: require(len(roots(11,11,datum)) == 2))
    datum = progression(1,2,3,22)
    require(prime_kind(3) == "O" and datum[1] % 3 == 0)
    rejected(lambda: require(len(roots(3,22,datum)) == 0))
    rejected(lambda: require(datum[1] % 3 != 0))
    return {
        "prefixes": list(PREFIXES), "d_samples": list(D_SAMPLES),
        "good_primes": list(GOOD_PRIMES),
        "all_combined_root_position_cases": root_cases,
        "additional_deep_target_cases": deep_cases,
        "inert_target_collision_cases": target_collision_cases,
        "O_prime_zero_second_slope_cases": slope_zero_cases,
        "incompatible_prefix_gcd_cases": bad_gcd_cases,
        "positive_progression_vs_direct_equation_cases": positive_cases,
        "negative_controls": [
            "wrong inert target collision retains two roots",
            "wrong omission of q|d removes the O-prime c-root",
            "wrong coefficient claim makes the second slope a unit at q|d",
        ],
    }


def squarefree_products(primes):
    items = [(1, ())]
    for p in primes:
        items += [(h*p, factors+(p,)) for h,factors in items]
    return sorted(items)


def mobius_squarefree(h, primes):
    count = 0
    for p in primes:
        if h % p == 0:
            count += 1
            h //= p
    require(h == 1)
    return -1 if count % 2 else 1


def multiplicative_value(h, prime_values):
    ans = Fraction(1)
    for p, value in prime_values.items():
        if h % p == 0:
            ans *= value
            h //= p
    require(h == 1)
    return ans


def selberg_data(target, z):
    pp = [p for p in TOY_PRIMES if expected_nu(p,target)]
    fprime = {p: Fraction(expected_nu(p,target),p) for p in pp}
    wprime = {p: Fraction(expected_nu(p,target),p-expected_nu(p,target)) for p in pp}
    sf = squarefree_products(pp)
    support = [h for h,_ in sf if h <= z]
    f = lambda h: multiplicative_value(h,fprime)
    w = lambda h: multiplicative_value(h,wprime)
    mu = lambda h: mobius_squarefree(h,pp)
    G = sum((w(h) for h in support),Fraction())
    weights = {}
    for h in support:
        Gh = sum((w(k) for k in support if k <= z//h and gcd(k,h)==1),Fraction())
        weights[h] = mu(h)*w(h)/f(h)*Gh/G
    require(weights[1] == 1)
    return pp,sf,support,f,w,mu,G,weights


def check_quadratic(support, f, w, mu, G, weights):
    require(weights[1] == 1)
    direct = sum((
        weights[h]*weights[k]*f(h*k//gcd(h,k))
        for h in support for k in support
    ),Fraction())
    yy = {
        r: sum((weights[h]*f(h) for h in support if h%r==0),Fraction())
        for r in support
    }
    diagonal = sum((yy[r]*yy[r]/w(r) for r in support),Fraction())
    require(direct == diagonal == 1/G)
    require(sum((mu(r)*yy[r] for r in support),Fraction()) == 1)
    for h in support:
        inverse = sum((
            mu(r//h)*yy[r] for r in support if r%h==0
        ),Fraction())
        require(inverse == weights[h]*f(h))
    return direct


def selberg_checks():
    quadratic_cases, interval_cases, crt_cases, lambda_mutants = 0,0,0,0
    targets = (1,11,17,121,187)
    z_values = (1,3,5,11,17,33,55)
    for target,z in product(targets,z_values):
        pp,sf,support,f,w,mu,G,weights = selberg_data(target,z)
        Q = check_quadratic(support,f,w,mu,G,weights)
        quadratic_cases += 1
        wrong = dict(weights)
        wrong[1] += Fraction(1,2)
        rejected(lambda: check_quadratic(support,f,w,mu,G,wrong))
        lambda_mutants += 1
        if len(support)>1:
            wrong = dict(weights)
            wrong[support[1]] += Fraction(1,2)
            rejected(lambda: check_quadratic(support,f,w,mu,G,wrong))
            lambda_mutants += 1
        for r,s,d in ((1,1,1),(1,2,3),(14,7,5)):
            actual_target = gcd(r,s)*target
            data = progression(r,s,d,actual_target)
            omega = {p: roots(p,actual_target,data) for p in pp}
            active = [p for p in pp if p <= z]
            nu = lambda h: prod(len(omega[p]) for p in pp if h%p==0)
            relevant = {h*k//gcd(h,k) for h in support for k in support}
            for start,length in ((-17,1),(-11,7),(0,31),(5,19),(23,43)):
                interval = range(start,start+length)
                def in_omega(k,h):
                    return all(k%p in omega[p] for p in pp if h%p==0)
                survivor = sum(all(k%p not in omega[p] for p in active) for k in interval)
                squared = sum((
                    sum((weights[h] for h in support if in_omega(k,h)),Fraction())**2
                    for k in interval
                ),Fraction())
                require(survivor <= squared)
                for h in relevant:
                    actual = sum(in_omega(k,h) for k in interval)
                    require(abs(Fraction(actual)-length*f(h)) <= nu(h))
                    crt_cases += 1
                endpoint = sum((
                    abs(weights[h]*weights[k])*nu(h*k//gcd(h,k))
                    for h in support for k in support
                ),Fraction())
                require(abs(squared-length*Q) <= endpoint)
                require(survivor <= Fraction(length,G)+endpoint)
                interval_cases += 1
    return {
        "toy_prime_universe": list(TOY_PRIMES),
        "target_samples": list(targets), "z_samples": list(z_values),
        "exact_quadratic_diagonal_and_inversion_cases": quadratic_cases,
        "integer_interval_sieve_upper_cases": interval_cases,
        "individual_CRT_endpoint_checks": crt_cases,
        "lambda_corruptions_rejected": lambda_mutants,
        "arithmetic": "fractions.Fraction throughout; direct quadratic matrix, independent diagonal sums, and literal interval enumeration",
    }


def euler_checks():
    cases, marginal_checks, markov_cases = 0,0,0
    for target in (1,11,17,121,187):
        pp = [p for p in TOY_PRIMES if expected_nu(p,target)]
        wp = {p:Fraction(expected_nu(p,target),p-expected_nu(p,target)) for p in pp}
        allh = squarefree_products(pp)
        weights = {h:multiplicative_value(h,wp) for h,_ in allh}
        total = sum(weights.values(),Fraction())
        euler = prod((1+wp[p] for p in pp), start=Fraction(1))
        require(total == euler)
        for p in pp:
            marginal = sum((w for h,w in weights.items() if h%p==0),Fraction())/total
            require(marginal == Fraction(expected_nu(p,target),p))
            marginal_checks += 1
        base = prod((
            Fraction(p,p-(2 if prime_kind(p)=="I" else 1))
            for p in pp
        ),start=Fraction(1))
        correction = prod((
            Fraction(p-2,p-1) for p in pp
            if prime_kind(p)=="I" and target%p==0
        ),start=Fraction(1))
        require(total == base*correction)
        # Exact rational majorant for the log-size Markov argument:
        # log_2(h) <= sum_{p|h} ceil(log_2 p), while h>2^k gives
        # log_2(h)>k. Choose k >= twice the expected upper bound.
        mean_upper = sum((
            Fraction(expected_nu(p,target),p)*p.bit_length() for p in pp
        ),Fraction())
        k = max(1,ceil_div(2*mean_upper.numerator,mean_upper.denominator))
        z = 2**k
        tail = sum((w for h,w in weights.items() if h>z),Fraction())/total
        require(tail <= mean_upper/k <= Fraction(1,2))
        truncated = sum((w for h,w in weights.items() if h<=z),Fraction())
        require(truncated >= total/2)
        markov_cases += 1
        cases += 1
    # A wrong target-correction direction is detected in a nontrivial case.
    require(Fraction(11-2,11-1) < 1)
    rejected(lambda: require(Fraction(11-2,11-1) == Fraction(11-1,11-2)))
    return {
        "Euler_normalization_cases": cases,
        "exact_prime_marginals_nu_over_q_checked": marginal_checks,
        "rational_log_size_Markov_half_mass_cases": markov_cases,
        "negative_controls": ["inverting the inert target-correction factor"],
    }


def run(bound):
    require(bound >= 1000)
    result = {
        "schema": "primitive-norm-class4-sieve-finite-check-v1",
        "scope": "Finite identities only. No finite extrapolation to the infinite pointwise class-error estimate, AP Mertens theorem, or eventual additive representability.",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "arithmetic_envelopes": arithmetic_checks(bound),
        "linear_forms_and_roots": progression_and_root_checks(),
        "rational_Selberg": selberg_checks(),
        "Euler_measure": euler_checks(),
        "status": "pass",
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bound",type=int,default=1_000_000)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    result = run(args.bound)
    text = json.dumps(result,indent=2,sort_keys=True)+"\n"
    if args.output:
        args.output.write_text(text,encoding="utf-8")
    else:
        print(text,end="")
    print(json.dumps({
        "status":result["status"],"bound":args.bound,
        "source_sha256":result["source_sha256"],
        "arithmetic_counts":result["arithmetic_envelopes"]["counts"],
        "quadratic_cases":result["rational_Selberg"]["exact_quadratic_diagonal_and_inversion_cases"],
        "root_cases":result["linear_forms_and_roots"]["all_combined_root_position_cases"],
        "output":str(args.output) if args.output else None,
    },sort_keys=True))


if __name__ == "__main__":
    main()
