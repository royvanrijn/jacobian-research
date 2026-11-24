#!/usr/bin/env python3
"""Exact support for refined bounds at discriminants -43, -67, -163.

This checks analytic certificates, local identities and a small initial sumset.
It does NOT cover the full finite domains needed for exact classification.
The written proof supplies the classical spectral input and infinite reductions.
"""
if not __debug__:
    raise RuntimeError('Active Python assertions are required; do not use -O or -OO.')

import argparse
import hashlib
import json
from bisect import bisect_right
from fractions import Fraction as F
from itertools import product
from math import gcd, isqrt, prod
from pathlib import Path

DS = (43, 67, 163)
CORE_EXPONENT = {43: 32, 67: 35, 163: 45}
M_CONSTANT = {43: 4, 67: 3, 163: 2}
DELTA_CONSTANT = {43: 7, 67: 9, 163: 14}
TAIL_EXPONENT = 96
R = 1 << 20
DEN = 1 << 80
EXCEPTIONS = {
    43: [4, 6, 8, 10, 16, 20, 38, 50, 56, 74],
    67: [4, 6, 8, 10, 12, 14, 16, 22, 26, 28, 32, 44, 50, 62, 68, 80, 86, 98, 114, 116, 138],
    163: [4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 46, 50, 52, 56, 58, 60, 64, 66, 68, 70, 74, 76, 78, 80, 92, 110, 116, 120, 128, 134, 146, 148, 162, 164, 170, 176, 182, 186, 188, 190, 206, 218, 230, 236, 254, 272, 278, 290, 300, 332, 338, 344, 362, 386, 392, 484, 568, 992],
}


def rational(x):
    x = F(x)
    return {'numerator': x.numerator, 'denominator': x.denominator}


def sha(x):
    return hashlib.sha256(x).hexdigest()


def primes(limit):
    flags = bytearray(b'\x01') * (limit + 1)
    flags[:2] = b'\0\0'
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p*p::p] = b'\0' * ((limit-p*p)//p+1)
    return [p for p, flag in enumerate(flags) if flag]


def sigma_power(p, a):
    return (p**(a+1)-1)//(p-1)


def root_floor(n, k):
    assert n >= 0 and k >= 1
    lo, hi = 0, 1 << ((n.bit_length()+k-1)//k)
    while lo+1 < hi:
        mid = (lo+hi)//2
        if mid**k <= n:
            lo = mid
        else:
            hi = mid
    assert lo**k <= n < (lo+1)**k
    return lo


def sqrt_fraction(x):
    x = F(x)
    assert x >= 0
    q = isqrt(x.numerator//x.denominator)
    assert q*q <= x < (q+1)**2
    return q


def upper_sixth_root(x):
    # Strict rational upper bound with denominator 10000.
    z = root_floor(x.numerator*10000**6//x.denominator, 6)+1
    value = F(z, 10000)
    assert (value-F(1,10000))**6 <= x < value**6
    return value


def log_upper_ratio(a, b):
    assert a >= b > 0
    k = ((a+b-1)//b-1).bit_length()
    assert a <= b*(1 << k)
    if k:
        assert a > b*(1 << (k-1))
    return F(7*k, 10)  # log(2) < 7/10


def divisor_frontier(split, upper):
    prefix, chosen = 1, []
    for p in split:
        prefix *= p
        chosen.append(p)
        if prefix > upper:
            break
    assert prefix > upper
    records = []
    def visit(i, last, n, d):
        records.append((n, d))
        if i == len(chosen):
            return
        power = 1
        for a in range(1, last+1):
            power *= chosen[i]
            if n*power > upper:
                break
            visit(i+1, a, n*power, d*(a+1))
    visit(0, upper.bit_length(), 1, 1)
    records.sort()
    frontier, best = [], 0
    for n, d in records:
        if d > best:
            best = d
            frontier.append((n,d))
    return frontier, len(records), chosen, prefix


def euler_factor(p, a):
    return F((a+1)**6*p**(5*a), sigma_power(p,a)**6)


def exponent_cap(p):
    a = 0
    while (a+2)**6 >= p*(a+1)**6:
        a += 1
    assert (a+2)**6 < p*(a+1)**6
    return a


def delta_frontiers(inert, chosen, upper):
    variable = []
    for p in inert:
        if p >= 64:
            break
        cap = max(exponent_cap(p), 2 if p in chosen else 0)
        factors = [euler_factor(p,a) for a in range(cap+1)]
        if p in chosen or max(factors) > 1:
            variable.append((p,cap,factors))
    records = {state: [] for state in product(range(3), repeat=len(chosen))}
    def visit(i, n, value, state):
        if i == len(variable):
            for t in range(1,7):
                if n*2**t <= upper:
                    records[state].append((n*2**t,value*F((t+3)**6,3**6*2**t)))
            return
        p, cap, factors = variable[i]
        for a in range(cap+1):
            nn = n*p**a
            if nn > upper:
                break
            ss = list(state)
            if p in chosen:
                ss[chosen.index(p)] = min(a,2)
            visit(i+1,nn,value*factors[a],tuple(ss))
    visit(0,1,F(1),(0,)*len(chosen))
    frontiers = {}
    for state, rr in records.items():
        rr.sort()
        best = F(0)
        nn, ff = [], []
        for n, value in rr:
            if value > best:
                nn.append(n)
                ff.append(value)
                best = value
        frontiers[state] = (nn,ff)
    return frontiers, {
        'variable_primes': [{'p': p, 'cap': cap} for p,cap,_ in variable],
        'tuple_count': sum(map(len,records.values())),
        'frontier_count': sum(len(x[0]) for x in frontiers.values()),
        'state_count': len(records),
        'upper': upper,
    }


def factor(n):
    out = []
    p = 2
    while p*p <= n:
        if n%p == 0:
            a = 0
            while n%p == 0:
                n //= p
                a += 1
            out.append((p,a))
        p += 1
    if n > 1:
        out.append((n,1))
    return out


def state_of(n, chosen):
    return tuple(0 if n%p else 1 if n%(p*p) else 2 for p in chosen)


def thin_constant(e):
    fac = factor(e)
    t = dict(fac).get(2,0)
    odd = [(p,a) for p,a in fac if p != 2]
    m = prod(p**a for p,a in odd)
    return F(2*(t+3)*prod(a+1 for p,a in odd)*m,
             3*prod(sigma_power(p,a) for p,a in odd))


def expect_rejected(name, action):
    try:
        action()
    except (AssertionError, ValueError):
        return name
    raise AssertionError('Negative control accepted: '+name)


def equal(a,b):
    assert a == b


class BoundData:
    def __init__(self, D, all_primes):
        self.D = D
        self.qr = frozenset(i*i%D for i in range(1,D))
        self.inert = [p for p in all_primes if p > 2 and p%D not in self.qr and p != D]
        self.split = [p for p in all_primes if p%D in self.qr]
        self.chosen = self.inert[:5]
        self.w, self.nu = F(6,D-1), F(12,D*D-1)
        self.r = (D-1)//2
        self.residues = [a for a in range(1,2*D,2) if a%D not in self.qr and a != D]
        assert len(self.residues) == self.r
        self.B = sum((F(2,p*(p-1)) for p in self.inert if p <= 1024),F(0))+3*self.ap_tail(1024)
        sqrt_lower = F(isqrt(D*10**12),10**6)
        assert sqrt_lower**2 < D < (sqrt_lower+F(1,10**6))**2
        self.area, self.boundary = F(33,7)/sqrt_lower, 3+6/sqrt_lower
        self.square_suffix = [0]*(len(self.inert)+1)
        self.harmonic_prefix = [0]*(len(self.inert)+1)
        for i,p in enumerate(self.inert):
            self.harmonic_prefix[i+1] = self.harmonic_prefix[i]+(DEN+p-1)//p
        for i in range(len(self.inert)-1,-1,-1):
            p = self.inert[i]
            self.square_suffix[i] = self.square_suffix[i+1]+(DEN+p*p-1)//(p*p)
        self.saves = {}
        for state in product(range(3),repeat=5):
            save = sum((F(2,p*(p-1))-(F(1,p*p) if s == 2 else 0)
                        for p,s in zip(self.chosen,state) if s),F(0))
            es = sum((1-((1-self.w)/p if s == 2 else 0)
                      for p,s in zip(self.chosen,state) if s),F(0))
            self.saves[state] = (save,es)

    def ap_tail(self,P):
        return F(self.r,2*self.D*P)+F(self.r,P*P)

    def prime_bounds(self,P,X):
        if P <= R:
            j = bisect_right(self.inert,P)
            pi = F(j)
            T = F(self.square_suffix[j],DEN)+self.ap_tail(R)
        else:
            j = None
            pi = F(self.r*P,2*self.D)+self.r-F(sum(self.residues),2*self.D)
            T = self.ap_tail(P)
        if X <= P:
            return pi,T,F(0),F(0)
        if P <= R:
            h = bisect_right(self.inert,min(R,X))
            U = F(self.harmonic_prefix[h]-self.harmonic_prefix[j],DEN)
            V = F(h-j)
            if X > R:
                U += F(self.r,2*self.D)*log_upper_ratio(X,R)+F(self.r,R)
                V += F(self.r*(X-R),2*self.D)+self.r
        else:
            U = F(self.r,2*self.D)*log_upper_ratio(X,P)+F(self.r,P)
            V = F(self.r*(X-P),2*self.D)+self.r
        return pi,T,U,V

    def interval_cost(self,L,U,d,delta,state):
        # The cutoff choice need not be optimal. All acceptance is rational.
        A = delta-F(4*d,3*L)
        B = F(2*d,3*isqrt(L))*self.boundary
        C = F(2*d,3)*self.area
        if A > 0:
            radical = F(sqrt_fraction((B*B+4*A*C)*10**12),10**6)
            P = int((B+radical)/(2*A))
        else:
            P = isqrt(U)
        P = max(P,self.chosen[-1])
        pi,T,harmonic,count = self.prime_bounds(P,isqrt(U))
        save,es = self.saves[state]
        coefficient = 2-self.w-2*self.nu+pi-es
        assert coefficient > 0
        cost = (F(4,self.D+1)+self.B-save
                +F(2*d,3*self.w)*(self.area*T+self.boundary*harmonic/isqrt(L)+2*count/L)
                +coefficient*delta/self.w)
        return cost,P


def exact_small_controls(data):
    D, qr = data.D,data.qr
    bound = 1024
    c = [n for n in range(1,bound+1) if all(p%D in qr for p,a in factor(n))]
    sums = {a+b for a in c for b in c if a+b <= bound}
    E = [n for n in range(2,bound+1,2) if n%D and n not in sums]
    assert E == EXCEPTIONS[D]
    primitive = set(c)|{D*a for a in c if D*a <= bound}
    psums = {a+b for a in primitive for b in primitive if a+b <= bound}
    residual = [n for n in E if n not in psums]
    full_candidates = sorted(residual+[D*e for e in E])
    coordinate = set()
    counts = [0]*(bound+1)
    for y in range(-isqrt(4*bound//D),isqrt(4*bound//D)+1):
        s = isqrt(4*bound-D*y*y)
        for x in range(-((s+y)//2),(s-y)//2+1):
            q = x*x+x*y+(D+1)//4*y*y
            assert 0 <= q <= bound
            if q:
                counts[q] += 1
                if gcd(x,y) == 1:
                    coordinate.add(q)
    assert coordinate == primitive
    odd_count = 0
    lattice_rows = 0
    for X in range(bound+1):
        if X%2:
            odd_count += counts[X]
        # The actual area constant is irrational; this rational upper is safe.
        sqrt_upper = isqrt(X)+1
        assert odd_count <= data.area*X+data.boundary*sqrt_upper+2
        lattice_rows += 1
    # Direct shell counts test the target-dividing inert-prime branches.
    shell = [sum(counts[a]*counts[n-a] for a in range(1,n,2))
             for n in range(bound+1)]
    local_cases = 0
    for p in data.chosen[:3]:
        for n in range(2,257,2):
            if n%D == 0:
                continue
            bad = sum(counts[a]*counts[n-a] for a in range(1,n,2)
                      if a%p == 0 or (n-a)%p == 0)
            if n%p == 0:
                assert bad == (shell[n//(p*p)] if n%(p*p) == 0 else 0)
                local_cases += 1
    # Check the exact local S and Delta ratio inequalities, including a=2.
    for p in data.chosen:
        for a in range(2,11):
            assert F(sigma_power(p,a-2),sigma_power(p,a)) <= F(1,p*p)
            assert F(a-1,(a+1)*p) <= F(1,p)
    for e in E:
        fac = dict(factor(e))
        for p in data.split[:10]:
            a = fac.get(p,0)
            assert F((a+2)*p*sigma_power(p,a),(a+1)*sigma_power(p,a+1)) <= 2
    # Minkowski cutoffs: use pi>3; the only admissible small ideal norms
    # are 1 or 4, and the ideal of norm 4 is (2).
    minkowski_ceiling = {43:5,67:6,163:9}[D]
    assert 4*D < 9*minkowski_ceiling*minkowski_ceiling
    assert all(p%D not in qr for p in primes(minkowski_ceiling-1))
    negatives = [
        expect_rejected('character-only membership admits 9',lambda: equal(9 in c,9%D in qr)),
        expect_rejected('excluded core 4 claimed represented',lambda: equal(4 in sums,True)),
        expect_rejected('forgetting a known exception',lambda: equal(E,EXCEPTIONS[D][1:])),
        expect_rejected('negative sixth-root input',lambda: root_floor(-1,6)),
    ]
    # Exact power bookkeeping for the thin tail. The boundary/count terms
    # keep M/sqrt(n)=n^(-1/3), unlike the area-secondary M/z^2=n^(-1/2).
    gamma = F(data.r,2*D)
    h = F(7)
    direct_tail = 4*(data.area*gamma/16+data.area*data.r/16**2
                    +(data.boundary*h+2)/64)
    correct_tail = (data.area*gamma/4+data.area*data.r/64
                    +(data.boundary*h+2)/16)
    wrong_tail = data.area*gamma/4+(data.area*data.r+data.boundary*h+2)/64
    assert direct_tail == correct_tail
    negatives.append(expect_rejected('thin-tail boundary wrongly decays as n^(-1/2)',
                                      lambda: equal(direct_tail,wrong_tail)))
    bad18 = sum(counts[a]*counts[18-a] for a in range(1,18,2)
                if a%3 == 0 or (18-a)%3 == 0)
    assert bad18 == shell[2] == 4 and shell[6] == 0
    negatives.append(expect_rejected('inert scaling uses p instead of p squared',
                                      lambda: equal(bad18,shell[6])))
    # Some failures are exact witnesses, not exceptions from the checker.
    assert 9%D in qr and 9 not in c
    return {
        'small_sumset_bound':bound,'core_exceptions_through_bound':E,
        'candidate_primitive_exceptions':full_candidates,
        'candidate_primitive_exception_count':len(full_candidates),
        'coordinate_primitive_values':len(coordinate),
        'odd_lattice_controls':lattice_rows,
        'target_dividing_local_shell_cases':local_cases,
        'thin_constants':[{'e':e,'A':rational(thin_constant(e))} for e in E],
        'negative_controls':negatives,
    }


def certify(D,all_primes):
    data = BoundData(D,all_primes)
    small = exact_small_controls(data)
    split_sixth = prod(max(F((a+1)**6,p**a) for a in range(11))
                      for p in data.split if p < 64)
    assert 12**6 < 2*11**6
    assert split_sixth < M_CONSTANT[D]**6
    delta_sixth = F(3,2)**6*prod(max(euler_factor(p,a) for a in range(5))
                              for p in data.inert if p < 64)
    assert 6**6 < 3*5**6 and 9**6 > 2*8**6 and 10**6 < 2*9**6
    assert delta_sixth < DELTA_CONSTANT[D]**6
    frontier,enum_count,chosen,prefix = divisor_frontier(data.split,1 << TAIL_EXPONENT)
    m_keys = [n for n,d in frontier]
    df,df_info = delta_frontiers(data.inert,data.chosen,1 << TAIL_EXPONENT)
    sixth_cache = {}
    intervals = []
    count_core = count_thin = 0
    trace = hashlib.sha256()
    for k in range(CORE_EXPONENT[D],TAIL_EXPONENT):
        L,U = 1 << k,(1 << (k+1))-1
        d = frontier[bisect_right(m_keys,U)-1][1]
        rootL = root_floor(L,3)
        worst_core = (F(-1),None,None,None)
        for state,(keys,values) in df.items():
            j = bisect_right(keys,U)-1
            if j < 0:
                continue
            value = values[j]
            if value not in sixth_cache:
                sixth_cache[value] = upper_sixth_root(value)
            delta = sixth_cache[value]/rootL
            cost,P = data.interval_cost(L,U,d,delta,state)
            assert cost < F(49,50), ('core interval failed',D,k,state,float(cost),P)
            if cost > worst_core[0]:
                worst_core = (cost,state,P,delta)
            trace.update(f'C:{k}:{state}:{cost.numerator}/{cost.denominator}:{P}\n'.encode())
            count_core += 1
        worst_thin = (F(-1),None,None,None)
        for e in EXCEPTIONS[D]:
            delta = thin_constant(e)/isqrt(L)
            cost,P = data.interval_cost(L,U,d,delta,state_of(e,data.chosen))
            assert cost < F(7,10), ('thin interval failed',D,k,e,float(cost),P)
            if cost > worst_thin[0]:
                worst_thin = (cost,e,P,delta)
            trace.update(f'T:{k}:{e}:{cost.numerator}/{cost.denominator}:{P}\n'.encode())
            count_thin += 1
        intervals.append({
            'k':k,'L':L,'U':U,'M_upper':d,
            'worst_core':{'cost':rational(worst_core[0]),'state':worst_core[1],
                          'P':worst_core[2],'delta_upper':rational(worst_core[3])},
            'worst_thin':{'cost':rational(worst_thin[0]),'e':worst_thin[1],
                          'P':worst_thin[2],'delta_upper':rational(worst_thin[3])},
        })
    # Uniform infinite tails. n^(-1/12), n^(-1/3), n^(-1/6),
    # n^(-1/2) at n=2^96 are exact powers of two.
    w,nu,r = data.w,data.nu,data.r
    cm,cd = M_CONSTANT[D],DELTA_CONSTANT[D]
    harmonic = 1+F(7*TAIL_EXPONENT,20)
    base = F(4,D+1)+data.B
    common = data.area*r+data.boundary*harmonic+2
    tail_core = (base+F(2*cm,3*w)*(data.area*F(r,2*D*2**8)+common/F(2**32))
                 +F(cd,w)*(F(r,2*D*2**8)+(2-w-2*nu+r)/2**32))
    ca = max(thin_constant(e) for e in EXCEPTIONS[D])
    tail_thin = (base+F(2*cm,3*w)*(data.area*F(r,2*D*2**16)+data.area*r/F(2**48)
                                    +(data.boundary*harmonic+2)/F(2**32))
                 +ca/w*(F(r,2*D*2**16)+(2-w-2*nu+r)/2**48))
    assert tail_core < F(19,20) and tail_thin < F(3,5)
    print(f'D={D}: PASS core and exceptional-quotient families from 2^{CORE_EXPONENT[D]}; {count_core} state inequalities, {count_thin} family inequalities',flush=True)
    return {
        'D':D,'analytic_onset':1 << CORE_EXPONENT[D],
        'large_finite_coverage_performed':False,
        'finite_classification_obligation':[2,1 << CORE_EXPONENT[D]],
        'mass':rational(F(D-1,12)),'principal_weight':rational(w),'ramified_point_weight':rational(nu),
        'odd_inert_residue_classes':data.residues,'branched_inert_primes':data.chosen,
        'B_upper':rational(data.B),'odd_lattice_area_upper':rational(data.area),
        'odd_lattice_boundary_upper':rational(data.boundary),
        'M_global_constant':cm,'M_sixth_product':rational(split_sixth),
        'delta_global_constant':cd,'delta_sixth_product':rational(delta_sixth),
        'M_frontier':{'enumerated':enum_count,'frontier':frontier,'prime_prefix':chosen,'prefix_product':prefix},
        'delta_frontier':df_info,'intervals':intervals,
        'core_state_inequalities':count_core,'thin_family_inequalities':count_thin,
        'inequality_trace_sha256':trace.hexdigest(),
        'infinite_tail_from':1 << TAIL_EXPONENT,'infinite_core_cost':rational(tail_core),
        'infinite_thin_cost':rational(tail_thin),'thin_global_constant':rational(ca),
        'small_exact_controls':small,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    all_primes = primes(R)
    rows = [certify(D,all_primes) for D in DS]
    result = {
        'status':'PASS',
        'scope':'Exact analytic bounds and finite-reduction theorem; no full classification or large finite coverage.',
        'checker_sha256':sha(Path(__file__).read_bytes()),
        'arithmetic':'Exact integer/rational arithmetic; no floating-point acceptance.',
        'prime_sieve_limit':R,'reciprocal_rounding_denominator':DEN,
        'tail_exponent':TAIL_EXPONENT,'results':rows,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
