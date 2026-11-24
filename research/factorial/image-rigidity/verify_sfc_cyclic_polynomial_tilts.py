#!/usr/bin/env python3
"""Exact finite controls for symmetric Gamma-tilted cyclic polynomial blocks."""
if not __debug__:
    raise RuntimeError("This checker requires assertions; do not run with -O.")

from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import permutations, product
from math import factorial, comb, gcd, lcm, prod
from pathlib import Path
import json
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
EXPECTED_NOTE_SHA256 = "ba7cc5fe27558993b93eb87d4a08e2cb29adb1c30d86ab8868b20cbfc5cfde0d"


def rising(x, n):
    return prod((F(x)+j for j in range(n)), start=F(1))


def exponent_sum(a, b):
    return tuple(x+y for x, y in zip(a, b))


def scalar_poly_mul(a, b):
    out = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = exponent_sum(ea, eb)
            out[e] = out.get(e, F(0)) + ca*cb
    return {e:c for e,c in out.items() if c}


def gamma_value(poly, beta):
    return sum(c*prod((rising(b, e) for b,e in zip(beta, alpha)), start=F(1))
               for alpha,c in poly.items())


def gamma_cumulants(poly, beta, bound):
    d = len(beta)
    assert poly and all(c > 0 for c in poly.values())
    assert all(len(e) == d and any(e) for e in poly)
    assert all(b > 0 for b in beta)
    moment, cumulant = [F(1)], [F(0)]
    power = {(0,)*d:F(1)}
    selected, coefficient = next(iter(poly.items()))
    active = next(i for i,e in enumerate(selected) if e)
    for n in range(1, bound+1):
        power = scalar_poly_mul(power, poly)
        moment.append(gamma_value(power, beta))
        cumulant.append(moment[n] - sum(comb(n-1,k-1)*cumulant[k]*moment[n-k]
                                       for k in range(1,n)))
        lower = F(factorial(n-1))*coefficient**n
        for i,(b,e) in enumerate(zip(beta,selected)):
            lower *= b**(n*(e-1)+1 if i == active else n*e)
        assert cumulant[n] >= lower > 0
    return moment, cumulant


@lru_cache(None)
def permutation_cycles(n):
    records = []
    for permutation in permutations(range(n)):
        seen, cycles = set(), []
        for i in range(n):
            if i in seen:
                continue
            cycle, j = [], i
            while j not in seen:
                seen.add(j)
                cycle.append(j)
                j = permutation[j]
            cycles.append(tuple(cycle))
        records.append(tuple(cycles))
    return tuple(records)


def connected_permutation_sum(poly, beta, n):
    terms = tuple(poly.items())
    total, tuple_count, connected_count = F(0), 0, 0
    for chosen in product(terms, repeat=n):
        coefficient = prod((c for _,c in chosen), start=F(1))
        owners = [[v for v,(exponents,_) in enumerate(chosen)
                   for _ in range(exponents[i])] for i in range(len(beta))]
        options = [permutation_cycles(len(labels)) for labels in owners]
        for selected in product(*options):
            tuple_count += 1
            parent = list(range(n))
            def root(i):
                while parent[i] != i:
                    parent[i] = parent[parent[i]]
                    i = parent[i]
                return i
            weight = coefficient
            for b,labels,cycles in zip(beta,owners,selected):
                weight *= b**len(cycles)
                for cycle in cycles:
                    first = labels[cycle[0]]
                    for leg in cycle[1:]:
                        parent[root(labels[leg])] = root(first)
            if len({root(i) for i in range(n)}) == 1:
                connected_count += 1
                total += weight
    return total, tuple_count, connected_count


class Ring:
    """Explicit exact cyclotomic quotient rings used by the finite cases."""
    moduli = {
        1:(-1,1), 2:(1,1), 3:(1,1,1), 4:(1,0,1),
        5:(1,1,1,1,1), 6:(1,-1,1), 10:(1,-1,1,-1,1),
    }
    def __init__(self, order):
        self.order, self.mod = order, self.moduli[order]
        self.degree = len(self.mod)-1
        self.zero = (F(0),)*self.degree
        self.one = (F(1),)+(F(0),)*(self.degree-1)
        self.powers = [self.reduce([0]*k+[1]) for k in range(order+1)]
        assert self.powers[-1] == self.one
        assert all(self.powers[k] != self.one for k in range(1,order))
    def reduce(self, values):
        values = list(map(F,values))
        while len(values)>self.degree:
            k, c = len(values)-1-self.degree, values[-1]
            for j,b in enumerate(self.mod):
                values[k+j] -= c*b
            values.pop()
        return tuple(values+[F(0)]*(self.degree-len(values)))
    def add(self,a,b):
        return tuple(x+y for x,y in zip(a,b))
    def scale(self,a,b):
        return tuple(x*b for x in a)
    def mul(self,a,b):
        out = [F(0)]*(2*self.degree-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):
                out[i+j] += x*y
        return self.reduce(out)


def ring_poly_mul(a,b,ring):
    out = {}
    for ea,ca in a.items():
        for eb,cb in b.items():
            e = exponent_sum(ea,eb)
            out[e] = ring.add(out.get(e,ring.zero),ring.mul(ca,cb))
    return {e:c for e,c in out.items() if c != ring.zero}


def ordinary_value(poly, shift, ring):
    total = ring.zero
    for e,c in poly.items():
        total = ring.add(total,ring.scale(c,prod(factorial(x+y) for x,y in zip(e,shift))))
    return total


def positive_exposing_weights(poly):
    d = len(next(iter(poly)))
    T = 2
    while True:
        w = [F(T**i) for i in range(d)]
        scores = {e:sum(x*y for x,y in zip(w,e)) for e in poly}
        maximum = max(scores.values())
        selected = [e for e,v in scores.items() if v == maximum]
        if len(selected) == 1:
            assert maximum > 0
            return [x/maximum for x in w],selected[0]
        T += 1


def reachable(generators, bound):
    out = [True]+[False]*bound
    for n in range(1,bound+1):
        out[n] = any(n >= q and out[n-q] for q in generators)
    return out


def check_orbit_case(specs, bound, label):
    order = lcm(*(q for _,q,_,_ in specs))
    ring = Ring(order)
    dimension = sum(q*len(a) for _,q,a,_ in specs)
    G, tilt, weights = {}, [], []
    cursor = 0
    for poly,q,a,rho in specs:
        assert rho > 0 and len(a) == len(next(iter(poly)))
        w,selected = positive_exposing_weights(poly)
        for j in range(q):
            phase = ring.powers[(j*order//q)%order]
            for alpha,c in poly.items():
                exponent = [0]*dimension
                exponent[cursor:cursor+len(a)] = alpha
                exponent = tuple(exponent)
                G[exponent] = ring.add(G.get(exponent,ring.zero),ring.scale(phase,rho*c))
            tilt.extend(a)
            weights.extend(w)
            cursor += len(a)
    assert all(c != ring.zero for c in G.values())
    assert len(G) == sum(q*len(poly) for poly,q,_,_ in specs)
    L = sum(q for _,q,_,_ in specs)
    assert sum(sum(x*y for x,y in zip(weights,e)) == 1 for e in G) == L
    members = reachable([q for _,q,_,_ in specs],bound)
    power = {(0,)*dimension:ring.one}
    records, power_supports = [], []
    for n in range(bound+1):
        if n:
            power = ring_poly_mul(power,G,ring)
        raw = ordinary_value(power,[n*a for a in tilt],ring)
        logseries = [F(0)]*(n+1)
        for poly,q,a,rho in specs:
            beta = [F(n*x+1) for x in a]
            _,cumulants = gamma_cumulants(poly,beta,n)
            for k in range(q,n+1,q):
                logseries[k] += q*cumulants[k]*rho**k/factorial(k)
        exponential = [F(1)]+[F(0)]*n
        for k in range(1,n+1):
            exponential[k] = sum(F(j,k)*logseries[j]*exponential[k-j]
                                 for j in range(1,k+1))
        scale = factorial(n)*prod(factorial(n*a) for a in tilt)
        assert raw == ring.scale(ring.one,scale*exponential[n])
        assert (exponential[n] > 0) == members[n]
        records.append(dict(n=n,normalized_moment=str(exponential[n]),
                            ordinary_moment=str(scale*exponential[n]),nonzero=members[n]))
        if 1 <= n <= 3:
            face_count = sum(sum(x*y for x,y in zip(weights,e)) == n for e in power)
            assert face_count == comb(L+n-1,n)
            assert len(power) >= face_count >= min(q for _,q,_,_ in specs)//gcd(
                min(q for _,q,_,_ in specs),n)
            power_supports.append(dict(power=n,actual_support=len(power),
                                       exposed_face_support=face_count,
                                       bound=comb(L+n-1,n)))
    return dict(label=label,dimension=dimension,cycles=[q for _,q,_,_ in specs],
                tilt_vectors=[list(a) for _,_,a,_ in specs],
                templates=[{str(e):str(c) for e,c in poly.items()} for poly,_,_,_ in specs],
                scales=[str(rho) for _,_,_,rho in specs],
                checked_orders=[0,bound],moments=records,power_supports=power_supports)


def verify_rejection_controls():
    optimized = subprocess.run([sys.executable, "-O", str(Path(__file__).resolve())],
                               capture_output=True, text=True)
    assert optimized.returncode != 0 and "do not run with -O" in optimized.stderr
    with tempfile.TemporaryDirectory(prefix="sfc_cyclic_tilt_guard_") as directory:
        temporary = Path(directory)
        script = temporary / Path(__file__).name
        script.write_bytes(Path(__file__).read_bytes())
        note = temporary / "SFC_CYCLIC_POLYNOMIAL_TILTS.md"
        note.write_bytes((HERE/note.name).read_bytes()+b"\nAltered-note rejection control.\n")
        altered = subprocess.run([sys.executable, str(script)], capture_output=True, text=True)
        assert altered.returncode != 0 and "Source SHA256 mismatch" in altered.stderr
    return dict(optimized_mode_rejected=True, altered_note_rejected=True)


def main():
    note = HERE / "SFC_CYCLIC_POLYNOMIAL_TILTS.md"
    if sha256(note.read_bytes()).hexdigest() != EXPECTED_NOTE_SHA256:
        raise RuntimeError("Source SHA256 mismatch: "+note.name)
    rejection_controls = verify_rejection_controls()
    diagram_records = []
    diagram_cases = [
        ({(1,):F(1),(2,):F(1)},[F(3,2)],3),
        ({(1,0):F(1),(0,1):F(1),(1,1):F(1,2)},[F(1,3),F(4)],3),
        ({(2,1):F(2)},[F(2),F(1,2)],2),
    ]
    for poly,beta,bound in diagram_cases:
        _,cumulants = gamma_cumulants(poly,beta,bound)
        for n in range(1,bound+1):
            value,tuples,connected = connected_permutation_sum(poly,beta,n)
            assert value == cumulants[n]
            diagram_records.append(dict(template={str(e):str(c) for e,c in poly.items()},
                                        beta=[str(b) for b in beta],n=n,
                                        cumulant=str(value),permutation_tuples=tuples,
                                        connected_tuples=connected))
    # Strict log positivity must not be justified by Kaluza in these shapes.
    moments,cumulants = gamma_cumulants({(1,):F(1)},[F(3)],8)
    a = [moments[n]/factorial(n) for n in range(9)]
    assert a[1]/a[0] > a[2]/a[1]
    assert all(cumulants[n] == 3*factorial(n-1) for n in range(1,9))
    cases = [
        ([({(1,):F(1),(2,):F(1)},2,(1,),F(1))],7,"two-cycle polynomial with common tilt"),
        ([({(1,0):F(1),(0,1):F(2),(1,1):F(1)},3,(1,2),F(1))],4,"three-cycle bivariate positive template"),
        ([({(2,0):F(1),(1,1):F(2),(0,2):F(1)},4,(2,0),F(1))],4,"four-cycle homogeneous template"),
        ([({(1,):F(1),(2,):F(1)},2,(1,),F(1,2)),
          ({(1,):F(1)},3,(2,),F(2))],6,"mixed semigroup generated by two and three"),
        ([({(1,0):F(1),(0,1):F(1),(1,1):F(1)},1,(2,1),F(1)),
          ({(1,):F(1),(3,):F(1)},2,(1,),F(1,3))],5,"cycle one gives no moment zeros"),
        ([({(2,0):F(1),(1,1):F(2),(0,2):F(1)},3,(1,2),F(1))],4,"before symmetric within-copy identification"),
        ([({(2,):F(4)},3,(3,),F(1))],4,"after symmetric within-copy identification"),
    ]
    orbit_records = []
    for specs,bound,label in cases:
        orbit_records.append(check_orbit_case(specs,bound,label))
        print("PASS",label,flush=True)
    assert [r['nonzero'] for r in orbit_records[-2]['moments']] == [
        r['nonzero'] for r in orbit_records[-1]['moments']]
    assert sum(c*prod(factorial(v) for v in e) for e,c in {
        (2,0):F(1),(1,1):F(-1)}.items()) == 1
    receipt = dict(status="PASS",evidence_type="exact finite controls, not a proof of universal quantifiers",
                   own_note_sha256=sha256((HERE/'SFC_CYCLIC_POLYNOMIAL_TILTS.md').read_bytes()).hexdigest(),
                   own_checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                   rejection_controls=rejection_controls,
                   expected_note_sha256=EXPECTED_NOTE_SHA256,
                   connected_permutation_controls=diagram_records,
                   orbit_tilt_controls=orbit_records,
                   gamma_three_kaluza_shortcut_refuted=True,
                   asymmetric_common_monomial_zero_pattern_refuted=True,
                   arbitrary_cross_copy_merging_proved=False,unrestricted_SFC_proved=False)
    output = HERE/'sfc_cyclic_polynomial_tilts_verification.json'
    output.write_text(json.dumps(receipt,indent=2)+'\n')
    print("PASS:",len(diagram_records),"connected-diagram controls and",len(orbit_records),"orbit cases")
    print(output)


if __name__ == '__main__':
    main()
