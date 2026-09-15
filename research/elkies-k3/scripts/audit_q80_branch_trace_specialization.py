#!/usr/bin/env python3
"""Exact generic-parity specialization gate at one prospective rational branch."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import resource
import time

ROOT = Path(__file__).resolve().parents[2]
SOURCE = 'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'


def prime(n):
    return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))


def evaluate(coeffs, t):
    value = Fraction(0)
    for c in reversed(coeffs):
        value = value*t + Fraction(c)
    return value


def point_at(section, t):
    return [evaluate(section[k]['numerator_coefficients_low_to_high'], t) /
            evaluate(section[k]['denominator_coefficients_low_to_high'], t) for k in ['X','Y']]


def reduction(value, p):
    return value.numerator * pow(value.denominator, -1, p) % p


def audit():
    resource.setrlimit(resource.RLIMIT_CPU, (15, 15))
    resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))
    start = time.monotonic()
    source = json.loads((ROOT / SOURCE).read_text())
    model = source['weierstrass_model']
    t = 0  # Fixed before arithmetic; no exceptional point or rank input.
    A = evaluate(model['A_coefficients_low_to_high'], t)
    B = evaluate(model['B_coefficients_low_to_high'], t)
    sections = source['sections']['records']
    assert [s['basis_index'] for s in sections] == list(range(17))
    points = [point_at(s, t) for s in sections]
    assert 4*A**3+27*B**2 != 0
    assert all(y*y == x*x*x+A*x+B for x,y in points)
    rows, selected, skipped = [], [], []
    pivots = {}
    pool = [p for p in range(5, 998) if prime(p)]
    for p in pool:
        if any(z.denominator % p == 0 for z in [A,B,*[c for point in points for c in point]]):
            skipped.append({'prime':p,'reason':'nonintegral literal coefficients or specialized point'})
            continue
        a,b = reduction(A,p),reduction(B,p)
        if (4*a**3+27*b*b)%p == 0:
            skipped.append({'prime':p,'reason':'singular literal reduction'})
            continue
        reduced = [[reduction(x,p),reduction(y,p)] for x,y in points]
        for e in range(p):
            if (e**3+a*e+b)%p:
                continue
            values = [(x-e)%p if x != e else (3*e*e+a)%p for x,y in reduced]
            assert all(values)
            code = sum((pow(v,(p-1)//2,p)==p-1)<<j for j,v in enumerate(values))
            row = {'prime':p,'root':e,'A':a,'B':b,'points':reduced,'values':values,'code':code}
            rows.append(row)
            remainder = code
            while remainder:
                lead = remainder.bit_length()-1
                if lead in pivots:
                    remainder ^= pivots[lead]
                else:
                    pivots[lead] = remainder
                    selected.append({'row':len(rows)-1,'pivot':lead,'reduced_code':remainder})
                    break
        if len(pivots) == 17:
            break
    return {'schema':'q80-branch-trace-specialization-v1','status':'COMPLETE_BOUNDED_GATE',
            'source':SOURCE,'source_sha256':sha256((ROOT/SOURCE).read_bytes()).hexdigest(),
            'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'branch_parameter':t,'prime_pool':pool,'limits':{'cpu_seconds':15,'memory_bytes':1024**3},
            'selection':'t=0 and ascending primes5..997, stop only at exact row rank17',
            'A':str(A),'B':str(B),'points':[[str(x),str(y)] for x,y in points],
            'rows':rows,'selected':selected,'skipped':skipped,'rank':len(pivots),
            'all_generic_trace_parities_excluded_at_branch':len(pivots)==17,
            'new_mw17_gain_constructed':False,'elapsed_seconds':time.monotonic()-start}


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',type=Path,required=True)
    args=parser.parse_args()
    assert not args.record.exists(),'Refusing to overwrite evidence'
    result=audit()
    args.record.parent.mkdir(parents=True,exist_ok=True)
    with args.record.open('x') as f:
        json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps({'rank':result['rank'],'rows':len(result['rows']),
                      'witness_primes':sorted({result['rows'][r['row']]['prime'] for r in result['selected']}),
                      'elapsed_seconds':result['elapsed_seconds']}))
