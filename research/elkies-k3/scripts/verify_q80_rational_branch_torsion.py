#!/usr/bin/env python3
"""Replay branch-field exclusions using polynomial gcd, not point enumeration.

The atlas' exact section construction and complete lattice census remain
inherited inputs. This checker independently verifies the new finite-field
obstructions, their full attachment, integrality and simple-root conditions.
"""
import argparse
from fractions import Fraction
import gzip
from hashlib import sha256
import json
from math import gcd, isqrt, lcm
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT/'artifacts/generated-results/elkies-k3-q80-rational-branch-torsion-v1'
ATLAS = 'artifacts/generated-results/elkies-k3-r17-smooth-character-replay-inputs-v1.json.gz'
PARENTS = 'artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/input.json'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def trim(a):
    while a and not a[-1]:
        a.pop()
    return a


def rem(a,b,p):
    a = trim(a[:])
    while len(a) >= len(b):
        c, shift = a[-1]*pow(b[-1],-1,p) % p, len(a)-len(b)
        for j,d in enumerate(b):
            a[j+shift] = (a[j+shift]-c*d) % p
        trim(a)
    return a


def mul(a,b,f,p):
    c = [0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            c[i+j] = (c[i+j]+x*y) % p
    return rem(c,f,p)


def root_free_cubic(a,b,p):
    f, acc, x, n = [b,a,0,1], [1], [0,1], p
    while n:
        if n & 1:
            acc = mul(acc,x,f,p)
        x, n = mul(x,x,f,p), n//2
    acc += [0]*max(0,2-len(acc))
    acc[1] = (acc[1]-1) % p
    acc = trim(acc)
    while acc:
        f, acc = acc, rem(f,acc,p)
    return len(f) == 1


def ev(coeff,t,p):
    return sum(c*pow(t,i,p) for i,c in enumerate(coeff)) % p


def check_witness(q,A,B,p,t):
    require(0 <= t < p, 'finite residue coordinate')
    require(q[2] % p and (q[1]**2-4*q[0]*q[2]) % p, 'branch leading unit and etale reduction')
    require(ev(q,t,p) == 0 and (q[1]+2*q[2]*t) % p, 'simple branch root')
    a,b = ev(A,t,p),ev(B,t,p)
    require((4*a**3+27*b*b) % p, 'singular parent fibre')
    require(root_free_cubic(a,b,p), 'cubic has a residue root')


def verify(out):
    start = time.monotonic()
    packet, result = read(out/'input.json'),read(out/'result.json')
    require(packet['checker_sha256'] == digest(Path(__file__)), 'checker binding')
    producer = ROOT/'elkies-k3/scripts/certify_q80_rational_branch_torsion.py'
    require(packet['producer_sha256'] == digest(producer), 'producer binding')
    require({ATLAS,PARENTS} <= set(packet['bindings']), 'missing source bindings')
    for path,expected in packet['bindings'].items():
        require(digest(ROOT/path) == expected, 'changed source '+path)
    require(result['input_sha256'] == digest(out/'input.json'), 'result binding')
    require(packet['chart'] == '11952' and packet['count'] == result['count'] == 39147, 'wrong atlas scope')
    require(result['positive_target_complete'] is False and result['whole_family_genus_bound'] == 'UNKNOWN', 'unsupported scope')
    parent = next(r for r in read(ROOT/PARENTS)['parents'] if r['name'] == 'alternate-q80')
    require(parent == packet['parent'], 'wrong parent equation')
    require(packet['bindings'][parent['source']] == parent['source_sha256'], 'parent source binding')
    atlas = read(ROOT/ATLAS)['atlases'][0]
    require(atlas['source']['chart'] == '11952' and len(atlas['rows']) == len(atlas['priority']) == 39147, 'atlas completeness')
    priority = dict(atlas['priority'])
    require(len(priority) == 39147, 'duplicate priority')
    branches, masks, labels = [],set(),set()
    for mask,label,word,num,den in atlas['rows']:
        require(mask not in masks and label not in labels and priority[mask] == word, 'branch attachment')
        require(len(num) == 3 and len(den) == 1 and Fraction(den[0]), 'branch chart')
        values = [Fraction(c)/Fraction(den[0]) for c in num]
        scale = lcm(*(c.denominator for c in values))
        q = [int(c*scale) for c in values]
        content = gcd(*q)
        q = [c//content for c in q]
        d = q[1]**2-4*q[0]*q[2]
        require(q[2] and d != 0 and (d < 0 or isqrt(d)**2 != d), 'irreducible quadratic branch field')
        branches.append(q)
        masks.add(mask); labels.add(label)
    coeff = [[Fraction(c) for c in parent[key]] for key in ('A','B')]
    remaining, count, summaries = set(range(39147)),0,[]
    require([s['prime'] for s in result['stages']] == packet['primes'][:len(result['stages'])], 'prime prefix')
    for meta in result['stages']:
        p = meta['prime']
        require(type(p) is int and p >= 101 and all(p % d for d in range(2,isqrt(p)+1)), 'prime field required')
        integral = all(c.denominator % p for line in coeff for c in line)
        if meta['status'] == 'DEFER_NONINTEGRAL_PARENT':
            require(not integral, 'unjustified deferral')
            continue
        require(integral, 'nonintegral parent')
        name = 'prime-%d.json' % p
        require(meta['path'] == name and meta['sha256'] == digest(out/name), 'stage binding')
        stage = read(out/name)
        require(stage['input_sha256'] == digest(out/'input.json') and stage['prime'] == p, 'stage attachment')
        A,B = [[c.numerator*pow(c.denominator,-1,p) % p for c in line] for line in coeff]
        fibres = [t for t in range(p) if (4*ev(A,t,p)**3+27*ev(B,t,p)**2) % p and root_free_cubic(ev(A,t,p),ev(B,t,p),p)]
        require(fibres == stage['root_free_smooth_fibres'], 'independent cubic gcd roster')
        indices = [i for i,t in stage['witnesses']]
        require(indices == sorted(set(indices)) and set(indices) <= remaining, 'duplicate or misplaced witness')
        for i,t in stage['witnesses']:
            check_witness(branches[i],A,B,p,t)
        remaining.difference_update(indices)
        require(stage['remaining_indices'] == sorted(remaining) and meta['remaining'] == len(remaining) and meta['excluded'] == len(indices), 'stage coverage')
        count += len(indices)
        summaries.append({'prime':p,'excluded':len(indices),'remaining':len(remaining)})
    require(result['remaining_indices'] == sorted(remaining) and result['excluded'] == count, 'final coverage')
    require(result['status'] == ('PASS' if not remaining else 'INCOMPLETE'), 'false completion')
    return {'status':result['status'], 'branch_fields_excluded':count,
            'remaining_indices':sorted(remaining), 'stages':summaries,
            'input_sha256':digest(out/'input.json'), 'result_sha256':digest(out/'result.json'),
            'new_branch_field_replay_independent':True,
            'inherited_section_and_lattice_replay':False, 'positive_target_complete':False,
            'elapsed_seconds':round(time.monotonic()-start,6)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=DEFAULT)
    parser.add_argument('--record',type=Path)
    args = parser.parse_args()
    result = verify(args.input)
    if args.record:
        with args.record.open('x') as stream:
            json.dump(result,stream,indent=2,sort_keys=True)
            stream.write('\n')
    print(json.dumps(result,sort_keys=True),flush=True)
