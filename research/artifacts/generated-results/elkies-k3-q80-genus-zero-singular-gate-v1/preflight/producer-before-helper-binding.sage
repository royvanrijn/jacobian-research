#!/usr/bin/env sage-python
"""Reuse the 49 full Q80 norm12 pencils to exclude rational singular members.

No new trace compilation, lattice enumeration or pair search is performed.
The arithmetic genus-one norm8 exclusion is an inherited theorem, with its
existing independent-replay gap unchanged. The genus-zero descent is written
in the canonical proof note, not certified by the finite-field calculation.
"""
import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import resource
import time

from sage.all import GF, QQ, PolynomialRing
from sage.version import version

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / 'artifacts/generated-results/elkies-k3-q80-genus-zero-singular-gate-v1'
OLD = 'artifacts/generated-results/elkies-k3-q80-norm12-moving-pencils-v1'
NORM8 = 'artifacts/generated-results/elkies-k3-r17-norm12-11952-singular-bisection-search-complete-v1.json'
PARITY = 'artifacts/generated-results/elkies-k3-r17-norm12-11952-product-tate-parity-v1.json'
CHECKER = 'elkies-k3/scripts/verify_q80_genus_zero_singular_gate.py'
INPUTS = [OLD + '/' + name for name in (
    'input.json.gz', 'frames-521.json.gz', 'frames-523.json.gz',
    'prime-521.json.gz', 'prime-523.json.gz', 'independent-replay.json')]
INPUTS += [NORM8, PARITY,
    'artifacts/generated-results/elkies-k3-common-quartic-branch-strata-v1/input.json',
    'artifacts/generated-results/elkies-k3-common-quartic-branch-strata-v1/result.json',
    'elkies-k3/scripts/verify_q80_norm12_moving_pencils.py',
    'elkies-k3/scripts/verify_q80_complete_genus_one_pairs.py',
    'elkies-k3/scripts/verify_r17_correlated_genus_one_pencils.py']


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def write_new(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(obj, stream, indent=2, sort_keys=True)
        stream.write('\n')


def freeze(out):
    write_new(out / 'input.json', {
        'schema': 'q80-genus-zero-singular-gate-input-v1',
        'bindings': {p: digest(ROOT / p) for p in INPUTS},
        'producer_sha256': digest(Path(__file__)),
        'checker_sha256': digest(ROOT / CHECKER),
        'primes': [521, 523], 'pencil_indices': list(range(49)),
        'cpu_seconds': 30, 'address_space_bytes': 4*1024**3,
        'sage_version': version,
        'selection': 'Reuse the two completed moving-pencil prime stages. A read-only discriminant preview closed all49; this packet records that discovered proof, not a prospective search.',
        'acceptance': 'For each pencil, the full homogeneous degree24 branch discriminant must have no projective root at one certified integral prime. No odd-part substitution.',
        'scope': 'New: full49 norm12 singular-member exclusion and12 genus-zero factor identities. Inherited: norm8 rational-normalization exclusion, parity census and Q80 unramified descent. The positive two-gain target and other factor strata remain open.',
    })
    print('Frozen 49 retained pencils at primes521,523', flush=True)


def genus_zero_strata():
    R = PolynomialRing(QQ, names='g,d,u,v,a,b,z,c,l')
    g,d,u,v,a,b,z,c,l = R.gens()
    h, q = c*g*u*v, c*l*g*d
    s1,s2 = (u*a+v*b)/2, (v*b-u*a)/2
    A = l*d*a*b - 3*z*z - 3*z*h - h*h
    B = q*s2*s2 + 2*z**3 + 3*z*z*h + z*h*h - l*z*d*a*b
    f = lambda x: x**3 + A*x + B
    assert f(z+h) == q*s1*s1 and f(z) == q*s2*s2
    rows = []
    for k in range(3):
        for j in range(5-k):
            weights = [k,2-k,j,4-k-j,5-j,1+k+j,4,0,0]
            degree = lambda p: max(sum(e*w for e,w in zip(ex,weights)) for ex in p.exponents())
            assert degree(q) == 2 and degree(h) == 4
            assert degree(s1) == degree(s2) == 5
            assert degree(A) <= 8 and degree(B) <= 12
            assert sum(weights[:4]) + (6-j) + (2+k+j) + 5 + 2 == 21
            cross = 4-k-2*j
            rows.append({'k':k, 'j':j, 'degree_bounds':weights,
                         'intersection':k+2*j, 'cross_height':cross,
                         'sum_height':16+2*cross, 'difference_height':16-2*cross})
    assert len(rows) == 12
    assert [r for r in rows if r['k'] == 2 and r['cross_height'] == 0][0]['j'] == 1
    assert (16-16)//2 == 0 and 1+(16-16)//4 == 1
    return rows


def run(out):
    packet = read(out / 'input.json')
    assert packet['producer_sha256'] == digest(Path(__file__))
    assert packet['checker_sha256'] == digest(ROOT / CHECKER)
    assert set(packet['bindings']) == set(INPUTS)
    for path, expected in packet['bindings'].items():
        assert digest(ROOT / path) == expected, path
    resource.setrlimit(resource.RLIMIT_CPU, (30,35))
    resource.setrlimit(resource.RLIMIT_AS, (4*1024**3,4*1024**3))
    start = time.process_time()
    remaining = set(packet['pencil_indices'])
    witnesses, stages = {}, []
    for p in packet['primes']:
        frame_path = ROOT / OLD / ('frames-%d.json.gz' % p)
        frame_packet = read(frame_path)
        assert frame_packet['input_sha256'] == digest(ROOT / OLD / 'input.json.gz')
        assert frame_packet['prime'] == p
        frames = {r['index']:r for r in frame_packet['records']}
        assert sorted(frames) == packet['pencil_indices']
        R = PolynomialRing(GF(p), 'u')
        rows, excluded = [], []
        selected = sorted(remaining)
        for i in selected:
            e,d,c,b,a = [R(row) for row in frames[i]['branch_matrix']]
            I = 12*a*e - 3*b*d + c*c
            J = 72*a*c*e + 9*b*c*d - 27*a*d*d - 27*b*b*e - 2*c**3
            disc = (4*I**3-J**2) * GF(p)(27)**-1
            assert disc and disc.degree() <= 24
            roots = sorted([[int(r),int(m)] for r,m in disc.roots()])
            infinity_order = 24-int(disc.degree())
            rows.append({'index':i, 'discriminant_coefficients':[int(disc[n]) for n in range(25)],
                         'finite_roots_with_multiplicity':roots,
                         'infinity_multiplicity':infinity_order})
            if not roots and infinity_order == 0:
                excluded.append(i)
                witnesses[str(i)] = p
        remaining.difference_update(excluded)
        name = 'prime-%d.json' % p
        write_new(out / name, {'prime':p, 'input_sha256':digest(out/'input.json'),
            'frames_sha256':digest(frame_path), 'selected_indices':selected,
            'records':rows, 'excluded_indices':excluded, 'remaining_indices':sorted(remaining)})
        stages.append({'prime':p, 'selected':len(selected), 'excluded':len(excluded),
                       'remaining':len(remaining), 'path':name, 'sha256':digest(out/name)})
        print(json.dumps(stages[-1]), flush=True)
    assert not remaining, 'unresolved pencils are not exclusions'
    result = {'schema':'q80-genus-zero-singular-gate-result-v1', 'status':'PASS',
        'input_sha256':digest(out/'input.json'), 'stages':stages, 'witnesses':witnesses,
        'genus_zero_strata':genus_zero_strata(), 'norm12_pencils_excluded':49,
        'inherited_norm8_pencils':63917,
        'q80_genus_zero_k2_j1':'EXCLUDED by the written descent and arithmetic-genus-one classification',
        'whole_family_genus_bound':'UNKNOWN', 'positive_target_complete':False,
        'cpu_seconds':round(time.process_time()-start,6)}
    write_new(out/'result.json', result)
    print(json.dumps({k:result[k] for k in ('status','norm12_pencils_excluded','cpu_seconds')}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['freeze','run'])
    parser.add_argument('--output', type=Path, default=DEFAULT)
    args = parser.parse_args()
    (freeze if args.mode == 'freeze' else run)(args.output)
