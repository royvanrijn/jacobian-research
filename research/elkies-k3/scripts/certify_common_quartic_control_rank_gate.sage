#!/usr/bin/env sage-python
"""Bound the rational Picard ranks of the two existing control K3 surfaces.

One good-prime count per control suffices to exclude every MW16/MW17
fibration on that surface. No section, lattice or specialization search runs.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import resource
import time

from sage.all import EllipticCurve, GF, PolynomialRing, ZZ
from sage.version import version

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT/'artifacts/generated-results/elkies-k3-common-quartic-control-rank-gate-v1'
SOURCE = 'artifacts/generated-results/elkies-k3-common-quartic-singularities-v1/input.json'
CHECKER = 'elkies-k3/scripts/verify_common_quartic_control_rank_gate.py'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def write_new(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as stream:
        json.dump(value,stream,indent=2,sort_keys=True)
        stream.write('\n')


def freeze(out):
    source = read(ROOT/SOURCE)
    write_new(out/'input.json', {
        'schema':'common-quartic-control-rank-input-v1',
        'source':SOURCE, 'source_sha256':digest(ROOT/SOURCE),
        'controls':source['controls'], 'primes':[5,7],
        'producer_sha256':digest(Path(__file__)), 'checker_sha256':digest(ROOT/CHECKER),
        'preview_sha256':digest(out/'preview.json'),
        'sage_version':version, 'cpu_seconds':10, 'address_space_bytes':4*1024**3,
        'selection':'First good primes for the two existing controls, discovered in the retained fixed-panel preview. This records that proof; it is not prospective model selection.',
        'bound':'rho_Q <= floor((22*p + #X(Fp)-1-p^2)/(2*p)); every Q-Jacobian fibration has arithmetic MW rank <= rho_Q-2.',
        'boundary':'Rational Picard-rank upper bounds, not exact or geometric ranks. The same bounds hold for every integral coefficient model with the identical good reduction.',
    })
    print('Frozen two controls at5 and7',flush=True)


def run(out):
    packet = read(out/'input.json')
    assert digest(Path(__file__)) == packet['producer_sha256']
    assert digest(ROOT/CHECKER) == packet['checker_sha256']
    assert packet['source'] == SOURCE and digest(ROOT/SOURCE) == packet['source_sha256']
    assert packet['controls'] == read(ROOT/SOURCE)['controls']
    assert digest(out/'preview.json') == packet['preview_sha256']
    start = time.process_time()
    R = PolynomialRing(ZZ,'t')
    records = []
    for control,p in zip(packet['controls'],packet['primes']):
        D,s = R(control['D']),R(control['s'])
        A,B = D*(2*s+1)-1,D*s*s
        assert A.degree()==8 and B.degree()==12
        delta = 4*A**3+27*B**2
        S = PolynomialRing(GF(p),'t')
        d = S(delta)
        assert d.degree()==24 and d.gcd(d.derivative()).degree()==0
        counts = []
        for t in list(range(p))+[None]:
            a = GF(p)(A[8] if t is None else A(t))
            b = GF(p)(B[12] if t is None else B(t))
            if 4*a**3+27*b*b:
                n = int(EllipticCurve(GF(p),[a,b]).cardinality())
            else:
                # The total K3 is smooth; retain the nodal plane cubic itself.
                n = 1+sum(1 if z==0 else 2 if z.is_square() else 0
                          for z in (x**3+a*x+b for x in GF(p)))
            counts.append(n)
        N = sum(counts)
        trace = N-1-p*p
        rho = (22*p+trace)//(2*p)
        record = {'name':control['name'],'prime':p,
                  'A':[int(A[i]) for i in range(9)], 'B':[int(B[i]) for i in range(13)],
                  'discriminant_mod_p':[int(d[i]) for i in range(25)],
                  'fibre_counts_finite_then_infinity':counts,
                  'surface_count':N, 'h2_trace':trace,
                  'rational_picard_rank_upper':rho, 'any_rational_fibration_mw_rank_upper':rho-2,
                  'exact_inherited_rank':'UNKNOWN', 'geometric_picard_rank':'UNKNOWN'}
        name = control['name']+'-checkpoint.json'
        write_new(out/name,record)
        records.append({'path':name,'sha256':digest(out/name)})
        print(json.dumps({k:record[k] for k in ('name','prime','surface_count','rational_picard_rank_upper','any_rational_fibration_mw_rank_upper')}),flush=True)
    write_new(out/'result.json', {'schema':'common-quartic-control-rank-result-v1','status':'PASS',
              'input_sha256':digest(out/'input.json'),'records':records,
              'any_mw16_or_mw17_fibration_on_either_control':'EXCLUDED',
              'positive_mw17_target_complete':False,'cpu_seconds':round(time.process_time()-start,6)})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['freeze','run'])
    parser.add_argument('--output',type=Path,default=DEFAULT)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(10,15))
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    (freeze if args.mode=='freeze' else run)(args.output)
