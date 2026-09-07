#!/usr/bin/env python3
"""Export and fully replay the exact maximum parity classes in all six frames."""
import argparse
from pathlib import Path
from tempfile import TemporaryDirectory
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from replay_r17_norm12_minima import replay
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL=ROOT/'artifacts/local/elliptic-curves'
OUTPUT=ART/'r17_exact_maximum_parity_classes_v1.json'


def expected():
    rows=[];inputs={str(spec.ATLAS.relative_to(ROOT)):cert.hashed(spec.ATLAS)}
    for f in cert.read(spec.ATLAS)['families']:
        family=f['family'];folder='r17-norm12-exact-minima-v1' if family in ('11952','08f72') else 'r17-norm12-exact-minima-remaining4-v1'
        path=LOCAL/folder/family/'result.json';proof_path=ART/f'r17_exact_parity_radius_{family}_v1.json'
        data,proof=cert.read(path),cert.read(proof_path)
        if proof['status']!='PASS' or proof['exact_maximum_parity_coset_minimum']!=12 or proof['parity_upper_bounds_checked']!=131072:raise ArithmeticError('exact parity proof missing')
        for name,h in proof['sources'].items():
            if cert.hashed(ROOT/name)!=h:raise ArithmeticError('parity proof input changed')
        if data['gram']!=f['generic_height_gram'] or any(r['exact_minimum']!=12 for r in data['rows']):raise ArithmeticError('maximum-class extraction differs')
        classes=[{'mask':r['mask'],'representative':r['minimum_witness']} for r in data['rows']]
        if len(classes)!=proof['classes_with_exact_minimum12']:raise ArithmeticError('maximum-class count differs')
        rows.append({'family':family,'gram':data['gram'],'exact_maximum_parity_minimum':12,
            'classes':classes,'initial43_masks':[r['mask'] for r in data['rows'] if r['used_in_initial43']],
            'omitted_maximum_masks':[r['mask'] for r in data['rows'] if not r['used_in_initial43']],
            'query_result_path':str(path.relative_to(ROOT)),'proof_path':str(proof_path.relative_to(ROOT))})
        inputs.update({str(p.relative_to(ROOT)):cert.hashed(p) for p in (path,proof_path)})
    return {'schema':'elliptic-curves.exact-r17-maximum-parity-classes.v1',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),ROOT/'elliptic-curves/cas/replay_r17_norm12_minima.py',ROOT/'elliptic-curves/cas/exact_parity_ellipsoid.py')},
        'inputs':inputs,'families':rows,'parity_upper_bounds_replayed':786432,'exact_maximum_classes':sum(len(r['classes']) for r in rows),
        'claim_boundary':'Exact maximum of the coset minima on each displayed integral Gram lattice modulo twice that lattice. Full exceptional class lists are43,49,43,43,43,49 in atlas order. This is not a continuous covering-radius theorem, full Mordell-Weil saturation statement, or specialized rank claim.'}


def check(path):
    data=cert.read(path)
    if data!=expected():raise ArithmeticError('exact parity dataset differs')
    with TemporaryDirectory(prefix='r17-parity-proof-') as temporary:
        for f in data['families']:
            output=Path(temporary)/(f['family']+'.json')
            replay(ROOT/f['query_result_path'],output)
            if cert.read(output)!=cert.read(ROOT/f['proof_path']):raise ArithmeticError('complete geometry replay differs')
    print('REPLAYED ALL SIX EXACT PARITY GEOMETRIES',data['exact_maximum_classes'],'maximum classes',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--check',type=Path);a=p.parse_args()
    if a.check:check(a.check)
    else:
        if a.output.exists():raise FileExistsError('preserve parity dataset')
        cert.write(a.output,expected())
