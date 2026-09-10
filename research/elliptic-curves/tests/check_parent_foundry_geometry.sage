#!/usr/bin/env sage-python
"""Exact invariant/section corruption regressions on constructed preflight parents."""
import argparse
import copy
import json
import sys
import tempfile
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import PolynomialRing,QQ

CAS=Path(__file__).resolve().parents[1]/'cas';sys.path.insert(0,str(CAS))
geo=SourceFileLoader('parent_geometry_regression',str(CAS/'parent_foundry_geometry.sage')).load_module()
p=argparse.ArgumentParser();p.add_argument('--runtime',type=Path,required=True);args=p.parse_args()
data=json.loads((args.runtime/'parent-inputs/generic.json').read_text())
keys={p['invariant']['key'] for p in data['known_a1']}
if len(keys)!=5:raise ArithmeticError('nine-to-five collision regression failed')
R=PolynomialRing(QQ,'t');t=R.gen()
for i in (0,1):
    path=args.runtime/'preflight'/f'construct-{i}'/'parent.json'
    parent=json.loads(path.read_text());a,b=R(parent['raw_A']),R(parent['raw_B'])
    invariant=geo.affine_invariant(a,b)
    if invariant['key'] in keys:raise ArithmeticError('new fibration is not separated')
    keys.add(invariant['key'])
    if geo.affine_invariant(a(3*t+7)/5**4,b(3*t+7)/5**6)!=invariant:
        raise ArithmeticError('affine/scaling invariance failed')
    if geo.affine_invariant(a+1,b)==invariant:
        raise ArithmeticError('perturbed equation was not separated')
    with tempfile.TemporaryDirectory() as name:
        bad=copy.deepcopy(parent);s=bad['sections'][0]['Y']['numerator_coefficients_low_to_high']
        s[0]=str(QQ(s[0])+1);out=Path(name)/'parent.json';out.write_text(json.dumps(bad))
        try:geo.verify(out)
        except ArithmeticError:pass
        else:raise ArithmeticError('corrupt section accepted')
    print('PASS_EXACT_PARENT_REGRESSION',parent['family'],flush=True)
