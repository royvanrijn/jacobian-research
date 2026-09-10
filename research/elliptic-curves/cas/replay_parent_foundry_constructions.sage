#!/usr/bin/env sage-python
"""Cold reconstruction and finite independence replay of the first two parents.

Rebuilds the complete equations/markings from generic source inputs in a new
temporary directory. No exceptional point enters construction. Then checks
the independent specialization certificates as a separate lower-bound proof.
"""
import argparse
import json
import sys
import tempfile
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import PolynomialRing,QQ

CAS=Path(__file__).resolve().parent;ROOT=CAS.parents[1];sys.path.insert(0,str(CAS))
geo=SourceFileLoader('parent_cold_reconstruction',str(CAS/'parent_foundry_geometry.sage')).load_module()
from verify_parent_foundry_certificate import verify


def run(directory):
    template=json.loads((ROOT/'elliptic-curves/data/a1_mw16_family_template_v1.json').read_text())
    R=PolynomialRing(QQ,'t');known=set();families={}
    for p in template['presentations']:
        a,b=(R(p['pencil'][k]) for k in ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
        key=geo.affine_invariant(a,b)['key'];known.add(key)
        families.setdefault(p['fibration_id'],set()).add(key)
    if len(known)!=5 or any(len(x)!=1 for x in families.values()):
        raise ArithmeticError('known five-fibration equivalence regression failed')
    for i in (0,1):
        path=directory/f'parent-{i}.json';parent=json.loads(path.read_text())
        with tempfile.TemporaryDirectory(prefix='parent-foundry-cold-') as name:
            rebuilt=geo.construct(parent['priority'],Path(name),parent['source_input'],parent['trace_word'])
            if rebuilt!=parent:raise ArithmeticError('cold generic reconstruction differs')
        key=parent['novelty_invariant']['key']
        if key in known:raise ArithmeticError('parent not separated from prior roster')
        known.add(key)
        verify(directory/f'point-{i}.json',path)
        print('PASS_COLD_CONSTRUCTION_AND_SPECIALIZATION',parent['family'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,default=ROOT/'artifacts/generated-results/elliptic-curves/parent-foundry-v2/commissioning')
    run(p.parse_args().directory)
