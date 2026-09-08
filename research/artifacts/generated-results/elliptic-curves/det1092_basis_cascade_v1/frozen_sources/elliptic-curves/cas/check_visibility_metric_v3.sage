#!/usr/bin/env sage-python
"""Independent recomputation of every V3 height decision form and LLL transport."""
from pathlib import Path
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
from sage.all import ZZ, matrix, pari
from research_runtime.store import checkpoint

import argparse
ap=argparse.ArgumentParser();ap.add_argument('--start',type=int,choices=[17,30],default=30)
args=ap.parse_args()
CAS = Path(__file__).resolve().parent
v2 = SourceFileLoader('metric_replay_v2',str(CAS/'adaptive_visibility_cascade_v3.sage')).load_module()
v2.v1.guard()
policy = v2.read(v2.D/'protocol.json')
assert policy['sources']==v2.sources()
assert all(v2.sha(v2.ROOT/k)==h for k,h in policy['inputs'].items())
geo = v2.load('prospective_half_lattice_v3.sage')
model, _ = v2.v1.seed({'parameter':'0','presentation':'normalized'})
stages = []
for path in sorted((v2.D/('terminal-M30' if args.start==30 else 'replay-M17')).glob('epoch-*/selection.json')):
    selection = v2.read(path)
    basis = tuple(tuple(map(F,q)) for q in selection['basis'])
    height, asym = geo.canonical_height_gram(model,basis)
    g = matrix(ZZ,geo.rounded_gram(height,1000000))
    assert [list(map(int,row)) for row in g.rows()]==selection['rounded_gram']
    assert str(asym)==selection['height_asymmetry']
    u = matrix(ZZ,pari(g).qflllgram()).transpose()
    assert [list(map(int,row)) for row in u.rows()]==selection['LLL']
    stages.append({'rank':len(basis),'selection':v2.rel(path),'sha256':v2.sha(path)})
    print('RECOMPUTED HEIGHT METRIC AND LLL',len(basis),flush=True)
checkpoint(v2.D/('metric-replay-M30.json' if args.start==30 else 'metric-replay-M17.json'),{'protocol_sha256':v2.sha(v2.D/'protocol.json'),
    'checker_sha256':v2.sha(Path(__file__)),'stages':stages,
    'claim':'Reproduced 384-bit height computation and exact rounding/LLL decisions; not an interval-certified exact canonical-height CVP.'})

