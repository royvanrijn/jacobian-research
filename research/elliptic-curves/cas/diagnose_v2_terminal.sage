#!/usr/bin/env sage-python
"""Retrospective terminal diagnostics; never writes into frozen V2."""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import ZZ, QQ, RealField, matrix, vector
from research_runtime.store import checkpoint

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT/'elliptic-curves/cas'
ART = ROOT/'artifacts/generated-results/elliptic-curves'
V2 = ROOT/'artifacts/local/elliptic-curves/adaptive-visibility-cascade-v2'
D = ROOT/'artifacts/local/elliptic-curves/v2-terminal-retrospective-v1'
def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p): return str(p.relative_to(ROOT))
def load(n): return SourceFileLoader('terminal_'+n.replace('.','_'),str(CAS/n)).load_module()


def freeze():
    if D.exists(): raise FileExistsError('preserve retrospective diagnostic')
    report = read(ART/'adaptive_visibility_cascade_v2.json')
    assert report['rank_lower_bound']==30 and report['status']=='COMPLETE_INDEPENDENTLY_REPLAYED_BOUNDED_EXPERIMENT'
    paths = [p for p in V2.rglob('*') if p.is_file()]
    paths += [ART/'adaptive_visibility_cascade_v2.json',CAS/'adaptive_visibility_cascade_v2.sage',
        CAS/'visibility_lattice_v2.py',CAS/'check_visibility_cascade_v2.sage',
        CAS/'check_visibility_metric_v2.sage',CAS/'finalize_visibility_cascade_v2.py',
        ROOT/'elliptic-curves/notes/ADAPTIVE_HALF_LATTICE_V2_2026-09-07.md']
    D.mkdir(parents=True)
    checkpoint(D/'V2-seal.json',{rel(p):sha(p) for p in paths})
    checkpoint(D/'fixed-M30.json',{'curve':report['final_rank_certificate']['curve'],
        'points':report['final_rank_certificate']['independent_points'],
        'certificate':report['final_rank_certificate'],'source_sha256':sha(ART/'adaptive_visibility_cascade_v2.json')})
    print('SEALED V2',len(paths),'files',flush=True)


def recognize():
    if (D/'coordinates.json').exists(): raise FileExistsError('preserve recognized coordinates')
    a2 = load('audit_curve302_exceptional_subgroup_landscape.sage')
    geo = load('prospective_half_lattice_v3.sage')
    group = load('half_lattice_pointed_sieve.py')
    fixed = read(D/'fixed-M30.json')
    model = tuple(map(F,fixed['curve']))
    basis = tuple(tuple(map(F,p)) for p in fixed['points'])
    vis = read(a2.VISIBILITY)
    targets = [a2.primary_target(row,basis[:17],model,group) for row in vis['directions']]
    ambient = (*basis[:17],*targets)
    field = RealField(192)
    heights = geo.canonical_height_gram(model,ambient)[0]
    h = matrix(field,[[field(str(x)) for x in row] for row in heights])
    words = [list(map(int,row)) for row in matrix.identity(ZZ,31).rows()[:17]]
    for point in basis[17:]:
        paired = geo.canonical_height_gram(model,(*ambient,point))[0]
        coeff = h.solve_right(vector(field,[field(str(paired[i][31])) for i in range(31)]))
        word = [int(x.round()) for x in coeff]
        assert max(abs(x-y) for x,y in zip(coeff,word))<field('1e-30')
        assert group.linear_combination(model,ambient,word)==point
        words.append(word)
    w = matrix(ZZ,words)
    choices = []
    for i in range(17,31):
        det = int(w.stack(matrix.identity(ZZ,31)[i:i+1,:]).det())
        if det: choices.append((abs(det),i,det))
    index,i,det = min(choices)
    result = {'model':fixed['curve'],'ambient_points':[list(map(str,p)) for p in ambient],
        'ambient_labels':['generic-%02d'%i for i in range(17)]+[r['id'] for r in vis['directions']],
        'M30_words':words,'missing_axis_index':i,'missing_axis_label':vis['directions'][i-17]['id'],
        'missing_point':list(map(str,ambient[i])),'completion_index_in_agent2_lattice':index,
        'noncontained_axis_completions':choices,
        'oracle_table_sha256':sha(a2.OUTPUT)}
    checkpoint(D/'coordinates.json',result)
    print('EXACT COORDINATES; missing lift',result['missing_axis_label'],'completion index',index,flush=True)


def certify():
    from memory_rank_certificate import checked_rank
    from mod2_reduction_independence import _primes_up_to
    from research_runtime.finite_reduction import ReductionCache
    from research_runtime.memory_store import MemoryFactStore
    fixed=read(D/'fixed-M30.json');data=read(D/'coordinates.json')
    points=tuple(tuple(map(F,p)) for p in fixed['points'])+(tuple(map(F,data['missing_point'])),)
    model=tuple(map(F,fixed['curve']));cache=ReductionCache(MemoryFactStore());primes=[]
    for p in _primes_up_to(1000):
        if p==2:continue
        try:cache.signature(model,points,p)
        except ValueError:continue
        primes.append(p)
    proof=checked_rank(model,points,primes,fixed['certificate']['rank_certificate']['no_rational_2_torsion_prime'])
    assert proof['rank_lower_bound']==31
    checkpoint(D/'B31-rank-certificate.json',{'points':[list(map(str,p)) for p in points],'proof':proof,
        'consequence':'All one- and two-generator omission subgroups in state-plan.json inherit exact independence from this 31-generator certificate.'})
    print('CERTIFIED FIXED B31 AND ALL OMISSION SUBGROUPS')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['freeze','recognize','certify','check-seal']);a=ap.parse_args()
    if a.action=='freeze': freeze()
    elif a.action=='recognize': recognize()
    elif a.action=='certify': certify()
    else:
        seal=read(D/'V2-seal.json');assert all(sha(ROOT/k)==h for k,h in seal.items());print('V2 SEAL PASSES',len(seal))
