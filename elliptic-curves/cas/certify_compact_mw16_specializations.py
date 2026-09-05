#!/usr/bin/env python3
"""Five fixed t=1 specialization checks for the new compact MW16 inputs.

These are usability and independence checks, not a discovery population.
All sixteen points come from the generic sections; no extra point is sought.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=ROOT/'artifacts/generated-results/elliptic-curves/compact_five_mw16_specializations_v1.json'


def sources():
    paths=(spec.ATLAS,Path(__file__).resolve(),Path(spec.__file__).resolve(),Path(cert.__file__).resolve(),
           ROOT/'elliptic-curves/cas/mod2_reduction_independence.py')
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}


def verify(record,family):
    model,points=spec.specialize(family,record['parameter'])
    if list(map(str,model))!=record['curve'] or [list(map(str,p)) for p in points]!=record['points']:
        raise ArithmeticError('specialized input differs')
    proof=record['rank_certificate']
    result=cert.checked_rank(model,points,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
    if result!=proof: raise ArithmeticError('independence certificate differs')


def build(output):
    if output.exists(): raise FileExistsError('preserve specialization certificate')
    rows=[]
    for family in cert.read(spec.ATLAS)['families']:
        model,points=spec.specialize(family,'1')
        row={'fibration_id':family['fibration_id'],'presentation_id':family['presentation_id'],
             'parameter':'1','curve':list(map(str,model)), 'points':[list(map(str,p)) for p in points],
             'rank_certificate':cert.checked_rank(model,points),
             'short_coefficient_bits':max(max(abs(q.numerator).bit_length(),q.denominator.bit_length()) for q in model)}
        verify(row,family);rows.append(row)
        print('CERTIFIED FIXED MW16',row['fibration_id'],'rank >=16','bits',row['short_coefficient_bits'],flush=True)
    cert.write(output,{'schema':'elliptic-curves.compact-five-mw16-specializations.v1','sources':sources(),'rows':rows,
        'claim_boundary':'Five fixed t=1 usability checks certify independence of their sixteen specialized generic points. No extra point search or novelty comparison. Any relation over Q(t) would specialize to a relation here, so these also independently witness generic independence of the exported point lists.'})


def check(path):
    data=cert.read(path)
    if data['sources']!=sources(): raise ArithmeticError('sources changed')
    families={f['fibration_id']:f for f in cert.read(spec.ATLAS)['families']}
    if len(data['rows'])!=5 or {r['fibration_id'] for r in data['rows']}!=set(families):
        raise ArithmeticError('fixed roster changed')
    for row in data['rows']:
        if row['parameter']!='1': raise ArithmeticError('fixed parameter changed')
        verify(row,families[row['fibration_id']])
        print('REPLAYED FIXED MW16',row['fibration_id'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--check',type=Path)
    a=p.parse_args();check(a.check) if a.check else build(a.output)
