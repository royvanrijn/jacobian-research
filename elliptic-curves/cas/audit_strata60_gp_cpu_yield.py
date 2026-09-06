#!/usr/bin/env python3
"""Secondary cost sensitivity from already audited timing markers; no search."""
import argparse
from fractions import Fraction
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
REPORT=ART/'retained_mw16_score_strata_experiment_v1.json'
AUDIT=ART/'retained_mw16_score_strata_accounting_replay_v1.json'
OUT=ART/'strata60_gp_cpu_yield_v1.json'


def expected():
    report=cert.read(REPORT);audit=cert.read(AUDIT)
    if report['status']!='COMPLETE_FIXED_COMPARISON' or audit['status']!='PASS':
        raise ArithmeticError('terminal independently audited comparison required')
    for data in (report,audit):
        if any(cert.hashed(ROOT/n)!=h for n,h in data['sources'].items()):
            raise ArithmeticError('audited timing evidence changed')
    rows=report['rows']
    if len(rows)!=60 or any(r['charts_without_cpu_time']!=0 or r['completed_boxes']!=43
                           or r['certified_gain'] is None for r in rows):
        raise ArithmeticError('all60 certified complete exposures and CPU markers required')
    arms={}
    for arm in ('top','moderate','lower'):
        selected=[r for r in rows if r['arm']==arm]
        if len(selected)!=20:raise ArithmeticError('equal20-curve arms required')
        arms[arm]={'certified_gain':sum(r['certified_gain'] for r in selected),
                   'gp_cpu_ms':sum(r['known_search_cpu_ms'] for r in selected)}
    gain=sum(a['certified_gain'] for a in arms.values())
    cpu=sum(a['gp_cpu_ms'] for a in arms.values())
    ratio=Fraction(arms['top']['certified_gain']*cpu,arms['top']['gp_cpu_ms']*gain)
    return {'schema':'elliptic-curves.strata60-gp-cpu-yield.v1','status':'PASS',
            'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in
                       (Path(__file__).resolve(),REPORT,AUDIT)},
            'arms':arms,'total_gp_cpu_ms':cpu,'top_to_pooled_yield_ratio':str(ratio),
            'claim_boundary':'Post-terminal secondary cost sensitivity, not a new policy criterion. All GP timing markers were independently audited. GP CPU measures point enumeration only; it excludes Python admission, map preparation and verification, so it is not total computation. No new search, rank claim, validation-prime input or change to the frozen matched comparison.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=d:raise ArithmeticError('CPU sensitivity replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve CPU sensitivity certificate')
        checkpoint(OUT,d)
    print('RECORDED GP CPU SENSITIVITY PASS',d['top_to_pooled_yield_ratio'])
