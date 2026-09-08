#!/usr/bin/env sage-python
"""Freeze a four-fibre adaptive test with any initial gain, rather than rank 24.

The earlier blind ladder recovered rank 29 from initial rank 20 and rank 28
from initial rank 21. A high initial-rank gate therefore discards demonstrated
adaptive opportunities. This finite cohort follows score order after any
positive exact gain; no public exceptional point is read.
"""
import argparse
from importlib.machinery import SourceFileLoader
from pathlib import Path
from fractions import Fraction as F
import sys

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path[:0]=[str(CAS),str(ROOT/'elliptic-curves')]
follow=SourceFileLoader('positive_followup_engine',str(CAS/'compact_r17_ladder_followup.sage')).load_module()
m=follow.m;certificate=follow.certificate
from research_runtime.store import checkpoint
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/compact-r17-positive-followup-v1'


def prepare(directory):
    if (directory/'protocol.json').exists():raise FileExistsError('positive-gain cohort already frozen')
    parent=m.read(follow.DIRECTORY/'protocol.json');used={r['parameter'] for r in parent['rows']};inventory={}
    for name in ('compact-r17-top64-v1','compact-r17-h16384-v1'):
        folder=ROOT/'artifacts/local/elliptic-curves'/name
        pop=m.read(folder/'population.json');ledger=m.read(folder/'ledger.json')
        for r in ledger['rows']:
            if r.get('rank_lower_bound',0)<=17 or r['status'] not in ('COMPLETE','REUSED_COMPLETE_MEASUREMENT') or r['parameter'] in used:continue
            candidate=pop['finalists'][r['index']]
            inventory[r['parameter']]={'parameter':r['parameter'],'score_units':candidate['score_units'],
                'rank_lower_bound':r['rank_lower_bound'],'path':r['path'],'sha256':m.hashed(ROOT/r['path'])}
    selected=sorted(inventory.values(),key=lambda r:(-r['score_units'],r['parameter']))[:4]
    if len(selected)!=4:raise ArithmeticError('positive-gain cohort incomplete')
    rows=[];database=certificate.read(certificate.DATABASE)
    for r in selected:
        path=ROOT/r['path'];d=m.read(path)
        if m.hashed(path)!=r['sha256'] or d['status']!='COMPLETE':raise ArithmeticError('selected measurement changed')
        model=tuple(map(F,d['curve']));points=[tuple(map(F,p)) for p in d['final_state']['state']['reductions']['points']]
        rank=certificate.checked_rank(model,points);certificate.family_check(d['parameter'],model,points)
        if any(certificate.isomorphic(model,x['ainvs']) for x in database['curves']):raise ArithmeticError('selected curve already in pinned catalogue')
        rows.append({'role':'prospective_positive_initial_gain','parameter':d['parameter'],'curve':d['curve'],
            'points':[list(map(str,p)) for p in points],'rank_certificate':rank,'input_path':r['path'],'input_sha256':r['sha256']})
    protocol={**parent,'schema':'elliptic-curves.compact-r17-positive-followup.v1','sources':follow.sources(),
        'selector_source_sha256':m.hashed(Path(__file__)),
        'selection_rule':'four highest full-prime scores among completed uncatalogued fibres with rank >17, excluding the already running rank-24 cohort',
        'selection_inventory':sorted(inventory.values(),key=lambda r:(-r['score_units'],r['parameter'])),
        'parent_adaptive_protocol_sha256':m.hashed(follow.DIRECTORY/'protocol.json'),
        'database_sha256':m.hashed(certificate.DATABASE),'rows':rows}
    checkpoint(directory/'protocol.json',protocol)
    checkpoint(directory/'population.json',{'protocol_hash':m.identity(protocol),'finalists':rows})
    print('POSITIVE FOLLOWUP FROZEN',[(r['parameter'],len(r['points'])) for r in rows],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run']);p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--index',type=int,default=0);a=p.parse_args()
    if a.stage=='prepare':prepare(a.directory)
    else:
        protocol=m.read(a.directory/'protocol.json')
        if protocol['selector_source_sha256']!=m.hashed(Path(__file__)):raise ArithmeticError('positive selector source changed')
        follow.run(a.directory,a.index)
