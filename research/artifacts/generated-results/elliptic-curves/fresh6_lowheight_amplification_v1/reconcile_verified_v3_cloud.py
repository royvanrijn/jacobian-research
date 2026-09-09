#!/usr/bin/env python3
"""Extend a replayed terminal basis using its saved audited clouds, without search."""
import argparse,json,hashlib
from pathlib import Path
from fractions import Fraction as F
from future_point_admission import FinitePointAdmission
from memory_rank_certificate import checked_rank
from v3_warm_engine import certified_state
from pointed_quartic_search import PointedQuarticSearch
import pari_pointed_backend as backend


def run(folder,output):
    read=lambda p:json.loads(p.read_text())
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    terminal=read(folder/'terminal.json');verified=read(folder/'verified.json')
    if verified['status']!='PASS_INDEPENDENT_PRODUCTIVE_V3_REPLAY' or verified['terminal_sha256']!=sha(folder/'terminal.json'):raise ArithmeticError('verified terminal required')
    model=tuple(map(F,terminal['curve']));basis=tuple(tuple(map(F,p)) for p in terminal['points'])
    admission=FinitePointAdmission(model,basis);bindings={};gains=[]
    for n in ['terminal.json','verified.json','protocol.json']:bindings[n]=sha(folder/n)
    for stage in terminal['stages']:
        wd=folder/f"epoch-{stage['epoch']:02d}";seed=read(wd/'seed.json')
        state=certified_state(model,tuple(tuple(map(F,p)) for p in seed['points']),seed['proof'])
        for i in range(stage['charts']):
            p=wd/f'chart-{i:04d}.json';chart=read(p);bindings[str(p.relative_to(folder))]=sha(p)
            mapping=chart['mapping'];search=PointedQuarticSearch(state=state,centre={'coefficients':mapping['centre']['representative']},coordinate_policy=mapping['coordinate_policy'])
            for point in backend.replay(search,mapping,chart['search']):
                before=len(admission.points);admission.consider(point)
                if len(admission.points)>before:gains.append({'chart':str(p.relative_to(folder)),'point':list(map(str,point)),'after':len(admission.points)})
    proof=checked_rank(model,admission.points,admission.primes,terminal['proof']['no_rational_2_torsion_prime'])
    result={'status':'PASS_REPLAYED_CLOUD_RECONCILIATION','initial_rank':len(basis),'rank_lower_bound':len(admission.points),'generic_rank':read(folder/'protocol.json')['generic_rank'],'curve':terminal['curve'],'points':[list(map(str,p)) for p in admission.points],'proof':proof,'gains':gains,'bindings':bindings,'source_sha256':sha(Path(__file__)),'point_searches':0,'claim_boundary':'Exact subgroup lower bound from saved point witnesses, not an exact rank or record.'}
    result=json.loads(json.dumps(result))
    if output.exists():
        if read(output)!=result:raise ArithmeticError('independent reconciliation differs')
    else:
        output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS_CLOUD_RECONCILIATION',len(basis),len(admission.points),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--folder',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();run(a.folder.resolve(),a.output.resolve())
