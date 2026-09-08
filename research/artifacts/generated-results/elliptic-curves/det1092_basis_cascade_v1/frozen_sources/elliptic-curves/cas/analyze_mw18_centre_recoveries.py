#!/usr/bin/env python3
"""Certify joint gains of the already recovered policy bases; no point search.

An increased joint lower bound shows complementary discoveries. Equal finite
ranks do not prove equality of rational spans or integral subgroups.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

from research_runtime.mw_state import MWState
from research_runtime.search_state import reduction_cache
from research_runtime.store import checkpoint
from run_mw18_centre_experiment import check_protocol


def analyze(protocol_path,directory,output,verify=False):
    protocol=json.loads(protocol_path.read_text());check_protocol(protocol)
    cache=reduction_cache();old=json.loads(output.read_text()) if verify else None
    if old:cache.store.import_snapshot(old['arithmetic_facts'])
    from mod2_reduction_independence import _is_prime
    primes=[p for p in range(3,protocol['limits']['prime_bound']+1) if _is_prime(p)]
    rows=[]
    for case in protocol['cases']:
        paths=[directory/'cells'/f'{case}--{policy}.json' for policy in protocol['policies']]
        if not all(p.exists() for p in paths):continue
        cells=[json.loads(p.read_text()) for p in paths]
        if any(c['status']!='COMPLETE' for c in cells):continue
        if any(c['protocol_hash']!=protocol['protocol_hash'] or c['case']!=case for c in cells):
            raise ArithmeticError('foreign recovery cell')
        for cell in cells:cache.store.import_snapshot(cell['arithmetic_facts'])
        initial=MWState.from_record(cells[0]['initial_state'],cache=cache)
        if any(c['initial_state']!=cells[0]['initial_state'] for c in cells):raise ArithmeticError('policy initial states differ')
        state=initial;individual={}
        for policy,cell in zip(protocol['policies'],cells):
            recovered=MWState.from_record(cell['final_state'],cache=cache)
            if recovered.basis[:18]!=initial.basis:raise ArithmeticError('policy changed inherited generators')
            individual[policy]=recovered.rank-18
            for point in recovered.basis[18:]:state=state.adjoin(point,cache=cache,extra_primes=primes)
        rows.append({'case':case,'individual_gains':individual,'joint_certified_gain':state.rank-18,
            'joint_gain_exceeds_each_individual_gain':state.rank-18>max(individual.values()),
            'final_state':state.record(),
            'input_cell_hashes':{p.name:sha256(p.read_bytes()).hexdigest() for p in paths}})
        print(f"MW18_JOINT|{case}|individual={list(individual.values())}|joint={state.rank-18}",flush=True)
    result={'schema':'elliptic-curves.mw18-joint-recoveries.v1','protocol_hash':protocol['protocol_hash'],
        'status':'COMPLETE' if len(rows)==len(protocol['cases']) else 'PARTIAL',
        'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'cases':rows,'arithmetic_facts':cache.store.snapshot(),
        'claim_boundary':'Certified lower bounds for unions of retained independent policy bases. Equal bounds do not prove equal rational spans. No extra chart was searched and this diagnostic does not alter the frozen policy ranking.'}
    if verify:
        if result['cases']!=old['cases'] or result['status']!=old['status']:raise ArithmeticError('joint recovery replay differs')
    else:checkpoint(output,result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--protocol',type=Path,required=True);p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--verify',action='store_true')
    a=p.parse_args();analyze(a.protocol,a.directory,a.output,a.verify)
