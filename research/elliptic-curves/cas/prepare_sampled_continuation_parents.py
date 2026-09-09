#!/usr/bin/env python3
"""Exact generic parents outside the span of the productive subset."""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
from sage.all import ZZ, matrix
from memory_rank_certificate import checked_rank
from visibility_generic_bank import enumerate_bank
from visibility_lattice_fast import IntegerExactParity
from visibility_lattice_v2 import ExactParity
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]


def run(preparation, output, continuation):
    read=lambda p:json.loads(p.read_text())
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    sealed=read(preparation/'prepared.json')
    assert all(sha(preparation/n)==h for n,h in sealed['files'].items())
    bank=read(preparation/'anchor-bank.json')
    assert bank['status']=='COMPLETE_FROZEN_PRODUCTIVE_SUBSET'
    seed_path=next(preparation.glob('seed-M*.json'));seed=read(seed_path)
    terminal=read(continuation/'terminal.json');verified=read(continuation/'verified.json')
    assert verified['status']=='PASS_INDEPENDENT_PRODUCTIVE_V3_REPLAY' and verified['terminal_sha256']==sha(continuation/'terminal.json')
    assert terminal['curve']==seed['curve'] and terminal['points'][:len(seed['points'])]==seed['points']
    seed={k:terminal[k] for k in ('curve','points','proof','rank_lower_bound')};seed['generic_rank']=bank['dimension']
    rank,dimension=seed['rank_lower_bound'],bank['dimension']
    model=tuple(map(F,seed['curve']));points=tuple(tuple(map(F,p)) for p in seed['points'])
    old=seed['proof'];fresh=checked_rank(model,points,[r['prime'] for r in old['signatures']],old['no_rational_2_torsion_prime'])
    assert json.loads(json.dumps(fresh))==old
    previous=read(preparation/'protocol.json')
    geometry=read(preparation/'generic-cvp-proofs.json')
    g,u=matrix(ZZ,geometry['gram']),matrix(ZZ,geometry['LLL'])
    assert g.nrows()==dimension and abs(u.det())==1
    inv=u.inverse();fast=IntegerExactParity((u*g*u.transpose()).rows())
    reference=ExactParity((u*g*u.transpose()).rows())
    old_masks={r['mask'] for r in bank['rows']}
    def reduce(mask,pivots):
        for k in sorted(pivots,reverse=True):
            if (mask>>k)&1:mask ^= pivots[k]
        return mask
    def insert(mask,pivots):
        mask=reduce(mask,pivots)
        if mask:pivots[mask.bit_length()-1]=mask
        return bool(mask)
    pivots={}
    for mask in sorted(old_masks):insert(mask,pivots)
    original_pivots=dict(pivots)
    bound=max(r['norm'] for r in bank['rows'])
    candidate_masks=[]
    i=0
    while len(candidate_masks)<128:
        raw=hashlib.sha256((previous['family']+'/outside-productive-span/'+str(i)).encode()).digest()
        i+=1
        mask=int.from_bytes(raw,'big') & ((1<<dimension)-1)
        if reduce(mask,original_pivots) and mask not in candidate_masks:candidate_masks.append(mask)
    candidate_rows=[];candidate_checks=[]
    for mask in candidate_masks:
        w=matrix(ZZ,1,dimension,[(mask>>j)&1 for j in range(dimension)])
        residue=tuple(int(x)%2 for x in (w*inv).row(0))
        starts,_=fast.babai(np.asarray([residue],dtype=np.int64));start=tuple(map(int,starts[0]))
        cert=fast.solve(residue,start,2000000)
        assert cert==reference.solve(residue,start,2000000)
        minima=[]
        for v in cert['minima']:
            z=tuple(map(int,(matrix(ZZ,1,dimension,v)*u).row(0)))
            if next(x for x in z if x)<0:z=tuple(-x for x in z)
            assert sum((x%2)<<j for j,x in enumerate(z))==mask
            minima.append(z)
        candidate_rows.append({'mask':mask,'word':list(min(minima)),'norm':cert['norm']})
        candidate_checks.append({'mask':mask,'reduced_seed':start,'proof':cert})
    complete={'status':'COMPLETE_FIXED_128_PARITY_MINIMA','rows':candidate_rows,'checks':candidate_checks,
        'claim_boundary':'All128 frozen sampled parities verified, not a complete shell enumeration.'}
    historical_bound=bound
    bound=max(r['norm'] for r in candidate_rows)
    candidates=sorted((r for r in candidate_rows if r['norm']==bound),key=lambda r:r['mask'])
    chosen=[]
    for row in candidates:
        if insert(row['mask'],pivots):chosen.append(row)
        if len(chosen)==16:break
    for row in candidates:
        if len(chosen)==16:break
        if row not in chosen:chosen.append(row)
    assert 0<len(chosen)<=16
    words={r['mask']:r['word'] for r in chosen}
    derivation={'original_masks':sorted(old_masks),'new_masks':sorted(words),
        'initial_span_rank':len(original_pivots),'extended_span_rank':len(pivots),
        'generic_scaled_shell':bound,'historical_productive_maximum':historical_bound,'candidate_count_outside_original_span':len(candidates),
        'rule':'Hash family/outside-productive-span/index with SHA256 into128 distinct parities outside the productive span. '
               'Certify all128 minima with both solvers, retain the largest sampled exact minimum, and in ascending mask order extend the productive F2 span; '
               'then fill to16 from remaining masks outside the original span. No specialized point labels.'}
    paths=[continuation/'terminal.json',continuation/'verified.json']+[preparation/n for n in sealed['files']]+[preparation/'prepared.json',Path(__file__),
        Path(__file__).with_name('visibility_generic_bank.py'), Path(__file__).with_name('visibility_lattice_fast.py'),Path(__file__).with_name('visibility_lattice_v2.py')]
    protocol={'family':previous['family'],'parameter':previous['parameter'],'initial_rank':rank,
        'generic_rank':dimension,'point_searches':0,'cvp_node_limit':2000000,
        'rule':derivation['rule'],
        'sampled_candidate_count':128, 'scaled_target_norm':bound,
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
    output.mkdir(exist_ok=False);checkpoint(output/'protocol.json',protocol)
    checkpoint(output/f'seed-M{rank}.json',seed)
    checkpoint(output/'sampled-parity-bank.json',complete)
    derivation['sampled_parities_sha256']=sha(output/'sampled-parity-bank.json')
    checkpoint(output/'derivation.json',derivation)
    rows,checks=[],[]
    for mask in sorted(words):
        w=matrix(ZZ,1,dimension,words[mask]);residue=tuple(int(x)%2 for x in (w*inv).row(0))
        starts,_=fast.babai(np.asarray([residue],dtype=np.int64));start=tuple(map(int,starts[0]))
        cert=fast.solve(residue,start,2000000)
        assert cert==reference.solve(residue,start,2000000)
        minima=[]
        for v in cert['minima']:
            z=tuple(map(int,(matrix(ZZ,1,dimension,v)*u).row(0)))
            if next(x for x in z if x)<0:z=tuple(-x for x in z)
            assert sum((x%2)<<j for j,x in enumerate(z))==mask
            minima.append(z)
        rows.append({'mask':mask,'word':list(min(minima)),'norm':cert['norm']})
        checks.append({'mask':mask,'reduced_seed':start,'proof':cert})
        checkpoint(output/f'cvp-{mask}.json',checks[-1])
    checkpoint(output/'generic-cvp-proofs.json',{'gram':geometry['gram'],'LLL':geometry['LLL'],'checks':checks})
    result={'status':'COMPLETE_FROZEN_COMPLEMENT_PARENT_SUBSET','dimension':dimension,
        'generic_gram_scale':2,'rows':rows,'shells':sorted({r['norm'] for r in rows}),
        'protocol_sha256':sha(output/'protocol.json'),'derivation_sha256':sha(output/'derivation.json'),
        'generic_cvp_sha256':sha(output/'generic-cvp-proofs.json'),
        'claim_boundary':'New exact parent classes outside the historical productive span; productivity '
                         'is untested. Not a complete generic shell bank or rank prediction.'}
    checkpoint(output/'anchor-bank.json',result)
    assert all(sha(ROOT/n)==h for n,h in protocol['inputs'].items())
    checkpoint(output/'prepared.json',{'status':'PASS_COMPLEMENT_PARENT_PREPARATION','point_searches':0,
        'anchor_count':len(rows),'files':{p.name:sha(p) for p in sorted(output.glob('*.json'))}})
    print('COMPLEMENT_PARENTS_PREPARED',len(rows),'scaled shells',result['shells'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--preparation',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--continuation',type=Path,required=True)
    a=p.parse_args();run(a.preparation.resolve(),a.output.resolve(),a.continuation.resolve())
