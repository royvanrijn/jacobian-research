#!/usr/bin/env python3
"""Exact new parent classes from pairwise XOR of a sealed productive subset."""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
from sage.all import ZZ, matrix
from memory_rank_certificate import checked_rank
from visibility_lattice_fast import IntegerExactParity
from visibility_lattice_v2 import ExactParity
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]


def run(preparation, output):
    read=lambda p:json.loads(p.read_text())
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    sealed=read(preparation/'prepared.json')
    assert all(sha(preparation/n)==h for n,h in sealed['files'].items())
    bank=read(preparation/'anchor-bank.json')
    assert bank['status']=='COMPLETE_FROZEN_PRODUCTIVE_SUBSET'
    seed_path=next(preparation.glob('seed-M*.json'));seed=read(seed_path)
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
    old_masks={r['mask'] for r in bank['rows']};derivation={};words={}
    for a,b in itertools.combinations(bank['rows'],2):
        mask=a['mask']^b['mask']
        if mask==0 or mask in old_masks:continue
        derivation.setdefault(mask,[]).append([a['mask'],b['mask']])
        words[mask]=[x+y for x,y in zip(a['word'],b['word'])]
    assert 0<len(words)<=64, 'bounded new-parent preparation required'
    paths=[preparation/n for n in sealed['files']]+[preparation/'prepared.json',Path(__file__),
        Path(__file__).with_name('visibility_lattice_fast.py'),Path(__file__).with_name('visibility_lattice_v2.py')]
    protocol={'family':previous['family'],'parameter':previous['parameter'],'initial_rank':rank,
        'generic_rank':dimension,'point_searches':0,'cvp_node_limit':2000000,
        'rule':'All distinct nonzero pairwise XORs of the sealed productive masks, excluding '
               'the original masks. Exact generic minima with independent rational CVP replay. '
               'No known higher point or specialized visibility label enters selection.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
    output.mkdir(exist_ok=False);checkpoint(output/'protocol.json',protocol)
    checkpoint(output/seed_path.name,seed)
    checkpoint(output/'derivation.json',{'original_masks':sorted(old_masks),
        'new_masks':sorted(words),'pairs':derivation})
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
    result={'status':'COMPLETE_FROZEN_PAIRWISE_PARENT_SUBSET','dimension':dimension,
        'generic_gram_scale':2,'rows':rows,'shells':sorted({r['norm'] for r in rows}),
        'protocol_sha256':sha(output/'protocol.json'),'derivation_sha256':sha(output/'derivation.json'),
        'generic_cvp_sha256':sha(output/'generic-cvp-proofs.json'),
        'claim_boundary':'New exact parent classes derived from historical ones; productivity '
                         'is untested. Not a complete generic shell bank or rank prediction.'}
    checkpoint(output/'anchor-bank.json',result)
    assert all(sha(ROOT/n)==h for n,h in protocol['inputs'].items())
    checkpoint(output/'prepared.json',{'status':'PASS_PAIRWISE_PARENT_PREPARATION','point_searches':0,
        'anchor_count':len(rows),'files':{p.name:sha(p) for p in sorted(output.glob('*.json'))}})
    print('PAIRWISE_PARENTS_PREPARED',len(rows),'scaled shells',result['shells'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--preparation',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.preparation.resolve(),a.output.resolve())
