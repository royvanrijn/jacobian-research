#!/usr/bin/env python3
"""Bounded native R17 seed confirmation; stop on the first certified M18."""
import argparse
import fcntl
from fractions import Fraction as F
import json
from pathlib import Path
import shutil
import numpy as np
from sage.all import ZZ, matrix, pari
from v3_warm_engine import load, certified_state
from run_complement_seed_v3 import sha, read, require, guard
from research_runtime.store import checkpoint
from visibility_lattice_fast import IntegerExactParity
from visibility_lattice_v2 import ExactParity
from half_lattice_pointed_sieve import linear_combination
from future_point_admission import FinitePointAdmission
from memory_rank_certificate import checked_rank
from pointed_quartic_search import PointedQuarticSearch
from pointed_box_equivalence import box_key
from lean_preconditioned_map_receipts import obtain
import pari_pointed_backend as backend

CAS=Path(__file__).resolve().parent
ROOT=CAS.parents[1]
PARITY=ROOT/'artifacts/generated-results/elliptic-curves/r17_exact_maximum_parity_classes_v1.json'
ORDER=('preconditioned_full','factor_free')


def geometry(seed,masks,reference=False):
    geo=load('fresh6_height_geometry',CAS/'prospective_half_lattice_v3.sage')
    model=tuple(map(F,seed['curve']));points=tuple(tuple(map(F,p)) for p in seed['points'])
    gram,asymmetry=geo.canonical_height_gram(model,points)
    g=matrix(ZZ,geo.rounded_gram(gram,1000000))
    u=matrix(ZZ,pari(g).qflllgram()).transpose();inv=u.inverse()
    require(inv.denominator()==1,'nonunimodular specialized LLL')
    solver=(ExactParity if reference else IntegerExactParity)((u*g*u.transpose()).rows())
    rows=[]
    for mask in masks:
        word=matrix(ZZ,1,17,[(mask>>j)&1 for j in range(17)])
        residue=tuple(int(x)%2 for x in (word*inv).row(0))
        starts,_=solver.babai(np.asarray([residue],dtype=np.int64))
        start=tuple(map(int,starts[0]));proof=solver.solve(residue,start,2000000)
        vectors=[]
        for v in proof['minima']:
            z=tuple(map(int,(matrix(ZZ,1,17,v)*u).row(0)))
            if next(x for x in z if x)<0:z=tuple(-x for x in z)
            require(sum((x%2)<<j for j,x in enumerate(z))==mask,'parity transport differs')
            vectors.append(z)
        rep=min(vectors);point=linear_combination(model,points,rep)
        require(point is not None,'infinite centre')
        rows.append(dict(mask=mask,representative=list(rep),point=list(map(str,point)),
                         metric_norm=proof['norm'],cvp=proof,reduced_start=start))
    rows.sort(key=lambda r:(-r['metric_norm'],r['mask']))
    return json.loads(json.dumps(dict(gram=[list(map(int,r)) for r in g.rows()],
        LLL=[list(map(int,r)) for r in u.rows()],centres=rows)))


def prepare(seed_path,folder):
    seed=read(seed_path);dataset=seed_path.parent.parent
    prepared=read(dataset/'prepared.json');verified=read(dataset/'verified.json')
    require(verified['status']=='PASS_INDEPENDENT_GENERIC17_REPLAY' and
            verified['prepared_sha256']==sha(dataset/'prepared.json') and
            any(r['sha256']==sha(seed_path) and ROOT/r['path']==seed_path for r in prepared['records']),
            'verified generic17 packet required')
    require(seed['rank_lower_bound']==seed['generic_rank']==len(seed['points'])==17,'native17 required')
    certified_state(seed['curve'],seed['points'],seed['proof'])
    family=next(r for r in read(PARITY)['families'] if r['family']==seed['family'])
    masks=sorted(r['mask'] for r in family['classes'])
    require(len(masks) in (43,49),'complete maximum-class roster required')
    paths=[seed_path,dataset/'prepared.json',dataset/'verified.json',dataset/'protocol.json',PARITY]
    names=['run_fresh6_seed_confirmation.py','prospective_half_lattice_v3.sage',
           'visibility_lattice_fast.py','visibility_lattice_v2.py','v3_warm_engine.py',
           'run_complement_seed_v3.py','future_point_admission.py','memory_rank_certificate.py',
           'half_lattice_pointed_sieve.py','pointed_box_equivalence.py','lean_preconditioned_map_receipts.py',
           'lean_preconditioned_map_worker.py','lean_preconditioned_full_pari_mapping.sage',
           'lean_factor_free_pari_mapping.sage']
    protocol=dict(schema='fresh6-first-M18.v1',family=seed['family'],parameter=seed['parameter'],
        masks=masks,mapping_order=list(ORDER),height=125000,seconds_per_chart=10,
        seconds_per_map=5,map_rss_bytes=1024**3,map_python=shutil.which('sage'),
        map_python_sha256=sha(Path(shutil.which('sage'))),gp_sha256=sha(Path('/usr/bin/gp')),
        max_point_invocations=2*len(masks),prime_bound=1000,initial_rank=17,target_rank=18,
        inputs={str(p.relative_to(ROOT)):sha(p) for p in paths},
        sources={**backend.sources(),**{str((CAS/n).relative_to(ROOT)):sha(CAS/n) for n in names}},
        rule='All exact generic maximum classes, with specialized exact CVP representatives '
             'ordered by descending rounded norm then mask. Bounded preconditioned and factor-free '
             'maps; exact signed-permutation box deduplication including infinity. Stop immediately '
             'after the first standalone certified additional direction. Preserve complete returned '
             'chart witnesses, including any points not processed after that first gain.',
        claim_boundary='Prospective bounded seed confirmation, not an exact rank or record. '
                       'CVP exactness concerns the recorded rounded height metric. Resource-limited '
                       'maps and bounded misses remain unresolved.')
    folder.mkdir(exist_ok=False);checkpoint(folder/'seed.json',seed);checkpoint(folder/'protocol.json',protocol)
    checkpoint(folder/'geometry.json',geometry(seed,masks))
    protocol['inputs'].update({str((folder/n).relative_to(ROOT)):sha(folder/n) for n in ('seed.json','geometry.json')})
    checkpoint(folder/'protocol.json',protocol);guard(protocol)


def execute(folder,replay=False):
    protocol=read(folder/'protocol.json');guard(protocol)
    if not replay and (folder/'terminal.json').exists():return
    seed=read(folder/'seed.json');selected=read(folder/'geometry.json')
    if replay:
        require(geometry(seed,protocol['masks'],reference=True)==selected,'independent geometry replay differs')
    model=tuple(map(F,seed['curve']));points=tuple(tuple(map(F,p)) for p in seed['points'])
    state=certified_state(model,points,seed['proof'])
    admission=FinitePointAdmission(model,points,prime_bound=protocol['prime_bound'])
    attempts=[];charts=[];found=False;proof=seed['proof'];previous=sha(folder/'geometry.json')
    for ci,centre in enumerate(selected['centres']):
        seen=set()
        for policy in ORDER:
            mapping,seal,limited=obtain(folder,ci,policy,model,points,centre,state,protocol,replay=replay)
            if mapping is not None and box_key(mapping['matrix']) in seen:
                limited='SAME_HEIGHT_PRESERVING_BOX'
            attempts.append(dict(centre_index=ci,policy=policy,receipt_sha256=seal,skip_reason=limited))
            if limited:continue
            seen.add(box_key(mapping['matrix']))
            search=PointedQuarticSearch(state=state,centre={'coefficients':centre['representative']},
                                       coordinate_policy=mapping['coordinate_policy'])
            path=folder/f'chart-{len(charts):04d}.json'
            if replay or path.exists():
                chart=read(path)
                require(chart['previous_sha256']==previous and chart['centre_index']==ci and
                        chart['policy']==policy and chart['mapping']==mapping,'chart schedule differs')
                transcript=chart['search'];returned=backend.replay(search,mapping,transcript)
            else:
                transcript,returned=backend.execute(search,mapping,protocol['height'],
                    protocol['seconds_per_chart'],protocol['gp_sha256'])
                require(backend.replay(search,mapping,transcript)==returned,'point replay differs')
                chart=dict(index=len(charts),previous_sha256=previous,centre_index=ci,
                           policy=policy,mapping=mapping,search=transcript)
                checkpoint(path,chart)
            require(transcript['height_bound']==protocol['height'] and
                    transcript['timeout_seconds']==protocol['seconds_per_chart'] and
                    transcript['gp_binary_sha256']==protocol['gp_sha256'],'point budget differs')
            charts.append(chart);previous=sha(path)
            for point in returned:
                admission.consider(point)
                if len(admission.points)>17:
                    proof=checked_rank(model,admission.points,admission.primes,
                                       seed['proof']['no_rational_2_torsion_prime'])
                    require(proof['rank_lower_bound']==18,'first gain must certify18')
                    found=True;break
            if not replay:
                checkpoint(folder/'progress.json',dict(charts=len(charts),rank_lower_bound=len(admission.points)))
                print('FRESH_SEED',seed['id'],len(charts),len(admission.points),flush=True)
            if found:break
        if found:break
    terminal=dict(status='FIRST_M18_CERTIFIED' if found else 'BOUNDED_SEED_ATTEMPT_NO_CERTIFIED_GAIN',
        protocol_sha256=sha(folder/'protocol.json'),charts=len(charts),map_attempts=attempts,
        point_timeouts=sum(c['search']['status']!='bounded_search_complete' for c in charts),
        last_chart_sha256=previous,curve=seed['curve'],points=[list(map(str,p)) for p in admission.points],
        proof=proof,rank_lower_bound=len(admission.points),claim_boundary=protocol['claim_boundary'])
    guard(protocol)
    if replay:
        require(read(folder/'terminal.json')==terminal,'independent seed replay differs')
        certified_state(model,admission.points,proof)
        checkpoint(folder/'verified.json',dict(status='PASS_INDEPENDENT_FIRST_SEED_REPLAY',
            terminal_sha256=sha(folder/'terminal.json'),rank_lower_bound=len(admission.points),charts=len(charts)))
    else:checkpoint(folder/'terminal.json',terminal)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('prepare','search','replay'))
    parser.add_argument('--folder',required=True,type=Path);parser.add_argument('--seed',type=Path)
    args=parser.parse_args();folder=args.folder.resolve()
    with (folder.parent/(folder.name+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if args.mode=='prepare':prepare(args.seed.resolve(),folder)
        else:execute(folder,args.mode=='replay')
