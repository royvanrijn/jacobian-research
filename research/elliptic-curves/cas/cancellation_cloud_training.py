#!/usr/bin/env sage -python
"""Fit cloud multiplicity and radii from all 545 retained development calls.

No new point searches. Full clouds are independently rank-certified, including
every additional witness that the historical first-gain worker did not admit.
"""
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import time

from cancellation_scheduler_prepare import coordinate_height,new_write,V2
from cancellation_cloud_policy import bucket
from finite_cancellation_corpus import ROOT, LOCAL, digest, write

OUT=ROOT/'artifacts/generated-results/elliptic-curves/cancellation_cloud_scheduler_v1'
RAW=LOCAL/'cancellation-cloud-scheduler-v1'
CAS=Path(__file__).resolve().parent


def train():
    from future_point_admission import FinitePointAdmission
    from memory_rank_certificate import checked_rank
    from v3_warm_engine import certified_state
    from pointed_quartic_search import PointedQuarticSearch
    from pari_pointed_backend import replay
    from half_lattice_pointed_sieve import linear_combination
    sage_rank=SourceFileLoader('cloud_training_independent',str(CAS/'verify_finite_cancellation_cpu.sage')).load_module().sage_rank
    start=time.process_time();base=json.loads((OUT/'fit.json').read_text())
    test={c['j_group'] for c in json.loads((OUT/'inputs.json').read_text())}
    development=json.loads((V2/'inputs.json').read_text());lookup={c['id']:c for c in development}
    if test&{c['j_group'] for c in development}:raise ArithmeticError('development/test overlap')
    records=[];radial=[];exposures=[0,0,0];gains=[0,0,0];cases={};packets=[]
    for source in base['development_rows']:
        case=lookup[source['case']];seed=case['seed'];path=ROOT/source['source']
        if digest(path.read_bytes())!=source['sha256']:raise ArithmeticError('development bytes changed')
        row=json.loads(path.read_text());curve=tuple(map(F,seed['curve']));basis=[tuple(map(F,p)) for p in seed['points']]
        if case['id'] not in cases:
            cases[case['id']]=(certified_state(curve,basis,seed['proof']),FinitePointAdmission(curve,basis,prime_bound=500))
        state,admission=cases[case['id']];before=len(admission.points)
        if before!=len(basis):raise ArithmeticError('development is not a first-gaining-call sequence')
        search=PointedQuarticSearch(state=state,centre={'coefficients':row['centre']},coordinate_policy=row['mapping']['coordinate_policy'])
        points=replay(search,row['mapping'],row['search']);new=[]
        for point in points:
            if admission.consider(point)['status']=='INDEPENDENT_FINITE_COLUMN':new.append(list(map(str,point)))
        if bool(new)!=bool(row.get('new_point')):raise ArithmeticError('original development outcome differs')
        ci=source['index'];b=bucket(ci);exposures[b]+=1;gains[b]+=len(new);heights=[]
        entry={'case':case['id'],'index':ci,'source':source['source'],'sha256':source['sha256'],'directions':len(new),'heights':heights}
        if new:
            anchor=linear_combination(curve,basis,row['centre'])
            heights.extend(coordinate_height(p,anchor,row['mapping']['matrix']) for p in new)
            if any(h>125000 or h<1 for h in heights):raise ArithmeticError('cloud radius outside actual box')
            radial.extend(heights)
            proof=checked_rank(curve,admission.points,admission.primes,seed['proof']['no_rational_2_torsion_prime'])
            packet={'curve':seed['curve'],'points':[list(map(str,p)) for p in admission.points],'proof':proof}
            independent=sage_rank(packet)
            if independent['rank']!=before+len(new):raise ArithmeticError('cloud rank failed')
            saved=RAW/'training'/case['id']/'cloud.json';new_write(saved,packet)
            entry.update(packet=str(saved.relative_to(ROOT)),packet_sha256=digest(saved.read_bytes()),independent=independent)
            packets.append(entry['packet'])
        records.append(entry)
        if len(records)%100==0:print(json.dumps({'development_calls':len(records),'cloud_directions':len(radial),'cpu':time.process_time()-start}),flush=True)
    fit={k:base[k] for k in ('local_weight','anchor_priors','backend_cpu_at_125000','per_call_overhead_cpu','independent_certificate_cpu_prior')}
    fit.update(radial_heights=sorted(radial),anchor_direction_rates=[(s+1)/(n+10) for s,n in zip(gains,exposures)],cloud_prior_exposure=4.)
    new_write(OUT/'cloud-fit.json',{'status':'PASS_RECONCILED_DEVELOPMENT_CLOUDS','fit':fit,'records':records,
        'calls_by_bucket':exposures,'directions_by_bucket':gains,'directions':len(radial),'positive_clouds':len(packets),
        'development_j_groups':sorted(c['j_group'] for c in development),'base_fit_sha256':digest((OUT/'fit.json').read_bytes()),
        'cpu_seconds':time.process_time()-start,
        'boundary':'Full returned clouds on retained V2 factor-free development calls, all reclassified at fixed primes through500 and every enlarged group independently Sage-certified. Radii describe a chosen independent subset, not every rational point. Counts remain dependent across anchors and bounded misses are censored. Gamma/Poisson intensity is a scheduling model, not a theorem.'})
    print(json.dumps({'status':'PASS_RECONCILED_DEVELOPMENT_CLOUDS','directions':len(radial),'positive_clouds':len(packets),'rates':fit['anchor_direction_rates'],'cpu':time.process_time()-start}),flush=True)


def verify():
    """Replay ranks and exact homogeneous radii, using the other formula."""
    from math import gcd,lcm
    from half_lattice_pointed_sieve import linear_combination
    checker=SourceFileLoader('cloud_training_replay',str(CAS/'verify_finite_cancellation_cpu.sage')).load_module().sage_rank
    fit=json.loads((OUT/'cloud-fit.json').read_text());radial=[];counts=[0,0,0];exposures=[0,0,0];start=time.process_time()
    cases={c['id']:c for c in json.loads((V2/'inputs.json').read_text())}
    for row in fit['records']:
        path=ROOT/row['source']
        if digest(path.read_bytes())!=row['sha256']:raise ArithmeticError('call binding differs')
        exposures[bucket(row['index'])]+=1;counts[bucket(row['index'])]+=row['directions']
        if not row['directions']:
            if json.loads(path.read_text()).get('new_point'):raise ArithmeticError('positive censored')
            continue
        packetpath=ROOT/row['packet'];packet=json.loads(packetpath.read_text());seed=cases[row['case']]['seed']
        if digest(packetpath.read_bytes())!=row['packet_sha256'] or packet['points'][:len(seed['points'])]!=seed['points']:raise ArithmeticError('cloud packet prefix differs')
        if checker(packet)!=row['independent'] or len(packet['points'])-len(seed['points'])!=row['directions']:raise ArithmeticError('cloud proof differs')
        source=json.loads(path.read_text());curve=tuple(map(F,seed['curve']));basis=[tuple(map(F,p)) for p in seed['points']]
        a,b=linear_combination(curve,basis,source['centre']);A,B,C,D=map(F,source['mapping']['matrix']);heights=[]
        returned={(F(p['x']),F(p['y'])) for p in source['search']['finite_curve_points']}
        for point in packet['points'][len(basis):]:
            x,y=map(F,point)
            if (x,y) not in returned:raise ArithmeticError('cloud point missing')
            hs=[]
            for sign in (1,-1):
                u=D*(sign*y+b)-B*(x-a);v=-C*(sign*y+b)+A*(x-a)
                scale=lcm(u.denominator,v.denominator);m,n=int(u*scale),int(v*scale)
                hs.append(max(abs(m),abs(n))//gcd(m,n))
            heights.append(min(hs))
        if heights!=row['heights']:raise ArithmeticError('independent homogeneous heights differ')
        radial.extend(heights)
    expected=[(g+1)/(n+10) for g,n in zip(counts,exposures)]
    if sorted(radial)!=fit['fit']['radial_heights'] or expected!=fit['fit']['anchor_direction_rates']:raise ArithmeticError('cloud fit differs')
    new_write(OUT/'cloud-training-replay.json',{'status':'PASS','directions':len(radial),'calls':sum(exposures),
        'cloud_fit_sha256':digest((OUT/'cloud-fit.json').read_bytes()),'cpu_seconds':time.process_time()-start,
        'boundary':'Independent Sage ranks and exact homogeneous coordinate calculations; statistical model remains a surrogate.'})
    print(json.dumps({'status':'PASS_CLOUD_TRAINING_REPLAY','directions':len(radial)}),flush=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['train','verify']);a=p.parse_args();globals()[a.command]()
