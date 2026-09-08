#!/usr/bin/env sage-python
"""Checkpointed complement-blind sensitivity calibration and gated replay.

Selection depends on exact generic data or previously recovered blind points.
Each setting is scored by exact recovered quotient rank. Historical public
point lists and per-curve rank labels are not inputs. Prospectives require a
separate control certificate with at least 54 total directions on five curves.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as Q
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
from types import SimpleNamespace
import subprocess
import time

import pointed_quartic_search as backend

ROOT=Path(__file__).resolve().parents[2]
CAS=Path(__file__).resolve().parent
FIXTURE=ROOT/'elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json'
PROSPECTIVE=ROOT/'artifacts/generated-results/elliptic-curves/icarm_mw16_nagao_finalist_specializations_h300_v1.json'
LOCAL=ROOT/'artifacts/local/elliptic-curves/mw16-sensitivity'


def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix('.tmp'); temporary.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n'); temporary.replace(path)


def readpoint(p):
    return Q(p['x']),Q(p['y'])


def pointrecord(p):
    return dict(zip(('x','y'),map(str,p)))


def classify(legacy,model,basis,discoveries):
    return legacy.classify_discovered_group(model=model,basis=basis,discoveries=discoveries,
        relation_chunk_size=64,relation_timeout_seconds=180,stack_bytes=1000000000)


def prepare(parent,seed,mode,centres,legacy,bits):
    model=tuple(map(Q,parent['target_short_model']))
    generic=tuple(map(readpoint,parent['specialized_generic_points']))
    if mode=='adaptive':
        initial=tuple(map(readpoint,seed['current_basis']))
        if initial[:16]!=generic:
            raise ArithmeticError('adaptive subset requires unchanged integral generic prefix')
        extra=tuple(sorted(initial[16:],key=legacy.point_sort_key))[:bits]
        if not extra:
            return model,generic,[],[]
        basis=generic+extra
        legacy.OLD_CLASS_COUNT=len(seed['deepest_masks'])
        args=SimpleNamespace(generic_points=generic,relation_chunk_size=64,
            relation_timeout_seconds=180,stack_bytes=1000000000)
        ranked,ranking=legacy.rank_lifts(model,basis,seed['deepest_masks'],args)
        jobs=[]
        for centre in centres:
            if centre not in ('cvp','residue'):
                raise ValueError('adaptive centre must be cvp or residue')
            for depth,old,word,residue,rep in ranked:
                jobs.append({'centre':centre,'mask':sum(int(v)<<i for i,v in enumerate(residue)),
                    'representative':list(map(int,rep if centre=='cvp' else residue)),
                    'old_mask':int(old),'quotient_word':int(word)})
        return model,basis,jobs,ranking
    jobs=[]
    for centre in centres:
        if centre not in ('specialized','generic'):
            raise ValueError('initial centre must be specialized or generic')
        for row in seed['cover_records']:
            rep=row[centre+'_representative']
            jobs.append({'centre':centre,'mask':int(row['mask']),'representative':rep})
    return model,generic,jobs,{}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=('initial','adaptive','prospective'),default='initial')
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--calibration',type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--curves',default='398,400,401,542,548')
    p.add_argument('--specifications',default='metric:16')
    p.add_argument('--heights',default='100000')
    p.add_argument('--centres',default='generic')
    p.add_argument('--adaptive-bits',type=int,default=5)
    p.add_argument('--seconds',type=float,default=20)
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--maximum-candidates',type=int,default=104)
    p.add_argument('--regression-controls', action='store_true', help='Explicitly use the preserved historical coordinate/search implementation')
    args=p.parse_args()
    global backend
    if args.regression_controls:
        import mw16_sensitivity_backend as backend
    specifications=args.specifications.split(';'); heights=list(map(int,args.heights.split(',')))
    centres=args.centres.split(',')
    if not 1<=args.adaptive_bits<=5 or not 1<=args.workers<=8 or not 0<args.seconds<=60 or any(not 1<=h<=1000000 for h in heights):
        raise ValueError('invalid declared finite campaign limits')
    inputs=[Path(__file__),args.seed,FIXTURE]
    seed=json.loads(args.seed.read_text())
    fixture=json.loads(FIXTURE.read_text())
    for key in ('public_point_lists_loaded','public_complement_coordinates_loaded','target_rank_lower_bounds_loaded'):
        if fixture['blindness_boundary'][key]:
            raise ArithmeticError('not complement-blind')
    selected={}
    for row in sorted(fixture['parents'],key=lambda x:int(x['priority_rank'])):
        selected.setdefault(int(row['curve_id']),row)
    if args.mode=='prospective':
        if not args.calibration:
            raise ValueError('prospective replay requires completed control calibration')
        gate=json.loads(args.calibration.read_text())
        if gate['status']!='COMPLETE' or {r['curve_id'] for r in gate['results']}!={398,400,401,542,548} or sum(r['exact_quotient_rank_recovered'] for r in gate['results'])<54:
            raise ArithmeticError('control sensitivity gate is not met')
        # The prospective settings must be among the settings evaluated in the
        # passing calibration; no uncalibrated replacement is silently selected.
        allowed={(s['specification'],s['height']) for r in gate['results'] for s in r['settings']}
        if any((s,h) not in allowed for s in specifications for h in heights):
            raise ArithmeticError('prospective setting was not calibrated')
        inputs += [args.calibration,PROSPECTIVE]
        candidates=json.loads(PROSPECTIVE.read_text())['candidates'][:args.maximum_candidates]
        parents=[{'parent_id':c['candidate_id'],'curve_id':int(c['parent_curve_id']),
            'target_short_model':c['raw_short_model'],'specialized_generic_points':c['raw_generic_points']}
            for c in candidates]
        lookup={r.get('candidate_id',r.get('parent_id')):r for r in seed['results']}
    else:
        curves=list(map(int,args.curves.split(',')))
        parents=[selected[c] for c in curves]
        lookup={r['parent_id']:r for r in seed['results']}
    hashes={**backend.sources(),**{str(f.resolve().relative_to(ROOT)):sha256(f.read_bytes()).hexdigest() for f in inputs}}
    for f in ('run_curve385_iterated_half_lattice_search.sage','mod2_reduction_independence.py','half_lattice_exact_relations.py'):
        path=CAS/f
        if path.exists(): hashes[str(path.relative_to(ROOT))]=sha256(path.read_bytes()).hexdigest()
    frozen={'mode':args.mode,'specifications':specifications,'heights':heights,'centres':centres,
        'adaptive_bits':args.adaptive_bits,'seconds_per_chart':args.seconds,'workers':args.workers,
        'parent_ids':[r['parent_id'] for r in parents]}
    fresh={'schema':'elliptic-curves.mw16-sensitivity.v1','status':'SEARCHING','inputs':hashes,
        'declared_budget':frozen,'results':[],
        'blindness_boundary':{'public_point_lists_loaded':False,'public_complements_loaded':False,'per_curve_target_ranks_loaded':False},
        'software':{'sage':subprocess.check_output(['sage','--version'],text=True).strip(),
                    'worker':subprocess.check_output([str(backend.base.compiled_worker()),'--version'],text=True).strip()},
        'claim_boundary':['Scores are certified quotient ranks; a bounded miss is not a rank upper bound.',
                          'Curve labels are calibration units; coordinates within a curve are nested trials.',
                          'Prospective replay is gated on at least 54 independently recovered control directions.']}
    payload=json.loads(args.output.read_text()) if args.output.exists() else fresh
    for key in ('schema','inputs','declared_budget'):
        if payload[key]!=fresh[key]: raise ArithmeticError('checkpoint configuration changed: '+key)
    write(args.output,payload)
    legacy=SourceFileLoader('sensitivity_legacy',str(CAS/'run_curve385_iterated_half_lattice_search.sage')).load_module()
    legacy.GENERIC_DIMENSION=16
    for parent in parents:
        prior=lookup[parent['parent_id']]
        result=next((r for r in payload['results'] if r['parent_id']==parent['parent_id']),None)
        if result and result['status']=='COMPLETE': continue
        model,basis,jobs,ranking=prepare(parent,prior,args.mode,centres,legacy,args.adaptive_bits)
        generic=tuple(map(readpoint,parent['specialized_generic_points']))
        if result is None:
            result={'parent_id':parent['parent_id'],'curve_id':parent['curve_id'],'status':'SEARCHING',
                'deepest_masks':prior.get('deepest_masks',prior.get('generic_half_lattice',{}).get('deepest_masks')),
                'search_basis':[pointrecord(p) for p in basis], 'ranking':ranking,
                'settings':[], 'current_basis':[pointrecord(p) for p in generic]}
            payload['results'].append(result)
        discoveries={}
        if args.mode=='adaptive':
            for pnt in map(readpoint,prior['current_basis']):
                discoveries.setdefault(legacy.canonical_point(pnt),set()).add('blind-initial-basis')
        for saved in result['settings']:
            for rec in saved['charts']:
                for pnt in map(readpoint,rec['finite_curve_points']):
                    discoveries.setdefault(legacy.canonical_point(pnt),set()).add(saved['key'])
        for centre in centres:
            for spec in specifications:
                for height in heights:
                    key=f'{centre}|{spec}|{height}'
                    if any(s['key']==key for s in result['settings']): continue
                    setting={'key':key,'centre':centre,'specification':spec,'height':height,'charts':[]}
                    selected_jobs=[j for j in jobs if j['centre']==centre]
                    def work(job):
                        rec=backend.checkpoint(LOCAL/'universal-charts',model=model,points=basis,
                            representative=job['representative'],mask=job['mask'],specification=spec,
                            height=height,seconds=args.seconds)
                        return {**rec,'centre_construction':job}
                    start=time.monotonic()
                    with ThreadPoolExecutor(max_workers=args.workers) as pool:
                        for index,rec in enumerate(pool.map(work,selected_jobs),1):
                            setting['charts'].append(rec)
                            if index%20==0:
                                print(f'SENSITIVITY|parent={parent["parent_id"]}|setting={key}|chart={index}/{len(selected_jobs)}',flush=True)
                    local={}
                    for rec in setting['charts']:
                        for point in map(readpoint,rec['finite_curve_points']):
                            cp=legacy.canonical_point(point)
                            local.setdefault(cp,set()).add(key); discoveries.setdefault(cp,set()).add(key)
                    local_basis,classification=classify(legacy,model,generic,local)
                    setting.update({'classification':classification,
                        'exact_quotient_rank_recovered':len(local_basis)-16,
                        'current_basis':[pointrecord(p) for p in local_basis],
                        'wall_seconds':time.monotonic()-start})
                    result['settings'].append(setting)
                    write(args.output,payload)
                    print(f'SENSITIVITY|parent={parent["parent_id"]}|setting={key}|gain={len(local_basis)-16}|status={classification["status"]}',flush=True)
        final,classification=classify(legacy,model,generic,discoveries)
        result.update({'classification':classification,'current_basis':[pointrecord(p) for p in final],
            'exact_quotient_rank_recovered':len(final)-16,
            'status':'COMPLETE' if classification['status']=='PASS_BASIS_EQUALS_DISCOVERED_GROUP' and all(c['status']=='bounded_search_complete' for s in result['settings'] for c in s['charts']) else 'INCOMPLETE'})
        write(args.output,payload)
        print(f'SENSITIVITY|parent={parent["parent_id"]}|union_gain={len(final)-16}|status={result["status"]}',flush=True)
    payload['status']='COMPLETE' if all(r['status']=='COMPLETE' for r in payload['results']) else 'INCOMPLETE'
    payload['total_control_quotient_rank']=sum(r['exact_quotient_rank_recovered'] for r in payload['results']) if args.mode!='prospective' else None
    payload['positive_candidate_count']=sum(r['exact_quotient_rank_recovered']>0 for r in payload['results']) if args.mode=='prospective' else None
    write(args.output,payload)


if __name__=='__main__': main()
