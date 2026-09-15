#!/usr/bin/env sage -python
"""Post-execution exact attribution of actual later gains; no point search."""
from collections import Counter
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
from pathlib import Path
import resource
import time

from cancellation_basis_certified import ROOT, CAS, OUT, RAW, read, sha, need, new_write, guard
from finite_cancellation_corpus import write

POSTMORTEM_CPU_LIMIT = 300


def audit():
    from sage.all import EllipticCurve, QQ
    import numpy as np
    from cancellation_basis_exploration import prepare
    from cancellation_basis_epoch import native_rank
    from pointed_quartic_search import PointedQuarticSearch
    from pari_pointed_backend import validate_map
    from search_observability import point_visibility
    plan, inputs = guard(); cases = {c['id']:c for c in inputs}
    supervision = read(OUT/'supervision.json')
    need(supervision['status']=='COMPLETE' and len(supervision['records'])==48, 'both validation blocks must finish')
    need(not (OUT/'mechanism.json').exists(), 'preserve prior mechanism audit')
    need(POSTMORTEM_CPU_LIMIT == plan['postmortem']['mechanism_cpu_limit'], 'postmortem cap differs')
    resource.setrlimit(resource.RLIMIT_CPU,(POSTMORTEM_CPU_LIMIT,POSTMORTEM_CPU_LIMIT+5))
    start = time.process_time(); rows=[]; counts=Counter(); blocks={str(i):Counter() for i in (0,1)}; map_hashes={}
    mapper = SourceFileLoader('early_attribution_mapper',str(CAS/'lean_factor_free_pari_mapping.sage')).load_module()
    mapper.pari.allocatemem(256000000,silent=True)
    geometry = SourceFileLoader('early_attribution_heights',str(CAS/'prospective_half_lattice_v3.sage')).load_module()
    growth = SourceFileLoader('early_attribution_growth',str(ROOT/'elkies-k3/scripts/rank_growth.py')).load_module()
    for receipt in supervision['records']:
        if receipt['arm']!='basis_refresh':continue
        dest=RAW/'arms'/receipt['case']/receipt['arm'];result=read(dest/'result.json');case=cases[receipt['case']]
        need(receipt['status']=='COMPLETE' and sha(dest/'result.json')==receipt['result_sha256'], 'unsealed case')
        need(sha(dest/'independent-verification.json')==result['independent_verification_sha256'], 'unsealed proof')
        initial=result['initial_rank'];block=blocks[str(case['validation_block'])]
        bank_events={e['epoch']:e for e in read(dest/'events.json') if e['kind']=='bank_ready'}
        curve=tuple(map(F,case['seed']['curve']));E=EllipticCurve(QQ,case['seed']['curve']);cache={}
        def models(q):
            key=tuple(map(str,q))
            if key not in cache:
                prepared=prepare(curve,[tuple(map(F,key))],[1],mapper)['models']
                for item in prepared:
                    search=PointedQuarticSearch(curve=curve,subgroup=[],centre={'point':key},coordinate_policy=item['mapping']['coordinate_policy'])
                    validate_map(search,item['mapping']);item['chart']=search.chart_record()
                cache[key]=prepared
            return cache[key]
        old_models=None
        for record in result['calls']:
            need(sha(dest/record['file'])==record['sha256'],'point call changed')
            call=read(dest/record['file'])
            if call['epoch']==0:continue
            counts['post_refresh_calls']+=1;block['post_refresh_calls']+=1
            uses_new=any(call['word'][initial:])
            counts['post_refresh_calls_using_new_generator']+=uses_new;block['post_refresh_calls_using_new_generator']+=uses_new
            if not call['gain']:continue
            counts['later_directions']+=call['gain'];block['later_directions']+=call['gain']
            if not uses_new:continue
            counts['later_directions_using_new_generator']+=call['gain'];block['later_directions_using_new_generator']+=call['gain']
            before=call['rank_after']-call['gain'];packet=read(dest/call['rank_file'])
            need(sha(dest/call['rank_file'])==call['rank_sha256'],'rank packet changed')
            need(native_rank(packet['packet'])==packet['independent_rank'],'actual gain proof failed')
            bank=read(dest/f'epoch-{call["epoch"]:02d}'/'bank.json');seed=bank['seed']
            need(seed['points']==packet['packet']['points'][:before],'pre-call subgroup differs')
            basis=[E(p) for p in seed['points']];word=call['word'];Q=sum((int(c)*p for c,p in zip(word,basis)),E(0))
            need(list(map(str,Q.xy()))==bank['centres'][call['centre']]['point'],'native anchor sum differs')
            search=PointedQuarticSearch(curve=curve,subgroup=[tuple(map(F,p)) for p in seed['points']],
                centre={'coefficients':word},coordinate_policy=call['mapping']['coordinate_policy'])
            validate_map(search,call['mapping'])
            if old_models is None:
                original=read(dest/'epoch-00/bank.json')
                old_models=[{'anchor_index':i,'point':c['point'],'models':models(c['point'])} for i,c in enumerate(original['centres'])]
            drops=[]
            for index in range(initial,len(word)):
                if not word[index]:continue
                q0=Q-int(word[index])*basis[index];dropped=word[:];dropped[index]=0
                need(q0==sum((int(c)*p for c,p in zip(dropped,basis)),E(0)),'coefficient deletion identity differs')
                drops.append({'generator_index':index,'coefficient':word[index],'generator':seed['points'][index],
                    'word':dropped,'anchor':None if q0.is_zero() else list(map(str,q0.xy())),
                    'models':[] if q0.is_zero() else models(list(map(str,q0.xy())))})
            # All counterfactual maps are fixed above, before evaluating any
            # newly returned point in this call. No target enters prepare().
            for P in packet['packet']['points'][before:]:
                R=E(P);reps=[list(map(str,R.xy())),list(map(str,(-R).xy()))]
                def minimum(bank_rows):
                    values=[]
                    for ai,anchor in enumerate(bank_rows):
                        for mi,item in enumerate(anchor['models']):
                            for si,rep in enumerate(reps):
                                witness=point_visibility(item['chart'],rep)
                                need(witness['status']=='OBSERVABLE_WITHOUT_TRANSCRIPT','independent target equals anchor')
                                h=max(abs(int(v)) for v in witness['coordinate'])
                                values.append((h,ai,mi,si,witness))
                    if not values:return None
                    h,ai,mi,si,witness=min(values,key=lambda v:v[:4])
                    return {'height':str(h),'anchor_index':ai,'model':mi,'sign_index':si,'witness':witness}
                actual=point_visibility(search.chart_record(),P);h=max(abs(int(v)) for v in actual['coordinate'])
                need(h<=plan['height'],'actual gained point outside executed box')
                pair=list(map(int,actual['coordinate']));root=int(actual['square_root_absolute'])
                need(tuple(map(F,P)) in {search.map_hit(*pair,root),search.map_hit(*pair,-root)},'actual square map differs')
                old=minimum(old_models);counterfactuals=[]
                for drop in drops:
                    found=minimum([drop]);counterfactuals.append({k:v for k,v in drop.items() if k!='models'} | {'minimum':found,'models':drop['models']})
                strict=any(c['minimum'] is not None and int(c['minimum']['height'])>plan['height'] for c in counterfactuals)
                old_exclusive=old is not None and int(old['height'])>plan['height']
                block['strict_deletion_witnesses']+=strict;counts['strict_deletion_witnesses']+=strict
                block['old_dictionary_exclusive_witnesses']+=old_exclusive;counts['old_dictionary_exclusive_witnesses']+=old_exclusive
                exact_basis=[tuple(map(F,p)) for p in seed['points']]
                gram,asym=geometry.canonical_height_gram(curve,exact_basis+[tuple(map(F,P))]);g=np.asarray(gram,dtype=float)
                row={'case':case['id'],'validation_block':case['validation_block'],'epoch':call['epoch'],'call':record['file'],
                    'call_sha256':record['sha256'],'bank_sha256':sha(dest/f'epoch-{call["epoch"]:02d}'/'bank.json'),
                    'pre_call_rank':before,'fresh_generator_start':bank_events[call['epoch']]['fresh_generator_start'],
                    'point':P,'anchor':list(map(str,Q.xy())),'anchor_word':word,'actual_chart':search.chart_record(),
                    'actual_witness':actual,'old_minimum':old,'counterfactuals':counterfactuals,
                    'strict_deletion_witness':strict,'old_dictionary_exclusive_witness':old_exclusive,
                    'numerical_geometry':{'gram':[[str(v) for v in r] for r in gram],'asymmetry':str(asym),
                        'before':growth.jsonable_metrics(growth.extension_metrics(g,initial,before)),
                        'after':growth.jsonable_metrics(growth.extension_metrics(g,before,before)),
                        'cascade':growth.jsonable_metrics(growth.cascade_metrics(g,initial,before,before))}}
                rows.append(row);write(RAW/'mechanism-prefix.json',{'status':'RUNNING','rows':rows,'cpu_seconds':time.process_time()-start})
        write(RAW/'mechanism-maps'/f'{case["id"]}.json',{'case':case['id'],'old_models':old_models,'map_cache':list(cache.values())})
        map_hashes[case['id']]=sha(RAW/'mechanism-maps'/f'{case["id"]}.json')
    for b in blocks.values():
        for key in ('strict_deletion_witnesses','later_directions_using_new_generator','old_dictionary_exclusive_witnesses'):b.setdefault(key,0)
    result={'status':'PASS','rows':rows,'counts':dict(counts),'blocks':{k:dict(v) for k,v in blocks.items()},
        'map_sha256':map_hashes,'cpu_seconds':time.process_time()-start,'cpu_limit':POSTMORTEM_CPU_LIMIT,'source_sha256':sha(Path(__file__)),
        'protocol_sha256':sha(OUT/'protocol.json'),'point_search_calls':0,
        'boundary':'After both fixed validation blocks finish, exact actual-gain and coefficient-deletion witnesses on the two representatives P,-P. Every old/deleted model is prepared without target input. These finite minima do not range over the full coset. Numerical Schur/cascade metrics are retrospective; no performance or arithmetic exclusion follows from them. Original point/rank calls remain independently certified inside their full CPU meters; this extra diagnostic CPU is reported separately.'}
    new_write(OUT/'mechanism.json',result)
    print({k:result[k] for k in ('status','counts','blocks','cpu_seconds','point_search_calls')},flush=True)


if __name__=='__main__':audit()
