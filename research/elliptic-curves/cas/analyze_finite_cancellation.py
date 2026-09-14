#!/usr/bin/env sage -python
"""Frozen grouped known-point tests; no map selection from target coordinates."""
from collections import Counter, defaultdict
import csv
import gzip
import json
from math import log2
from pathlib import Path
import time
import numpy as np

from finite_cancellation_corpus import OUT, canonical, digest, write

REAL=['real_logS','coefficient_bits']
FINITE=REAL+['local_g_uniform','local_g_differential','log_soluble_mass','unknown_mass_sum','missing_local_expectations']

def summary(z):
    z=np.array(list(z),float)
    if len(z)==0:return {'count':0}
    return dict(count=len(z),mean=float(z.mean()),median=float(np.median(z)),
        q10=float(np.quantile(z,.1)),q90=float(np.quantile(z,.9)),minimum=float(z.min()),maximum=float(z.max()))
def clustered(values,ratio=False):
    a=np.array([np.mean(v) for _,v in sorted(values.items())],float)
    if not len(a):return {'groups':0}
    rng=np.random.default_rng(20260914)
    bs=np.array([rng.choice(a,len(a),replace=True).mean() for _ in range(2000)])
    f=(lambda x:float(2**x)) if ratio else float
    return {'groups':len(a),'estimate':f(a.mean()),'ci95':[f(np.quantile(bs,.025)),f(np.quantile(bs,.975))],
        'scope':'Equal weight per j group;2000 deterministic cluster bootstrap replicates'}
def ridge(train,test,names):
    x=np.array([[r['x'][k] for k in names] for r in train]);y=np.array([r['y'] for r in train])
    z=np.array([[r['x'][k] for k in names] for r in test]);mu=x.mean(axis=0);sd=x.std(axis=0);sd[sd<1e-12]=1
    x=np.column_stack((np.ones(len(x)),(x-mu)/sd));z=np.column_stack((np.ones(len(z)),(z-mu)/sd))
    penalty=np.eye(x.shape[1])*10;penalty[0,0]=0
    w=np.linalg.solve(x.T@x+penalty,x.T@y)
    return z@w,{'features':names,'coefficients':w.tolist(),'means':mu.tolist(),'scales':sd.tolist()}

def analyze():
    start=time.process_time();packets=[];gaps=[];hashes={}
    for p in sorted((OUT/'cases').glob('*.json.gz')):
        b=p.read_bytes();hashes[p.name]=digest(b);d=json.loads(gzip.decompress(b))
        if d['status']!='PASS_EXACT_ACCESSIBILITY':gaps.append({k:v for k,v in d.items() if k!='source_sha256'});continue
        assert digest(canonical(d['prepared']))==d['prepared_sha256']
        packets.append(d)
    pairs=[];units=[];residue=defaultdict(lambda:defaultdict(list));logs=defaultdict(list)
    byprime=defaultdict(lambda:defaultdict(list));unknown=Counter();model_counts=Counter()
    for packet in packets:
        prep=packet['prepared'];models=prep['models'];group=packet['j_group'];family=packet['family']
        model_counts[len(models)]+=1
        development=family=='Curve302-development'
        obs=packet['observations'];base=models[0]['features']
        for r in obs:
            M=models[r['model_index']]
            logs['log2_g'].append(r['logg']);logs['unprocessed_g_bits'].append(r['cofactor_log2'])
            if r['model_index']==0:
                logs['base_log2_H'].append(r['logH']);logs['base_log2_g'].append(r['logg'])
                logs['panel_fraction_of_logg'].append((r['logg']-r['cofactor_log2'])/r['logg'] if r['logg'] else 1.)
            for location,tree in zip(r['local'],M['local']):
                p=tree['prime'];leaf=tree['leaves'][location['leaf']];unknown['locations']+=1
                if leaf['status']!='square':unknown['censored_locations']+=1;continue
                good=[x for x in tree['leaves'] if x['status']=='square']
                def percentile(differential):
                    from fractions import Fraction as F
                    weights=[float(F(x['mass']))*(p**(x['vq']/2) if differential else 1) for x in good]
                    return sum(w*(int(x['vg']<location['vg'])+.5*int(x['vg']==location['vg'])) for w,x in zip(weights,good))/sum(weights)
                if not development and r['model_index']==0:
                    for diff in [False,True]:
                        label='differential' if diff else 'uniform_soluble'
                        value=percentile(diff)
                        residue[label][group].append(value)
                        byprime[p][label].append(value)
                    byprime[p]['vg'].append(location['vg'])
        for ti in sorted({r['target_index'] for r in obs}):
            selected=[min((r for r in obs if r['target_index']==ti and r['model_index']==mi),key=lambda r:r['logH']) for mi in range(len(models))]
            ident=f'{packet["case_id"]}-{prep["anchor_index"]}-{ti}'
            unit={'id':ident,'group':group,'family':family,'development':development,
                  'logH':[r['logH'] for r in selected],'models':[m['name'] for m in models],
                  'selected':{'factor_free':0},'f':[m['features'] for m in models]}
            unit['selected']['oracle_literal']=min(range(len(models)),key=lambda i:selected[i]['logH'])
            unit['selected']['real_only']=max(range(len(models)),key=lambda i:models[i]['features']['real_logS'])
            for label,feat in [('uniform_local_g_over_S','local_g_uniform'),('differential_local_g_over_S','local_g_differential')]:
                unit['selected'][label]=min(range(len(models)),key=lambda i:models[i]['features'][feat]-models[i]['features']['real_logS'])
            units.append(unit)
            for mi,M in enumerate(models):
                pairs.append({'unit':ident,'group':group,'family':family,'development':development,'model':mi,
                    'x':{k:M['features'][k]-base[k] for k in FINITE},'y':selected[mi]['logH']-selected[0]['logH']})
    eligible=[r for r in pairs if not r['development']];fit_records=[]
    for scheme in ['five_j_folds','leave_family_out']:
        group_key=lambda r: int(r['group'][:8],16)%5 if scheme=='five_j_folds' else r['family']
        predictions={name:{} for name in ['ridge_real','ridge_real_and_finite']}
        for fold in sorted({group_key(r) for r in eligible},key=str):
            test=[r for r in eligible if group_key(r)==fold];excluded={r['group'] for r in test}
            train=[r for r in eligible if r['group'] not in excluded]
            if not train:continue
            for name,features in [('ridge_real',REAL),('ridge_real_and_finite',FINITE)]:
                predicted,record=ridge(train,test,features)
                fit_records.append({'scheme':scheme,'fold':str(fold),'predictor':name,'train':len(train),'test':len(test),**record})
                for r,pred in zip(test,predicted):predictions[name][(r['unit'],r['model'])]=float(pred)
        for unit in units:
            if unit['development']:continue
            for name,values in predictions.items():
                unit['selected'][name+('' if scheme=='five_j_folds' else '_family')]=min(range(len(unit['models'])),key=lambda i:values[(unit['id'],i)])
    outcomes={};families={};comparisons={}
    for name in sorted({k for u in units if not u['development'] for k in u['selected']}):
        diff=defaultdict(list);hits=defaultdict(list);same=0;cheap0=cheap1=0;byfamily=defaultdict(list)
        for u in units:
            if u['development']:continue
            i=u['selected'][name];delta=2*(u['logH'][i]-u['logH'][0]);diff[u['group']].append(delta)
            oracle=min(u['logH']);hits[u['group']].append(float(abs(u['logH'][i]-oracle)<1e-10))
            same+=i==0;cheap0+=u['logH'][0]<=log2(125000);cheap1+=u['logH'][i]<=log2(125000)
            byfamily[u['family']].append(delta)
        outcomes[name]={'H_squared_ratio':clustered(diff,True),'oracle_model_match_fraction':clustered(hits),
            'baseline_selected':same,'literal_directions_under_H125000':cheap1,'baseline_under_H125000':cheap0}
        families[name]={f:{'units':len(v),'H_squared_ratio':2**float(np.mean(v))} for f,v in byfamily.items()}
    for name,ref in [('uniform_local_g_over_S','real_only'),('differential_local_g_over_S','real_only'),
                     ('ridge_real_and_finite','ridge_real'),('ridge_real_and_finite_family','ridge_real_family')]:
        diff=defaultdict(list)
        for u in units:
            if not u['development']:diff[u['group']].append(2*(u['logH'][u['selected'][name]]-u['logH'][u['selected'][ref]]))
        comparisons[name+'_vs_'+ref]=clustered(diff,True)
    passed=[]
    for name,ref in [('uniform_local_g_over_S','real_only'),('differential_local_g_over_S','real_only'),('ridge_real_and_finite','ridge_real')]:
        a=outcomes[name]['H_squared_ratio'];b=comparisons[name+'_vs_'+ref]
        if a['estimate']<=.9 and b['estimate']<=.9 and a['ci95'][1]<1 and b['ci95'][1]<1:passed.append(name)
    rows=[{'unit':u['id'],'group':u['group'],'family':u['family'],'models':u['models'],
           'logH':u['logH'],'selected':u['selected'],'development':u['development']} for u in units]
    with gzip.GzipFile(str(OUT/'evaluation.json.gz'),'wb',mtime=0) as f:f.write(canonical(rows))
    result={'status':'ACCESSIBILITY_GATE_PASSED' if passed else 'NO_FROZEN_FINITE_PREDICTOR_PASSED_ACCESSIBILITY_GATE',
        'passed_predictors':passed,'anchors':len(packets),'censored_anchors':len(gaps),
        'j_groups':len({x['j_group'] for x in packets}),'literal_direction_anchor_units':len(units),
        'models_per_anchor':dict(model_counts),'exact_signed_model_evaluations':sum(len(p['observations']) for p in packets),
        'distributions':{k:summary(v) for k,v in logs.items()},
        'residue_percentile_null_half':{k:clustered(v) for k,v in residue.items()},
        'by_prime':{str(p):{k:summary(v) for k,v in values.items()} for p,values in byprime.items()},
        'local_censoring':dict(unknown),'outcomes':outcomes,'finite_ablation_comparisons':comparisons,
        'family_sensitivity':families,'gaps':gaps,'fit_records':fit_records,
        'preparation_cpu_seconds':sum(p['prepared']['preparation_cpu_seconds'] for p in packets),
        'baseline_map_cpu_seconds':sum(p['prepared']['base_map_cpu_seconds'] for p in packets),
        'case_files_sha256':hashes,'analysis_source_sha256':digest(Path(__file__).read_bytes()),
        'protocol_sha256':digest((OUT/'protocol.json').read_bytes()),'analysis_cpu_seconds':time.process_time()-start,
        'point_search_calls':0,
        'boundary':'Conditional literal-point accessibility, both signs; no claim about all translated representatives of a direction. H squared is not CPU. Small-prime panel only; unprocessed g cofactor retained exactly. Historic selection effects persist despite j/family holdouts.'}
    write(OUT/'analysis.json',result)
    print(json.dumps({k:result[k] for k in ['status','passed_predictors','anchors','j_groups','exact_signed_model_evaluations','residue_percentile_null_half','finite_ablation_comparisons','preparation_cpu_seconds']}))

if __name__=='__main__':analyze()
