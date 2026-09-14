#!/usr/bin/env sage -python
"""Post-primary attribution diagnostic; never changes the frozen CPU selector.

Separate cancellation expectations from q-solubility/censoring features. Keep
the primary penalty, folds, labels and j-cluster weighting unchanged.
"""
from collections import defaultdict
import gzip
import json
from pathlib import Path

from analyze_finite_cancellation import REAL, FINITE, ridge, clustered
from finite_cancellation_corpus import OUT, canonical, digest, write

def main():
    names=REAL+['log_soluble_mass','unknown_mass_sum','missing_local_expectations']
    pairs=[];units={r['unit']:r for r in json.loads(gzip.decompress((OUT/'evaluation.json.gz').read_bytes())) if not r['development']}
    for p in sorted((OUT/'cases').glob('*.json.gz')):
        d=json.loads(gzip.decompress(p.read_bytes()))
        if d['status']!='PASS_EXACT_ACCESSIBILITY' or d['family']=='Curve302-development':continue
        prep=d['prepared'];base=prep['models'][0]['features']
        for ti in {r['target_index'] for r in d['observations']}:
            uid=f'{d["case_id"]}-{prep["anchor_index"]}-{ti}';u=units[uid]
            for mi,m in enumerate(prep['models']):pairs.append({'unit':uid,'group':d['j_group'],'family':d['family'],'model':mi,
                'x':{k:m['features'][k]-base[k] for k in FINITE},'y':u['logH'][mi]-u['logH'][0]})
    results={}
    for scheme in ['five_j_folds','leave_family_out']:
        key=lambda r:int(r['group'][:8],16)%5 if scheme=='five_j_folds' else r['family']
        preds={}
        for fold in sorted({key(r) for r in pairs},key=str):
            test=[r for r in pairs if key(r)==fold];excluded={r['group'] for r in test}
            train=[r for r in pairs if r['group'] not in excluded];values,_=ridge(train,test,names)
            preds.update({(r['unit'],r['model']):float(v) for r,v in zip(test,values)})
        vsbase=defaultdict(list);fullvsq=defaultdict(list);same=0
        full='ridge_real_and_finite'+('_family' if scheme=='leave_family_out' else '')
        for uid,u in units.items():
            i=min(range(len(u['models'])),key=lambda j:preds[uid,j]);j=u['selected'][full]
            vsbase[u['group']].append(2*(u['logH'][i]-u['logH'][0]))
            fullvsq[u['group']].append(2*(u['logH'][j]-u['logH'][i]));same+=i==j
        results[scheme]={'real_plus_q_vs_factor_free':clustered(vsbase,True),
            'full_finite_vs_real_plus_q':clustered(fullvsq,True),'identical_model_choices':same,'units':len(units)}
    write(OUT/'q_ablation.json',{'status':'POST_PRIMARY_ATTRIBUTION_ONLY','features':names,'results':results,
        'primary_analysis_sha256':digest((OUT/'analysis.json').read_bytes()),
        'source_sha256':digest(Path(__file__).read_bytes()),
        'boundary':'Added after the primary outcome to distinguish q information from explicit gcd-expectation information. No hyperparameter search, CPU-policy changes, new gate or new search. These are exploratory attribution results.'})
    print(json.dumps(results))

if __name__=='__main__':main()
