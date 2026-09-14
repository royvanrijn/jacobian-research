#!/usr/bin/env sage -python
"""Check retained controls, fit on excluded-curve data, and seal the experiment."""
import argparse
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import time

from finite_cancellation_corpus import ROOT, OUT as V1, canonical, digest, pointkey, short, write
from finite_cancellation_validation_audit import OUT


def preflight():
    from memory_rank_certificate import checked_rank
    from importlib.machinery import SourceFileLoader
    sage_check=SourceFileLoader('validation_sage_rank',str(Path(__file__).with_name('verify_finite_cancellation_cpu.sage'))).load_module().sage_rank
    start=time.process_time();inputs=json.loads((OUT/'inputs.json').read_text())
    oracle={o['id']:o for o in json.loads(gzip.decompress((OUT/'oracle.json.gz').read_bytes()))}
    roster=json.loads((OUT/'roster.json').read_text());rows=[]
    assert digest((OUT/'inputs.json').read_bytes())==roster['inputs_sha256']
    assert digest((OUT/'oracle.json.gz').read_bytes())==roster['oracle_sha256']
    for case,source in zip(inputs,roster['chosen']):
        assert case['id']==source['id']
        for key,name in [('source_seed','seed_sha256'),('source_landscape','landscape_sha256')]:
            assert digest((ROOT/source[key]).read_bytes())==source[name]
        seed=case['seed'];curve,points=short(seed['curve'],seed['points']);o=oracle[case['id']]
        assert curve==o['curve'] and {pointkey(p) for p in points}<={pointkey(p) for p in o['points']}
        assert len(seed['points'])<o['endpoint_rank']
        actual=checked_rank(tuple(map(F,seed['curve'])),[tuple(map(F,p)) for p in seed['points']],
            [s['prime'] for s in seed['proof']['signatures']],seed['proof']['no_rational_2_torsion_prime'])
        assert json.loads(json.dumps(actual))==seed['proof']
        first=sage_check(seed)
        path=o['source'];i=path.find('.json');p=ROOT/path[:i+5];packet=json.loads(p.read_text())
        for k in path[i+6:].split('/'):
            if k:packet=packet[int(k)] if isinstance(packet,list) else packet[k]
        proof=packet.get('rank_certificate',packet.get('proof'))
        assert digest(canonical(proof))==o['proof_sha256']
        model=packet.get('curve',packet.get('ainvs',packet.get('model')))
        cm,cp=short(model,packet['points'])
        assert cm==curve and {pointkey(p) for p in cp}=={pointkey(p) for p in o['points']}
        # The source order is preserved for its finite signature certificate.
        original={'curve':model,'points':packet['points'],'proof':proof}
        final=sage_check(original)
        rows.append({'id':case['id'],'stratum':case['stratum'],'seed_rank':first,'endpoint_rank':final,
            'source_packet_file':str(p.relative_to(ROOT)),'source_packet_file_sha256':digest(p.read_bytes()),
            'literal_withheld_count':len(cp)-len(points)})
        if len(rows)%8==0:print(json.dumps({'preflight_cases':len(rows),'cpu_seconds':time.process_time()-start}),flush=True)
    result={'status':'PASS_ALL_SEEDS_AND_LARGER_ENDPOINTS','rows':rows,'cpu_seconds':time.process_time()-start,
        'inputs_sha256':roster['inputs_sha256'],'oracle_sha256':roster['oracle_sha256'],
        'source_sha256':digest(Path(__file__).read_bytes()),
        'boundary':'Portable exact seed signatures and independent Sage quotient arithmetic for every seed and larger endpoint. Initial points form an explicit subset of the independently certified endpoint. No new point search.'}
    write(OUT/'preflight.json',result);print(json.dumps({'status':result['status'],'cases':len(rows),'cpu_seconds':result['cpu_seconds']}),flush=True)


def fit():
    from analyze_finite_cancellation import ridge
    from finite_cancellation_validation_features import FEATURES, model_features
    start=time.process_time();roster=json.loads((OUT/'roster.json').read_text());excluded=set(roster['excluded_training_j_groups'])
    rows=[];records=[];hashes={}
    for path in sorted((V1/'cases').glob('*.json.gz')):
        raw=path.read_bytes();packet=json.loads(gzip.decompress(raw))
        if packet['status']!='PASS_EXACT_ACCESSIBILITY' or packet['j_group'] in excluded or packet['family']=='Curve302-development':continue
        hashes[path.name]=digest(raw);prep=packet['prepared'];fs=[]
        for m in prep['models']:
            n,d,q=([int(v) for v in m[k]] for k in ('N','D','q'))
            f,_=model_features(n,d,q,True);fs.append(f)
        records.append({'source':path.name,'case_id':packet['case_id'],'anchor':prep['anchor_index'],
                        'j_group':packet['j_group'],'features':fs})
        for ti in sorted({o['target_index'] for o in packet['observations']}):
            hs=[min(o['logH'] for o in packet['observations'] if o['target_index']==ti and o['model_index']==mi) for mi in range(len(fs))]
            for mi,f in enumerate(fs):
                rows.append({'group':packet['j_group'],'x':{k:f[k]-fs[0][k] for k in FEATURES['real_q_g']},'y':hs[mi]-hs[0]})
        if len(records)%200==0:print(json.dumps({'training_anchors':len(records),'cpu_seconds':time.process_time()-start}),flush=True)
        if time.process_time()-start>1800:raise TimeoutError('Declared offline preparation cap; no CPU execution.')
    assert rows and not {r['group'] for r in rows}&excluded
    fits={arm:ridge(rows,[rows[0]],names)[1] for arm,names in FEATURES.items()}
    with gzip.GzipFile(str(OUT/'training-features.json.gz'),'wb',mtime=0) as f:f.write(canonical(records))
    result={'status':'FROZEN_EXCLUDED_CURVE_FITS','fits':fits,'training_model_rows':len(rows),'training_anchors':len(records),
            'training_j_groups':sorted({r['group'] for r in rows}),'excluded_training_j_groups':sorted(excluded),
            'ridge_penalty':10,'source_cases_sha256':hashes,'cpu_seconds':time.process_time()-start,
            'training_features_sha256':digest((OUT/'training-features.json.gz').read_bytes()),
            'feature_source_sha256':digest(Path(__file__).with_name('finite_cancellation_validation_features.py').read_bytes()),
            'source_sha256':digest(Path(__file__).read_bytes())}
    write(OUT/'fit.json',result);print(json.dumps({'status':result['status'],'anchors':len(records),'rows':len(rows),'cpu_seconds':result['cpu_seconds']}),flush=True)


def seal():
    import platform
    if (OUT/'protocol.json').exists():raise FileExistsError('Do not retag a frozen protocol.')
    design=json.loads((OUT/'design.json').read_text());roster=json.loads((OUT/'roster.json').read_text())
    pre=json.loads((OUT/'preflight.json').read_text());fitted=json.loads((OUT/'fit.json').read_text())
    assert pre['status']=='PASS_ALL_SEEDS_AND_LARGER_ENDPOINTS'
    assert fitted['status']=='FROZEN_EXCLUDED_CURVE_FITS'
    assert len(pre['rows'])==len(roster['chosen'])
    assert not set(fitted['training_j_groups'])&set(roster['excluded_training_j_groups'])
    from finite_cancellation_validation_features import FEATURES
    assert {k:v['features'] for k,v in fitted['fits'].items()}==FEATURES
    sources=['finite_cancellation_validation_audit.py','finite_cancellation_validation_prepare.py',
        'finite_cancellation_validation_features.py','finite_cancellation_validation_cpu.py',
        'run_finite_cancellation_validation.py','finite_cancellation_features.py','finite_cancellation_corpus.py',
        'lean_factor_free_pari_mapping.sage','pari_pointed_backend.py','pointed_quartic_search.py',
        'half_lattice_pointed_sieve.py','future_point_admission.py','memory_rank_certificate.py',
        'v3_warm_engine.py','search_observability.py','analyze_finite_cancellation.py']
    result={**design,'status':'FROZEN_BEFORE_VALIDATION_CPU','cases':len(roster['chosen']),
        'fits':fitted['fits'],'source_sha256':{n:digest(Path(__file__).with_name(n).read_bytes()) for n in sources},
        'input_sha256':{n:digest((OUT/n).read_bytes()) for n in ['design.json','roster.json','inputs.json','preflight.json','fit.json','training-features.json.gz']},
        'gp_sha256':digest(Path('/usr/bin/gp').read_bytes()),'python':platform.python_version(),
        'order':'All six arm permutations cycle by case ordinal; one process at a time. No outcome-based timing repeats.',
        'cost':'Each outer receipt charges interpreter, all map candidates and discarded features, point backend, admission and producer rank proof. Training and shared input verification are reported separately. Cached retained centre orders are common inputs.',
        'q_feature_definition':'q square status alone determines the common q tree. Soluble q leaves are separately refined for gcd constancy only in the full arm. All q features and both real features agree exactly between fitted arms on the same model.'}
    write(OUT/'protocol.json',result);print(json.dumps({'status':result['status'],'cases':result['cases'],'sha256':digest((OUT/'protocol.json').read_bytes())}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=['preflight','fit','seal']);a=parser.parse_args()
    globals()[a.command]()
