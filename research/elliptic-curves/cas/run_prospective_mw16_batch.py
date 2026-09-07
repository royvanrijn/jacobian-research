#!/usr/bin/env python3
"""Bounded fixed20-address batch, with no catalogue read or rank-based batch stop."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=Path(__file__).resolve().parents[2]
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/prospective-mw16-h1024-v1'

def batch(directory):
    if (directory/'point-launch.json').exists():raise FileExistsError('preserve the previous batch')
    p=cert.read(directory/'point-protocol.json');rows=[]
    for name,h in p['sources'].items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('worker source changed')
    for i in range(4):
        for f in p['families']:
            path=directory/f/'population.json'
            if cert.hashed(path)!=p['population_hashes'][f]:raise ArithmeticError('population changed')
            rows.append({'family':f,'index':i,'parameter':cert.read(path)['finalists'][i]['parameter'],'status':'PENDING'})
    launch={'schema':'elliptic-curves.prospective-mw16-point-launch.v1','protocol_hash':digest(p),
        'batch_source_sha256':cert.hashed(Path(__file__).resolve()),'maximum_workers':4,'roster':rows,
        'scope':'Exactly20 frozen addresses in round-robin index then family order. No known equation, parameter, rank, point or novelty inputs.'}
    checkpoint(directory/'point-launch.json',launch)
    ledger={'launch_sha256':cert.hashed(directory/'point-launch.json'),'status':'RUNNING','rows':[dict(r) for r in rows]}
    checkpoint(directory/'point-ledger.json',ledger)
    def run(i):
        row=rows[i];f=row['family'];j=row['index'];path=directory/f/f'candidate-{j:02}'/'result.json'
        try:
            r=capture(['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',str(ROOT/'elliptic-curves/cas/search_prospective_mw16_atlas.sage'),
                'run','--directory',str(directory),'--family',f,'--index',str(j)],
                limits=Limits(300,1610612736),log_path=directory/f/f'worker-{j:02}.log')
            d=cert.read(path)
            if d['status'] not in ('COMPLETE_DECLARED_PILOT','INCOMPLETE_GENERIC_MOD2_CERTIFICATE','TARGET_REACHED_PENDING_INDEPENDENT_REPLAY'):raise ArithmeticError('missing terminal checkpoint')
            return {**row,'status':d['status'],'rank_lower_bound':d['rank_lower_bound'],'charts':len(d['charts']),
                    'result_path':str(path.relative_to(ROOT)),'result_sha256':cert.hashed(path),'supervision':r.supervision}
        except Exception as e:
            return {**row,'status':'FAILED_OR_CENSORED','error':str(e),'result_path':str(path.relative_to(ROOT)),
                    'result_sha256':cert.hashed(path) if path.exists() else None}
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(run,i):i for i in range(20)}
        for future in as_completed(futures):
            i=futures[future];ledger['rows'][i]=future.result();checkpoint(directory/'point-ledger.json',ledger)
            r=ledger['rows'][i];print('PROSPECTIVE MW16 BATCH',r['family'],r['parameter'],r['status'],r.get('rank_lower_bound'),flush=True)
    ledger['status']='COMPLETE_FIXED_BATCH_ATTEMPTS';checkpoint(directory/'point-ledger.json',ledger)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--directory',type=Path,default=DIRECTORY);batch(p.parse_args().directory.resolve())
