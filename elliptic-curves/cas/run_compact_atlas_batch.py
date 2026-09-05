#!/usr/bin/env python3
"""Run the frozen 24-address balanced compact-atlas pilot with exact post-freeze deduplication."""
import argparse
from concurrent.futures import ThreadPoolExecutor,wait,FIRST_COMPLETED
from fractions import Fraction as F
from pathlib import Path
import time
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/local/elliptic-curves/compact-six-r17-h1024-v2'


def batch(directory):
    if (directory/'ledger.json').exists():raise FileExistsError('preserve the prior batch; any continuation needs its own declared limits')
    protocol=cert.read(directory/'protocol.json');atlas=cert.read(spec.ATLAS)
    families={r['family']:r for r in atlas['families']};populations={f:cert.read(directory/f/'population.json') for f in protocol['families']}
    for pop in populations.values():
        if pop['protocol_hash']!=digest(protocol) or len(pop['finalists'])!=4:raise ArithmeticError('population incomplete')
    for name,h in protocol['sources'].items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('source changed before launch')
    catalogue=cert.read(cert.DATABASE)['curves'];own=[];own_bindings={}
    for name in ('compact_r17_new_curves_v1.json','compact_r17_wide_new_curves_v1.json',
                 'compact_r17_top64_interim_curves_v1.json','compact_r17_largest_gain_curve_v1.json'):
        path=ROOT/'artifacts/generated-results/elliptic-curves'/name;own_bindings[str(path.relative_to(ROOT))]=cert.hashed(path)
        own.extend(cert.read(path)['curves'])
    rows=[];seen=[]
    for index in range(4):
        for f in protocol['families']:
            candidate=populations[f]['finalists'][index];t=F(candidate['parameter']);d=t.denominator
            A=spec.polynomial(families[f]['A_coefficients_low_to_high'],t)*d**8
            B=spec.polynomial(families[f]['B_coefficients_low_to_high'],t)*d**12
            model=(F(0),F(0),F(0),A,B);row={'family':f,'index':index,'parameter':str(t),'status':'PENDING'}
            if 4*A**3+27*B**2==0:row['status']='EXACT_SINGULAR_SPECIALIZATION'
            else:
                matches=[r['id'] for r in catalogue if cert.isomorphic(model,r['ainvs'])]
                own_matches=[r['parameter'] for r in own if cert.isomorphic(model,r['curve'])]
                duplicates=[address for previous,address in seen if cert.isomorphic(model,previous)]
                if matches:row.update(status='SKIPPED_CATALOGUED',matches=matches)
                elif own_matches:row.update(status='SKIPPED_PREVIOUSLY_CERTIFIED',previous_parameters=own_matches)
                elif duplicates:row.update(status='SKIPPED_ROSTER_DUPLICATE',earlier_addresses=duplicates)
                seen.append((model,f+':'+str(t)))
            rows.append(row)
    launch={'schema':'elliptic-curves.compact-atlas-launch.v1','protocol_hash':digest(protocol),'batch_source_sha256':cert.hashed(Path(__file__).resolve()),
        'population_hashes':{f:cert.hashed(directory/f/'population.json') for f in protocol['families']},
        'own_certificate_hashes':own_bindings,'order':'round-robin finalist index, then frozen family order',
        'point_workers':4,'maximum_addresses':24,'submission_wall_seconds':3600,'roster':rows}
    checkpoint(directory/'launch.json',launch)
    ledger={'schema':'elliptic-curves.compact-atlas-ledger.v1','launch_sha256':cert.hashed(directory/'launch.json'),
        'status':'RUNNING','rows':[dict(r) for r in rows]};checkpoint(directory/'ledger.json',ledger)
    started=time.monotonic();pending=[i for i,r in enumerate(rows) if r['status']=='PENDING'];stop=False
    def run(index):
        r=rows[index];f=r['family'];j=r['index'];path=directory/f/f'candidate-{j:02}'/'result.json'
        try:
            result=capture(['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',str(ROOT/'elliptic-curves/cas/search_compact_r17_atlas.sage'),
                'run','--directory',str(directory),'--family',f,'--index',str(j)],
                limits=Limits(wall_seconds=300,rss_bytes=1610612736),log_path=directory/f/f'worker-{j:02}.log')
            data=cert.read(path)
            if data['status'] not in ('COMPLETE_DECLARED_PILOT','TARGET_REACHED_PENDING_INDEPENDENT_REPLAY','SKIPPED_CATALOGUED'):
                raise ArithmeticError('worker returned without a terminal checkpoint')
            return {**r,'status':data['status'],'rank_lower_bound':data.get('rank_lower_bound'),
                'path':str(path.relative_to(ROOT)),'result_sha256':cert.hashed(path),'supervision':result.supervision}
        except Exception as e:
            return {**r,'status':'FAILED_OR_CENSORED','path':str(path.relative_to(ROOT)),'error':str(e)}
    with ThreadPoolExecutor(max_workers=4) as pool:
        active={}
        while pending or active:
            while pending and len(active)<4 and not stop and time.monotonic()-started<3600:
                i=pending.pop(0);active[pool.submit(run,i)]=i
            if not active:break
            ready,_=wait(active,return_when=FIRST_COMPLETED)
            for future in ready:
                i=active.pop(future);r=future.result();ledger['rows'][i]=r
                if (r.get('rank_lower_bound') or 0)>=28:stop=True
                checkpoint(directory/'ledger.json',ledger);print('ATLAS BATCH',r['family'],r['parameter'],r['status'],r.get('rank_lower_bound'),flush=True)
        for i in pending:ledger['rows'][i]['status']='UNRUN_TARGET_STOP' if stop else 'UNRUN_SUBMISSION_LIMIT'
    ledger['status']='TARGET_PENDING_INDEPENDENT_REPLAY' if stop else 'COMPLETE_DECLARED_BATCH' if not pending else 'SUBMISSION_LIMIT'
    checkpoint(directory/'ledger.json',ledger)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--directory',type=Path,default=DEFAULT)
    batch(p.parse_args().directory.resolve())
