#!/usr/bin/env python3
"""Exact terminal-prefix proofs and independent exposure accounting for all arms."""
import argparse
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import nearcut60v2_mw16_pari_batch as batch
import certify_compact_r17_candidates as cert
import replay_corrected60_mw16_geometry as geometry
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import run, Limits
ROOT=batch.ROOT; CAS=batch.CAS; D=batch.D; ART=batch.ART


def paths(row):
    folder=D/row['id']
    stem='nearcut60v2_mw16_'+row['id'].replace('-','_')
    return folder, ART/(stem+'_mod2_v1.json'), ART/(stem+'_modl_v1.json')


def exposure_geometry(index, check=False):
    p=batch.protocol(); row=p['rows'][index]; folder,cloud,odd=paths(row)
    terminal=cert.read(folder/'worker.supervisor.json')
    if terminal['outcome']=='running':raise ArithmeticError('terminal worker required')
    data=cert.read(folder/'result.json'); maps=cert.read(folder/'maps.json')
    exposure=cert.read(folder/'exposure.json'); baseline=cert.read(folder/'baseline.json')
    if any(x['protocol_hash']!=digest(p) for x in (data,maps,exposure,baseline)):
        raise ArithmeticError('exposure protocol differs')
    if any(x['maps_sha256']!=cert.hashed(folder/'maps.json') for x in (data,exposure,baseline)):
        raise ArithmeticError('exposure maps differ')
    if baseline['state']!=data['initial_state'] or baseline['rank']!=16:
        raise ArithmeticError('pre-search independent16 baseline differs')
    rows=exposure['charts']
    if len(rows)>43 or len(rows)<len(data['charts']) or len(rows)>len(data['charts'])+1:
        raise ArithmeticError('exposure/admission prefix differs')
    returned=[]
    for i,r in enumerate(rows):
        if r['index']!=i or r['centre']!=maps['rows'][i]['centre']:
            raise ArithmeticError('attempt roster differs')
        if r['status']=='RETURNED':returned.append(r)
        elif r['status']!='STARTED' or i!=len(rows)-1 or 'search' in r:
            raise ArithmeticError('only terminal in-flight chart may lack transcript')
    for a,b in zip(data['charts'],returned):
        if any(a[k]!=b[k] for k in ('index','centre','search')):
            raise ArithmeticError('exposure/admission transcript differs')
    if len(returned)<len(data['charts']):raise ArithmeticError('admitted chart has no transcript')
    if terminal['outcome']=='completed' and terminal['returncode']==0:
        if data['status']!='COMPLETE_DECLARED_POINT_ATTEMPT' or len(data['charts'])!=43 or len(returned)!=43:
            raise ArithmeticError('successful worker omitted declared exposure')
    proof={**data,'charts':returned,'status':'TERMINAL_RETAINED_EXPOSURE_PREFIX',
           'admission_result_sha256':cert.hashed(folder/'result.json'),
           'exposure_sha256':cert.hashed(folder/'exposure.json'),
           'worker_supervisor_sha256':cert.hashed(folder/'worker.supervisor.json')}
    initial=geometry.tuples(data['generic_points'])
    geometry.geometry(proof,maps,initial,{**p,'maps_path':folder/'maps.json'})
    target=folder/'proof-input.json'
    if check:
        if cert.read(target)!=proof:raise ArithmeticError('proof input changed')
        geometry.cloud_check(proof,returned,cloud,target)
    else:
        if target.exists():raise FileExistsError('preserve terminal proof input')
        checkpoint(target,proof)


def verify(index):
    p=batch.protocol(); row=p['rows'][index]; folder,cloud,odd=paths(row)
    target=folder/'verification.json'
    if target.exists():raise FileExistsError('preserve verification attempt')
    ledger=cert.read(D/'ledger.json')
    if ledger['status']!='COMPLETE_FIXED_BATCH_ATTEMPTS':raise ArithmeticError('all search attempts terminal first')
    binding=ledger['rows'][index]
    if binding['id']!=row['id'] or binding['point_supervision']!=cert.read(folder/'worker.supervisor.json'):raise ArithmeticError('terminal roster or supervision differs')
    steps=[]
    record={'id':row['id'],'status':'RUNNING','steps':steps,'inputs':{}}
    for name in ('result.json','exposure.json','worker.supervisor.json','maps.json','baseline.json'):
        path=folder/name
        if path.exists():record['inputs'][name]=cert.hashed(path)
    checkpoint(target,record)
    if 'result.json' not in record['inputs'] or 'exposure.json' not in record['inputs']:
        record['status']='NO_RETAINED_POINT_RESULT';checkpoint(target,record);return record
    if binding['result_sha256']!=record['inputs']['result.json']:
        raise ArithmeticError('terminal point result changed')
    proof=folder/'proof-input.json'
    commands=[
        ('history',[str(CAS/'nearcut60v2_mw16_pari_batch.py'),'replay','--index',str(index)],300),
        ('geometry',[str(Path(__file__).resolve()),'geometry','--index',str(index)],120),
        ('cloud-build',[str(CAS/'audit_recorded_point_mod2_rank_v3.py'),'--input',str(proof),'--input-sha256','PROOF_HASH','--output',str(cloud),'--prime-bound','997'],120),
        ('cloud-check',[str(CAS/'audit_recorded_point_mod2_rank_v3.py'),'--check',str(cloud)],120),
        ('odd-build',[str(CAS/'audit_retained_cloud_modl.py'),'--input',str(cloud),'--output',str(odd)],300),
        ('odd-check',[str(CAS/'audit_retained_cloud_modl.py'),'--check',str(odd)],300),
        ('provenance',[str(Path(__file__).resolve()),'geometry-check','--index',str(index)],120)]
    for label,args,seconds in commands:
        args=[cert.hashed(proof) if arg=='PROOF_HASH' else arg for arg in args]
        s=run([sys.executable,*args],limits=Limits(seconds,1610612736),cwd=ROOT,
              log_path=folder/(label+'.verification.log'),checkpoint_path=folder/(label+'.verification.supervisor.json'))
        ok=s['outcome']=='completed' and s['returncode']==0
        steps.append({'name':label,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s})
        checkpoint(target,record)
        if not ok:
            record['status']='FAILED_OR_CENSORED';checkpoint(target,record);return record
    lower=max(cert.read(cloud)['rank_lower_bound'],*(a['finite_column_rank'] for a in cert.read(odd)['audits']))
    if lower<16:raise ArithmeticError('cloud failed certified baseline')
    record.update(status='PASS',rank_lower_bound=lower,certified_gain=lower-16,
                  proof_input_sha256=cert.hashed(proof),cloud_path=str(cloud.relative_to(ROOT)),
                  cloud_sha256=cert.hashed(cloud),odd_path=str(odd.relative_to(ROOT)),odd_sha256=cert.hashed(odd))
    checkpoint(target,record);print('VERIFIED STRATA',row['id'],lower,flush=True)
    return record


def all_rows():
    p=batch.protocol(); out=D/'verification-ledger.json'
    if out.exists():raise FileExistsError('preserve all-row verification')
    ledger={'status':'RUNNING','rows':[]};got={};checkpoint(out,ledger)
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(verify,i):i for i in range(len(p['rows']))}
        for f in as_completed(pending):
            i=pending[f]
            try:record=f.result()
            except Exception as exc:
                folder,_,_=paths(p['rows'][i]);v=folder/'verification.json'
                record=cert.read(v) if v.exists() else {'id':p['rows'][i]['id'],'steps':[],'inputs':{}}
                record.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(v,record)
            got[i]=record;ledger['rows']=[got[j] for j in sorted(got)];checkpoint(out,ledger)
    ledger['status']='COMPLETE_ALL_ALLOCATED_VERIFICATIONS';checkpoint(out,ledger)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['all','geometry','geometry-check'])
    parser.add_argument('--index',type=int);a=parser.parse_args()
    if a.stage=='all':all_rows()
    else:exposure_geometry(a.index,check=a.stage=='geometry-check')
