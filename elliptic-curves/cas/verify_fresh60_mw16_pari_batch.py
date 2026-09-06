#!/usr/bin/env python3
"""Bounded exact history and full-cloud checks for terminal fresh outer attempts."""
import argparse,sys
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import fresh60_mw16_pari_batch as batch
import certify_compact_r17_candidates as cert
from research_runtime.supervisor import run,Limits
from research_runtime.store import checkpoint
ROOT=batch.ROOT;ART=batch.ART;D=batch.D;CAS=batch.CAS

def verify(index):
    p=batch.protocol();gate=cert.read(ROOT/'artifacts/local/elliptic-curves/memory-cloud-audit-regression-v3/agreement.json');
    if gate['status']!='PASS':raise ArithmeticError('memory cloud regression incomplete')
    row=p['rows'][index];folder=D/row['id'];path=folder/'result.json';terminal=cert.read(folder/'worker.supervisor.json')
    if terminal['outcome']=='running' or not path.exists():raise ArithmeticError('terminal result required')
    name='fresh60_mw16_'+row['id'].replace('-','_')+'_mod2_v1.json';out=ART/name;steps=[]
    for label,args,seconds in [('history',[str(CAS/'fresh60_mw16_pari_batch.py'),'replay','--index',str(index)],300),('cloud-build',[str(CAS/'audit_recorded_point_mod2_rank_v3.py'),'--input',str(path),'--input-sha256',cert.hashed(path),'--output',str(out),'--prime-bound','997'],120),('cloud-check',[str(CAS/'audit_recorded_point_mod2_rank_v3.py'),'--check',str(out)],120)]:
        record=folder/(label+'.verification.json')
        if record.exists():raise FileExistsError('preserve paired verification')
        old=folder/'replay.supervisor.json'
        if label=='history' and old.exists():
            s=cert.read(old);bound=next(r for r in cert.read(D/'ledger.json')['rows'] if r['id']==row['id'])
            if s['outcome']!='completed' or s['returncode']!=0 or Path(s['command'][1]).resolve()!=CAS/'fresh60_mw16_pari_batch.py' or s['command'][2:]!=['replay','--index',str(index)] or cert.hashed(Path(s['log']))!=s['log_sha256'] or bound['result_sha256']!=cert.hashed(path):raise ArithmeticError('prior exact replay cannot be reused')
            checkpoint(record,s)
        else:s=run([sys.executable,*args],limits=Limits(seconds,1610612736),log_path=folder/(label+'.verification.log'),checkpoint_path=record,cwd=ROOT)
        steps.append({'name':label,'status':'PASS' if s['outcome']=='completed' and s['returncode']==0 else 'FAILED_OR_CENSORED','supervision':s});checkpoint(folder/'verification.json',{'status':'RUNNING','input_sha256':cert.hashed(path),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'steps':steps})
        if steps[-1]['status']!='PASS':raise ArithmeticError('paired verification failed/censored')
    checkpoint(folder/'verification.json',{'status':'PASS','input_sha256':cert.hashed(path),'cloud_certificate':str(out.relative_to(ROOT)),'cloud_sha256':cert.hashed(out),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'steps':steps});print('VERIFIED FRESH60 MW16',row['id'],'rank >=',cert.read(out)['rank_lower_bound'],flush=True)

def all_rows():
    p=batch.protocol();ledger=cert.read(D/'ledger.json')
    if ledger['status']!='COMPLETE_FIXED_BATCH_ATTEMPTS':raise ArithmeticError('terminal batch required')
    out=D/'verification-ledger.json'
    if out.exists():raise FileExistsError('preserve paired verification ledger')
    result={'status':'RUNNING','rows':[]};checkpoint(out,result)
    def one(i):
        row=p['rows'][i];v=D/row['id']/'verification.json'
        if not v.exists():verify(i)
        r=cert.read(v)
        if r['status']!='PASS' or r['verifier_sha256']!=cert.hashed(Path(__file__).resolve()) or r['input_sha256']!=cert.hashed(D/row['id']/'result.json') or r['cloud_sha256']!=cert.hashed(ROOT/r['cloud_certificate']):raise ArithmeticError('retained verification differs')
        for step in r['steps']:
            if step['status']!='PASS' or cert.hashed(Path(step['supervision']['log']))!=step['supervision']['log_sha256']:raise ArithmeticError('retained verification log differs')
        return {'id':row['id'],**r}
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(one,i):i for i in range(len(p['rows']))}
        collected={}
        for f in as_completed(pending):
            collected[pending[f]]=f.result();result['rows']=[collected[i] for i in sorted(collected)];checkpoint(out,result)
    result['status']='PASS';checkpoint(out,result)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int);a=p.parse_args();all_rows() if a.index is None else verify(a.index)
