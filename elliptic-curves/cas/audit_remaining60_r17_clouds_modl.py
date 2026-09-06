#!/usr/bin/env python3
"""Bounded mod3/5 audit of every complete cloud in the fixed60 fresh outer cohort."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import certify_compact_r17_candidates as cert
import remaining60_r17_pari_batch as batch
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=batch.ROOT;ART=batch.ART;D=batch.D/'odd-cloud-audit'
def main():
    if (D/'protocol.json').exists():raise FileExistsError('preserve odd-cloud audit')
    p=batch.protocol();v=cert.read(batch.D/'verification-ledger.json')
    if v['status']!='PASS' or len(p['rows'])!=60:raise ArithmeticError('all60 terminal exact histories and mod2 clouds required')
    rows=[]
    for r in p['rows']:
        path=ART/('remaining60_r17_'+r['id'].replace('-','_')+'_mod2_v1.json')
        rows.append({'id':r['id'],'input':str(path.relative_to(ROOT)),'input_sha256':cert.hashed(path),'output':str((ART/('remaining60_r17_'+r['id'].replace('-','_')+'_modl_v1.json')).relative_to(ROOT))})
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.remaining60_r17-odd-cloud-audit.v1','sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in (Path(__file__).resolve(),ROOT/'elliptic-curves/cas/audit_retained_cloud_modl.py',batch.D/'protocol.json',batch.D/'verification-ledger.json')},'rows':rows,'prime_bound':997,'moduli':[3,5],'maximum_workers':2,'seconds_per_stage':300,'rss_bytes_per_worker':1073741824,'gate':'Every point attempt completed and replayed. A finite mod2 dependence need not be rational dependence; check all recorded points modulo3 and5 before interpreting the low admission counts. No new point, trace or selection.','boundaries':'Exact finite-quotient ranks are lower bounds only; matching mod2/3/5 ranks do not prove saturation, exact rank or absence of more directions.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(D/'ledger.json',ledger)
    def one(row):
        ip=ROOT/row['input'];op=ROOT/row['output']
        if cert.hashed(ip)!=row['input_sha256']:raise ArithmeticError('cloud changed')
        for label,args in [('build',['--input',str(ip),'--output',str(op)]),('check',['--check',str(op)])]:
            folder=D/row['id'];s=run(['/usr/bin/python3',str(ROOT/'elliptic-curves/cas/audit_retained_cloud_modl.py'),*args],limits=Limits(300,1073741824),log_path=folder/(label+'.log'),checkpoint_path=folder/(label+'.supervisor.json'),cwd=ROOT)
            if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('odd-cloud check failed/censored')
        d=cert.read(op);return {**row,'status':'PASS','output_sha256':cert.hashed(op),'mod2_lower_bound':d['original_lower_bound'],'odd_lower_bounds':{str(a['modulus']):a['finite_column_rank'] for a in d['audits']}}
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures={pool.submit(one,r):i for i,r in enumerate(rows)};done={}
        for f in as_completed(futures):
            done[futures[f]]=f.result();ledger['rows']=[done[k] for k in sorted(done)];checkpoint(D/'ledger.json',ledger);print('PAIRED CLOUD',done[futures[f]]['id'],done[futures[f]]['mod2_lower_bound'],done[futures[f]]['odd_lower_bounds'],flush=True)
    ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
if __name__=='__main__':main()
