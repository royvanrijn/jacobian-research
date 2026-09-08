#!/usr/bin/env python3
"""Complete exact equation and experiment certificates after the late64 proofs."""
import argparse,sys,time
from pathlib import Path
import certify_compact_r17_candidates as cert
import full11952_late64_r17_pari_batch as batch
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=batch.ROOT;CAS=batch.CAS;D=batch.D/'post-batch';CONTROL=ROOT/'artifacts/local/elliptic-curves/full11952-late64-controller-v1'
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve late64 finalization')
    paths=[Path(__file__).resolve(),CAS/'certify_full11952_late64_r17_results.py',CAS/'report_full11952_late64_experiment.py',CONTROL/'protocol.json',batch.D/'protocol.json']
    checkpoint(D/'protocol.json',{'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'wait_seconds':24000,'rss_bytes':4294967296,
        'jobs':[{'name':n,'script':s,'args':a,'seconds':1200} for n,s,a in [('certify','certify_full11952_late64_r17_results.py',[]),('proof-replay','certify_full11952_late64_r17_results.py',['--check']),('aggregate-build','report_full11952_late64_experiment.py',[]),('aggregate-check','report_full11952_late64_experiment.py',['--check'])]],
        'scope':'After every fixed point, history, modulo2/3/5 cloud and rational geometry stage succeeds, certify all64 independent point subgroups and recheck rational-isomorphism novelty against the pinned catalogue and all853 prior equations. No search, inventory promotion or retries.'})
def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen finalization sources differ')
    out=D/'ledger.json'
    if out.exists():raise FileExistsError('preserve late64 finalization execution')
    ledger={'status':'WAITING_FOR_POINT_PROOFS','rows':[]};checkpoint(out,ledger);deadline=time.monotonic()+p['wait_seconds']
    while True:
        upstream=cert.read(CONTROL/'ledger.json')
        if upstream['status']=='PASS':break
        if upstream['status']!='RUNNING' or time.monotonic()>deadline:raise ArithmeticError('point pipeline failed/censored or wait deadline reached')
        time.sleep(5)
    ledger['status']='RUNNING';checkpoint(out,ledger)
    for job in p['jobs']:
        path=D/(job['name']+'.supervisor.json')
        if path.exists():raise FileExistsError('preserve finalization stage')
        result=run([sys.executable,str(CAS/job['script']),*job['args']],limits=Limits(job['seconds'],p['rss_bytes']),log_path=D/(job['name']+'.log'),checkpoint_path=path,cwd=ROOT)
        ok=result['outcome']=='completed' and result['returncode']==0
        ledger['rows'].append({'name':job['name'],'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':result});checkpoint(out,ledger)
        print('LATE64 FINALIZATION',job['name'],ledger['rows'][-1]['status'],flush=True)
        if not ok:raise ArithmeticError('fixed finalization failed/censored; no retry')
    summary=cert.read(batch.ART/'full11952_late64_experiment_v1.json');ledger['stronger_odd_prime_bounds']=summary['stronger_odd_prime_bounds'];ledger['status']=summary['status'];checkpoint(out,ledger)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args()
    try:globals()[a.stage]()
    except Exception as exc:
        if a.stage=='launch' and (D/'ledger.json').exists():
            d=cert.read(D/'ledger.json');d.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'ledger.json',d)
        raise
