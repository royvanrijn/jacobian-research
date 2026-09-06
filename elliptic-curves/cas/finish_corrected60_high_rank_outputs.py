#!/usr/bin/env python3
"""Terminal corrected-trial exports and stable inventory promotion; no search."""
import argparse,sys,time
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/corrected60-high-rank-outputs-v1'
UPSTREAM=ROOT/'artifacts/local/elliptic-curves/corrected60-mw16-point-portable-controller-v1'
SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'


def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve corrected output protocol')
    names=['finish_corrected60_high_rank_outputs.py','certify_corrected60_high_rank_minimal.py',
           'export_corrected60_high_rank_sage.py','export_new_high_rank_curve_index_v19.py',
           'replay_inventory_v19_memory.py','certify_discarded_rank26_minimal.py']
    paths=[*(CAS/n for n in names),ART/'new_high_rank_curve_index_v18.json',ART/'corrected60_mw16_results_v1.json',UPSTREAM/'protocol.json']
    checkpoint(D/'protocol.json',{'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'wait_seconds':12000,'rss_bytes':2147483648,
        'stages':[
            ['minimal','certify_corrected60_high_rank_minimal.py',[],300],
            ['minimal-check','certify_corrected60_high_rank_minimal.py',['--check'],300],
            ['export','export_corrected60_high_rank_sage.py',[],300],
            ['export-check','export_corrected60_high_rank_sage.py',['--check'],300],
            ['sage',str(ART/'new_corrected60_high_rank_curves.sage'),[],180],
            ['inventory','export_new_high_rank_curve_index_v19.py',[],300],
            ['inventory-replay','replay_inventory_v19_memory.py',['--output',str(ART/'new_high_rank_curve_index_v19_memory_replay_v1.json')],1800]],
        'scope':'After all182 corrected standalone checks pass, certify global minimal integral models and transported independent points for its unmatched>=22 discoveries, export executable Sage equations, extend the198-curve stable-ID inventory and replay every point certificate and CSV entry. No point search, score or selection changes; no exact rank, conductor or universal novelty claim.'})


def launch():
    p=cert.read(D/'protocol.json')
    if (D/'ledger.json').exists():raise FileExistsError('preserve outputs attempt')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('output source changed')
    ledger={'status':'WAITING_FOR_CORRECTED_STANDALONE_PROOFS','rows':[]};checkpoint(D/'ledger.json',ledger);deadline=time.monotonic()+p['wait_seconds']
    try:
        while True:
            state=cert.read(UPSTREAM/'ledger.json')['status']
            if state=='PASS':break
            if state not in ('RUNNING','WAITING_FOR_POINT_FINALIZATION') or time.monotonic()>deadline:raise ArithmeticError('portable proof gate failed or wait elapsed')
            time.sleep(5)
        if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('output source changed during wait')
        ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
        for name,script,args,seconds in p['stages']:
            command=[SAGE if name=='sage' else sys.executable,str(Path(script) if Path(script).is_absolute() else CAS/script),*args]
            s=run(command,limits=Limits(seconds,p['rss_bytes']),cwd=ROOT,log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'))
            ok=s['outcome']=='completed' and s['returncode']==0
            ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger)
            print('CORRECTED HIGH-RANK OUTPUTS',name,ledger['rows'][-1]['status'],flush=True)
            if not ok:raise ArithmeticError('output verification failed/censored; no retry')
        ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
    except Exception as exc:
        ledger.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'ledger.json',ledger);raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()
