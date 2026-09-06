#!/usr/bin/env python3
"""Promote only the completed and isolated full11952 point proofs, preserving IDs."""
import argparse,sys,time
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/full11952-inventory-v14-controller-v1';UP=ROOT/'artifacts/local/elliptic-curves/full11952-64-point-portable-controller-v1'
SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve full11952 inventory protocol')
    paths=[Path(__file__).resolve(),CAS/'export_new_high_rank_curve_index_v14.py',CAS/'replay_inventory_v14_memory.py',CAS/'export_full11952_high_rank_models.py',UP/'protocol.json',ART/'new_high_rank_curve_index_v13.json']
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.full11952-inventory-v14-controller.v1','sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in paths},'wait_seconds':22000,'rss_bytes':2147483648,'jobs':[{'name':'inventory-build','script':'export_new_high_rank_curve_index_v14.py','args':[],'seconds':600},{'name':'inventory-replay','script':'replay_inventory_v14_memory.py','args':['--output',str(ART/'new_high_rank_curve_index_v14_memory_replay_v1.json')],'seconds':600},{'name':'models-build','script':'export_full11952_high_rank_models.py','args':[],'seconds':120},{'name':'models-check','script':'export_full11952_high_rank_models.py','args':['--check'],'seconds':120},{'name':'sage-export','script':str(ART/'new_full11952_high_rank_curves.sage'),'args':[],'seconds':60,'executable':SAGE}],'scope':'After all194 isolated full64 point-proof stages pass, extract every pinned-catalogue-unmatched lower bound at least22 into inventoryV14 with existing IDs preserved. Independently replay every inventory point proof, equation comparison and CSV. Export integral models and transported points for unmatched bounds at least27, with global minimality claimed only when the separate cheap certificate passes. This changes no point search or selection and asserts no exact ranks, records or universal novelty. Any failure stops this fixed chain without retry.'})

def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen inventory pipeline changed')
    if (D/'ledger.json').exists():raise FileExistsError('preserve one inventory pipeline')
    ledger={'status':'WAITING_FOR_ISOLATED_POINT_PROOFS','rows':[]};checkpoint(D/'ledger.json',ledger);deadline=time.monotonic()+p['wait_seconds']
    while True:
        q=UP/'controller.supervisor.json'
        if q.exists():
            s=cert.read(q)
            if s['outcome']!='running':
                if s['outcome']!='completed' or s['returncode']!=0 or cert.read(UP/'ledger.json')['status']!='PASS':raise ArithmeticError('isolated point proof failed/censored; no promotion')
                break
        if time.monotonic()>deadline:raise TimeoutError('declared point-proof wait expired')
        time.sleep(5)
    ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
    for job in p['jobs']:
        n=job['name'];script=Path(job['script']);script=script if script.is_absolute() else CAS/script
        s=run([job.get('executable',sys.executable),str(script),*job['args']],limits=Limits(job['seconds'],p['rss_bytes']),log_path=D/(n+'.log'),checkpoint_path=D/(n+'.supervisor.json'),cwd=ROOT)
        ok=s['outcome']=='completed' and s['returncode']==0;ledger['rows'].append({'name':n,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger);print('FULL11952 INVENTORY',n,ledger['rows'][-1]['status'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(D/'ledger.json',ledger);raise ArithmeticError('inventory/model stage failed; no retry')
    ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()
