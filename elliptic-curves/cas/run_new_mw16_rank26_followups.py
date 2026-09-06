#!/usr/bin/env python3
"""Run the fixed three MW16 follow-ups with at most two independent workers."""
import argparse,sys
from concurrent.futures import ThreadPoolExecutor,as_completed
import followup_new_mw16_rank26 as follow
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
D=follow.LOCAL/'mw16-new26-followups-v1'
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve three-curve aggregate protocol')
    rows=[]
    for identifier in follow.IDS:
        follow.configure(identifier);follow.prepare();rows.append({'id':identifier,'protocol_path':str((follow.D/'protocol.json').relative_to(follow.ROOT)),'protocol_sha256':cert.hashed(follow.D/'protocol.json')})
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.three-new-mw16-rank26-followups.v1','sources':{str(Path(__file__).resolve().relative_to(follow.ROOT)):cert.hashed(Path(__file__).resolve())},'rows':rows,'maximum_workers':2,'declared_charts':903,'scope':'Exactly three previously certified new26 curves. One fixed301-centre attempt each, unchanged100000 height,10 seconds per chart, own26-point seed and no oracle or new parameter search. Maximum two subprocess workers; all maps precede each curve point search; stop separately at28 pending independent replay. No automatic retries or continuation.'})
def launch():
    p=cert.read(D/'protocol.json')
    if p['sources']!={str(Path(__file__).resolve().relative_to(follow.ROOT)):cert.hashed(Path(__file__).resolve())}:raise ArithmeticError('aggregate runner differs')
    if (D/'ledger.json').exists():raise FileExistsError('preserve follow-up ledger')
    ledger={'status':'RUNNING','rows':[]};checkpoint(D/'ledger.json',ledger)
    def one(row):
        path=follow.ROOT/row['protocol_path'];q=cert.read(path);folder=path.parent
        if cert.hashed(path)!=row['protocol_sha256']:raise ArithmeticError('individual protocol differs')
        stages=[]
        for name,command,seconds in [('maps',[follow.SAGE if hasattr(follow,'SAGE') else '/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',str(follow.CAS/'prepare_new_mw16_rank26_adaptive.sage'),'--id',row['id']],q['geometry_wall_seconds']),('worker',[sys.executable,str(follow.CAS/'followup_new_mw16_rank26.py'),'worker','--id',row['id']],q['worker_wall_seconds']),('replay',[sys.executable,str(follow.CAS/'followup_new_mw16_rank26.py'),'replay','--id',row['id']],600)]:
            s=run(command,limits=Limits(seconds,q['rss_bytes']),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=follow.ROOT);ok=s['outcome']=='completed' and s['returncode']==0;stages.append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});print('MW16 NEW26 FOLLOWUP',row['id'],name,s['outcome'],s['returncode'],flush=True)
            if not ok:return {**row,'status':'FAILED_OR_CENSORED','stages':stages}
        data=cert.read(folder/'result.json');return {**row,'status':'PASS','stages':stages,'rank_lower_bound':data['rank_lower_bound'],'charts':len(data['charts']),'result_sha256':cert.hashed(folder/'result.json')}
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(one,r):i for i,r in enumerate(p['rows'])};collected={}
        for f in as_completed(pending):collected[pending[f]]=f.result();ledger['rows']=[collected[i] for i in sorted(collected)];checkpoint(D/'ledger.json',ledger)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(D/'ledger.json',ledger)
from pathlib import Path
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()
