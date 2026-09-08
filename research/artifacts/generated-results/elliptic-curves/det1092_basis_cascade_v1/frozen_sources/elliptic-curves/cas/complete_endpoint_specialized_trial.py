#!/usr/bin/env python3
"""Continue only the frozen endpoint roster; retain the first replay timeout."""
import sys,argparse
from pathlib import Path
import endpoint_specialized_parity_trial as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=batch.ROOT;CAS=batch.CAS;D=batch.BATCH/'completion-v1'

def prepare():
    p=batch.protocol();out=D/'protocol.json'
    if out.exists():raise FileExistsError('preserve endpoint completion protocol')
    failed=batch.BATCH/'103b2-infinity/replay.supervisor.json';s=cert.read(failed)
    if s['outcome']!='strict_wall_timeout' or s['timeout_seconds']!=180:raise ArithmeticError('retained original180-second replay timeout required')
    paths=[Path(__file__).resolve(),batch.BATCH/'protocol.json',failed,batch.BATCH/'launch.supervisor.json']
    for row in p['rows']:
        folder=batch.BATCH/row['id'];paths.extend([folder/'maps.json',folder/'maps.supervisor.json'])
        if cert.read(folder/'maps.json')['status']!='COMPLETE_DECLARED_MAPS':raise ArithmeticError('all21 original maps must already be frozen')
        for name in ('worker.supervisor.json','replay.supervisor.json','result.json'):
            if (folder/name).exists():paths.append(folder/name)
    checkpoint(out,{'schema':'elliptic-curves.endpoint-specialized-completion.v1','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'rows':p['rows'],'worker_wall_seconds':300,'replay_wall_seconds':600,'rss_bytes':2147483648,'maximum_workers':1,'gate':'All21 original chart sets were frozen before searching. The first two workers complete12 boxes each; the second returns18966 raw points and its exact history replay exceeds180 seconds. Retain that timeout. Permit one600-second replay per remaining or previously timed-out history. Continue only previously unattempted workers from the original21-curve roster. No point worker is repeated and no result is overwritten.','scope':'No new parameter, parity mask, chart, height or ten-second PARI budget. Every unattempted original worker retains its300-second supervisor cap. A worker failure/censor is recorded and the remaining roster continues; it is not automatically retried. Histories with completed workers receive one600-second exact replay. Failed/censored results do not count as completed point coverage.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen endpoint completion inputs differ')
    batch.protocol();return p

def launch():
    p=protocol();out=D/'ledger.json'
    if out.exists():raise FileExistsError('preserve endpoint completion ledger')
    ledger={'status':'RUNNING','rows':[]};checkpoint(out,ledger)
    for index,row in enumerate(p['rows']):
        folder=batch.BATCH/row['id'];wp=folder/'worker.supervisor.json'
        if wp.exists():worker=cert.read(wp)
        else:worker=run([sys.executable,str(CAS/'endpoint_specialized_parity_trial.py'),'worker','--index',str(index)],limits=Limits(p['worker_wall_seconds'],p['rss_bytes']),log_path=folder/'worker.log',checkpoint_path=wp,cwd=ROOT)
        item={'id':row['id'],'worker_supervision':str(wp.relative_to(ROOT)),'worker_sha256':cert.hashed(wp),'worker_outcome':worker['outcome'],'status':'FAILED_OR_CENSORED'}
        if worker['outcome']=='completed' and worker['returncode']==0:
            old=folder/'replay.supervisor.json';prior=cert.read(old) if old.exists() else None
            if prior and prior['outcome']=='completed' and prior['returncode']==0:rp=old;replay=prior
            else:
                rp=D/row['id']/'replay.supervisor.json'
                if rp.exists():raise FileExistsError('preserve single600-second replay attempt')
                replay=run([sys.executable,str(CAS/'endpoint_specialized_parity_trial.py'),'replay','--index',str(index)],limits=Limits(p['replay_wall_seconds'],p['rss_bytes']),log_path=rp.with_name('replay.log'),checkpoint_path=rp,cwd=ROOT)
            item.update(replay_supervision=str(rp.relative_to(ROOT)),replay_sha256=cert.hashed(rp),replay_outcome=replay['outcome'])
            if replay['outcome']=='completed' and replay['returncode']==0:item['status']='PASS'
        ledger['rows'].append(item);checkpoint(out,ledger);print('ENDPOINT COMPLETION',row['id'],item['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES_OR_CENSORING';checkpoint(out,ledger)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()
