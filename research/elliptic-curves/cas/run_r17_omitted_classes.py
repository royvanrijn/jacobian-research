#!/usr/bin/env python3
"""Supervise every member of the fixed eight-address omitted-class experiment."""
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2]
D=ROOT/'artifacts/local/elliptic-curves/r17-omitted-generic-classes-v1'


def batch():
    protocol=cert.read(D/'protocol.json')
    if (D/'ledger.json').exists():raise FileExistsError('preserve fixed batch')
    for path,h in protocol['sources'].items():
        if cert.hashed(ROOT/path)!=h:raise ArithmeticError('worker source changed')
    checkpoint(D/'launch.json',{'protocol_sha256':cert.hashed(D/'protocol.json'),'runner_sha256':cert.hashed(Path(__file__).resolve()),'maximum_workers':2})
    ledger={'status':'RUNNING','rows':[]};checkpoint(D/'ledger.json',ledger)
    def worker(entry):
        family,index=entry['family'],entry['index'];folder=D/family/f'candidate-{index:02}'
        r=run(['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',str(ROOT/'elliptic-curves/cas/search_r17_omitted_classes.sage'),'run','--family',family,'--index',str(index)],limits=Limits(90,1073741824),log_path=folder/'worker.log',checkpoint_path=folder/'supervisor.json',cwd=ROOT)
        row={**entry,'supervision':r,'status':'PASS' if r['outcome']=='completed' and r['returncode']==0 else 'FAILED_OR_CENSORED'}
        path=folder/'result.json'
        if path.exists():
            d=cert.read(path);row.update(result_path=str(path.relative_to(ROOT)),result_sha256=cert.hashed(path),charts=len(d['charts']),initial_rank=d['initial_dimension'],rank_lower_bound=d['rank_lower_bound'],search_status=d['status'])
        return row
    roster=sorted(protocol['roster'],key=lambda r:(r['index'],r['family']))
    with ThreadPoolExecutor(max_workers=2) as pool:
        for f in as_completed([pool.submit(worker,e) for e in roster]):
            r=f.result();ledger['rows'].append(r);checkpoint(D/'ledger.json',ledger);print('OMITTED-CLASS ATTEMPT',r['family'],r['index'],r['status'],r.get('rank_lower_bound'),flush=True)
    ledger['status']='COMPLETE_DECLARED_ATTEMPTS';checkpoint(D/'ledger.json',ledger)


if __name__=='__main__':batch()
