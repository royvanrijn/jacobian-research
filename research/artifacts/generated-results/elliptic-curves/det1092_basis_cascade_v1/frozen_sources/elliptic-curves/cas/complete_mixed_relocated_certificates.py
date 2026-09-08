#!/usr/bin/env python3
"""Finish unchanged search evidence after the live repository directory relocation.

The interrupted ledgers remain immutable. New directories hold certificate-only
replays for the five unfinished rows; no point worker is invoked.
"""
import sys, shutil
from pathlib import Path
import mixed_reduced_parent_exposure as batch
from research_runtime.store import checkpoint
from research_runtime.supervisor import run, Limits

ROOT, CAS, OLD = batch.ROOT, batch.CAS, batch.BATCH
D = ROOT/'artifacts/local/elliptic-curves/mixed-reduced-parent-relocation-replay-v1'

def main():
    p = batch.protocol()
    old = batch.cert.read(OLD/'ledger.json')
    assert old['status'] == 'PASS' and len(old['rows']) == 9
    assert not (D/'protocol.json').exists()
    interrupted = batch.cert.read(OLD/'post-ledger.json')
    assert [r['id'] for r in interrupted['rows']] == [r['id'] for r in p['rows'][:4]]
    # Preserve both the moved running snapshots and the terminal records written
    # to the original directory by processes which were alive during the move.
    originals = [OLD/n for n in ['post-ledger.json','run-ledger.json']]
    leftovers = ROOT.parent/'artifacts/local/elliptic-curves/mixed-reduced-parent-exposure-v1/d1092-fibre-000'
    originals += list(leftovers.glob('*.json'))
    evidence = D/'interruption'; evidence.mkdir(parents=True, exist_ok=False)
    for i, src in enumerate(originals): shutil.copy2(src,evidence/(str(i)+'-'+src.name))
    checkpoint(D/'protocol.json',dict(source_sha256=batch.cert.hashed(Path(__file__)),
        original_protocol_sha256=batch.cert.hashed(OLD/'protocol.json'),
        inputs={str((OLD/r['id']/n).relative_to(ROOT)):batch.cert.hashed(OLD/r['id']/n)
                for r in p['rows'] for n in ['seed.json','maps.json','result.json']},
        scope='Certificate-only continuation after repository relocation to research/. All441 point boxes and9 histories already completed. Preserve original interrupted ledgers; copy the same frozen inputs and replay six certificate stages for the five unfinished rows. No search or selection change.',
        maximum_workers=1,rss_bytes=2147483648,seconds_per_certificate_stage=180,
        seconds_per_curve_driver=1800,interruption_files={x.name:batch.cert.hashed(x) for x in evidence.iterdir()}))
    ledger=dict(status='RUNNING',rows=[]);checkpoint(D/'ledger.json',ledger)
    for i,row in enumerate(p['rows']):
        ident=row['id'];folder=OLD/ident
        if i < 4:
            c=batch.cert.read(folder/'certification-ledger.json');assert c['status']=='PASS'
            ledger['rows'].append(dict(id=ident,status='PASS',reused=True,folder=str(folder.relative_to(ROOT)),rank_lower_bound=c['rank_lower_bound']))
        else:
            folder=D/ident;folder.mkdir(exist_ok=False)
            for n in ['seed.json','maps.json','result.json']:shutil.copy2(OLD/ident/n,folder/n)
            s=run([sys.executable,str(CAS/'certify_factor_free_exposure_v3.py'),
                '--run',str(folder),'--protocol',str(OLD/'protocol.json'),
                '--prefix','mixed_relocated_'+ident.replace('-','_')],
                limits=Limits(1800,2147483648),log_path=folder/'driver.log',
                checkpoint_path=folder/'driver.supervisor.json',cwd=ROOT)
            ok=s['outcome']=='completed' and s['returncode']==0
            ledger['rows'].append(dict(id=ident,status='PASS' if ok else 'FAILED_OR_CENSORED',reused=False,folder=str(folder.relative_to(ROOT)),supervision=s))
            if not ok:
                ledger['status']='FAILED_OR_CENSORED';checkpoint(D/'ledger.json',ledger)
                raise ArithmeticError('preserve certificate continuation failure')
            ledger['rows'][-1]['rank_lower_bound']=batch.cert.read(folder/'certification-ledger.json')['rank_lower_bound']
        checkpoint(D/'ledger.json',ledger);print('CERTIFICATE CONTINUATION',ident,ledger['rows'][-1]['rank_lower_bound'],flush=True)
    assert all(batch.cert.hashed(ROOT/n)==h for n,h in batch.cert.read(D/'protocol.json')['inputs'].items())
    ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)

if __name__=='__main__':main()
