#!/usr/bin/env python3
"""Isolated history, point-cloud, geometry and small-prime proof replay."""
from pathlib import Path
import sys,zipfile
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/compact192-followup-supplement-portable-v2'
OUT=ART/'compact192_followup_supplement_portable_replay_v2.json'
SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'

def main():
    path=ART/'compact192_followup_supplement_evidence_v2.json';manifest=cert.read(path);workspace=D/'workspace';folder=D/'verification'
    if workspace.exists() or OUT.exists():raise FileExistsError('preserve isolated supplement replay')
    archive=ROOT/manifest['archive']
    if manifest['required_base_archives'] or cert.hashed(archive)!=manifest['archive_sha256']:raise ArithmeticError('standalone supplement binding differs')
    workspace.mkdir(parents=True)
    with zipfile.ZipFile(archive) as z:
        if any(Path(name).is_absolute() or '..' in Path(name).parts for name in z.namelist()):raise ArithmeticError('unsafe archive path')
        z.extractall(workspace)
    for row in manifest['files']:
        if cert.hashed(workspace/row['path'])!=row['sha256']:raise ArithmeticError('isolated member differs')
    cas=workspace/'elliptic-curves/cas';local=workspace/'artifacts/local/elliptic-curves';jobs=[]
    p=cert.read(local/'compact192-specialized-followup-v1/protocol.json');proof=cert.read(local/'compact192-specialized-followup-v1/verification-ledger.json')
    if len(p['rows'])!=5 or [r['id'] for r in p['rows']]!=[r['id'] for r in proof['rows']]:raise ArithmeticError('fixed five-curve proof roster differs')
    for i,row in enumerate(proof['rows']):
        jobs.append((row['id']+'-history',sys.executable,[str(cas/'compact192_specialized_followup.py'),'replay','--index',str(i)],600))
        for key,script in [('mod2_certificate','audit_recorded_point_mod2_rank_v3.py'),('modl_certificate','audit_retained_cloud_modl.py')]:
            jobs.append((row['id']+'-'+key,sys.executable,[str(cas/script),'--check',str(workspace/row[key])],300))
    for name,executable,script,args,seconds in [
        ('geometry',SAGE,'replay_compact192_specialized_geometry.sage',[],1200),
        ('followup-summary',sys.executable,'report_compact192_specialized_followup.py',['--check'],300),
        ('resultants',sys.executable,'audit_r17_scaling_prime_support.py',['check'],120),
        ('independent-resultants',SAGE,'verify_r17_scaling_prime_support.sage',[],120),
        ('residue-trees',sys.executable,'classify_r17_other_small_prime_scalings.py',['check'],300),
        ('independent-residues',sys.executable,'verify_r17_other_small_prime_scalings.py',[],300),
        ('universal-small-prime',sys.executable,'report_r17_small_prime_minimality.py',['--check'],300)]:
        jobs.append((name,executable,[str(cas/script),*args],seconds))
    checkpoint(folder/'protocol.json',{'manifest_sha256':cert.hashed(path),'verifier_sha256':cert.hashed(Path(__file__).resolve()),
        'jobs':[{'name':n,'executable':e,'args':a,'wall_seconds':s} for n,e,a,s in jobs],'rss_bytes':2147483648,'maximum_workers':1,'scope':manifest['claim_boundary']})
    ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'ledger.json',ledger)
    for name,executable,args,seconds in jobs:
        s=run([executable,*args],limits=Limits(seconds,2147483648),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=workspace)
        ok=s['outcome']=='completed' and s['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(folder/'ledger.json',ledger)
        print('ISOLATED SUPPLEMENT',name,ledger['rows'][-1]['status'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(folder/'ledger.json',ledger);raise ArithmeticError('isolated supplement failed; no retry')
    for row in manifest['files']:
        if cert.hashed(workspace/row['path'])!=row['sha256']:raise ArithmeticError('proof member changed during replay')
    ledger['status']='PASS';checkpoint(folder/'ledger.json',ledger)
    paths=[Path(__file__).resolve(),path,folder/'protocol.json',folder/'ledger.json']
    checkpoint(OUT,{'schema':'elliptic-curves.compact192-followup-supplement-portable-replay.v2','status':'PASS','logical_stages':len(jobs),
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'archive_sha256':manifest['archive_sha256'],'ledger':ledger,'claim_boundary':manifest['claim_boundary']})
    print('ISOLATED FOLLOWUP AND SMALL-PRIME SUPPLEMENT',len(jobs),'STAGES PASS',flush=True)
if __name__=='__main__':main()
