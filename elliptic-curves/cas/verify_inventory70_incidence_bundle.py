#!/usr/bin/env python3
"""Six bounded isolated replays from the self-contained incidence archive."""
from pathlib import Path
import zipfile,sys
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/inventory70-incidence-portable-v1';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'
def main():
    manifest=cert.read(ART/'inventory70_incidence_evidence_v1.json');workspace=D/'workspace'
    if workspace.exists():raise FileExistsError('preserve isolated incidence replay')
    archive=ROOT/manifest['archive']
    if cert.hashed(archive)!=manifest['archive_sha256']:raise ArithmeticError('archive changed')
    with zipfile.ZipFile(archive) as z:
        if any(Path(n).is_absolute() or '..' in Path(n).parts for n in z.namelist()):raise ArithmeticError('unsafe archive path')
        z.extractall(workspace)
    for row in manifest['files']:
        if cert.hashed(workspace/row['path'])!=row['sha256']:raise ArithmeticError('member changed')
    C=workspace/'elliptic-curves/cas';A=workspace/'artifacts/generated-results/elliptic-curves';stages=[]
    for name in ['compact','latest7','latest8','latest23']:
        stages.append((name,[sys.executable,str(C/('replay_'+name+'_cross_family_incidence.py')),'--input',str(A/(name+'_cross_family_j_incidence_v1.json')),'--output',str(workspace/(name+'-replayed.json'))]))
    stages += [('generic-transport',[SAGE,str(C/'audit_compact_published_r17_transport_v3.sage'),'--check',str(A/'compact_published_r17_generic_transport_v1.json')]),('aggregate',[sys.executable,str(C/'certify_inventory70_incidence.py'),'--check',str(A/'inventory70_cross_family_incidence_v1.json')])]
    p={'manifest_sha256':cert.hashed(ART/'inventory70_incidence_evidence_v1.json'),'verifier_sha256':cert.hashed(Path(__file__).resolve()),'wall_seconds_per_stage':120,'rss_bytes':1610612736,'stages':[n for n,_ in stages],'scope':'All840 exact incidence checks, generic17-section transport and fixed70 cohort binding; no rank replay or point search.'};checkpoint(D/'protocol.json',p);rows=[]
    for name,command in stages:
        r=run(command,limits=Limits(120,1610612736),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=workspace)
        rows.append({'name':name,'status':'PASS' if r['outcome']=='completed' and r['returncode']==0 else 'FAILED_OR_CENSORED','supervision':r,'log':(D/(name+'.log')).read_text()});checkpoint(D/'ledger.json',{'status':'RUNNING','rows':rows})
        if rows[-1]['status']!='PASS':raise ArithmeticError('isolated incidence stage failed/censored')
    checkpoint(D/'ledger.json',{'status':'PASS','rows':rows});checkpoint(ART/'inventory70_incidence_portable_replay_v1.json',{'schema':'elliptic-curves.inventory70-incidence-portable-replay.v1','status':'PASS','protocol':p,'rows':rows});print('ISOLATED INVENTORY70 INCIDENCE SIX STAGES PASS',flush=True)
if __name__=='__main__':main()
