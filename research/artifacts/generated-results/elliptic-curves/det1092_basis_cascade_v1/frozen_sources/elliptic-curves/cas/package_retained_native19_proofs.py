#!/usr/bin/env python3
"""Four-file, fresh-directory independent Sage replay of the retained rank proofs."""
import argparse,json,tempfile,zipfile
from hashlib import sha256
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/retained-native19-point-trial-v3'
ZIP=ART/'retained_native19_rank_standalone_v1.zip';OUT=ZIP.with_suffix('.json')
VERIFIER=CAS/'verify_retained_native19_rank.sage'
INPUTS=[ART/('retained_native19_trial_v3_'+n+'_mod2_v1.json') for n in ('08','01')]+[ART/'factor_free_known28_control_v1.json']
REFERENCE=ART/'retained_native19_sage_replay_v1.json'

def payload():
    source=VERIFIER.read_text().replace("ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'", "ROOT=Path(__file__).resolve().parent;ART=ROOT")
    return {VERIFIER.name:source.encode(),**{p.name:p.read_bytes() for p in INPUTS}}

def main(check):
    files=payload();hashes={n:sha256(v).hexdigest() for n,v in files.items()}
    sources={str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),VERIFIER,*INPUTS,REFERENCE]}
    if check:
        m=cert.read(OUT);assert m['sources']==sources and m['files']==hashes and m['zip_sha256']==cert.hashed(ZIP)
        with zipfile.ZipFile(ZIP) as z:assert {n:z.read(n) for n in z.namelist()}==files
        s=cert.read(D/'standalone.supervisor.json');assert s['outcome']=='completed' and s['returncode']==0 and s['log_sha256']==cert.hashed(D/'standalone.log')
        assert m['supervision']==s
        print('PASS standalone bundle bindings');return
    if ZIP.exists() or OUT.exists() or (D/'standalone.supervisor.json').exists():raise FileExistsError('preserve standalone evidence')
    with zipfile.ZipFile(ZIP,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for name,data in sorted(files.items()):
            info=zipfile.ZipInfo(name,date_time=(2026,9,7,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o100644<<16;z.writestr(info,data)
    with tempfile.TemporaryDirectory(prefix='native19-standalone-') as temp:
        folder=Path(temp)
        with zipfile.ZipFile(ZIP) as z:z.extractall(folder)
        s=run(['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',str(folder/VERIFIER.name)],limits=Limits(120,1610612736),log_path=D/'standalone.log',checkpoint_path=D/'standalone.supervisor.json',cwd=folder)
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('standalone verification failed')
        independent=json.loads((folder/REFERENCE.name).read_text());reference=cert.read(REFERENCE)
        for a,b in zip(independent['rows'],reference['rows']):
            assert Path(b['input']).name==a['input']
            assert {k:v for k,v in a.items() if k!='input'}=={k:v for k,v in b.items() if k!='input'}
        assert len(independent['rows'])==len(reference['rows'])==3 and independent['status']==reference['status']=='PASS'
    checkpoint(OUT,{'status':'PASS','sources':sources,'files':hashes,'zip_sha256':cert.hashed(ZIP),'zip_bytes':ZIP.stat().st_size,'supervision':s,
        'rank_lower_bounds':[19,19,28],'boundary':'Fresh-directory Sage full finite-group enumeration on two retained native19 clouds and one known28 control, with no repository imports. Proves lower bounds only. Search exposure and exact maps are separately replayed.'})
    print('PASS standalone19,19,28',ZIP.stat().st_size)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');main(p.parse_args().check)
