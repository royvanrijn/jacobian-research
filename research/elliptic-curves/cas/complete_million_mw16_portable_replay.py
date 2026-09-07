#!/usr/bin/env python3
"""Repair the isolated missing-input check without rerunning successful arithmetic."""
import sys,zipfile
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/million-height-mw16-portable-v1';OUT=ART/'million_height_mw16_portable_replay_v1.json'
def main():
    if OUT.exists():raise FileExistsError('preserve final portable report')
    ledger=cert.read(D/'verification/ledger.json');manifest=cert.read(ART/'million_height_mw16_evidence_v1.json');supp=cert.read(ART/'million_height_mw16_population_supplement_v1.json');workspace=D/'workspace'
    if len(ledger['rows'])!=47 or [(r['name'],r['status']) for r in ledger['rows'] if r['status']!='PASS']!=[('diagnostic-traces','FAILED_OR_CENSORED')]:raise ArithmeticError('only declared omitted-population failure may be repaired')
    for row in ledger['rows']:
        s=row['supervision']
        if cert.hashed(Path(s['log']))!=s['log_sha256']:raise ArithmeticError('original isolated log differs')
    if supp['base_manifest_sha256']!=cert.hashed(ART/'million_height_mw16_evidence_v1.json') or cert.hashed(ROOT/supp['archive'])!=supp['archive_sha256']:raise ArithmeticError('supplement differs')
    for r in [*manifest['files'],*manifest['inherited_exact_members']]:
        if cert.hashed(workspace/r['path'])!=r['sha256']:raise ArithmeticError('successful-stage source changed')
    with zipfile.ZipFile(ROOT/supp['archive']) as z:
        if any(Path(n).is_absolute() or '..' in Path(n).parts for n in z.namelist()):raise ArithmeticError('unsafe archive path')
        if any((workspace/n).exists() for n in z.namelist()):raise FileExistsError('supplement must add missing files only')
        z.extractall(workspace)
    for r in supp['files']:
        if cert.hashed(workspace/r['path'])!=r['sha256']:raise ArithmeticError('new input hash differs')
    command=[sys.executable,str(workspace/'elliptic-curves/cas/replay_million_and_mw16_diagnostics.py'),'traces'];s=run(command,limits=Limits(120,1610612736),log_path=D/'verification/diagnostic-traces-supplement.log',checkpoint_path=D/'verification/diagnostic-traces-supplement.supervisor.json',cwd=workspace)
    if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('supplemented exact benchmark replay failed')
    paths=[Path(__file__).resolve(),ART/'million_height_mw16_evidence_v1.json',ART/'million_height_mw16_population_supplement_v1.json',D/'verification/protocol.json',D/'verification/ledger.json',D/'verification/diagnostic-traces-supplement.supervisor.json']
    checkpoint(OUT,{'schema':'elliptic-curves.million-height-mw16-portable-replay.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'logical_stages':47,'original_passes':46,'repaired_stage':'diagnostic-traces','supplemented_replay':s,'original_ledger':ledger,'base_archives':manifest['required_base_archives'],'main_archive':{'path':manifest['archive'],'sha256':manifest['archive_sha256']},'supplement_archive':{'path':supp['archive'],'sha256':supp['archive_sha256']},'claim_boundary':'All47 declared isolated logical stages pass after adding five originally pinned population files and rerunning only the failed trace benchmark check. The original46 successful stages and their unchanged source hashes are retained; the original missing-input failure remains visible. No new point search or trace computation. The supplied version2 extractor applies both archives before running the same47 checks.'})
    print('ISOLATED MILLION MW16 PASS:46 ORIGINAL STAGES PLUS REPAIRED TRACE CHECK',flush=True)
if __name__=='__main__':main()
