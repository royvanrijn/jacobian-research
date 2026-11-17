"""Check retained replay evidence; optionally run fresh bounded CAS workers.

Default mode is an integrity and coverage-reference check, not a fresh CAS run.
The mathematical assembly and normalization review are written proofs.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
EXPECTED={'prelim':19,'binary':148,'unary1':12,'unary2':21,'unary3':47,'ternary':29}

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def author_receipt(job,directory):
    receipt=json.loads((directory/f'author-{job}.json').read_text())
    log=directory/f'author-{job}.log'
    assert sha(log)==receipt['log_sha256']
    assert sha(HERE/receipt['source'])==receipt['source_sha256']
    assert receipt['terminal']=='completed' and receipt['exit_code']==0 and receipt['outcome']=='PASS'
    text=log.read_text()
    good=len(re.findall(r'^\s*(?:\[OK\]|PASS\s)',text,re.M))
    bad=len(re.findall(r'^\s*(?:\[ECHEC\]|FAIL\s)',text,re.M))
    assert good==receipt['passed_checks']==receipt['expected_checks']==EXPECTED[job]
    assert bad==receipt['failed_checks']==0
    assert 1<=receipt['cap_seconds']<=240
    return good

def bounded(command,log,cap):
    with log.open('w') as stream:
        process=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        try:code=process.wait(timeout=cap)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGTERM)
            try:process.wait(timeout=3)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL);process.wait(timeout=5)
            raise RuntimeError(f'Bounded replay reached its {cap}s cap; retained {log}')
    assert code==0,f'Worker failed; inspect {log}'

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--fresh-author',action='store_true')
    parser.add_argument('--fresh-zquad',action='store_true')
    args=parser.parse_args()
    for row in json.loads((HERE/'sources.json').read_text()):
        if 'file' in row:assert sha(HERE/row['file'])==row['sha256']
    for row in json.loads((HERE/'source-files.json').read_text()):
        assert sha(HERE/row['file'])==row['sha256']
    total=sum(author_receipt(job,HERE) for job in EXPECTED)
    assert total==276
    z=json.loads((HERE/'zquad-receipt.json').read_text())
    for name,digest in z['hashes'].items():assert sha(HERE/name)==digest
    result=json.loads((HERE/'zquad-result.json').read_text())
    assert result['unit_ideal'] and result['basis_size']==1 and result['basis']==['1']
    inp=json.loads((HERE/'zquad-input.json').read_text())
    assert len(inp['monomials'])==len(set(map(tuple,inp['monomials'])))==35
    assert len(inp['equations'])==167
    fresh={}
    if args.fresh_author or args.fresh_zquad:
        local=ROOT/'artifacts/local';local.mkdir(parents=True,exist_ok=True)
        output=Path(tempfile.mkdtemp(prefix='hc4-degree-five-fresh-',dir=local))
        print('Retained fresh output:',output,flush=True)
        if args.fresh_author:
            for job in EXPECTED:
                cap=240 if job=='ternary' else 180
                bounded([sys.executable,str(HERE/'replay_author.py'),job,'--cap',str(cap),'--output',str(output)],output/f'wrapper-{job}.log',cap+15)
                fresh[job]=author_receipt(job,output)
        if args.fresh_zquad:
            (output/'independent_zquad.py').write_bytes((HERE/'independent_zquad.py').read_bytes())
            bounded(['/home/royvanrijn/.local/bin/sage','-python',str(output/'independent_zquad.py')],output/'zquad.log',60)
            zr=json.loads((output/'zquad-result.json').read_text())
            assert zr['unit_ideal'] and zr['basis']==['1']
            fresh['zquad']='unit ideal reproduced'
    print(json.dumps(dict(retained_author_checks=total,source_and_log_hashes='PASS',complete_35_coefficient_Z_input='PASS',retained_Z_unit_ideal='PASS',fresh_replays=fresh,normalization_exhaustiveness='written review in NORMALIZATION_REVIEW.md',full_independent_implementation=False,formal_verification=False,external_review=False,JC2_or_unrestricted_HC4_resolved=False),indent=2))

if __name__=='__main__':main()
