"""Selected, checkpointed author-source replay; no full-proof claim from counts."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time

HERE=Path(__file__).resolve().parent
JOBS={
    'prelim': ('2611.00161','verify_hc4_prelim.py',19),
    'binary': ('2611.00161','verify_hc4_binary.py',148),
    'unary1': ('2611.00161','verify_hc4_unary_1.py',12),
    'unary2': ('2611.00161','verify_hc4_unary_2.py',21),
    'unary3': ('2611.00161','verify_hc4_unary_3.py',47),
    'ternary': ('2611.00134','verify_hc4_ternary.py',29),
}

def run(job, cap, output_dir=None):
    for row in json.loads((HERE/'source-files.json').read_text()):
        assert hashlib.sha256((HERE/row['file']).read_bytes()).hexdigest()==row['sha256']
    ident,filename,expected=JOBS[job]
    source=HERE/f'source-{ident}/anc'/filename
    output_dir=HERE if output_dir is None else Path(output_dir).resolve()
    output_dir.mkdir(parents=True,exist_ok=True)
    log=output_dir/f'author-{job}.log'
    receipt=output_dir/f'author-{job}.json'
    assert not receipt.exists() and not log.exists(), 'Preserve existing checkpoints; use a new output directory for a fresh replay.'
    command=['/home/royvanrijn/.local/bin/sage','-python',str(HERE/'author_adapter.py'),str(source)]
    start=time.monotonic()
    with log.open('w') as output:
        process=subprocess.Popen(command,stdout=output,stderr=subprocess.STDOUT,
                                 cwd=source.parent,start_new_session=True)
        try:
            code=process.wait(timeout=cap)
            terminal='completed'
        except subprocess.TimeoutExpired:
            terminal='time_cap'
            os.killpg(process.pid,signal.SIGTERM)
            try:code=process.wait(timeout=3)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL);code=process.wait(timeout=5)
    text=log.read_text()
    passed=len(re.findall(r'^\s*(?:\[OK\]|PASS\s)',text,re.M))
    failed=len(re.findall(r'^\s*(?:\[ECHEC\]|FAIL\s)',text,re.M))
    outcome='PASS' if terminal=='completed' and code==0 and failed==0 and passed==expected else 'INCOMPLETE' if terminal=='time_cap' else 'FAIL'
    result=dict(job=job,source=str(source.relative_to(HERE)),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),cap_seconds=cap,elapsed_seconds=round(time.monotonic()-start,3),terminal=terminal,exit_code=code,expected_checks=expected,passed_checks=passed,failed_checks=failed,outcome=outcome,independent_implementation=False,full_case_exhaustiveness_verified=False,log_sha256=hashlib.sha256(log.read_bytes()).hexdigest())
    receipt.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    print('\n'.join(text.splitlines()[-5:]))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('job',choices=JOBS);parser.add_argument('--cap',type=int,default=180);parser.add_argument('--output')
    args=parser.parse_args();assert 1<=args.cap<=240
    run(args.job,args.cap,args.output)
