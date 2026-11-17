"""Run this scoped audit only; each worker is capped at 60 seconds."""
from pathlib import Path
import argparse
import hashlib
import json
from importlib.metadata import version
import subprocess
import sys
import tarfile

ROOT=Path(__file__).resolve().parent
SOURCE_SHA256='0e48fed1f0804ad263c24e7f4b146eb8d0ff618bf22a8f975907829b5cfe49ba'


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--bounded-search',action='store_true',help='repeat the specified 512 grouping and 16 basis trials')
    args=parser.parse_args()
    assert hashlib.sha256((ROOT/'source-v1.tar').read_bytes()).hexdigest()==SOURCE_SHA256
    source_files=0
    with tarfile.open(ROOT/'source-v1.tar') as archive:
        for member in archive.getmembers():
            if member.isfile():
                assert (ROOT/'source'/member.name).read_bytes()==archive.extractfile(member).read()
                source_files+=1
    for name in ['verify.py','equivalence.py','polarize.py','druzkowski.py','check_independent.py','graph_audit.py']+(['bounded_search.py'] if args.bounded_search else []):
        subprocess.run([sys.executable,str(ROOT/name)],cwd=ROOT,check=True,timeout=60)
    graphs=json.loads((ROOT/'graph_audit.json').read_text())
    assert len(graphs)==18
    assert all(r['augmented_rank']==r['rank']+1 for r in graphs)
    files=['source-v1.tar','verify.py','equivalence.py','polarize.py','druzkowski.py',
           'check_independent.py','graph_audit.py','bounded_search.py',
           'audit.json','equivalence.json','polarization.json','druzkowski90.json','independent.json','graph_audit.json']
    if (ROOT/'bounded_search.json').is_file():files.append('bounded_search.json')
    result={'source_url':'https://arxiv.org/src/2610.24889v1','source_sha256':SOURCE_SHA256,
            'status':'PASS','coverage':'12/18/90 constructions, exact six-cubic conjugacy and nine literal graph obstructions',
            'source_files_checked':source_files,
            'software':{'python':sys.version.split()[0],'sympy':version('sympy'),'python-flint':version('python-flint')},
            'not_claimed':['complete audit of arbitrary-ring theorems','independent Vanishing Conjecture operator replay',
                           'new dimension record','global minimality','formal verification','external review'],
            'sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files}}
    (ROOT/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS SCOPED DE BONDT AUDIT')


if __name__=='__main__':main()
