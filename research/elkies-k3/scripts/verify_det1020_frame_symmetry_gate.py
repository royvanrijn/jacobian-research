#!/usr/bin/env python3
"""Bounded exact automorphism computation on one retained rank17 frame."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import resource
import subprocess
import time
ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'artifacts/generated-results/elkies-k3-det1020-root-gate-v1/certificate.json'
PACKET=ROOT/'artifacts/generated-results/elkies-k3-det1020-frame-symmetry-gate-v1'

def run():
    started=time.monotonic();d=json.loads(SOURCE.read_text());H=d['frame_gram'];n=17
    assert d['rank']==n and d['determinant']==1020 and len(H)==n
    assert all(len(row)==n for row in H) and all(H[i][j]==H[j][i] for i in range(n) for j in range(n))
    literal='['+';'.join(','.join(map(str,row)) for row in H)+']'
    code='default(nbthreads,1);\ndefault(parisize,134217728);\nH='+literal+';a=qfauto(H);if(a[1]!=2,error("unexpected group order"));if(#a[2]!=1,error("unexpected generators"));if(a[2][1]!=-matid(17),error("unexpected generator"));print("PASS_ORDER_2_GENERATOR_MINUS_IDENTITY");quit;\n'
    def limits():
        resource.setrlimit(resource.RLIMIT_CPU,(20,20))
        resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
    p=subprocess.run(['gp','-q','-f'],input=code,text=True,capture_output=True,timeout=25,preexec_fn=limits,check=True)
    assert p.stdout.strip()=='PASS_ORDER_2_GENERATOR_MINUS_IDENTITY',p.stderr
    return {'status':'PASS','source_sha256':sha256(SOURCE.read_bytes()).hexdigest(),
            'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'pari_version':subprocess.run(['gp','--version-short'],capture_output=True,text=True,check=True).stdout.strip(),
            'rank':n,'frame_determinant':1020,'automorphism_order':2,'generator':'-I_17',
            'involution_traces':[-17],'required_pointed_symplectic_involution_trace':1,
            'compatible_involutions':0,'independent_implementation':False,
            'scope':'Only the frame in the retained4E6 determinant1020 witness; other frames and nonsymmetry correlated-gain constructions are not excluded.',
            'elapsed_seconds':time.monotonic()-started}
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--write',action='store_true');args=a.parse_args();out=run()
    if args.write:
        PACKET.mkdir(exist_ok=True);(PACKET/'result.json').write_text(json.dumps(out,indent=2)+'\n')
    else:
        old=json.loads((PACKET/'result.json').read_text())
        assert {k:v for k,v in out.items() if k!='elapsed_seconds'}=={k:v for k,v in old.items() if k!='elapsed_seconds'}
    print(json.dumps(out,sort_keys=True))
