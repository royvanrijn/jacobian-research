#!/usr/bin/env sage-python
"""Classify the retained projective matches for the sole unexcluded pair."""
from hashlib import sha256
import importlib.util
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import resource
import time

from sage.all import GF, PolynomialRing

ROOT=Path(__file__).resolve().parents[2]
PATH=ROOT/'artifacts/generated-results/elkies-k3-r17-one-node-complete-pairs-v1'


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(60,65));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    started=time.process_time()
    loader=SourceFileLoader('complete_pair_producer',str(ROOT/'elkies-k3/scripts/compare_r17_one_node_complete_families.sage'))
    spec=importlib.util.spec_from_loader(loader.name,loader);M=importlib.util.module_from_spec(spec);loader.exec_module(M)
    packet=json.loads((PATH/'input.json').read_text());result=json.loads((PATH/'result.json').read_text())
    assert result['survivors']==[[27,3]]
    records=[]
    for p in packet['limits']['primes']:
        C=M.mod_coefficients(packet['source_families'][27]['channels'],p)
        B=M.mod_coefficients(packet['target_matrices'][3],p)
        assert C is not None and B is not None
        targets={}
        for k,q in enumerate(M.target_images(B,p)):
            lam,rr=divmod(k,p+1);lam=None if lam==p else lam;rr=None if rr==p else rr
            targets.setdefault(q,[]).append((lam,rr))
        matches=[];R=PolynomialRing(GF(p),'t')
        for k,q in enumerate(M.source_images(C,p)):
            if q not in targets:continue
            if k<p*p:a,b=divmod(k,p);c=1
            elif k<p*p+p:a,b,c=k-p*p,1,0
            else:a,b,c=1,0,0
            for lam,r in targets[q]:
                d=([sum(row[j]*pow(lam,j,p) for j in range(5))%p for row in B] if lam is not None else [row[4] for row in B])
                f=R(d)
                matches.append({'source_parameter':[a,b,c],'lambda':lam,'node':r,
                                'target_quartic':d,'target_affine_degree':int(f.degree()),
                                'target_repeated_factor_degree':int(f.gcd(f.derivative()).degree())})
        records.append({'prime':p,'matches':matches})
    out={'schema':'r17-one-node-survivor-diagnostic-v1','source_index':27,'target_index':3,
         'input_sha256':sha256((PATH/'input.json').read_bytes()).hexdigest(),'records':records,
         'cpu_seconds':time.process_time()-started,
         'scope':'These finite-field intersections are necessary local conditions, not rational solutions or a rank construction.'}
    M.write_new(PATH/'survivor-diagnostic.json',out)
    print(json.dumps(out,sort_keys=True),flush=True)


if __name__=='__main__':main()
