#!/usr/bin/env sage-python
"""Three-prime bounded continuation on the sole retained family pair."""
import argparse
from hashlib import sha256
import importlib.util
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import resource
import time

ROOT=Path(__file__).resolve().parents[2]
PATH=ROOT/'artifacts/generated-results/elkies-k3-r17-one-node-complete-pairs-v1'


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=('prepare','run'));args=parser.parse_args()
    loader=SourceFileLoader('complete_pair_continuation',str(ROOT/'elkies-k3/scripts/compare_r17_one_node_complete_families.sage'))
    spec=importlib.util.spec_from_loader(loader.name,loader);M=importlib.util.module_from_spec(spec);loader.exec_module(M)
    if args.mode=='prepare':
        result=json.loads((PATH/'result.json').read_text());assert result['survivors']==[[27,3]]
        M.write_new(PATH/'continuation-input.json',{'schema':'r17-one-node-pair-continuation-input-v1',
                    'prior_input_sha256':sha256((PATH/'input.json').read_bytes()).hexdigest(),
                    'prior_result_sha256':sha256((PATH/'result.json').read_bytes()).hexdigest(),
                    'pairs':[[27,3]],'primes':[127,131,137],'cpu_seconds':20,'address_space_gib':4,
                    'reason':'The sole surviving pair has only singular target quartics at the five original primes. A projective image miss at another good prime would exclude even those possible rational lifts; singular reduction alone is not an exclusion.'})
        print('Prepared the sole retained pair at127,131,137.',flush=True);return
    resource.setrlimit(resource.RLIMIT_CPU,(20,25));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    started=time.process_time();plan=json.loads((PATH/'continuation-input.json').read_text());packet=json.loads((PATH/'input.json').read_text())
    assert plan['prior_input_sha256']==sha256((PATH/'input.json').read_bytes()).hexdigest()
    assert plan['prior_result_sha256']==sha256((PATH/'result.json').read_bytes()).hexdigest()
    records=[];excluded=False
    for p in plan['primes']:
        C=M.mod_coefficients(packet['source_families'][27]['channels'],p);B=M.mod_coefficients(packet['target_matrices'][3],p)
        if C is None or B is None or not int(M.matrix(M.QQ,B).det())%p:
            records.append({'prime':p,'status':'DEFERRED_BAD_COEFFICIENT_REDUCTION'});continue
        targets=set(M.target_images(B,p));assert None not in targets
        zeros=0;matches=0;count=0
        for q in M.source_images(C,p):
            count+=1
            if q is None:zeros+=1
            elif q in targets:matches+=1
        row={'prime':p,'source_parameters':count,'source_base_points':zeros,'source_parameters_meeting_target':matches,
             'status':'EXCLUDED' if not zeros and not matches else 'SURVIVES'}
        records.append(row)
        if row['status']=='EXCLUDED':excluded=True;break
    out={'schema':'r17-one-node-pair-continuation-result-v1',
         'continuation_input_sha256':sha256((PATH/'continuation-input.json').read_bytes()).hexdigest(),
         'records':records,'additional_excluded_pairs':[[27,3]] if excluded else [],
         'remaining_pairs':[] if excluded else [[27,3]],'cpu_seconds':time.process_time()-started}
    M.write_new(PATH/'continuation-result.json',out);print(json.dumps(out,sort_keys=True),flush=True)


if __name__=='__main__':main()
