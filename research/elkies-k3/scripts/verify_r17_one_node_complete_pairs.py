#!/usr/bin/env python3
"""Independent Fraction/integer proof replay of the complete family gate.

The producer expands slopes and convolves polynomials at each parameter.
This checker reconstructs the exact coefficient channels, expands their
binomial tensor once, and evaluates the projective source by row Horner
evaluation. It checks every reduction point, including all boundaries.
"""
import argparse
from hashlib import sha256
import importlib.util
import json
from math import comb
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/"artifacts/generated-results/elkies-k3-r17-one-node-complete-pairs-v1"
SPEC=importlib.util.spec_from_file_location('one_node_base_replay',ROOT/'elkies-k3/scripts/verify_r17_one_node_correlated.py')
N=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(N)
V=N.V;require=V.require


def normalize(v,p):
    for x in v:
        if x%p:
            inv=pow(x%p,-1,p);return tuple(c*inv%p for c in v)
    return None


def reduce_channels(rows,p):
    rows=[V.poly(row) for row in rows]
    if any(c.denominator%p==0 for row in rows for c in row):return None
    return [[c.numerator*pow(c.denominator,-1,p)%p for c in row] for row in rows]


def source_image(C,p):
    tensor={}
    for j,f in enumerate(C):
        for v in range(j+1):
            vec=[0]*7
            for k,c in enumerate(f):vec[k+v]=comb(j,v)*c%p
            tensor[j-v,v]=vec
    for a in range(p):
        ap=[pow(a,j,p) for j in range(5)]
        row=[[sum(ap[u]*tensor[u,v][k] for u in range(5-v))%p for k in range(7)] for v in range(5)]
        for b in range(p):
            out=[0]*7
            for vec in reversed(row):out=[(x*b+y)%p for x,y in zip(out,vec)]
            yield normalize(out,p)
    # c=0 uses exactly the homogeneous degree-four tensor.
    for a in range(p):
        yield normalize([sum(pow(a,4-v,p)*tensor[4-v,v][k] for v in range(5))%p for k in range(7)],p)
    yield normalize(tensor[4,0],p)


def target_image(B,p):
    for lam in list(range(p))+[None]:
        d=([sum(row[j]*pow(lam,j,p) for j in range(5))%p for row in B]
           if lam is not None else [row[4] for row in B])
        for r in range(p):
            q=[((r*r*d[k] if k<5 else 0)-(2*r*d[k-1] if 0<=k-1<5 else 0)+(d[k-2] if 0<=k-2<5 else 0))%p for k in range(7)]
            yield normalize(q,p)
        yield normalize(d+[0,0],p)


def verify_continuation(path,packet,remaining):
    input_path=path/'continuation-input.json';output_path=path/'continuation-result.json'
    require(input_path.is_file() and output_path.is_file(),'missing completion certificate: the final pair remains UNKNOWN')
    plan=json.loads(input_path.read_text());result=json.loads(output_path.read_text())
    require(plan['schema']=='r17-one-node-pair-continuation-input-v1','continuation schema')
    require(plan['prior_input_sha256']==sha256((path/'input.json').read_bytes()).hexdigest()
            and plan['prior_result_sha256']==sha256((path/'result.json').read_bytes()).hexdigest(),'continuation prior binding')
    require(plan['pairs']==sorted(map(list,remaining))==[[27,3]] and plan['primes']==[127,131,137],'continuation selection')
    require(result['continuation_input_sha256']==sha256(input_path.read_bytes()).hexdigest(),'continuation input binding')
    counts=[];parameters=0;excluded=[]
    for k,row in enumerate(result['records']):
        p=row['prime'];require(p==plan['primes'][k] and remaining,'continuation prime order')
        C=reduce_channels(packet['source_families'][27]['channels'],p)
        B=V.reduce_matrix([[V.Q(c) for c in r] for r in packet['target_matrices'][3]],p)
        if C is None or B is None:
            require(row=={'prime':p,'status':'DEFERRED_BAD_COEFFICIENT_REDUCTION'},'continuation bad reduction');continue
        target=set(target_image(B,p));require(None not in target,'continuation target base point')
        size=zeros=hits=0
        for q in source_image(C,p):
            size+=1
            if q is None:zeros+=1
            elif q in target:hits+=1
        status='EXCLUDED' if not zeros and not hits else 'SURVIVES'
        require(size==p*p+p+1 and row=={'prime':p,'source_parameters':size,'source_base_points':zeros,
                'source_parameters_meeting_target':hits,'status':status},'continuation image audit')
        parameters+=size
        if status=='EXCLUDED':remaining.remove((27,3));excluded.append([27,3])
        counts.append({'prime':p,'new_exclusions':int(status=='EXCLUDED'),'remaining':len(remaining)})
    require(result['additional_excluded_pairs']==excluded and result['remaining_pairs']==sorted(map(list,remaining)),'continuation conclusion')
    require(not remaining or len(result['records'])==len(plan['primes']),'incomplete continuation')
    return counts,parameters


def verify(path,check_dependencies=True):
    packet=json.loads((path/'input.json').read_text());result=json.loads((path/'result.json').read_text())
    require(packet['schema']=='r17-one-node-complete-pairs-input-v1','input schema')
    require(result['input_sha256']==sha256((path/'input.json').read_bytes()).hexdigest(),'input binding')
    src_path=N.DEFAULT/'input.json';prior=json.loads(src_path.read_text())
    pencil_path=V.DEFAULT/'pencils.json';pencils=json.loads(pencil_path.read_text())
    require(packet['source_input_sha256']==sha256(src_path.read_bytes()).hexdigest() and packet['target_pencils_sha256']==sha256(pencil_path.read_bytes()).hexdigest(),'prior input binding')
    if check_dependencies:
        N.verify(N.DEFAULT)
        V.verify(json.loads((V.DEFAULT/'input.json').read_text()),pencils,json.loads((V.DEFAULT/'result.json').read_text()))
    A=V.poly(prior['generic_source']['A'])
    require(len(packet['source_families'])==67 and len(packet['target_matrices'])==25,'family counts')
    for k,(trace,row) in enumerate(zip(prior['traces'],packet['source_families'])):
        require(row['index']==k and row['word']==trace['word'],'source trace attachment')
        require([V.poly(f) for f in row['channels']]==N.q_channels(trace,A),'source coefficient identity')
    require(packet['target_matrices']==[r['branch_matrix'] for r in pencils],'target branch matrices')
    require(packet['limits']['primes']==[101,103,107,109,113],'prime budget')
    remaining={(i,j) for i in range(67) for j in range(25)};counts=[];total_parameters=0
    for file in result['prime_certificates']:
        raw=(path/file['path']).read_bytes();require(sha256(raw).hexdigest()==file['sha256'],'prime checkpoint hash')
        certificate=json.loads(raw);p=certificate['prime']
        require(p in packet['limits']['primes'] and file['path']==f'prime-{p}.json','prime identity')
        targets={};good=[]
        for j,rows in enumerate(packet['target_matrices']):
            B=V.reduce_matrix([[V.Q(c) for c in r] for r in rows],p)
            if B is None:continue
            good.append(j)
            for q in target_image(B,p):
                require(q is not None,'target projective base point')
                targets[q]=targets.get(q,0)|(1<<j)
        require(good==certificate['good_targets'],'good target roster')
        excluded=[]
        for i,(row,audit) in enumerate(zip(packet['source_families'],certificate['source_audits'])):
            C=reduce_channels(row['channels'],p)
            if C is None:
                require(audit=={'source':i,'status':'DEFERRED_NONINTEGRAL_COEFFICIENTS'},'source coefficient deferral');continue
            seen=0;zeros=0;size=0
            for q in source_image(C,p):
                size+=1
                if q is None:zeros+=1
                else:seen|=targets.get(q,0)
            total_parameters+=size
            require(size==p*p+p+1,'source projective coverage')
            require(audit=={'source':i,'projective_parameters':size,'base_points':zeros,'intersecting_targets':[j for j in good if (seen>>j)&1]},'source image audit')
            if not zeros:
                for j in good:
                    if (i,j) in remaining and not (seen>>j)&1:
                        remaining.remove((i,j));excluded.append([i,j])
        require(len(certificate['source_audits'])==67 and excluded==certificate['excluded_pairs'],'pair exclusion coverage')
        require(sorted(map(list,remaining))==certificate['remaining_pairs'],'prime survivor reconciliation')
        counts.append({'prime':p,'new_exclusions':len(excluded),'remaining':len(remaining)})
    require([r['prime'] for r in counts]==packet['limits']['primes'][:len(counts)],'prime execution order')
    require(result['pairs']==1675 and result['excluded_pairs']==1675-len(remaining) and result['survivors']==sorted(map(list,remaining)),'result reconciliation')
    continuation,parameters=verify_continuation(path,packet,remaining)
    counts.extend(continuation);total_parameters+=parameters
    require(not remaining,'complete family exclusion remains UNKNOWN')
    return {'status':'PASS_EXACT_PROJECTIVE_FAMILY_EXCLUSIONS','pairs':1675,'excluded_pairs':1675-len(remaining),
            'survivors':sorted(map(list,remaining)),'prime_counts':counts,'source_parameter_evaluations':total_parameters,
            'scope':'Each excluded pair has no rational coefficient match anywhere on its P2-by-P1-by-P1 parameter space. Surviving pairs remain UNKNOWN; other families and quadratic covers are not excluded.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-dir',type=Path,default=DEFAULT);p.add_argument('--output',type=Path)
    a=p.parse_args();started=time.monotonic();result=verify(a.input_dir);result['elapsed_seconds']=round(time.monotonic()-started,3)
    if a.output:
        with a.output.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(result,sort_keys=True))


if __name__=='__main__':main()
