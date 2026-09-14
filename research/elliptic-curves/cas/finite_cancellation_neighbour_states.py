#!/usr/bin/env python3
"""Exact two-state finite-height law and target-blind q residue-ball priors."""
from collections import Counter, defaultdict
from fractions import Fraction as F
import gzip
import json
from math import gcd
from pathlib import Path

from finite_cancellation_corpus import OUT, digest, write
from half_lattice_pointed_sieve import binary_transform
from analyze_finite_cancellation import clustered

def v(n,p):
    k=0
    while n%p==0:n//=p;k+=1
    assert n==1
    return k

def main():
    laws=[];constants=Counter();prediction_errors=defaultdict(lambda:defaultdict(list));state=defaultdict(lambda:defaultdict(list))
    for path in sorted((OUT/'cases').glob('*.json.gz')):
        d=json.loads(gzip.decompress(path.read_bytes()))
        if d['status']!='PASS_EXACT_ACCESSIBILITY':continue
        models=d['prepared']['models'];base=models[0];N,D=[list(map(int,base[k])) for k in ['N','D']]
        obs={(r['target_index'],r['sign'],r['model_index']):r for r in d['observations']}
        for mi,M in enumerate(models[1:],1):
            p,r=M['prime'],M['residue'];T=tuple(map(F,M['mapping']['second_matrix']))
            assert all(a.denominator==1 for a in T);T=tuple(map(int,T));a,b,c,z=T;assert abs(a*z-b*c)==p
            nn=binary_transform(N,T);dd=binary_transform(D,T);content=gcd(*(nn+dd));e=v(content,p);assert 0<=e<=4
            assert [x//content for x in nn]==list(map(int,M['N'])) and [x//content for x in dd]==list(map(int,M['D']))
            constants[e]+=1;tree=next(t for t in base['local'] if t['prime']==p);good=[l for l in tree['leaves'] if l['status']=='square']
            probabilities={}
            for label,weighted in [('uniform_soluble',False),('differential',True)]:
                ws=[float(F(l['mass']))*(p**(l['vq']/2) if weighted else 1) for l in good]
                probabilities[label]=sum(w for l,w in zip(good,ws) if l['kind']=='affine' and l['residue']%p==r)/sum(ws) if sum(ws) else None
            counts=Counter()
            for (ti,sign,j),row in obs.items():
                if j!=mi:continue
                old=obs[ti,sign,0];m,n=map(int,old['coordinate']);root=(m-r*n)%p==0
                delta=int(not root);ratio=F(p)**(4*delta-e)
                assert F(int(row['g']),int(old['g']))==ratio
                assert F(int(row['H'])**4,int(old['H'])**4)==ratio*F(old['S'])/F(row['S'])
                counts['root' if root else 'outside']+=1
                if d['family']!='Curve302-development':
                    key='root_ball' if root else 'outside_ball';group=d['j_group']
                    state[key][group].append(float(int(row['H'])<int(old['H'])))
                    for label,prob in probabilities.items():
                        if prob is not None:prediction_errors[label][group].append(float(root)-prob)
            laws.append({'case':d['case_id'],'anchor':d['prepared']['anchor_index'],'model':mi,'p':p,'r':r,
                'content_valuation':e,'predicted_root_probabilities':probabilities,'observed_state_counts':dict(counts)})
    write(OUT/'neighbour_states.json',{'status':'PASS_ALL_EXACT_TWO_STATE_LAWS','neighbours':len(laws),'content_valuations':dict(constants),
        'root_probability_calibration_error':{k:clustered(z) for k,z in prediction_errors.items()},
        'probability_neighbour_has_smaller_literal_H':{k:clustered(z) for k,z in state.items()},'laws':laws,
        'source_sha256':digest(Path(__file__).read_bytes()),
        'boundary':'General integral degree4 change-of-model identity, checked on every retained neighbour and target. Ball priors use only the base q residue tree; calibration conditions on resolved soluble mass. Explanatory post-primary diagnostic, no selector retuning.'})
    print(json.dumps({'status':'PASS_ALL_EXACT_TWO_STATE_LAWS','neighbours':len(laws),'content_valuations':dict(constants),
        'calibration':{k:clustered(z) for k,z in prediction_errors.items()},'cheap_by_state':{k:clustered(z) for k,z in state.items()}}))

if __name__=='__main__':main()
