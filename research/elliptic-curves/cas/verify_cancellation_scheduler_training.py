#!/usr/bin/env python3
"""Independent arithmetic/partition accounting for the frozen scheduler fit."""
from collections import defaultdict
from fractions import Fraction as F
import gzip
import json
from math import gcd, lcm
import statistics
import time

from cancellation_scheduler_prepare import OUT, CORPUS
from finite_cancellation_corpus import ROOT, digest, write


def main():
    from finite_cancellation_validation_features import q_tree
    from finite_cancellation_validation_replay import verify_tree
    from half_lattice_pointed_sieve import linear_combination
    fit=json.loads((OUT/'fit.json').read_text());roster=json.loads((OUT/'roster.json').read_text())
    train=set(fit['training_j_groups']);test=set(fit['withheld_j_groups'])
    if train & set(roster['excluded_training_j_groups']):raise ArithmeticError('fit leakage')
    start=time.process_time();data=defaultdict(list);leaves=0;anchors=0
    for name,h in fit['corpus_cases_sha256'].items():
        path=CORPUS/'cases'/name
        if digest(path.read_bytes())!=h:raise ArithmeticError('corpus packet changed')
        packet=json.loads(gzip.decompress(path.read_bytes()));prep=packet['prepared']
        q=list(map(int,prep['models'][0]['q']));distributions={}
        for p in sorted({m['prime'] for m in prep['models'][1:]}):
            tree=q_tree(q,p,depth=4,cap=512)
            leaves+=verify_tree(tree,[1,0,0,0,0],[1,0,0,0,0],q)
            mass=defaultdict(F)
            for leaf in tree['leaves']:
                if leaf['status']=='square':
                    key=leaf['residue']%p if leaf['kind']=='affine' else None
                    mass[key]+=F(leaf['mass'])*p**(leaf['vq']//2)
            total=sum(mass.values())
            distributions[p]=(mass,total)
        for model in prep['models'][1:]:
            p,r=model['prime'],model['residue'];mass,total=distributions[p]
            for observation in packet['observations']:
                if observation['model_index']!=0:continue
                m,n=map(int,observation['coordinate'])
                y=int(n%p!=0 and m%p==r*(n%p)%p)
                data[packet['j_group']].append((1/(p+1),float(mass[r]/total) if total else 1/(p+1),y))
        anchors+=1
    numerator=sum(sum((d-u)*(y-u) for u,d,y in data[g])/len(data[g]) for g in train)
    denominator=sum(sum((d-u)**2 for u,d,y in data[g])/len(data[g]) for g in train)
    beta=min(.98,max(0.,numerator/denominator))
    if abs(beta-fit['local_weight'])>1e-12:raise ArithmeticError('independent probability fit differs')
    for name,b in [('uniform',0),('local',beta)]:
        brier=statistics.mean(sum((u+b*(d-u)-y)**2 for u,d,y in data[g])/len(data[g]) for g in test)
        if abs(brier-fit['withheld_brier'][name])>1e-12:raise ArithmeticError('withheld calibration differs')
    radial=[];positives=0
    for record in fit['development_rows']:
        path=ROOT/record['source']
        if digest(path.read_bytes())!=record['sha256']:raise ArithmeticError('development call changed')
        row=json.loads(path.read_text())
        if record['independent_height'] is None:
            if row.get('new_point'):raise ArithmeticError('positive witness censored')
            continue
        seedpath=path.parent/'rank-input.json';packet=json.loads(seedpath.read_text())
        # The independently verified development endpoint must be exactly the
        # baseline prefix followed by this witness, up to its recorded sign.
        if packet['points'][-1]!=row['new_point']:raise ArithmeticError('development certificate point differs')
        curve=tuple(map(F,packet['curve']));basis=[tuple(map(F,p)) for p in packet['points'][:-1]]
        a,b=linear_combination(curve,basis,row['centre']);x,y=map(F,row['new_point'])
        A,B,C,D=map(F,row['mapping']['matrix']);hs=[]
        for sign in (1,-1):
            u=D*(sign*y+b)-B*(x-a);v=-C*(sign*y+b)+A*(x-a)
            den=lcm(u.denominator,v.denominator);m,n=int(u*den),int(v*den);g=gcd(m,n)
            hs.append(max(abs(m),abs(n))//g)
        if min(hs)!=record['independent_height']:raise ArithmeticError('independent radius replay differs')
        radial.append(min(hs));positives+=1
    if sorted(radial)!=fit['radial_heights']:raise ArithmeticError('radius distribution differs')
    result={'status':'PASS','corpus_anchors':anchors,'training_j_groups':len(train),'withheld_j_groups':len(test),
        'residue_leaves':leaves,'development_calls':len(fit['development_rows']),
        'certified_development_witnesses':positives,'cpu_seconds':time.process_time()-start,
        'fit_sha256':digest((OUT/'fit.json').read_bytes()),'source_sha256':digest(__import__('pathlib').Path(__file__).read_bytes()),
        'boundary':'Independently checked q partitions and exact rational differential weights; reconstructed mixture, withheld calibration and development witness radii. Uses the retained V2 independent rank replay, with exact endpoint membership checks; no point search.'}
    if (OUT/'training-replay.json').exists():raise FileExistsError('preserve replay')
    write(OUT/'training-replay.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
