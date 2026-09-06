#!/usr/bin/env python3
"""Simultaneous local square conditions; no rational-point or parameter search."""
import argparse
from collections import Counter
import gzip
import json
from math import gcd, isqrt
from pathlib import Path
import retrospective as r
import solubility_first as s

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'LOCAL_SOLUBILITY_BLOCK_PROTOCOL.json'
INPUT=r.OUT/'rank_jump_local_solubility_block_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_local_solubility_blocks_v1.json'


def qp_square(n,p):
    if not n:return True
    v=0
    while n%p==0:n//=p;v+=1
    return v%2==0 and (n%8==1 if p==2 else pow(n%p,(p-1)//2,p)==1)


def mask_of(indices):return sum(1<<i for i in indices)
def members(mask,n):return [i for i in range(n) if mask>>i&1]
def value(q,T,Z):return q[0]*Z*Z+q[1]*T*Z+q[2]*T*T


def capture():
    old=r.read(s.OUTPUT);lifts=r.read(s.INPUT)['split_lift_maps']
    geometry=json.loads(gzip.decompress(s.GEOMETRY.read_bytes()))
    covers=[]
    for c in geometry:
        if c['label'] not in lifts:continue
        content=gcd(*c['integer_quadratic']);root=isqrt(content)
        assert root*root==content
        covers.append({'label':c['label'],'form':[x//content for x in c['integer_quadratic']],
                       'removed_rational_square_root':str(r.F(root,c['denominator_scale']))})
    assert len(covers)==14
    labels=[c['label'] for c in covers]
    groups=[{'source_id':row['source_id'], 'mask':mask_of(labels.index(h['label']) for h in row['nonzero_square_hits']),
             'published_parameter':row['published_parameter']}
            for row in old['rows'] if row['split_cover_count']]
    r.write_new(INPUT,{'schema':'rank-jump.local-solubility-block-inputs.v1','covers':covers,'observed_groups':groups,
                       'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,s.OUTPUT,s.INPUT,s.GEOMETRY)},
                       'boundary':'Successful-cover labels are retrospective. No ranks or exceptional point coordinates enter local arithmetic.'})


def downward(masks,n):
    yes=bytearray(1<<n)
    for mask in masks:yes[mask]=1
    for bit in range(n):
        b=1<<bit
        for mask in range(1<<n):
            if mask&b and yes[mask]:yes[mask^b]=1
    return yes


def maximal(masks):
    selected=[]
    for m in sorted(set(masks),key=lambda x:(-x.bit_count(),x)):
        if not any(m&z==m for z in selected):selected.append(m)
    return sorted(selected)


def place(covers,p):
    modulus=32 if p==2 else p*p
    squares={u*u%modulus for u in range(modulus)}
    projective=[(t,1) for t in range(modulus)]+[(1,p*z) for z in range(modulus//p)]
    possible=set();witnesses={}
    for T,Z in projective:
        vals=[value(c['form'],T,Z) for c in covers]
        maybe=mask_of(i for i,v in enumerate(vals) if v%modulus in squares)
        proven=mask_of(i for i,v in enumerate(vals) if qp_square(v,p))
        assert proven&maybe==proven
        possible.add(maybe)
        witnesses.setdefault(proven,[T,Z])
    poss=maximal(possible);prov=maximal(witnesses)
    return {'prime':p,'modulus':modulus,'projective_residue_classes':len(projective),
            'maximal_possible_masks':poss,
            'maximal_proved_masks':[{'mask':m,'rational_base_point':witnesses[m]} for m in prov]}


def compute():
    inp=r.read(INPUT);protocol=r.read(PROTOCOL)
    for path,sha in inp['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    covers=inp['covers'];n=len(covers);places=[]
    allowed=bytearray([1])*(1<<n);proven=bytearray([1])*(1<<n)
    first_obstruction={}
    for p in protocol['local_places']:
        row=place(covers,p);places.append(row)
        a=downward(row['maximal_possible_masks'],n)
        b=downward([x['mask'] for x in row['maximal_proved_masks']],n)
        for mask in range(1<<n):
            if not a[mask]:first_obstruction.setdefault(mask,p)
            allowed[mask]&=a[mask];proven[mask]&=b[mask]
        print('completed p',p,flush=True)
    groups=inp['observed_groups']
    for group in groups:
        t=r.F(group['published_parameter'])
        for i in members(group['mask'],n):
            v=value(covers[i]['form'],t.numerator,t.denominator)
            assert v>0 and isqrt(v)**2==v
        # Enumeration must not declare a known rational group locally impossible.
        assert allowed[group['mask']]
    pairs=[]
    for i in range(n):
        for j in range(i+1,n):
            m=(1<<i)|(1<<j)
            pairs.append({'indices':[i,j],'mask':m,
                          'observed_together':any(g['mask']&m==m for g in groups),
                          'status':'OBSTRUCTED' if not allowed[m] else 'PROVED_AT_TESTED_PLACES' if proven[m] else 'UNKNOWN',
                          'first_obstruction_prime':first_obstruction.get(m)})
    minimal=[]
    for m in range(1,1<<n):
        if not allowed[m] and all(allowed[m^(1<<i)] for i in members(m,n)):
            minimal.append({'mask':m,'indices':members(m,n),'size':m.bit_count(),
                            'first_obstruction_prime':first_obstruction[m]})
    summary=[]
    for size in range(n+1):
        masks=[m for m in range(1<<n) if m.bit_count()==size]
        summary.append({'size':size,'total':len(masks),'obstructed':sum(not allowed[m] for m in masks),
                        'proved_at_tested_places':sum(proven[m] for m in masks),
                        'unknown':sum(allowed[m] and not proven[m] for m in masks)})
    return {'schema':'rank-jump.local-solubility-blocks.v1','status':'PASS',
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
            'labels':[c['label'] for c in covers], 'places':places,'observed_groups':groups,'pairs':pairs,
            'pair_counts':{role:dict(Counter(x['status'] for x in pairs if x['observed_together']==flag)) for role,flag in [('observed',True),('cross_group',False)]},
            'subset_counts':summary,'minimal_obstructed_subsets':minimal,
            'boundary':protocol['boundary']}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','build','check']);a=p.parse_args()
    if a.mode=='capture':capture()
    else:
        result=compute()
        if a.mode=='build':r.write_new(OUTPUT,result)
        else:assert r.read(OUTPUT)==result
        print(result['pair_counts']);print(result['subset_counts'])
        print('minimal obstruction counts',Counter(x['size'] for x in result['minimal_obstructed_subsets']))
