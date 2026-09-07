#!/usr/bin/env python3
"""Check residue obstructions, exact local witnesses and subset accounting."""
import argparse
from collections import Counter
from pathlib import Path
import retrospective as r

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_local_solubility_block_inputs_v1.json'
RESULT=r.OUT/'rank_jump_local_solubility_blocks_v1.json'
OUTPUT=r.OUT/'rank_jump_local_solubility_block_verification_v1.json'


def local_square(v,p):
    if not v:return True
    exponent=0
    while v%p==0:v//=p;exponent+=1
    modulus=8 if p==2 else p
    return exponent%2==0 and v%modulus in {x*x%modulus for x in range(1,modulus) if x%p}


def verify():
    inp=r.read(INPUT);out=r.read(RESULT);covers=inp['covers'];n=len(covers)
    for path,sha in out['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    allowed=[True]*(1<<n);proved=[True]*(1<<n);first={}
    witness_checks=0;residue_checks=0
    for row in out['places']:
        p,m=row['prime'],row['modulus']
        roots={x*x%m for x in range(m)}
        raw=set()
        for a,b in [(a,1) for a in range(m)]+[(1,b) for b in range(0,m,p)]:
            mask=0
            for k,c in enumerate(covers):
                e,d,f=c['form']
                if (e*b*b+d*a*b+f*a*a)%m in roots:mask|=1<<k
            raw.add(mask);residue_checks+=1
        maximal=row['maximal_possible_masks']
        assert all(x in raw for x in maximal)
        assert all(any(x&y==x for y in maximal) for x in raw)
        assert all(not any(x!=y and x&y==x for y in maximal) for x in maximal)
        positive=[]
        for w in row['maximal_proved_masks']:
            a,b=w['rational_base_point'];mask=w['mask'];positive.append(mask)
            for k,c in enumerate(covers):
                if mask>>k&1:
                    e,d,f=c['form'];assert local_square(e*b*b+d*a*b+f*a*a,p)
                    witness_checks+=1
        for mask in range(1<<n):
            yes=any(mask&x==mask for x in maximal)
            if not yes:first.setdefault(mask,p)
            allowed[mask]&=yes
            proved[mask]&=any(mask&x==mask for x in positive)
    for counts in out['subset_counts']:
        masks=[m for m in range(1<<n) if m.bit_count()==counts['size']]
        assert counts=={'size':counts['size'],'total':len(masks),
                        'obstructed':sum(not allowed[m] for m in masks),
                        'proved_at_tested_places':sum(proved[m] for m in masks),
                        'unknown':sum(allowed[m] and not proved[m] for m in masks)}
    minimal=[]
    for m in range(1,1<<n):
        indices=[i for i in range(n) if m>>i&1]
        if not allowed[m] and all(allowed[m^(1<<i)] for i in indices):
            minimal.append({'mask':m,'indices':indices,'size':len(indices),'first_obstruction_prime':first[m]})
    assert minimal==out['minimal_obstructed_subsets']
    for pair in out['pairs']:
        m=pair['mask'];assert m==sum(1<<i for i in pair['indices'])
        assert pair['observed_together']==any(g['mask']&m==m for g in inp['observed_groups'])
        assert pair['status']==('OBSTRUCTED' if not allowed[m] else 'PROVED_AT_TESTED_PLACES' if proved[m] else 'UNKNOWN')
        assert pair['first_obstruction_prime']==first.get(m)
    assert out['pair_counts']=={role:dict(Counter(x['status'] for x in out['pairs'] if x['observed_together']==flag)) for role,flag in [('observed',True),('cross_group',False)]}
    # A compact independently readable four-way obstruction over Q_23.
    quartet=[1,2,3,8];prime=23
    supports=[]
    for i in quartet:
        e,d,f=covers[i]['form'];squares={x*x%prime for x in range(prime)}
        support=[t for t in range(prime) if (e+d*t+f*t*t)%prime in squares]
        if f%prime in squares:support.append('infinity')
        supports.append(support)
    assert not set.intersection(*(set(x) for x in supports))
    triple_witnesses=[]
    for omitted,t in zip(quartet,[13,18,17,1]):
        roots=[]
        for i in quartet:
            if i==omitted:continue
            e,d,f=covers[i]['form'];v=(e+d*t+f*t*t)%prime
            root=next(u for u in range(1,prime) if u*u%prime==v)
            roots.append({'index':i,'root_mod_23':root})
        triple_witnesses.append({'omitted_index':omitted,'parameter_mod_23':t,'unit_roots':roots})
    return {'schema':'rank-jump.local-solubility-block-verification.v1','status':'PASS',
            'result_sha256':r.digest(RESULT.read_bytes()),'script_sha256':r.digest(Path(__file__).read_bytes()),
            'projective_residue_classes_checked':residue_checks,'individual_exact_local_witness_checks':witness_checks,
            'four_way_obstruction':{'indices':quartet,'prime':prime,'square_supports_on_P1':supports,'triple_unit_witnesses':triple_witnesses},
            'boundary':'Complete replay at the declared places and depths only. The displayed Q23 obstruction alone proves no rational simultaneous lift of those four covers.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args()
    result=verify()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert r.read(OUTPUT)==result
    print(result['status'],result['projective_residue_classes_checked'],result['individual_exact_local_witness_checks'])
