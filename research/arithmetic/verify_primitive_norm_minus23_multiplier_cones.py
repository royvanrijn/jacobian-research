#!/usr/bin/env python3
"""Exact C3 phase covers for the universal multiplier seeds16,19,39.

The retained C3 individual primitive criterion is the mathematical
input. This finite certificate checks every possible bank phase and
outside-bank active-count case; it is not a bounded integer scan.
"""
if not __debug__:
    raise RuntimeError('This certificate requires Python assertions enabled.')
from itertools import product
from pathlib import Path
from math import isqrt
import argparse,hashlib,json

BANK=(2,3,13)
SEEDS={16:((8,8),(3,13),(4,12)),
       19:((1,18),(3,16),(6,13)),
       39:((3,36),(12,27),(13,26))}

def factor_vector(n):
    result=[]
    for p in BANK:
        e=0
        while n%p==0:
            n//=p;e+=1
        result.append(e%3)
    if n!=1:
        raise ValueError('A seed summand has a factor outside the fixed bank.')
    return tuple(result)

def active_after(phase,vector,outside):
    return outside+sum((x+y)%3!=0 for x,y in zip(phase,vector))

def find_witness(pairs,phase,outside):
    for index,(u,v)in enumerate(pairs):
        ru=active_after(phase,factor_vector(u),outside)
        rv=active_after(phase,factor_vector(v),outside)
        if ru!=1 and rv!=1:
            return index,ru,rv
    return None

def prime_controls():
    witnesses={2:(1,0),3:(0,1),13:(2,1),29:(-2,3)}
    controls=[]
    for p,(u,v)in witnesses.items():
        assert pow(p,11,23)==1
        assert 2*u*u+u*v+3*v*v==p
        # Q=(x+y/2)^2+23y^2/4 bounds every principal representation.
        ybound=isqrt(4*p//23)
        represented=[]
        for y in range(-ybound,ybound+1):
            xbound=isqrt(p)+abs(y)+2
            for x in range(-xbound,xbound+1):
                if x*x+x*y+6*y*y==p:represented.append((x,y))
        assert not represented
        controls.append({'p':p,'split':True,'nonprincipal_form_witness':[u,v],
                         'principal_representations':represented})
    assert pow(19,11,23)==22
    return controls

def build_packet():
    seeds=[];table={tuple(e):{}for e in product(range(3),repeat=3)}
    for seed,pairs in SEEDS.items():
        assert all(u>0 and v>0 and u+v==seed for u,v in pairs)
        vectors=[(factor_vector(u),factor_vector(v))for u,v in pairs]
        rows=[]
        for phase in product(range(3),repeat=3):
            for outside in(0,1,2):
                witness=find_witness(pairs,phase,outside)
                assert witness is not None,(seed,phase,outside)
                index,ru,rv=witness
                assert ru!=1 and rv!=1
                if outside==2:assert ru>=2 and rv>=2
                label=chr(ord('A')+index)
                rows.append({'q_bank_phase':phase,
                             'outside_active_case':'>=2'if outside==2 else str(outside),
                             'pair_label':label,'u':pairs[index][0],'v':pairs[index][1],
                             'active_counts_at_case_minimum':[ru,rv]})
                if outside<2:table[phase][f'{seed}_R{outside}']=label
        deletion_controls=[]
        for removed in range(3):
            smaller=[pair for k,pair in enumerate(pairs)if k!=removed]
            missing=next(((phase,outside)for phase in product(range(3),repeat=3)
                          for outside in(0,1)if find_witness(smaller,phase,outside)is None),None)
            assert missing is not None
            deletion_controls.append({'removed_label':chr(ord('A')+removed),'uncovered_phase_case':missing})
        seeds.append({'seed':seed,'pairs':[{'label':chr(ord('A')+k),'u':u,'v':v,
             'u_bank_phase':vectors[k][0],'v_bank_phase':vectors[k][1]}
             for k,(u,v)in enumerate(pairs)],'complete_phase_cases':len(rows),
             'witnesses':rows,'pair_deletion_controls':deletion_controls})
    reverse=[]
    for phase in product(range(3),repeat=3):
        q0=1
        for p,e in zip(BANK,phase):q0*=p**e
        assert factor_vector(q0)==phase
        reverse.extend([{'bank_phase':phase,'outside_active_count':0,'q':q0},
                        {'bank_phase':phase,'outside_active_count':1,'q':29*q0}])
    return {'bank':BANK,'prime_controls':prime_controls(),'seeds':seeds,
        'compact_table':[{'phase':phase,**values}for phase,values in table.items()],
        'reverse_phase_realizations':reverse,
        'finite_cover_equivalence':'For the specified bank-smooth pair list, the complete R0/R1 phase cover is equivalent to one of those displayed pairs giving two primitive values after every multiplier q in C. Reverse realizes phase e by the bank product and R1 by its product with nonprincipal29. This does not claim all universal P+P membership must use the displayed list.',
        'outside_case_completeness':'All seed factors lie in the bank, so outside active primes are unchanged. R>=2 makes every final count at least2; exact R0 andR1 cover every remaining integer.',
        'ramified_scope':'Arbitrary q with no inert factor and v23(q)<=1. Every seed summand is prime to23, so both resulting summands retain that ramified valuation.',
        'scope':'Universal unbounded multiplier cones proved by a complete finite C3 phase reduction. No unrestricted exact -23 exception completion, external review or formal verification is asserted.'}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('primitive_norm_minus23_multiplier_cone_checks.json'))
    args=parser.parse_args()
    result=build_packet()
    result['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'seeds':list(SEEDS),
                      'phase_cases_per_seed':81,'output':str(args.output)}))

if __name__=='__main__':main()
