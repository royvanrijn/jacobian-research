#!/usr/bin/env python3
"""Certify an elementary ordinary class-group direct factor containing the new block."""
import argparse
from pathlib import Path
import retrospective as r
import constructed_class_half_ideal as half
ARTIN=r.OUT/'rank_jump_constructed_class_direct_artin_v1.json'
COMPACT=r.OUT/'rank_jump_constructed_class_compaction_v1.json'
OUTPUT=r.OUT/'rank_jump_constructed_half_block_verification_v1.json'

def compute():
    from sage.all import QQ,GF,PolynomialRing,pari,matrix,vector
    h=r.read(half.OUTPUT);d=r.read(ARTIN);c=r.read(COMPACT);ref=r.read(half.REFERENCE)
    for source in [h,d,c]:
        for path,sha in source['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'z');f=R(ref['cubic_ascending']);nf=pari.nfinit([pari(f),ref['S_finite']])
    gammas=[pari.Mod(pari(R(g['beta_ascending'])),pari(f)) for g in ref['generic_classes']]
    squares=[]
    for i,record in enumerate(h['columns']):
        if i<6:
            beta=pari.Mod(1,pari(f))
            for label in record['factor_labels']:assert label['kind']=='generic';beta*=gammas[label['index']]
            for step in record['half_ideal_reduction_steps']:
                a=pari.nfbasistoalg(nf,pari(step['principal_multiplier_GP']));beta/=a*a
        else:beta=pari.Mod(pari(R(c['cases'][i-6]['beta_ascending'])),pari(f))
        J=pari(record['final_reduced_half_ideal_hnf'])
        assert pari.idealpow(nf,J,2)==pari.idealhnf(nf,beta)
        squares.append({'column':i,'ideal_hnf':str(J),'square_generator_ascending':[str(pari.lift(beta).polcoef(j)) for j in range(3)]})
    A=matrix(GF(2),[row['artin_bits'] for row in d['columns']]);assert A.rank()==8 and A[:6].rank()==6
    duals=[]
    for i in range(8):
        target=vector(GF(2),[int(i==j) for j in range(8)]);x=A.solve_right(target);assert A*x==target
        duals.append(list(map(int,x)))
    new_self=[[d['columns'][j]['artin_bits'][i] for j in [6,7]] for i in [10,11]]
    assert matrix(GF(2),new_self).det()==1
    return {'schema':'rank-jump.constructed-half-block-verification.v1','status':'PASS','square_ideal_certificates':squares,
        'ordinary_class_group_elementary_direct_factor_rank':8,'generic_strict_factor_rank':6,'additional_factor_rank':2,
        'dual_character_coefficient_rows':duals,'new_two_by_two_self_artin_matrix':new_self,
        'proof':'The eight displayed ideals have principal squares. The eight displayed combinations of certified unramified characters evaluate as the identity matrix on them. Thus their order-two subgroup is (Z/2)^8 and the Artin map supplies a retraction from Cl(K). The first six ideals are generic strict images; the last two form an additional elementary direct factor relative to them.',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),half.OUTPUT,ARTIN,COMPACT,half.REFERENCE]},
        'boundary':'An explicit ordinary ideal-class direct factor, not the whole class group. The self-Artin matrix is not Cassels-Tate; its nonzero diagonal is not a solubility obstruction.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS elementary class-group factor: generic6 plus constructed2')
