#!/usr/bin/env python3
"""Independent polynomial reconstruction and dyadic lifting checks."""
import argparse
from itertools import product
from pathlib import Path
import retrospective as r
import residual_double_covers as original
import residual_alignment as aligned

OUTPUT=r.OUT/'rank_jump_residual_alignment_verification_v1.json'


def convolution(a,b):
    out=[r.F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def polynomial_form(coeff,coords):
    out=[r.F(0)]*5
    for (i,j),c in zip(original.PAIRS,coeff):
        for k,v in enumerate(convolution(coords[i],coords[j])):out[k]+=r.F(c)*v
    return out


def compute():
    from sage.all import QQ,Qp,PolynomialRing
    old=r.read(original.OUTPUT);data=r.read(aligned.OUTPUT);source=r.read(original.old.OUTPUT)
    for record in [old,data]:
        for name,digest in record['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==digest
    for saved,src in zip(old['quartics'],source['records'],strict=True):
        z=[list(map(r.F,v)) for v in saved['parametrization_coefficients']]
        assert polynomial_form(src['Q2'],z)==[0]*5
        assert polynomial_form(src['Q1'],z)==list(map(r.F,saved['quartic_coefficients']))
        base=list(map(r.F,saved['qfsolve_point']))
        assert sum(r.F(c)*base[i]*base[j] for (i,j),c in zip(original.PAIRS,src['Q2']))==0
    R=PolynomialRing(QQ,'x');x=R.gen();field=Qp(2,40);rows=[]
    for row in data['rows']:
        coeff=[list(map(r.F,q)) for q in row['quartic_coefficients']]
        expected=[[r.F(c)/d for c in q['quartic_coefficients']] for q,d in zip(old['quartics'],[100,25])]
        if row['alignment']=='swap_second':expected[1].reverse()
        assert coeff==expected
        pols=[sum(QQ(str(c))*x**(4-i) for i,c in enumerate(q)) for q in coeff]
        assert pols[0].gcd(pols[1])==1
        assert str(pols[0].resultant(pols[1]))==row['resultant']
        assert all(f.degree()==4 and f.gcd(f.derivative())==1 for f in pols)
        signs=[]
        for result in row['signs']:
            e=result['epsilon'];square_residues={y*y%8 for y in range(8)}
            survivors=[]
            # Independent exhaustive primitive pairs, without projective normalization.
            for s,t in product(range(8),repeat=2):
                if s%2==0 and t%2==0:continue
                if all((e*original.evaluate(q,s,t))%8 in square_residues for q in coeff):survivors.append([s,t])
            if result['Q2_solubility']=='NO':assert not survivors
            witness=result['Q2_witness']
            if witness:
                s,t=witness['parameter'];values=[e*original.evaluate(q,s,t) for q in coeff]
                assert list(map(str,values))==witness['values']
                assert all(field(QQ(str(v))).is_square() for v in values)
                assert result['Q2_solubility']=='YES'
            signs.append({'epsilon':e,'Q2_solubility':result['Q2_solubility'],
                'primitive_mod8_survivor_count':len(survivors),'p_adic_witness_verified':witness is not None})
        rows.append({'alignment':row['alignment'],'status':'PASS','signs':signs,
            'branch_points':8,'geometric_character_rank':2,
            'genus_by_Riemann_Hurwitz':1+2**0*(8-4),
            'quotient_character_genera':[1,1,3],'Jacobian_dimension':5,
            'signed_Jacobian_rank_difference_from_retained_E_ranks':2*(3-1)})
    return {'schema':'rank-jump.residual-alignment-verification.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),original.OUTPUT,aligned.OUTPUT)},
        'rows':rows,'boundary':'The elliptic ranks 3 and 1 are retained exact control results. No genus-three or genus-five Jacobian rank computation, full local-solubility claim, or rational-point existence claim is made.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();data=compute()
    if a.mode=='check':assert r.read(OUTPUT)==data;print('PASS independent residual polynomials and dyadic witnesses')
    else:r.write_new(OUTPUT,data)
