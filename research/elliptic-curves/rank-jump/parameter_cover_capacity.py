#!/usr/bin/env python3
"""Retrospective parameter-cover geometry on the original published K3."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r
import branch_blocks

PROTOCOL=Path(__file__).with_name('PARAMETER_COVER_CAPACITY_PROTOCOL.json')
MODEL=r.ROOT/'elkies-k3/data/fibrations/elkies_2026_published_r17_model.json'
QUARTICS=r.ROOT/'artifacts/generated-results/elkies-k3-r17-rank28-genus-one-bisection-pilot-v1.json'
INPUT=r.OUT/'rank_jump_parameter_cover_capacity_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_parameter_cover_capacity_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-parameter-cover-capacity-v1'


def capacity(root_rank,branch,bad):
    assert branch>0 and branch%2==0 and 0<=bad<=branch
    chi=2+branch//2;roots=root_rank+4*branch+bad;h11=10*chi
    return {'base_curve_genus':branch//2-1,'twist_chi':chi,'twist_root_rank':roots,
        'twist_h11':h11,'twist_geometric_MW_rank_upper_bound':h11-2-roots}


def export():
    model=r.read(MODEL);covers=[]
    for row in r.read(branch_blocks.INPUT)['covers']:
        covers.append({'label':row['label'],'kind':'retained_quadratic',
            'q':row['residual_chord']['q_coefficients']})
    for row in r.read(QUARTICS)['traces'][0]['targets']:
        covers.append({'label':row['target_label'],'kind':'oracle_quartic',
            'q':row['branch_polynomial_q_coefficients_low_to_high']})
    r.write_new(INPUT,{'schema':'rank-jump.parameter-cover-capacity-inputs.v1',
        'source_hashes':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (MODEL,branch_blocks.INPUT,QUARTICS)},
        'A':model['A_coefficients_low_to_high'],'B':model['B_coefficients_low_to_high'],'covers':covers})


def compute():
    from sage.all import QQ,PolynomialRing
    from sage.version import version
    data=r.read(INPUT);R=PolynomialRing(QQ,'u');u=R.gen()
    A=R(list(map(QQ,data['A'])));B=R(list(map(QQ,data['B'])))
    delta=-16*(4*A**3+27*B**2);c4=-48*A
    assert [A.degree(),B.degree(),delta.degree()]==[8,12,24]
    assert delta.gcd(delta.derivative())==1 and delta.gcd(c4)==1
    factor=delta.factor();assert all(e==1 for f,e in factor)
    factors=[{'degree':int(f.degree()),'coefficients_ascending':list(map(str,f.monic()))} for f,e in factor]
    covers=[]
    for old in data['covers']:
        q=R(list(map(QQ,old['q'])));b=int(q.degree())
        assert b in [2,4] and q.gcd(q.derivative())==1 and q.gcd(delta)==1
        # Explicit twisting invariants in the finite chart.
        At=A*q**2;Bt=B*q**3;Dt=-16*(4*At**3+27*Bt**2)
        assert Dt==delta*q**6 and -48*At==c4*q**2
        chi=2+b//2
        assert At.degree()<=4*chi and Bt.degree()<=6*chi and Dt.degree()==12*chi
        covers.append({'label':old['label'],'kind':old['kind'],'branch_points':b,
            'bad_branch_points':0,'branch_separable':True,'branch_avoids_singular_fibres':True,
            'twist_fibres':{'I1':24,'I0*':b},**capacity(0,b,0)})
    table=[]
    for family,root,generic in [('R17',0,17),('A1_MW16',1,16)]:
        for b in [2,4]:
            for bad in range(b+1):
                row=capacity(root,b,bad)
                table.append({'family':family,'root_rank':root,'original_generic_rank':generic,
                    'branch_points':b,'multiplicative_branch_points':bad,**row,
                    'pullback_geometric_rank_upper_bound':generic+row['twist_geometric_MW_rank_upper_bound']})
    return {'schema':'rank-jump.parameter-cover-capacity.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,INPUT)},
        'software':{'sage':version},'original_surface':{'finite_fibres':'24 I1','infinity':'smooth',
            'chi':2,'geometric_MW_rank':17,'Picard_rank':19,
            'rank_source':'retained published R17 lattice/K3 proof; not recomputed here'},
        'discriminant_coefficients_ascending':list(map(str,delta)),
        'discriminant_factorization_monic':factors,'smallest_singular_closed_point_degree':min(x['degree'] for x in factors),
        'retained_covers':covers,'capacity_table':table,
        'boundary':'Geometric generic rank upper bounds for specified quadratic twists. No new section, original-fibre rank bound, rationality or specialization-independence certificate.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker','check']);a=p.parse_args()
    if a.mode=='export':export()
    elif a.mode=='worker':r.write_new(WORK/'checkpoint-complete.json',compute())
    elif a.mode=='check':assert r.read(OUTPUT)==compute();print('PASS parameter-cover capacity')
    else:
        WORK.mkdir(parents=True,exist_ok=True);path=WORK/'checkpoint-complete.json'
        if not path.exists():
            with (WORK/'worker-complete.log').open('x') as log:
                try:
                    proc=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker'],
                        cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                    reason=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:reason='30-second timeout'
                if reason and not path.exists():r.write_new(path,{'status':'UNKNOWN','reason':reason})
        data=r.read(path);r.write_new(OUTPUT,data)
        print(data.get('status','PASS'),data.get('smallest_singular_closed_point_degree'),len(data.get('retained_covers',[])))
