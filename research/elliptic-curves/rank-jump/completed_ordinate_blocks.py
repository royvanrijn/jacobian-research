#!/usr/bin/env python3
"""Exact invariant horizontal triples in six retained paired controls."""
import argparse
from itertools import combinations
from pathlib import Path
import retrospective as r
import bad_prime_support as bad

PROTOCOL=Path(__file__).with_name('COMPLETED_ORDINATE_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_completed_ordinate_blocks_v1.json'


def on_curve(model,P):
    a1,a2,a3,a4,a6=map(r.F,model);x,y=map(r.F,P)
    return y*y+a1*x*y+a3*y==x**3+a2*x*x+a4*x+a6


def change(model,points,u=2,shift=3,s=5,t=7):
    a1,a2,a3,a4,a6=map(r.F,model);u,shift,s,t=map(r.F,(u,shift,s,t))
    assert u
    transformed=[(a1+2*s)/u,(a2-s*a1+3*shift-s*s)/u**2,
        (a3+shift*a1+2*t)/u**3,
        (a4-s*a3+2*shift*a2-(t+shift*s)*a1+3*shift**2-2*s*t)/u**4,
        (a6+shift*a4+shift**2*a2+shift**3-t*a3-shift*t*a1-t*t)/u**6]
    mapped=[]
    for P in points:
        x,y=map(r.F,P);xp=(x-shift)/u**2;yp=(y-u*u*s*xp-t)/u**3
        assert on_curve(model,P) and on_curve(transformed,(xp,yp));mapped.append((xp,yp))
        assert 2*y+a1*x+a3==u**3*(2*yp+transformed[0]*xp+transformed[2])
    return transformed,mapped


def grouping(model,points):
    a1,_,a3,_,_=map(r.F,model);groups={}
    for i,P in enumerate(points):
        x,y=map(r.F,P);assert on_curve(model,P)
        groups.setdefault(abs(2*y+a1*x+a3),{}).setdefault(x,[]).append(i)
    return sorted([sorted(sum(group.values(),[])) for group in groups.values() if len(group)>1])


def horizontal(model,points,signatures=None,generic_count=0):
    short,ps=r.short(model,points);A,B=map(r.F,short[3:]);groups={}
    for i,P in enumerate(ps):
        x,y=map(r.F,P);groups.setdefault(abs(y),{}).setdefault(x,[]).append(i)
    rows=[];allresidual=[]
    for y,xs in sorted(groups.items()):
        if len(xs)<2:continue
        assert len(xs)<=3
        x1,x2=sorted(xs)[:2];x3=-x1-x2
        assert all(x**3+A*x+B==y*y for x in (x1,x2,x3))
        # Distinct x1,x2 give slope zero, so their sum is (x3,-y).
        indices=sorted(sum(xs.values(),[]));residual=None;lower=None
        if signatures is not None:
            generic=r.basis(signatures[:generic_count])
            residual=[r.reduce(signatures[i],generic) for i in indices]
            lower=r.rank(residual);allresidual+=residual
            assert lower<=2
        rows.append({'absolute_short_ordinate':str(y),'input_indices':indices,
            'distinct_x_count':len(xs),'short_x_coordinates':list(map(str,sorted(xs))),
            'third_short_x':str(x3),'third_point_already_listed':x3 in xs,
            'relation':'(x1,y)+(x2,y)+(x3,y)=O, counting a repeated intersection twice',
            'generic_indices_in_group':[i for i in indices if i<generic_count],
            'finite_quotient_signatures':residual,'quotient_rank_lower_bound':lower})
    return rows,r.rank(allresidual)


def calculate():
    protocol=r.read(PROTOCOL);rows=[]
    for index in protocol['production_cases']:
        row=bad.cases()[index];points=row['generic_points']+row['points'];g=len(row['generic_points'])
        assert len(points)<=protocol['limits']['max_input_points_per_case']
        profile,_,signatures=r.characterize(row)
        groups,lower=horizontal(row['model'],points,signatures,g)
        changed,mapped=change(row['model'],points)
        assert grouping(row['model'],points)==grouping(changed,mapped)
        # Manufacturing raw-y equality must not manufacture completed-ordinate equality.
        a1,_,a3,_,_=map(r.F,row['model'])
        i,j=next((i,j) for i,j in combinations(range(g),2)
            if r.F(points[i][0])!=r.F(points[j][0]) and
               abs(2*r.F(points[i][1])+a1*r.F(points[i][0])+a3)!=abs(2*r.F(points[j][1])+a1*r.F(points[j][0])+a3))
        P,Q=[tuple(map(r.F,points[k])) for k in (i,j)]
        slope=(Q[1]-P[1])/(Q[0]-P[0])
        shear,sheared=change(row['model'],[P,Q],u=1,shift=0,s=slope,t=0)
        assert sheared[0][1]==sheared[1][1]
        assert grouping(shear,sheared)==[]
        rows.append({'case_index':index,'id':row['id'],'generic_rank':g,
            'known_independent_subgroup_rank':profile['certified_independent_subgroup_rank_exact'],
            'observed_quotient_rank':profile['certified_independent_quotient_rank_exact'],
            'input_point_count':len(points),'horizontal_groups':groups,
            'group_quotient_union_rank_lower_bound':lower,
            'fixed_Weierstrass_change_invariance':'PASS',
            'raw_y_false_positive_control':{'generic_indices':[i,j],'slope':str(slope),
                'transformed_ainvariants':list(map(str,shear)),
                'transformed_points':[[str(a) for a in P] for P in sheared],
                'raw_ordinates_equal':True,'completed_ordinates_equal_up_to_sign':False},
            'boundary':'Only the supplied point list is tested. Finite quotient signatures give lower bounds, not exact global relations for arbitrary points.'})
    controls=[]
    for c,n in protocol['small_controls']:
        m=r.F(n*n-c*c-2,2);model=[0,m,0,-m-3,c*c];points=[(0,c),(-1,n),(2,n)]
        groups,_=horizontal(model,points,generic_count=1)
        assert len(groups)==1
        controls.append({'c':c,'n':n,'m':str(m),'groups':groups,
            'retained_subgroup_rank_from_prior_certificate':1 if (c,n)==(1,1) else 3,
            'retained_rank_modulo_R_from_prior_certificate':0 if (c,n)==(1,1) else 2})
    return {'schema':'rank-jump.completed-ordinate-blocks.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,r.INPUT,r.RESULT,
             r.OUT/'rank_jump_shared_value_soluble_block_completion_v1.json')},
        'production_rows':rows,'small_controls':controls,
        'classification':'retrospective exact arithmetic pattern; no prospective selection or rank upper bounds'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();data=calculate()
    if args.mode=='check':assert r.read(OUTPUT)==data;print('PASS completed-ordinate blocks and shear controls')
    else:r.write_new(OUTPUT,data)
    for row in data['production_rows']:print(row['id'],len(row['horizontal_groups']),row['group_quotient_union_rank_lower_bound'])
