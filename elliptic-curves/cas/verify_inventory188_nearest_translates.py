#!/usr/bin/env python3
"""Independent rational group-law, metric transport and coordinate replay."""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
import inventory188_nearest_translate_control as control
from half_lattice_pointed_sieve import linear_combination_python
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint, digest

OUT = control.ART/'inventory188_nearest_translate_replay_v1.json'


def matmul(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def transpose(a):
    return list(map(list,zip(*a)))


def determinant(a):
    a=[list(map(F,r)) for r in a];det=F(1)
    for j in range(len(a)):
        pivot=next((k for k in range(j,len(a)) if a[k][j]),None)
        if pivot is None:return F(0)
        if pivot!=j:a[pivot],a[j]=a[j],a[pivot];det=-det
        q=a[j][j];det*=q
        for k in range(j+1,len(a)):
            ratio=a[k][j]/q
            for t in range(j+1,len(a)):a[k][t]-=ratio*a[j][t]
    return det


def positive(a):
    a=[list(map(F,r)) for r in a]
    assert a==transpose(a)
    for j in range(len(a)):
        assert a[j][j]>0
        for k in range(j+1,len(a)):
            for t in range(j+1,len(a)):
                a[k][t]-=a[k][j]*a[j][t]/a[j][j]


def expected():
    p=control.protocol();result=control.read(control.OUT);public=control.read(control.PUBLIC)
    assert result['status']=='PASS_EXACT_TRANSLATE_COORDINATES'
    assert result['protocol_sha256']==control.hashed(control.D/'protocol.json')
    assert result['charts']==control.charts()
    model=tuple(map(F,public['curve']));old=[tuple(map(F,P)) for P in public['points'][:27]]
    proofs=control.read(control.PREVIOUS)['witness_rank_certificates']
    observations=points_checked=0;visible=0;minima=[]
    assert [r['public_index'] for r in result['rows']]==p['public_indices']
    for row in result['rows']:
        index=row['public_index'];P=tuple(map(F,public['transported_public_points'][index]))
        proof=proofs[str(index)]
        replay=checked_rank(model,old+[P],[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
        assert digest(replay)==digest(proof) and replay['rank_lower_bound']==28
        H=row['rounded_full_gram'];assert len(H)==28 and all(len(r)==28 for r in H)
        with localcontext() as context:
            context.prec=110
            assert [[int((Decimal(x)*10**6).to_integral_value()) for x in r] for r in row['metric_gram']]==H
        positive(H)
        G=[r[:27] for r in H[:27]];U=row['change_of_basis']
        assert abs(determinant(U))==1
        C=matmul(matmul(U,G),transpose(U));assert C==row['reduced_gram']
        t=list(map(F,row['target_reduced_coefficients']));b=matmul(U,[[r[27]] for r in H[:27]])
        assert all(sum(C[i][j]*t[j] for j in range(27))==-b[i][0] for i in range(27))
        q=row['nearest_reduced_word'];offset=[F(x)-y for x,y in zip(q,t)]
        distance=sum(offset[i]*C[i][j]*offset[j] for i in range(27) for j in range(27))
        assert distance==F(row['exact_rounded_distance'])
        choices=[tuple([0]*27),tuple(q)]
        for i in range(27):
            for sign in (-1,1):
                choice=q.copy();choice[i]+=sign;choices.append(tuple(choice))
        choices=list(dict.fromkeys(choices))
        assert [tuple(r['reduced_word']) for r in row['representatives']]==choices
        best=[]
        for j,r in enumerate(row['representatives']):
            word=matmul([r['reduced_word']],U)[0];assert word==r['old_word']
            Q=linear_combination_python(model,old+[P],word+[1])
            assert Q and list(map(str,Q))==r['point'];points_checked+=1
            vv=[]
            for c in result['charts']:
                x,y=map(F,c['base_point']);a,b,c0,d=map(F,c['matrix'])
                assert a*d-b*c0 and Q[0]!=x
                for sign in (-1,1):
                    slope=(sign*Q[1]+y)/(Q[0]-x);num=d*slope-b;den=a-c0*slope
                    if den:
                        s=num/den;coordinate=[str(s.numerator),str(s.denominator)]
                        height=max(abs(s.numerator),s.denominator);infinity=False
                    else:coordinate=['1','0'];height=1;infinity=True
                    vv.append({'arm':c['arm'],'chart':c['index'],'sign':sign,'coordinate':coordinate,
                               'height':height,'infinity':infinity,'within_height_or_infinity':infinity or height<=p['height']})
            assert vv==r['observations'];observations+=len(vv);visible+=sum(v['within_height_or_infinity'] for v in vv)
            minimum=min(vv,key=lambda x:(x['height'],x['arm'],x['chart'],x['sign']));assert minimum==r['best']
            best.append({'representative_index':j,**minimum})
        minimum=min(best,key=lambda x:(x['height'],x['representative_index'],x['arm'],x['chart'],x['sign']))
        assert minimum==row['best'];minima.append({'public_index':index,**minimum})
    assert observations==result['observation_count'] and visible==result['visible_count']
    return {'schema':'elliptic-curves.inventory188-nearest-translates-replay.v1','status':'PASS',
            'sources':{str(path.relative_to(control.ROOT)):control.hashed(path) for path in
                       [control.OUT,control.D/'protocol.json',control.CAS/'verify_inventory188_nearest_translates.py']},
            'independent_rational_group_sums':points_checked,'independent_coordinates':observations,
            'within_height_or_infinity':visible,'minima':minima,
            'scope':'Exact metric transports, finite independence, group sums and coordinates. Floating CVP optimality is not certified; oracle-assisted coverage is not a prospective recovery.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
    r=expected()
    if args.check:assert r==control.read(OUT)
    else:
        if OUT.exists():raise FileExistsError('preserve independent replay')
        checkpoint(OUT,r)
    print(r['status'],r['independent_rational_group_sums'],r['independent_coordinates'],r['within_height_or_infinity'],r['minima'])
