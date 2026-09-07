#!/usr/bin/env python3
"""Sylvester determinants and Sturm certificates for the fixed lift criterion."""
import argparse
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,matrix,prod
import retrospective as r

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_soluble_quartet_compression_inputs_v1.json'
SOURCE=r.OUT/'rank_jump_collision_prime_lift_v1.json'
OUTPUT=r.OUT/'rank_jump_collision_prime_lift_verification_v1.json'


def compute():
    inp=r.read(INPUT);src=r.read(SOURCE)
    for path,sha in src['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'t');t=R.gen();rows=[]
    for case,record in zip(inp['cases'],src['rows'],strict=True):
        assert record['execution']['status']=='COMPLETE';old=record['result'];qs=[R(c['form']) for c in case['covers']]
        primes=set();identities=[]
        for pair in old['pair_resultants']:
            i,j=pair['indices'];a,b,c=qs[i][2],qs[i][1],qs[i][0];d,e,f=qs[j][2],qs[j][1],qs[j][0]
            determinant=matrix(QQ,[[a,b,c,0],[0,a,b,c],[d,e,f,0],[0,d,e,f]]).det()
            assert determinant==QQ(pair['resultant'])
            factors=[(ZZ(p),e) for p,e in pair['prime_factorization']]
            assert all(p.is_prime(proof=True) for p,e in factors)
            assert prod(p**e for p,e in factors)==abs(determinant)
            primes.update(int(p) for p,e in factors)
            charts=[]
            for pols in ([qs[i],qs[j]],[R(qs[i].list()[::-1]),R(qs[j].list()[::-1])]):
                g,s,u=pols[0].xgcd(pols[1]);assert g==1
                s*=determinant;u*=determinant
                assert s*pols[0]+u*pols[1]==determinant
                assert all(x in ZZ for x in s.list()+u.list())
                charts.append({'bezout_left':list(map(str,s.list())),'bezout_right':list(map(str,u.list()))})
            identities.append({'indices':[i,j],'two_chart_Bezout_identities':charts})
        assert list(map(str,sorted(primes)))==old['collision_primes']
        F=prod(qs);sturm=[F,F.derivative()]
        while True:
            remainder=-(sturm[-2]%sturm[-1])
            if not remainder:break
            sturm.append(remainder)
        def sign(v):return 1 if v>0 else -1 if v<0 else 0
        def variations(signs):
            ss=[s for s in signs if s];return sum(ss[i]!=ss[i-1] for i in range(1,len(ss)))
        def variation_at(a):return variations([sign(pol(a)) for pol in sturm])
        nreal=variations([sign(pol.leading_coefficient())*(-1)**pol.degree() for pol in sturm])-variations([sign(pol.leading_coefficient()) for pol in sturm])
        intervals=F.real_root_intervals();assert len(intervals)==nreal
        verified_roots=[];seen=[0]*4;last=None
        for (lo,hi),multiplicity in intervals:
            assert multiplicity==1 and lo<hi and F(lo)!=0 and F(hi)!=0
            assert last is None or last<=lo;last=hi
            assert variation_at(lo)-variation_at(hi)==1
            labels=[i for i,q in enumerate(qs) if q(lo)*q(hi)<0];assert len(labels)==1
            i=labels[0];verified_roots.append({'cover_index':i,'root_index':seen[i],'isolating_interval':[str(lo),str(hi)]});seen[i]+=1
        assert [(x['cover_index'],x['root_index']) for x in verified_roots]==[(x['cover_index'],x['root_index']) for x in old['real_branch_order']]
        signs=[sign(q[2]) for q in qs];cells=[];branches=[]
        def cell():return {'signs':list(signs),'product_positive':prod(signs)>0,'all_native_positive':all(s>0 for s in signs)}
        cells.append(cell())
        for root in verified_roots:
            i=root['cover_index'];branches.append(all(signs[j]>0 for j in range(4) if j!=i));signs[i]*=-1;cells.append(cell())
        assert cells==old['real_open_cells'] and branches==[x['other_native_values_positive'] for x in old['real_branch_order']]
        onto=all(not c['product_positive'] or c['all_native_positive'] for c in cells) and all(branches)
        onto=onto and (not old['real_product_points_at_infinity'] or old['real_native_lifts_at_infinity'])
        assert onto==old['real_lift_map_surjective']
        if case['observed_parameter'] is not None:
            param=QQ(case['observed_parameter']);a,b=param.numerator(),param.denominator()
            values=[ZZ(q[0]*b*b+q[1]*a*b+q[2]*a*a) for q in qs]
            assert all(v>0 and v.is_square() for v in values)
            assert all(v.valuation(p)%2==0 for p in primes for v in values)
        rows.append({'id':case['id'],'collision_prime_count':len(primes),'collision_primes':old['collision_primes'],
            'real_lift_map_surjective':bool(onto),'Sturm_isolating_intervals':verified_roots,
            'resultant_Bezout_certificates':identities,'finite_valuation_criterion_needs_sign_test':not onto})
    assert [r['collision_prime_count'] for r in rows]==[18,23,25]
    assert [r['real_lift_map_surjective'] for r in rows]==[True,True,False]
    return {'schema':'rank-jump.collision-prime-lift-verification.v1','status':'PASS','rows':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,SOURCE,Path(__file__),HERE/'retrospective.py')},
        'boundary':'The note proves the valuation-support implication for arbitrary rational product points. Sturm and Bezout identities verify the specific input geometry; no rational-point existence or new rank follows from the finite criterion alone.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    print('PASS: collision supports18,23,25; real surjectivity true,true,false')
