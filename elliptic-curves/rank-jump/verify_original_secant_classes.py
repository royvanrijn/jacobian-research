#!/usr/bin/env python3
"""Independent Sage rational tests, including every pair without local filtering."""
import argparse
from itertools import combinations
from pathlib import Path
import retrospective as r
import original_secant_classes as secants

OUTPUT=r.OUT/'rank_jump_original_secant_class_verification_v1.json'


def compute():
    from sage.all import QQ,GF,matrix
    import bad_prime_support as bad
    inp=r.read(secants.INPUT);out=r.read(secants.OUTPUT)
    assert inp['source_sha256']==r.digest(r.INPUT.read_bytes()) and out['bindings']==secants.bindings()
    rows=[]
    for data,result,source in zip(inp['rows'],out['rows'],bad.cases(),strict=True):
        assert result['status']=='PASS' and data['case_index']==result['case_index']
        assert r.short(source['model'],source['generic_points'])==(data['short_model'],data['generic_points'])
        A,B=map(QQ,data['short_model'][3:]);points=[tuple(map(QQ,P)) for P in data['generic_points']]
        expected=[([i,j],sign) for i,j in combinations(range(len(points)),2) for sign in [1,-1]]
        assert [(x['pair'],x['relative_sign']) for x in result['rows']]==expected
        values=[];sigs=[]
        primes=r.read(secants.PROTOCOL)['limits']['fingerprint_primes']
        for row in result['rows']:
            assert row['status']=='PASS'
            i,j=row['pair'];a,p=points[i];b,q=points[j];q*=row['relative_sign']
            assert p*p==a**3+A*a+B and q*q==b**3+A*b+B
            # Determinant formula, independently of the producer's slope formula.
            x0=(a*q-b*p)/(q-p);C=x0**3+A*x0+B
            assert x0==QQ(row['x0']) and C==QQ(row['C']) and C
            assert bool(C.is_square())==row['rationally_soluble']==False
            bits=[int(C<0)]
            for prime in primes:
                v=C.valuation(prime);unit=C/QQ(prime)**v
                bits.extend([int(v%2),int(not GF(prime)(unit).is_square())])
            sig=sum(bit<<i for i,bit in enumerate(bits))
            assert sig==row['local_signature'];sigs.append(sig);values.append(C)
        # Independent exhaustive class separation: no local-signature prefilter.
        comparisons=0
        for a,b in combinations(values,2):
            assert not (a/b).is_square();comparisons+=1
        assert len(result['classes'])==len(values)
        for i,group in enumerate(result['classes']):
            assert group=={'class_index':i,'representative_C':str(values[i]),
                'local_signature':sigs[i],'members':[i]}
        selected=[];basis=[]
        for i,s in enumerate(sigs):
            if r.rank(basis+[s])>len(basis):selected.append(i);basis.append(s)
        M=matrix(GF(2),[[(s>>i)&1 for i in range(1+2*len(primes))] for s in basis])
        assert M.rank()==len(basis)
        profile,_,_=r.characterize(source)
        rows.append({'case_index':data['case_index'],'id':source['id'],'status':'PASS',
            'original_generic_rank':len(points),'original_known_independent_rank':profile['certified_independent_subgroup_rank_exact'],
            'original_observed_quotient_rank':profile['certified_independent_quotient_rank_exact'],
            'secants':len(values),'rationally_soluble':0,'repeated_quadratic_classes':0,
            'independent_pairwise_square_tests':comparisons,
            'quadratic_character_rank_lower_bound':len(basis),
            'independent_local_signature_row_indices':selected,
            'compositum_degree_lower_bound':2**len(basis)})
        print('verified',data['case_index'],len(values),'classes; character rank >=',len(basis),flush=True)
    assert sum(x['secants'] for x in rows)==1504
    return {'schema':'rank-jump.original-secant-class-verification.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),secants.INPUT,secants.OUTPUT,secants.PROTOCOL)},
        'rows':rows,'boundary':'Character rank is a certified lower bound, not an exact rank. No twist Mordell-Weil rank or absence theorem outside the fixed secant dictionary.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();data=compute()
    if args.mode=='check':assert r.read(OUTPUT)==data;print('PASS independent original secant classes')
    else:r.write_new(OUTPUT,data)
