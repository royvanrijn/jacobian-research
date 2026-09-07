#!/usr/bin/env python3
"""Independent resultant pairing and permutation-module verification."""
import argparse
from itertools import permutations
from pathlib import Path
import retrospective as r
import residual_quartic_incidence as run

OUTPUT=r.OUT/'rank_jump_residual_quartic_incidence_verification_v1.json'


def compute():
    import sympy as s
    data=r.read(run.INPUT);out=r.read(run.OUTPUT)
    for obj in (data,out):
        for name,sha in obj['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    t=s.Symbol('t')
    poly=lambda cs:s.Poly.from_list([s.Rational(c) for c in cs[::-1]],t,domain=s.QQ)
    A=poly(data['A']);B=poly(data['B']);pts=[(poly(z['x']),poly(z['y'])) for z in data['sections']]
    rows=[0]*17;pairdata={(d['i'],d['j']):d for d in out['pairs']}
    for i,(x,y) in enumerate(pts):
        assert y*y==x**3+A*x+B and x.degree()==4 and y.degree()==6
        assert s.gcd(y,y.diff()).degree()==0 and s.gcd(y,3*x*x+A).degree()==0
        for j in range(i):
            xx,yy=pts[j];dx=x-xx
            assert dx.degree()==4 and s.gcd(y,yy).degree()==0
            # Common infinity terms cancel the leading-y normalization since deg(dx)=4.
            numerator=s.resultant(yy,dx);denominator=s.resultant(y,-dx)
            assert denominator!=0 and numerator in (denominator,-denominator)
            bit=int(numerator==-denominator)
            assert bit==pairdata[i,j]['weil_pairing_bit']==pairdata[i,j]['theta_translate_h0']%2
            rows[i]|=bit<<j;rows[j]|=bit<<i
    assert rows==out['weil_pairing_rows'] and r.rank(rows)==16
    radical=out['generic_radical_mask']
    def b(a,c):return sum((rows[i]&c).bit_count() for i in range(17) if (a>>i)&1)%2
    def q(a):return sum((rows[i]&a&((1<<i)-1)).bit_count() for i in range(17) if (a>>i)&1)%2
    assert radical!=0 and q(radical)==0 and all(b(radical,1<<i)==0 for i in range(17))
    # Independent symplectic elimination, rather than enumeration of 65536 vectors.
    basis=[1<<i for i in range(1,17)];arf=0;hyperbolic=[]
    while basis:
        a=basis.pop(0);k=next(i for i,c in enumerate(basis) if b(a,c));c=basis.pop(k)
        arf^=q(a)*q(c);hyperbolic.append([a,c,q(a),q(c)])
        basis=[v^(a if b(v,c) else 0)^(c if b(v,a) else 0) for v in basis]
    assert arf==1 and out['nondegenerate_generic_zero_count']==2**15-2**7
    group=out['residual_group'];letters=group['quartic_letters'];rad4=group['chosen_radical']
    assert r.rank(letters)==4 and letters[0]^letters[1]^letters[2]^letters[3]==rad4
    coordinate={0:0}
    for i,v in enumerate(letters):
        for x,mask in list(coordinate.items()):coordinate[x^v]=mask^(1<<i)
    assert len(coordinate)==16
    models=[]
    for row in group['elements']:
        p=row['permutation'];cols=row['columns'];fixed=[]
        for v,mask in coordinate.items():
            image=sum(((mask>>i)&1)<<p[i] for i in range(4))
            actual=0
            for i,col in enumerate(cols):
                if (v>>i)&1:actual^=col
            assert coordinate[actual]==image
            if image==mask:fixed.append(mask)
            assert run.q4(v)==(mask.bit_count()*(mask.bit_count()-1)//2)%2
        assert r.rank(fixed)==row['fixed_dimension']
        assert row['trace_mod_two']==sum(p[i]==i for i in range(4))%2
        models.append(p)
    assert sorted(map(tuple,models))==list(permutations(range(4)))
    # Rebuild every subgroup in the permutation model independently.
    ident=(0,1,2,3);allp=list(permutations(range(4)))
    def close(h,p):
        gens=list(h)+[p];result={ident};changed=True
        while changed:
            changed=False
            for a in list(result):
                for g in gens:
                    v=tuple(a[g[i]] for i in range(4))
                    if v not in result:result.add(v);changed=True
        return frozenset(result)
    groups={frozenset([ident])};pending=list(groups)
    while pending:
        h=pending.pop()
        for p in allp:
            if p in h:continue
            k=close(h,p)
            if k not in groups:groups.add(k);pending.append(k)
    assert len(groups)==30
    histogram={}
    for h in groups:
        oddtrace=any(sum(p[i]==i for i in range(4))%2 for p in h)
        transposition=any(sum(p[i]!=i for i in range(4))==2 for p in h)
        if not(oddtrace and transposition):continue
        fixedletters=[i for i in range(4) if all(p[i]==i for p in h)]
        fixedvectors=[m for m in range(16) if all(sum(((m>>i)&1)<<p[i] for i in range(4))==m for p in h)]
        assert len(h) in (6,24) and r.rank(fixedvectors)==1+len(fixedletters)
        histogram[str(len(h))]=histogram.get(str(len(h)),0)+1
    assert histogram=={'6':4,'24':1}
    return {'schema':'rank-jump.residual-quartic-incidence-verification.v1','status':'PASS',
        'generic_identities_and_simple_divisors_verified':17,'independent_resultant_pairings':136,
        'symplectic_elimination_pairs':hyperbolic,'generic_complement_arf':arf,
        'residual_module_is_quartic_permutation_module':True,'permutations_verified':24,
        'subgroups_verified':30,'admissible_subgroup_orders':histogram,
        'extra_global_class_constructed':False,'coefficient_quartic_constructed':False,
        'bindings':run.bindings([Path(__file__),run.INPUT,run.OUTPUT,Path(run.__file__),Path(r.__file__)])}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS independent resultant pairings, Arf invariant and quartic permutation representation')
