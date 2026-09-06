#!/usr/bin/env python3
"""Independent multiplication norms and complete fourth-root component fields."""
import argparse
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,matrix,prime_range
import retrospective as r

HERE=Path(__file__).resolve().parent
SOURCE=r.OUT/'rank_jump_fourth_lift_on_relation_norms_v1.json'
OUTPUT=r.OUT/'rank_jump_fourth_lift_on_relation_verification_v1.json'


def compute():
    src=r.read(SOURCE)
    for path,sha in src['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'t');t=R.gen();q=R(src['fourth_cover']['form'])
    Z=PolynomialRing(QQ,'z');z=Z.gen();T=PolynomialRing(Z,'t');tt=T.gen()
    rows=[]
    for row in src['rows']:
        components=[];total=Z(1)
        for c in row['components']:
            f=R(c['parameter_factor']);n=int(f.degree());assert f.is_irreducible() and f.gcd(q)==1
            # Multiplication by q in the residue algebra, independent of resultant.
            cols=[[(q*t**j % f)[i] for i in range(n)] for j in range(n)]
            M=matrix(QQ,cols).transpose();norm=M.det();assert norm==QQ(c['fourth_value_norm'])
            norm_square=bool(norm.is_square());assert norm_square==c['norm_is_square']
            cp=Z(T(f.list()).resultant((z-tt)**2-T(q.list()))).monic()
            assert cp.degree()==2*n and cp.is_squarefree();total*=cp
            witnesses={}
            if not norm_square:
                if norm<0:witnesses={'place':'real','norm_sign':-1}
                else:
                    for p in prime_range(3,1000):
                        if norm.numerator()%p and norm.denominator()%p and not GF(p)(norm).is_square():
                            witnesses={'place':int(p),'norm_unit_residue':int(GF(p)(norm))};break
                    assert witnesses,'No bounded rational norm obstruction witness'
                # The quadratic extension of the degree-n field is a field.
                # Squarefree cp means t+sqrt(q) has 2n distinct conjugates.
                assert cp.is_irreducible();degrees=[2*n]
            else:
                assert n==1
                root=-f[0];a=q(root).sqrt()
                assert a*a==q(root) and cp==(z-(root+a))*(z-(root-a))
                degrees=[1,1];witnesses={'rational_parameter':str(root),'fourth_roots':[str(a),str(-a)]}
            components.append({'parameter_degree':n,'lifted_component_degrees':degrees,
                'multiplication_norm':str(norm),'norm_obstruction_or_lift':witnesses,
                'primitive_element':'t+sqrt(q_13109(t))','primitive_polynomial':list(map(str,cp.list()))})
        assert total.degree()==24 and total.is_squarefree()
        factors=sorted(int(f.degree()) for f,e in total.factor() for _ in range(int(e)))
        assert factors==sorted(d for c in components for d in c['lifted_component_degrees'])
        rows.append({'relation_index':row['relation_index'],'components':components,
            'total_lift_degree':24,'factor_degrees':factors,
            'total_primitive_polynomial':list(map(str,total.list())),
            'rational_lift_count':factors.count(1)})
    assert [row['factor_degrees'] for row in rows]==[[24],[1,1,22],[24]]
    return {'schema':'rank-jump.fourth-lift-on-relation-verification.v1','status':'PASS','rows':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (SOURCE,Path(__file__),HERE/'retrospective.py')},
        'boundary':'Complete pullbacks of the three already complete finite relation schemes, not complete rational points of the genus-17 carrier. The two rational components are opposite fourth-root signs at one parameter.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    print('PASS: fourth-lift component degrees',[row['factor_degrees'] for row in d['rows']])
