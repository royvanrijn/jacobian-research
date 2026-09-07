#!/usr/bin/env python3
"""Replay rational bisection maps and certify the reducible residuals as generic sections."""
import argparse
from pathlib import Path
import retrospective as r
import curve302_conic_triple_bisection as source
OUTPUT=r.OUT/'rank_jump_curve302_conic_triple_bisection_verification_v1.json'

def compute():
    from sage.all import QQ,PolynomialRing,matrix,vector,prod,gcd
    d=r.read(source.INPUT);out=r.read(source.OUTPUT)
    for path,sha in out['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,['X','Y','Z','W']);X,Y,Z,W=R.gens()
    def decode(rows):return sum(QQ(a['coefficient'])*prod(v**e for v,e in zip(R.gens(),a['powers'])) for a in rows)
    F=decode(d['quartic']);L=decode(d['residual_line'])
    B=PolynomialRing(QQ,['s0','s1']);s0,s1=B.gens();P=PolynomialRing(QQ,'t');t=P.gen()
    classes=[];groups=[];verified=[];sections=[]
    def parameter_map(curve,degree):
        common=gcd(curve);curve=[q//common for q in curve]
        assert F(*curve)==0 and all(q.total_degree()==degree for q in curve)
        numerator=curve[3];denominator=B(L(*curve));base=gcd(numerator,denominator)
        return B(numerator//base),B(denominator//base)
    for row in out['rows']:
        if row['status']=='EXPLICIT_RATIONAL_BISECTION':
            curve=[B(p) for p in row['curve_parameterization']];q=decode(row['quadric'])
            assert q(*curve)==0
            assert matrix(QQ,[[p.monomial_coefficient(s0**(3-i)*s1**i) for i in range(4)] for p in curve]).rank()==4
            n,den=parameter_map(curve,3);assert n.total_degree()==den.total_degree()==2 and gcd(n,den)==1
            ns=list(map(QQ,row['parameter_numerator']));ds=list(map(QQ,row['parameter_denominator']))
            assert n==ns[0]*s0*s0+ns[1]*s0*s1+ns[2]*s1*s1
            assert den==ds[0]*s0*s0+ds[1]*s0*s1+ds[2]*s1*s1
            delta=(ns[1]-t*ds[1])**2-4*(ns[0]-t*ds[0])*(ns[2]-t*ds[2])
            assert delta==P(row['discriminant_polynomial_ascending']) and delta.degree()==2 and delta.discriminant()!=0
            value=delta(0);assert value and not value.is_square() and not row['rational_zero_fibre_points']
            match=next((i for i,c in enumerate(classes) if (value/c).is_square()),None)
            if match is None:match=len(classes);classes.append(value);groups.append([])
            groups[match].append(row['configuration_index']);verified.append(row['configuration_index'])
        elif row['status']=='REDUCIBLE_RESIDUAL':
            ca,cb=map(B,row['residual_biruling_coefficients']);g=gcd(ca,cb);assert g.total_degree()==1
            triple=row['line_triple'];outer=[i for i in triple if sum(j in d['adjacency'][i] for j in triple)==1]
            q=decode(row['quadric']);gram=matrix(QQ,4,4,lambda i,j:q.derivative(R.gen(i)).derivative(R.gen(j))/2)
            ps=[vector(QQ,p+[1,0]) for p in d['plane_points']];vs=[vector(QQ,p+[0,1]) for p in d['line_directions']]
            i,j=outer;A=[s0*ps[i][k]+s1*vs[i][k] for k in range(4)]
            def polar(a,b):return sum(a[k]*gram[k,l]*b[l] for k in range(4) for l in range(4))
            r0,r1=-polar(A,vs[j]),polar(A,ps[j]);C=[r0*ps[j][k]+r1*vs[j][k] for k in range(4)]
            u,v=-g.monomial_coefficient(s1),g.monomial_coefficient(s0)
            line=[s0*a(u,v)+s1*b(u,v) for a,b in zip(A,C)]
            conic=[-(cb//g)*a+(ca//g)*b for a,b in zip(A,C)]
            maps=[]
            for curve,degree in [(line,1),(conic,2)]:
                n,den=parameter_map(curve,degree)
                assert n.total_degree()==den.total_degree()==1 and gcd(n,den)==1
                assert matrix(QQ,[[p.monomial_coefficient(s0),p.monomial_coefficient(s1)] for p in [n,den]]).det()!=0
                maps.append({'curve_degree':degree,'numerator':str(n),'denominator':str(den)})
            sections.append({'configuration_index':row['configuration_index'],'degree_one_parameter_maps':maps})
    assert len(verified)==178 and len(sections)==14 and len(classes)==114
    old=r.read(source.prior.OUTPUT)
    oldclasses=[QQ(row['residual_discriminant']) for row in old['rows'] if 'residual_discriminant' in row and not row['zero_fibre_split']]
    assert all(not (x/y).is_square() for x in classes for y in oldclasses)
    return {'schema':'rank-jump.curve302-conic-triple-bisection-verification.v1','status':'PASS',
        'verified_twisted_cubic_bisections':len(verified),'distinct_zero_fibre_quadratic_fields':len(classes),
        'zero_fibre_squareclass_groups':groups,'overlap_with_prior_four_line_fields':0,
        'reducible_residual_section_certificates':sections,'new_rational_directions':0,
        'bindings':source.prior.binding([Path(__file__),source.INPUT,source.OUTPUT,source.prior.OUTPUT,source.prior.PARENT]),
        'boundary':'All178 explicit degree-two maps are nonsplit at zero. Each of14 reducible residuals is a rational line plus a rational conic, both degree-one over t and hence generic sections in the already full saturated basis. This closes only the266 fixed configurations, not all bisections, covers or Selmer classes.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS178 nonsplit bisections and14 reducible generic residuals')
