#!/usr/bin/env python3
"""Independent elliptic/nodal group verification of projective local gates."""
import argparse
from itertools import product
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,EllipticCurve,prod
import retrospective as r

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_direct_relation_local_gate_inputs_v1.json'
FIRST=r.OUT/'rank_jump_direct_relation_local_gate_v1.json'
SECOND=r.OUT/'rank_jump_direct_relation_good_primes_v1.json'
OUTPUT=r.OUT/'rank_jump_direct_relation_local_gate_verification_v1.json'
PRIOR=r.OUT/'rank_jump_triple_translate_controls_verification_v1.json'


def verify():
    inp=r.read(INPUT);sources=[r.read(FIRST),r.read(SECOND)]
    for data in (inp,*sources):
        for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    results=[]
    for old in [row for source in sources for row in source['rows']]:
        p=old['prime'];Fp=GF(p);R=PolynomialRing(Fp,'t')
        A,B=[R([Fp(QQ(x)) for x in inp[k]]) for k in ('A','B')]
        lifts=[{k:R([Fp(QQ(x)) for x in c[k]]) for k in ('q','x0','x1','y0','y1')} for c in inp['lifts']]
        assert A.degree()<=8 and B.degree()<=12
        for c in lifts:
            assert all(c[k].degree()<=bound for k,bound in [('q',2),('x0',4),('x1',3),('y0',6),('y1',5)])
        delta=-16*(4*A**3+27*B**2);qs=[c['q'] for c in lifts]
        good=delta.degree()==24 and delta.is_squarefree() and 4*A[8]**3+27*B[12]**2!=0
        good=good and all(q.degree()==2 and q.is_squarefree() and q.gcd(delta)==1 for q in qs)
        good=good and all(qs[i].gcd(qs[j])==1 for i in range(3) for j in range(i))
        if not good:
            assert old['status']=='UNKNOWN_BAD_GEOMETRY' and 'relations' not in old
            results.append({'prime':p,'status':'VERIFIED_UNKNOWN_BAD_GEOMETRY'});continue
        assert old['status']=='PASS'
        basis=[]
        for c in inp['sections']:
            x=R([Fp(QQ(z)) for z in c['x_coefficients_low_to_high']])
            if 'y_coefficients_low_to_high' in c:y=R([Fp(QQ(z)) for z in c['y_coefficients_low_to_high']])
            else:
                ch=c['chord'];ref=basis[ch['reference_basis_index']]
                y=ref[1]+R([Fp(QQ(z)) for z in ch['slope_coefficients_low_to_high']])*(x-ref[0])
            assert x.degree()<=4 and y.degree()<=6 and y*y==x**3+A*x+B
            basis.append((x,y))
        for c in lifts:
            q,x0,x1,y0,y1=[c[k] for k in ('q','x0','x1','y0','y1')]
            assert y0*y0+q*y1*y1==x0**3+3*q*x0*x1*x1+A*x0+B
            assert 2*y0*y1==3*x0*x0*x1+q*x1**3+A*x1
        hits=[set() for _ in inp['generic_words']];carrier_count=0;nodal=[]
        for t in list(Fp)+[None]:
            def ev(poly,weight):return poly[weight] if t is None else poly(t)
            a,b=ev(A,8),ev(B,12);bs=[(ev(x,4),ev(y,6)) for x,y in basis]
            roots=[]
            for c in lifts:
                value=ev(c['q'],2)
                roots.append([Fp(0)] if not value else list({value.sqrt(),-value.sqrt()}) if value.is_square() else [])
            choices=list(product(*roots));carrier_count+=len(choices);smooth=4*a**3+27*b*b!=0
            tag='infinity' if t is None else int(t)
            if smooth:
                E=EllipticCurve(Fp,[a,b]);gen=[E(x,y) for x,y in bs]
                targets=[sum((n*P for n,P in zip(word,gen,strict=True)),E(0)) for word in inp['generic_words']]
                def encode(P):return E(*P)
                def relation(ps):return ps[0]-ps[1]+ps[2]
            else:
                # Independent normalization of the smooth nodal cubic to Gm.
                X=PolynomialRing(Fp,'x');x=X.gen();cubic=x**3+a*x+b;double=cubic.gcd(cubic.derivative())
                assert double.degree()==1
                rr=-double[0]/double[1];assert rr and cubic==(x-rr)**2*(x+2*rr)
                field=Fp if (3*rr).is_square() else GF(p**2,'alpha');k=field(3*rr).sqrt()
                def encode(P):
                    x,y=P;assert y*y==x**3+a*x+b and x!=rr
                    z=field(y)/field(x-rr);return (z-k)/(z+k)
                gen=[encode(P) for P in bs]
                targets=[prod(P**n for n,P in zip(word,gen,strict=True)) for word in inp['generic_words']]
                def relation(ps):return ps[0]/ps[1]*ps[2]
                nodal.append({'base_t':tag,'double_root':int(rr),'split':field==Fp,'carrier_points_above':len(choices)})
            for us in choices:
                ps=[encode((ev(c['x0'],4)+u*ev(c['x1'],3),ev(c['y0'],6)+u*ev(c['y1'],5))) for c,u in zip(lifts,us,strict=True)]
                value=relation(ps)
                for i,S in enumerate(targets):
                    if value==S:hits[i].add((tag,tuple(map(int,us)),bool(smooth)))
        assert carrier_count==old['carrier_Fp_point_count']
        for i,points in enumerate(hits):
            actual={(row['base_t'],tuple(row['roots']),row['smooth_elliptic_fibre']) for row in old['relations'][i]['Fp_points']}
            assert points==actual and len(points)==old['relations'][i]['Fp_point_count']
            assert old['relations'][i]['local_obstruction_proved']==(not points)
        assert len(old['base_fibre_counts'])==p+1 and sum(x['root_choice_count'] for x in old['base_fibre_counts'])==carrier_count
        results.append({'prime':p,'status':'PASS','carrier_Fp_points':carrier_count,'relation_Fp_points':list(map(len,hits)),
                        'nodal_fibres':nodal,'infinity_included':True})
    # Only after completing the independent equation-only replay compare with
    # the earlier characteristic-zero finite schemes, including Hensel lifting.
    comparisons=[]
    for row in results:
        if row['status']!='PASS':continue
        p=row['prime'];Rp=PolynomialRing(GF(p),'t')
        for i,count in enumerate(row['relation_Fp_points']):
            path=r.OUT/f'rank_jump_triple_translate_control_{i}_v1.json';data=r.read(path)
            Q=PolynomialRing(QQ,'t');f=Q(data['intersection_polynomial']);f*=f.denominator();fp=Rp(f)
            comparison={'prime':p,'index':i,'direct_Fp_points':count}
            if fp.degree()==12 and fp.is_squarefree():
                nr=len(fp.roots(multiplicities=False));assert nr==count
                comparison.update({'finite_scheme_simple_Fp_roots':nr,'relation_Qp_solubility':'YES' if nr else 'NO'})
            else:comparison['finite_scheme_comparison']='UNKNOWN_BAD_POLYNOMIAL_REDUCTION'
            comparisons.append(comparison)
    assert next(x for x in results if x['prime']==131)['relation_Fp_points']==[1,1,0]
    assert next(x for x in results if x['prime']==137)['relation_Fp_points']==[1,1,1]
    return {'schema':'rank-jump.direct-relation-local-gate-verification.v1','status':'PASS','rows':results,
        'after_detector_comparison_with_finite_schemes':comparisons,
        'bindings':{str(path.relative_to(r.ROOT)):r.digest(path.read_bytes()) for path in
                    (INPUT,FIRST,SECOND,PRIOR,Path(__file__),HERE/'retrospective.py',
                     *(r.OUT/f'rank_jump_triple_translate_control_{i}_v1.json' for i in range(3)))},
        'boundary':'Sage elliptic groups and independent nodal normalization replay every projective carrier point; polynomial-root comparisons enter only afterwards. Bad surface reductions remain UNKNOWN.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();data=verify()
    if a.mode=='build':r.write_new(OUTPUT,data)
    else:assert r.read(OUTPUT)==data
    for row in data['rows']:print(row['prime'],row['status'],row.get('relation_Fp_points'),row.get('nodal_fibres'))
    print(data['after_detector_comparison_with_finite_schemes'])
