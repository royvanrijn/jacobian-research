#!/usr/bin/env python3
"""Explicit rational bisections from one generic conic and three generic lines."""
import argparse
from itertools import combinations
from pathlib import Path
import subprocess
import sys
import retrospective as r
import curve302_section_star_v2 as prior

PROTOCOL=Path(__file__).with_name('CURVE302_CONIC_TRIPLE_BISECTION_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_curve302_conic_triple_bisection_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_curve302_conic_triple_bisection_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-curve302-conic-triple-bisection-v1'

def prepare():
    d=r.read(prior.INPUT);adj=d['adjacency'];configs=[]
    for a,b in combinations(range(15),2):
        if b not in adj[a]:continue
        eligible=[i for i in range(15) if i not in [a,b] and i not in adj[a] and i not in adj[b]]
        for triple in combinations(eligible,3):
            if sum(j in adj[i] for i,j in combinations(triple,2))==2:configs.append({'conic_pair':[a,b],'line_triple':list(triple)})
    assert len(configs)==266
    keys=['quartic','plane_cubic','residual_line','quartic_W_quotient','plane_points','line_directions','generic_line_words','adjacency']
    r.write_new(INPUT,{'schema':'rank-jump.curve302-conic-triple-bisection-inputs.v1',**{k:d[k] for k in keys},'configurations':configs,
        'bindings':prior.binding([Path(__file__),PROTOCOL,prior.INPUT,prior.PARENT])})
    print('PREPARED266 generic configurations',flush=True)

def compute():
    from sage.all import QQ,PolynomialRing,matrix,vector,prod,gcd
    d=r.read(INPUT)
    for path,sha in d['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,['X','Y','Z','W']);X,Y,Z,W=R.gens()
    def decode(rows):return sum(QQ(a['coefficient'])*prod(v**e for v,e in zip(R.gens(),a['powers'])) for a in rows)
    def encode(q):return [{'powers':list(map(int,e)),'coefficient':str(c)} for e,c in sorted(q.dict().items())]
    F=decode(d['quartic']);L=decode(d['residual_line']);node=vector(QQ,[1,1,0,0])
    P=[vector(QQ,p+[1,0]) for p in d['plane_points']];V=[vector(QQ,v+[0,1]) for v in d['line_directions']]
    mons=[u*v for i,u in enumerate(R.gens()) for v in R.gens()[i:]]
    Plane=PolynomialRing(QQ,['a','b','c']);a,b,c=Plane.gens();pmons=[a*a,a*b,a*c,b*b,b*c,c*c]
    T=PolynomialRing(QQ,['s0','s1','la','mu']);s0,s1,la,mu=T.gens()
    B=PolynomialRing(QQ,['s0','s1']);z0,z1=B.gens();U=PolynomialRing(QQ,'z');z=U.gen()
    line_rows=[]
    for p,v in zip(P,V):
        values=[m(*[a*x+b*y for x,y in zip(p,v)]) for m in mons]
        line_rows.append([[q.monomial_coefficient(term) for q in values] for term in [a*a,a*b,b*b]])
    conics={};seen={};rows=[]
    for number,config in enumerate(d['configurations']):
        pair=tuple(config['conic_pair']);triple=config['line_triple'];record={'configuration_index':number,**config}
        if pair not in conics:
            i,j=pair;basis=matrix(QQ,[P[i],V[i],P[j],V[j]]).row_space().basis_matrix();assert basis.nrows()==3
            param=[a*basis[0,k]+b*basis[1,k]+c*basis[2,k] for k in range(4)]
            ell=matrix(QQ,[P[i],V[i],P[j],V[j]]).right_kernel().basis()[0]
            lines=[]
            for k in pair:
                cp=basis.transpose().solve_right(P[k]);cv=basis.transpose().solve_right(V[k])
                coeff=matrix(QQ,[cp,cv]).right_kernel().basis()[0];lines.append(sum(x*y for x,y in zip(coeff,Plane.gens())))
            conic=Plane(F(*param)//(lines[0]*lines[1]));assert F(*param)==conic*lines[0]*lines[1]
            smooth=matrix(QQ,3,3,lambda i,j:conic.derivative(Plane.gen(i)).derivative(Plane.gen(j))/2).det()!=0
            conics[pair]=(basis,param,ell,conic,smooth)
        basis,param,ell,conic,smooth=conics[pair]
        if not smooth:record['status']='REDUCIBLE_GENERIC_CONIC';rows.append(record);continue
        restricted=[m(*param) for m in mons]
        constraints=[[q.monomial_coefficient(term) for q in restricted]+[-conic.monomial_coefficient(term)] for term in pmons]
        constraints += [row+[0] for i in triple for row in line_rows[i]]
        kernel=matrix(QQ,constraints).right_kernel();record['quadric_space_dimension']=int(kernel.dimension())
        if kernel.dimension()!=1:record['status']='NONUNIQUE_QUADRIC';rows.append(record);continue
        coefficients=kernel.basis()[0];q=sum(x*m for x,m in zip(coefficients[:10],mons));q/=q.coefficients()[0]
        key=str(q)
        if key in seen:record.update(status='DUPLICATE_QUADRIC',duplicate_of=seen[key]);rows.append(record);continue
        seen[key]=number;record['quadric']=encode(q)
        gram=matrix(QQ,4,4,lambda i,j:q.derivative(R.gen(i)).derivative(R.gen(j))/2)
        if not gram.det():record['status']='SINGULAR_QUADRIC';rows.append(record);continue
        if q(*node)==0:record['status']='QUADRIC_THROUGH_NODE';rows.append(record);continue
        outer=[i for i in triple if sum(j in d['adjacency'][i] for j in triple)==1]
        central=next(i for i in triple if i not in outer);assert len(outer)==2
        i,j=outer;A=[s0*P[i][k]+s1*V[i][k] for k in range(4)]
        def polar(a,b):return sum(a[k]*gram[k,l]*b[l] for k in range(4) for l in range(4))
        t0=-polar(A,V[j]);t1=polar(A,P[j]);C=[t0*P[j][k]+t1*V[j][k] for k in range(4)]
        qparam=[la*x+mu*y for x,y in zip(A,C)];assert q(*qparam)==0
        ann=matrix(QQ,[P[central],V[central]]).right_kernel().basis()
        central_form=next(sum(x*y for x,y in zip(v,A)) for v in ann if sum(x*y for x,y in zip(v,A))!=0)
        plane_form=sum(x*y for x,y in zip(ell,qparam));divisor=la*mu*central_form*plane_form
        pulled=T(F(*qparam));residual=T(pulled//divisor);assert residual*divisor==pulled
        assert all(e[0]+e[1]==2 and e[2]+e[3]==1 for e in residual.dict())
        ca=B(residual(s0,s1,1,0));cb=B(residual(s0,s1,0,1));common=gcd(ca,cb)
        record['residual_biruling_coefficients']=[str(ca),str(cb)]
        if common.total_degree()>0:record['status']='REDUCIBLE_RESIDUAL';rows.append(record);continue
        curve=[B(q0(s0,s1,-cb,ca)) for q0 in qparam]
        assert all(p.total_degree()==3 for p in curve) and gcd(curve)==1
        assert F(*curve)==0 and q(*curve)==0
        numerator=curve[3];denominator=B(L(*curve));base=gcd(numerator,denominator)
        n=B(numerator//base);den=B(denominator//base)
        assert base.total_degree()==1 and n.total_degree()==den.total_degree()==2
        ns=[n.monomial_coefficient(t) for t in [z0*z0,z0*z1,z1*z1]]
        ds=[den.monomial_coefficient(t) for t in [z0*z0,z0*z1,z1*z1]]
        Parameter=PolynomialRing(QQ,'t');t=Parameter.gen();delta=(ns[1]-t*ds[1])**2-4*(ns[0]-t*ds[0])*(ns[2]-t*ds[2])
        record.update(status='EXPLICIT_RATIONAL_BISECTION',curve_parameterization=list(map(str,curve)),
            parameter_numerator=list(map(str,ns)),parameter_denominator=list(map(str,ds)),
            discriminant_polynomial_ascending=list(map(str,delta.list())),zero_fibre_discriminant=str(delta(0)),
            zero_fibre_split=bool(delta(0).is_square()),zero_fibre_etale=bool(delta(0)!=0))
        rational=[];rootpairs=[(r0,QQ(1)) for r0,m in U(n(z,1)).roots()]
        if n(1,0)==0:rootpairs.append((QQ(1),QQ(0)))
        for u,v in rootpairs:
            if den(u,v)==0:continue
            point=[p(u,v) for p in curve];assert point[3]==0 and any(point)
            scale=next(x for x in point if x);point=[x/scale for x in point]
            assert F(*point)==0
            rational.append({'bisection_parameter':[str(u),str(v)],'plane_point':list(map(str,point[:3]))})
        record['rational_zero_fibre_points']=rational;rows.append(record)
        r.write_new(WORK/('configuration_%03d.json'%number),record)
        print('BISECTION',number,'RATIONAL',len(rational),'DELTA_BITS',max(abs(x.numerator()).nbits() for x in delta.list()),flush=True)
    return {'schema':'rank-jump.curve302-conic-triple-bisection.v1','status':'PASS','rows':rows,
        'candidate_count':len(rows),'distinct_quadrics':len(seen),
        'explicit_bisections':sum(row['status']=='EXPLICIT_RATIONAL_BISECTION' for row in rows),
        'split_bisections':sum(bool(row.get('rational_zero_fibre_points')) for row in rows),
        'bindings':prior.binding([Path(__file__),PROTOCOL,INPUT,prior.PARENT]),
        'boundary':'Explicit geometric bisections and pre-point quadratic splitting conditions at the fixed zero fibre. New quotient directions and strict classes await arithmetic verification.'}

def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    with (WORK/'worker.log').open('x') as log:
        try:
            p=subprocess.run([sys.executable,__file__,'worker'],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['bounds']['geometry_seconds'])
            error='worker failure' if p.returncode else None
        except subprocess.TimeoutExpired:error='bounded timeout'
    if error and not OUTPUT.exists():r.write_new(OUTPUT,{'status':'UNKNOWN','reason':error})
    print(r.read(OUTPUT)['status'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','capture','worker']);a=p.parse_args()
    if a.mode=='worker':r.write_new(OUTPUT,compute())
    else:globals()[a.mode]()
