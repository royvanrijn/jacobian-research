#!/usr/bin/env python3
"""Exact fixed section-star intersections on the recovered 302 parent."""
import argparse
from itertools import combinations
from pathlib import Path
import subprocess
import sys
import retrospective as r

PROTOCOL=Path(__file__).with_name('CURVE302_SECTION_STAR_PROTOCOL.json')
SOURCE=r.OUT/'curve302_recovered_mw17_construction_inputs_v1.json'
PARENT=r.OUT/'rank_jump_curve302_new_parent_inputs_v1.json'
INPUT=r.OUT/'rank_jump_curve302_section_star_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_curve302_section_star_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-curve302-section-star-v1'

def binding(paths):return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}

def prepare():
    from sage.all import QQ,matrix,vector
    d=r.read(SOURCE)['quartic']
    keys=['quartic','plane_cubic','residual_line','quartic_W_quotient','plane_points','line_directions']
    data={k:d[k] for k in keys};data['generic_line_words']=d['polarization']['line_section_words']
    P=[vector(QQ,xy+[1,0]) for xy in d['plane_points']]
    V=[vector(QQ,uv+[0,1]) for uv in d['line_directions']]
    adj=[[] for _ in P]
    for i,j in combinations(range(15),2):
        if matrix([P[i],V[i],P[j],V[j]]).det()==0:adj[i].append(j);adj[j].append(i)
    stars=[(i,*js) for i in range(15) for js in combinations(adj[i],3)
           if all(b not in adj[a] for a,b in combinations(js,2))]
    assert len(stars)==60 and sum(map(len,adj))==54
    data.update(schema='rank-jump.curve302-section-star-inputs.v1',adjacency=adj,stars=stars,
                bindings=binding([Path(__file__),PROTOCOL,SOURCE,PARENT]))
    r.write_new(INPUT,data);print('Prepared60 stars from generic lines',flush=True)

def compute():
    from sage.all import QQ,ZZ,PolynomialRing,matrix,vector,prod
    d=r.read(INPUT)
    R=PolynomialRing(QQ,names=('X','Y','Z','W'));X,Y,Z,W=R.gens()
    def decode(rows):return sum(QQ(a['coefficient'])*prod(v**e for v,e in zip(R.gens(),a['powers'])) for a in rows)
    def encode(q):return [{'powers':list(e),'coefficient':str(c)} for e,c in sorted(q.dict().items())]
    quartic,C,L,Q3=map(decode,[d['quartic'],d['plane_cubic'],d['residual_line'],d['quartic_W_quotient']])
    assert quartic==C*L+W*Q3
    P=[vector(QQ,xy+[1,0]) for xy in d['plane_points']];V=[vector(QQ,uv+[0,1]) for uv in d['line_directions']]
    T=PolynomialRing(QQ,names=('z','w'));z,w=T.gens()
    for A,B in zip(P,V):assert quartic(*[z*a+w*b for a,b in zip(A,B)])==0
    node=vector(QQ,[1,1,0,0]);assert quartic(*node)==0 and all(quartic.derivative(v)(*node)==0 for v in R.gens())
    mons=[a*b for i,a in enumerate(R.gens()) for b in R.gens()[i:]]
    restrictions=[]
    for A,B in zip(P,V):
        pols=[m(*[z*a+w*b for a,b in zip(A,B)]) for m in mons]
        restrictions.append([[q.monomial_coefficient(v) for q in pols] for v in [z*z,z*w,w*w]])
    affine=PolynomialRing(QQ,names=('x','y'));x,y=affine.gens();cf=affine(C(x,y,1,0));Ux=PolynomialRing(QQ,'x')
    seen={};rows=[]
    for n,star in enumerate(d['stars']):
        mat=matrix(QQ,[row for i in star for row in restrictions[i]])
        ker=mat.right_kernel();record={'star_index':n,'star':star,'quadric_space_dimension':ker.dimension()}
        if ker.dimension()!=1:rows.append(record);continue
        v=ker.basis()[0];q=sum(a*m for a,m in zip(v,mons));q/=q.coefficients()[0]
        key=str(q)
        if key in seen:record['duplicate_of']=seen[key];rows.append(record);continue
        seen[key]=n
        gram=matrix(QQ,4,4,lambda i,j:q.derivative(R.gen(i)).derivative(R.gen(j))/2)
        smooth=gram.det()!=0;record.update(quadric=encode(q),quadric_smooth=bool(smooth),
            contains_node=bool(q(*node)==0),additional_known_lines=[i for i in range(15) if i not in star
                and all(sum(v[j]*restrictions[i][k][j] for j in range(10))==0 for k in range(3))])
        if not smooth:record['status']='DEGENERATE_QUADRIC';rows.append(record);continue
        qf=affine(q(x,y,1,0));res=cf.resultant(qf,y);assert res.degree(y)==0
        res=Ux([QQ(res.monomial_coefficient(x**i)) for i in range(res.degree(x)+1)])
        prescribed=Ux.prod(Ux.gen()-P[i][0] for i in star)
        residual,rem=res.quo_rem(prescribed);assert rem==0 and residual
        residual=residual.monic();record.update(residual_degree=int(residual.degree()),
            residual_quadratic_ascending=list(map(str,residual.list())))
        if residual.degree()==2:
            disc=residual.discriminant();record['residual_discriminant']=str(disc)
            record['zero_fibre_split']=bool(disc.is_square());record['zero_fibre_etale']=bool(disc!=0)
        elif residual.degree()<2:record['zero_fibre_at_infinity']=True
        else:raise ArithmeticError('unexpected residual degree')
        rational=[]
        for x0,mult in residual.roots():
            Uy=PolynomialRing(QQ,'y');yy=Uy.gen()
            g=Uy(cf(x0,yy)).gcd(Uy(qf(x0,yy)))
            for y0,ymult in g.roots():
                assert cf(x0,y0)==qf(x0,y0)==0
                rational.append({'plane_point':[str(x0),str(y0),'1'],'x_multiplicity':int(mult),'y_multiplicity':int(ymult)})
        record.update(status='EXACT_ZERO_FIBRE_RESIDUAL',rational_residual_points=rational)
        rows.append(record)
        r.write_new(WORK/('star_%02d.json'%n),record)
        print(n,star,'extra lines',len(record['additional_known_lines']),'node',record['contains_node'],
            'degree',record['residual_degree'],'rational',len(rational),flush=True)
    ell=L.monomial_coefficient(Z);at=C(X,X+ell*Z,Z,0);bt=Q3(X,X+ell*Z,Z,0)
    assert at%Z==bt%Z==0;at=R(at//Z);bt=R(bt//Z)
    base0=Ux([at.monomial_coefficient(Z*Z),at.monomial_coefficient(X*Z),at.monomial_coefficient(X*X)])
    return {'schema':'rank-jump.curve302-section-star.v1','status':'PASS','bindings':binding([Path(__file__),PROTOCOL,INPUT,PARENT]),
        'star_count':len(rows),'distinct_quadric_count':len(seen),'rows':rows,
        'base_line_control':{'constant_quadratic':encode(at),'linear_parameter_quadratic':encode(bt),
                             'zero_fibre_discriminant':str(base0.discriminant()),'zero_fibre_split':bool(base0.discriminant().is_square())},
        'boundary':'Geometric intersection construction only. No new class, strictness, independence or irreducibility of the generic residual asserted before the separate arithmetic stage.'}

def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    with (WORK/'worker.log').open('x') as log:
        try:
            p=subprocess.run([sys.executable,__file__,'worker'],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['bounds']['geometry_seconds'])
            error='worker failure' if p.returncode else None
        except subprocess.TimeoutExpired:error='bounded timeout'
    if error and not OUTPUT.exists():r.write_new(OUTPUT,{'status':'UNKNOWN','reason':error})
    out=r.read(OUTPUT);print(out['status'],out.get('distinct_quadric_count'),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','capture','worker']);a=p.parse_args()
    if a.mode=='worker':r.write_new(OUTPUT,compute())
    else:globals()[a.mode]()
