#!/usr/bin/env python3
"""Verify star intersections, quadratic fields and exact generic membership."""
import argparse
from itertools import combinations
from pathlib import Path
import retrospective as r
import curve302_section_star_v2 as source

ENDPOINT=r.OUT/'curve302_recovered_mw17_parent_v1.json'
OUTPUT=r.OUT/'rank_jump_curve302_section_star_verification_v1.json'

def compute():
    from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector,prod
    from sage.schemes.elliptic_curves.jacobian import Jacobian
    d=r.read(source.INPUT);out=r.read(source.OUTPUT);parent=r.read(source.PARENT);endpoint=r.read(ENDPOINT)
    for name,sha in out['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha
    R=PolynomialRing(QQ,names=('X','Y','Z'));X,Y,Z=R.gens()
    C=sum(QQ(a['coefficient'])*prod(v**e for v,e in zip(R.gens(),a['powers'][:3])) for a in d['plane_cubic'])
    O=vector(QQ,[1,1,0]);T=PolynomialRing(QQ,'v');v=T.gen()
    dy,dz=[C.derivative(z)(*O) for z in [Y,Z]];direction=vector(QQ,[0,dz,-dy])
    restriction=T(C(*(O+v*direction)));assert restriction[0]==restriction[1]==0 and restriction[3]
    P2=O-restriction[2]/restriction[3]*direction
    dy,dz=[C.derivative(z)(*P2) for z in [Y,Z]];P3=vector(QQ,[0,dz,-dy])
    M=matrix(QQ,[O,P2,P3]).transpose();assert M.det()
    F2=R(M.act_on_polynomial(C));pulled=R(F2(X*X,Y*Z,X*Z))
    assert all(e[0]>=2 and e[2]>=1 for e in pulled.dict())
    F3=R({(e[0]-2,e[1],e[2]-1):c for e,c in pulled.dict().items()})
    aa,bb=map(QQ,[F3.monomial_coefficient(X**3),F3.monomial_coefficient(Y*Y*Z)])
    WW=R(F3(-X,Y/bb,aa*bb*Z)/aa);inter=EllipticCurve(WW(X,Y,1));J=Jacobian(C)
    iso=inter.isomorphism_to(J);u,shift,s,tau=map(QQ,[-6,15,-3,-108])
    E=J.change_weierstrass_model([u,shift,s,tau]);assert list(map(str,E.a_invariants()))==parent['specialized_model']
    MI=M.inverse()
    def convert(point):
        v,w,z=MI*vector(QQ,point);assert z
        raw=iso(inter(-aa*bb*v/z,aa*bb*bb*v*w/z**2))
        return E((raw[0]-shift)/u**2,(raw[1]-s*(raw[0]-shift)-tau)/u**3)
    basis=[E(list(map(QQ,P))) for P in parent['specialized_generic_points']]
    first=convert(d['plane_points'][0]+[1]);assert first in [basis[0],-basis[0]]
    sign=1 if first==basis[0] else -1
    words=matrix(ZZ,endpoint['basis_words_in_recovered_core']);assert abs(words.det())==1
    coords=matrix(ZZ,matrix(QQ,d['generic_line_words'])*words.inverse())
    lines=[]
    for row,point in zip(coords.rows(),d['plane_points']):
        P=sum((int(n)*g for n,g in zip(row,basis)),E(0))
        assert P==sign*convert(point+[1]);lines.append(P)
    # In the plane embedding, H cuts out O+P+Q; these conic sections are
    # precisely basis positions13 and14 in the certified endpoint.
    assert sign*convert([1,0,0])==basis[13] and sign*convert([0,1,0])==basis[14]
    Q4=PolynomialRing(QQ,names=('X','Y','Z','W'));A,B,D,W=Q4.gens()
    P4=[vector(QQ,p+[1,0]) for p in d['plane_points']];V4=[vector(QQ,p+[0,1]) for p in d['line_directions']]
    pairs=[];classes=[];known=[];membership=[];identities=0
    for row in out['rows']:
        if 'duplicate_of' in row:continue
        q=sum(QQ(a['coefficient'])*prod(v**e for v,e in zip(Q4.gens(),a['powers'])) for a in row['quadric'])
        assert matrix(QQ,4,4,lambda i,j:q.derivative(Q4.gen(i)).derivative(Q4.gen(j))/2).det()!=0
        for i in row['star']:
            assert q(*P4[i])==q(*V4[i])==q(*(P4[i]+V4[i]))==0;identities+=3
        residual=PolynomialRing(QQ,'x')(row['residual_quadratic_ascending'])
        assert residual.degree()==2 and str(residual.discriminant())==row['residual_discriminant']
        rational=row['rational_residual_points']
        if rational:
            assert len(rational)==2 and len(row['additional_known_lines'])==1
            j=row['additional_known_lines'][0]
            trace=2*(basis[13]+basis[14])-sum((lines[i] for i in row['star']),E(0))
            actual=[sign*convert(z['plane_point']) for z in rational]
            expected=[lines[j],trace-lines[j]];assert set(actual)==set(expected)
            tracecoords=[2*int(k in [13,14])-sum(int(coords[i,k]) for i in row['star']) for k in range(17)]
            extra=list(map(int,coords[j]));other=[a-b for a,b in zip(tracecoords,extra)]
            membership.append({'star_index':row['star_index'],'generic_coordinate_rows':[extra,other]})
        else:
            disc=residual.discriminant();assert disc and not disc.is_square()
            for old in classes:assert not (disc/old).is_square();pairs.append([str(disc),str(old)])
            classes.append(disc);known.append(row['star_index'])
    assert len(classes)==54 and len(membership)==3 and len(pairs)==1431
    return {'schema':'rank-jump.curve302-section-star-verification.v1','status':'PASS',
        'generic_line_transports_verified':15,'line_quadric_identities':identities,
        'distinct_nonsplit_quadratic_fields':54,'pairwise_squareclass_checks':1431,
        'rational_residual_generic_membership':membership,'new_rational_Kummer_directions_from_catalogue':0,
        'plane_to_normalized_sign':sign,'cubic_to_intermediate_matrix':[[str(c) for c in row] for row in M.rows()],
        'bindings':source.binding([Path(__file__),source.INPUT,source.OUTPUT,source.PARENT,ENDPOINT]),
        'boundary':'Exact membership and field comparisons for this fixed generic-section construction. No exclusion of other classes or larger carriers, and no extra strict class constructed.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS exact generic membership;54 distinct quadratic fields')
