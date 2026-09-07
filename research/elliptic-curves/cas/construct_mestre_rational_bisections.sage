#!/usr/bin/env sage-python
"""Exact bounded section additions on the retained degree-two fibration."""
import argparse,hashlib,json
from itertools import combinations,product
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/mestre-rational-bisections-v1'
def write(path,obj):
    with path.open('x') as f:json.dump(obj,f,indent=2);f.write('\n')
def inputs():
    p=json.loads((D/'protocol.json').read_text())
    for name,h in p['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    ns=json.loads((ART/'mestre_rational_ns_gram_v2.json').read_text())['rows'][0]
    heights=json.loads((ART/'mestre_468_replay_bundle_v1.json').read_text())['rows'][0]['generic_heights']
    return p,ns,heights
def admission():
    p,ns,h=inputs();G=matrix(QQ,ns['rational_NS_Gram']);H=matrix(QQ,h['seed_height_gram'])
    std=matrix.identity(QQ,18).rows();F,O=std[:2];C=G[2:7,2:7];fibre=O+std[7];phi=[]
    for i,r in enumerate(ns['section_profiles']):
        v=std[7+i]-O-(2+r['zero_section_intersection'])*F
        v[2:7]-=C.inverse()*vector(QQ,G[2:7,7+i].column(0));phi.append(v)
    domain=json.loads((ROOT/'artifacts/local/elliptic-curves/mestre-all-section-triangles-v1/height-domain.json').read_text())
    rows=[]
    for word in domain['section_words_up_to_sign']:
        for sign in (-1,1):
            w=sign*vector(ZZ,word)
            comp=[sum(w[i]*ns['section_profiles'][i]['components'][j] for i in range(11))%n for j,n in enumerate((2,2,4))]
            correction=QQ(comp[0]+comp[1])/2+QQ(comp[2]*(4-comp[2]))/4
            oo=(w*H*w-4+correction)/2;assert oo in ZZ and oo>=0
            v=O+(2+oo)*F+sum((w[i]*phi[i] for i in range(11)),vector(QQ,18))
            cross=vector(QQ,comp[:2]+[int(comp[2]==k) for k in (1,2,3)])
            v[2:7]+=C.inverse()*cross
            assert v*G*v==-2 and all(z in ZZ for z in v)
            if fibre*G*v==1:
                rows.append({'word':list(map(int,w)),'class':list(map(int,v)),'old_O_intersection':int(oo),'height':str(w*H*w)})
    rows.sort(key=lambda r:(sum(abs(c) for c in r['word']),r['word']))
    write(D/'admission.json',{'schema':'mestre-rational-bisections.admission.v1','status':'PASS','eligible_old_sections':rows,'selected_old_sections':rows[:3],'sources':p['sources']})
    print('PASS eligible old sections',len(rows),'selected',rows[:3],flush=True)
def curves():
    p,ns,h=inputs();ad=json.loads((D/'admission.json').read_text())
    pencil=json.loads((ART/'mestre_u11_visible_two_neighbor_v1.json').read_text())
    R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field()
    Z=PolynomialRing(QQ,'z');z=Z.gen();L=Z.fraction_field();U=PolynomialRing(L,'u');u=U.gen()
    A=R(pencil['old_A']);B=R(pencil['old_B']);old=EllipticCurve(K,[A,B])
    pole=old([K(c) for c in pencil['pole_point']]);den=R(pencil['q']);offset=R(pencil['c'])
    quartic=U([L(Z(c)) for c in pencil['quartic_T_coefficients']])
    points=[(L(t),L(w)) for t,w in pencil['quartic_sections']]
    t0,q0=points[0];shift=quartic(u+t0);ee,d,c,b,a=shift.list();assert ee==q0*q0
    E=EllipticCurve(L,[0,c,0,d*b-4*q0*q0*a,q0*q0*b*b+a*d*d-4*q0*q0*a*c])
    def forward(t,w):
        assert w*w==quartic(t)
        v=t-t0
        if v==0:
            if w==q0:return E(0)
            assert w==-q0
            X=d*d/(4*q0*q0)-c;return E([X,-d*X/(2*q0)-q0*b])
        X=(2*q0*(w+q0)+d*v)/(v*v)
        Y=(2*(X*X-4*q0*q0*a)*v-2*d*X-4*q0*q0*b)/(4*q0)
        return E([X,Y])
    newpoints=[forward(*pt) for pt in points]
    generators=[newpoints[1],newpoints[2]]
    oldseed=[old([K(v) for v in pt]) for pt in ns['section_points']]
    for row in ad['selected_old_sections']:
        Q=sum((k*P for k,P in zip(row['word'],oldseed)),old(0));x,y=Q.xy()
        zmap=K((den*(y+pole[1])/(x-pole[0])+offset)/(den*den))
        n,dd=zmap.numerator(),zmap.denominator();assert max(n.degree(),dd.degree())==1
        t=L((n[0]-z*dd[0])/(z*dd[1]-n[1]));assert zmap(t)==z
        m=den(t)*z-offset(t)/den(t)
        w=(2*x(t)+pole[0](t)-m*m)/den(t)
        assert w*w==quartic(t)
        generators.append(forward(t,w))
    words=[]
    for size in (1,2):
        for indices in combinations(range(len(generators)),size):
            for signs in product((-1,1),repeat=size):
                word=[0]*len(generators)
                for i,s in zip(indices,signs):word[i]=s
                words.append(word)
    write(D/'maps-and-words.json',{'generators':[[str(v) for v in P] for P in generators],'new_ainvs':list(map(str,E.a_invariants())),'quartic_origin':[str(t0),str(q0)],'words':words,'old_section_maps':ad['selected_old_sections']})
    rows=[]
    for i,word in enumerate(words):
        Q=sum((k*P for k,P in zip(word,generators)),E(0));r={'index':i,'word':word}
        if Q.is_zero():r['status']='ZERO_ON_ALTERNATE_FIBRATION'
        else:
            X,Y=Q.xy();v=(2*q0*Y+d*X+2*q0*q0*b)/(X*X-4*q0*q0*a)
            t=v+t0;w=(X*v*v-d*v)/(2*q0)-q0
            assert w*w==quartic(t) and forward(t,w)==Q
            degree=max(t.numerator().degree(),t.denominator().degree());r['old_base_degree']=int(degree)
            if degree<=0:r['status']='OLD_VERTICAL_CURVE'
            else:
                m=den(t)*z-offset(t)/den(t);x=(den(t)*w-pole[0](t)+m*m)/2;y=m*(x-pole[0](t))-pole[1](t)
                assert y*y==x*x*x+A(t)*x+B(t)
                assert (den(t)*(y+pole[1](t))/(x-pole[0](t))+offset(t))/(den(t)**2)==z
                dx=x.denominator();dt=t.denominator()
                finite=dx.degree()-dx.gcd(dt**4).degree()
                infinite=max(0,x.numerator().degree()-dx.degree()-4*max(0,t.numerator().degree()-dt.degree()))
                assert (finite+infinite)%2==0
                oo=(finite+infinite)//2
                r.update(status='EXACT_RATIONAL_CURVE',old_O_intersection=int(oo),old_T=str(t),old_x=str(x),old_y=str(y))
                if degree==2:
                    n,dd=t.numerator(),t.denominator();deck=-(n[1]-t*dd[1])/(n[2]-t*dd[2])-z
                    assert t(deck)==t and deck!=z and (x(deck)!=x or y(deck)!=y)
                    r['deck']=str(deck);r['section_plus_bisection_pencil_candidate']=oo==2
        rows.append(r);write(D/('curve'+str(i)+'.json'),r)
        print('CURVE',i,r['status'],r.get('old_base_degree'),r.get('old_O_intersection'),flush=True)
    write(D/'result.json',{'schema':'mestre-rational-bisections.v1','status':'PASS','rows':rows,'sources':p['sources'],'scope':'Exact section additions on one retained alternate fibration. Birational inverse proves rational-curve parametrizations. Degree2 with old O intersection2 gives effective nef D=O+R of old degree3, but new pencil equation, rational section, rank and search admission require further proofs. No elliptic point search or new parent.'})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['admission','curves']);args=p.parse_args()
    {'admission':admission,'curves':curves}[args.stage]()
