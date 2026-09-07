#!/usr/bin/env sage-python
"""Translate two exact bisections by zero and +-the eleven full old sections."""
import hashlib,json
from itertools import product
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/mestre-translated-bisections-v1'
def write(p,r):
    with p.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
def main():
    protocol=json.loads((D/'protocol.json').read_text())
    for name,h in protocol['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h
    source=json.loads((ROOT/'artifacts/local/elliptic-curves/mestre-rational-bisections-v1/result.json').read_text())
    candidates=[r for r in source['rows'] if r.get('old_base_degree')==2];assert [r['index'] for r in candidates]==[27,35]
    ns=json.loads((ART/'mestre_rational_ns_gram_v2.json').read_text())['rows'][0]
    h=json.loads((ART/'mestre_468_replay_bundle_v1.json').read_text())['rows'][0]['generic_heights']
    G=matrix(QQ,ns['rational_NS_Gram']);inv=G.inverse();std=matrix.identity(QQ,18).rows();F,O=std[:2]
    R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();A=R(h['curve_A_coefficients']);B=R(h['curve_B_coefficients'])
    Z=PolynomialRing(QQ,'z');z=Z.gen();L=Z.fraction_field()
    words=[[0]*11]
    for i in range(11):
        for s in (-1,1):w=[0]*11;w[i]=s;words.append(w)
    def intersection(x,t):
        dx=x.denominator();dt=t.denominator()
        n=dx.degree()-dx.gcd(dt**4).degree()+max(0,x.numerator().degree()-dx.degree()-4*max(0,t.numerator().degree()-dt.degree()))
        assert n%2==0;return n//2
    # The complete old component intersection possibilities for a bisection:
    # each I2 component count <=2; three nonidentity I4 counts sum<=2.
    component_options=[(i,j,a,b,c) for i,j in product(range(3),repeat=2) for a,b,c in product(range(3),repeat=3) if a+b+c<=2]
    assert len(component_options)==90
    old_pool=json.loads((ROOT/'artifacts/local/elliptic-curves/mestre-degree-three-pencils-v1/roster.json').read_text())
    visible=[vector(QQ,v) for v in old_pool['visible_curve_classes']]
    rows=[]
    for source_row in candidates:
        t=L(source_row['old_T']);deck=L(source_row['deck']);E=EllipticCurve(L,[A(t),B(t)])
        base=E([L(source_row['old_x']),L(source_row['old_y'])])
        seed=[E([K(v)(t) for v in coords]) for coords in ns['section_points']]
        for word in words:
            point=base+sum((a*P for a,P in zip(word,seed)),E(0));x,y=point.xy()
            assert x(deck)!=x or y(deck)!=y
            oo=intersection(x,t)
            row={'source_curve_index':source_row['index'],'translation_word':word,'old_base_degree':2,'old_O_intersection':int(oo),'old_T':str(t),'old_x':str(x),'old_y':str(y),'deck':str(deck)}
            if oo==2:
                section_intersections=[intersection((point-P)[0],t) for P in seed]
                classes=[]
                for comp in component_options:
                    rhs=vector(QQ,[2,oo,*comp,*section_intersections]);v=inv*rhs
                    if all(c in ZZ for c in v) and v*G*v==-2:classes.append(list(map(int,v)))
                row['old_basis_section_intersections']=list(map(int,section_intersections));row['possible_NS_classes']=classes
                if len(classes)==1:
                    v=vector(QQ,classes[0]);fibre=O+v;assert fibre*G*fibre==0 and fibre*G*F==3
                    curves=visible+[v];cross=[fibre*G*c for c in curves];assert min(cross)>=0
                    one=[i for i,c in enumerate(cross) if c==1];zero=[i for i,c in enumerate(cross) if c==0]
                    vertical=matrix(QQ,[curves[i] for i in zero]);bound=17-vertical.rank()
                    row.update(fibre_class=list(map(int,fibre)),visible_section_indices=one,vertical_curve_indices=zero,visible_MW_upper_bound=int(bound))
            rows.append(row);write(D/('curve'+str(len(rows)-1)+'.json'),row)
            print('TRANSLATE',len(rows)-1,'old O intersection',oo,'NS possibilities',len(row.get('possible_NS_classes',[])),'upper',row.get('visible_MW_upper_bound'),flush=True)
    write(D/'result.json',{'schema':'mestre-translated-bisections.v1','status':'PASS','rows':rows,'sources':protocol['sources'],
        'scope':'Two rational bisections translated by zero and+-the eleven known full basis sections,46 exact curves. O intersection2 gives a nef old-degree3 pencil. All ninety old component-intersection possibilities are checked; a unique integral square-minus-two class identifies the curve in the full rational NS lattice. Visible degree-one curves prove Jacobian sections. No new pencil equation, generic MW lower bound, parameter population, point search or new parent.'})
if __name__=='__main__':main()
