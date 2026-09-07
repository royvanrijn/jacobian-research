#!/usr/bin/env sage-python
"""Fixed off-path Kihara parents, exact sections and surface fingerprints."""
import json,hashlib,sys
from pathlib import Path
from sage.all import QQ,PolynomialRing,EllipticCurve,prod,matrix
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
from audit_kihara_parent_distinctness import count
D=ROOT/'artifacts/local/elliptic-curves/kihara-positive-parents-v1'
def write(p,r):
    with p.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
def construct(v):
    R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();S=PolynomialRing(R,'x');x=S.gen()
    roots=[QQ(0),(2*v*v+v+2)**2,2*(v+1)**2*(2*v*v+v+1),4*v*v-v+4,v*(2*v-1)*(2*v*v+4*v+5),4*v**4+8*v**3+9*v*v-2*v+2]
    normalizer=roots[1];roots=[c/normalizer for c in roots];shifted=[c+T for c in roots]+[c-T for c in roots]
    product=prod(x-r for r in shifted);square=x**6
    for j in range(5,-1,-1):square+=(product[6+j]-(square*square)[6+j])/2*x**j
    remainder=square*square-product;assert remainder.degree()==4
    coeff=[]
    for c in remainder.list():
        q,r=c.quo_rem(T*T);assert not r;coeff.append(q)
    quartic=S(coeff);points=[(K(r),K(square(r)/T)) for r in shifted]
    den=2*v*v+2*v+3;constant=(8*v**6+28*v**5+58*v**4+69*v**3+76*v*v+40*v+22)/(den*normalizer)
    xx=K(constant+(2*v*v+4*v+5)*T/den);yy=K(quartic(xx)).sqrt()
    if yy.numerator().leading_coefficient()/yy.denominator().leading_coefficient()<0:yy=-yy
    points.append((xx,yy));assert all(Y*Y==quartic(X) for X,Y in points)
    e,d,c,b,a=coeff;A=-27*(12*a*e-3*b*d+c*c);B=-27*(72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3)
    disc=-16*(4*A**3+27*B**2)
    row={'ratio':str(v),'normalized_roots':list(map(str,roots)),'normalizer':str(normalizer),'raw_A':list(map(str,A.list())),'raw_B':list(map(str,B.list())),'quartic_coefficients':[list(map(str,c.list())) for c in coeff],'quartic_points':[[str(c) for c in p] for p in points]}
    if (A.degree(),B.degree(),disc.degree())!=(8,12,20) or disc.gcd(disc.derivative())!=1 or disc.gcd(A)!=1:
        row.update(status='OUTSIDE_FIXED_SMOOTH_FIBRE_GATE',discriminant_factors=str(disc.factor()));return row
    node=-3*B[12]/(2*A[8]);tangent=QQ(3*node).sqrt();assert node**3+A[8]*node+B[12]==0
    X0,Y0=points[0];trans=quartic(x+X0);ee,dd,cc,bb,aa=[K(trans[i]) for i in range(5)];E=EllipticCurve(K,[A,B]);images=[]
    for X,Y in points[1:]:
        xx=X-X0;wx=(2*Y0*(Y+Y0)+dd*xx)/xx**2;wy=(2*(wx*wx-4*Y0*Y0*aa)*xx-2*dd*wx-4*Y0*Y0*bb)/(4*Y0)
        images.append(E([9*wx+3*cc,27*wy]))
    wx=dd*dd/(4*Y0*Y0)-cc;wy=-dd*wx/(2*Y0)-Y0*bb;involution=E([9*wx+3*cc,27*wy])
    assert 6*involution==sum(images[:11],E(0));original=list(images);images[0]=involution
    def infinity(f,weight):
        if not f:return QQ(0)
        n,d=f.numerator(),f.denominator();valuation=d.degree()-n.degree()+weight
        return None if valuation<0 else QQ(0) if valuation>0 else n.leading_coefficient()/d.leading_coefficient()
    def height(P):
        if not P:return QQ(0)
        xx,yy=P.xy();n,d=xx.numerator(),xx.denominator();twice=max(d.degree(),n.degree()-4);comp=0
        if infinity(xx,4)==node and infinity(yy,6)==0:
            dx=xx-node*T**4;vx=10**6 if not dx else 4+dx.denominator().degree()-dx.numerator().degree();vy=10**6 if not yy else 6+yy.denominator().degree()-yy.numerator().degree()
            if min(vx,vy)>=2:comp=2
            else:
                slope=infinity(yy/dx,2);assert slope in (tangent,-tangent);comp=1 if slope==tangent else 3
        return 4+twice-QQ(comp*(4-comp))/4
    H=matrix(QQ,12)
    for i,P in enumerate(images):
        H[i,i]=height(P)
        for j in range(i):H[i,j]=H[j,i]=(height(P+images[j])-H[i,i]-H[j,j])/2
    row.update(status='EXACT_GEOMETRY_PENDING_REPLAY',generic_sections=[[str(c) for c in P.xy()] for P in images],original_sections=[[str(c) for c in P.xy()] for P in original],height_Gram=[list(map(str,r)) for r in H.rows()],height_determinant=str(H.det()),generic_Q_rank_lower_bound=int(H.rank()) if H.is_positive_definite() else 'UNKNOWN',surface_counts=[count(A,B,p) for p in (131,239,251)])
    return row
def main():
    protocol=json.loads((D/'protocol.json').read_text())
    for n,h in protocol['sources'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h
    rows=[]
    for i,v in enumerate(protocol['parent_ratios']):
        row=construct(QQ(v));write(D/('parent'+str(i)+'.json'),row);rows.append(row)
        print('parent',v,row['status'],'rank',row.get('generic_Q_rank_lower_bound'),'det',row.get('height_determinant'),flush=True)
    write(D/'result.json',{'schema':'kihara-positive-parent-intake.v1','status':'PENDING_INDEPENDENT_REPLAY','rows':rows,'sources':protocol['sources'],'scope':'Three fixed positive ratios outside the retained rank14 path. Exact generic sections and finite surface fingerprints; no new parent inequivalence claim until independent replay and comparison, no full generic rank or full NS claim, no fibre point-search exposure.'})
if __name__=='__main__':main()
