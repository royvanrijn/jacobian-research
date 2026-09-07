#!/usr/bin/env sage-python
"""Standalone quadratic section and full geometric lattice calculation."""
import argparse,json
from pathlib import Path
from sage.all import QQ,ZZ,QuadraticField,PolynomialRing,EllipticCurve,matrix,vector
def verify(bundle):
    data=json.loads(bundle.read_text());row=data['parent'];old=data['rational_rank'];found=data['section']
    K=QuadraticField(-3,'r');r=K.gen();bar=K.hom([-r]);R=PolynomialRing(K,'T');T=R.gen();L=R.fraction_field()
    a=R(row['raw_A']);b=R(row['raw_B']);disc=-16*(4*a**3+27*b**2)
    assert (a.degree(),b.degree(),disc.degree())==(8,12,20) and disc.gcd(disc.derivative())==1 and disc.gcd(a)==1
    scale=QQ(found['coordinate_scale']);shift=QQ(found['base_shift']);X=R(found['X']);Y=R(found['Y']);s=T*T-shift
    Qx=scale**2*X(s);Qy=scale**3*Y(s)
    E=EllipticCurve(L,[a,b]);P=E([Qx,Qy]);sigma=E([R([bar(c) for c in Qx.list()]),R([bar(c) for c in Qy.list()])])
    assert sigma!=P and X.degree()==2 and Y.degree()==2
    basis=[E([L(x),L(y)]) for x,y in old['basis']]
    H0=matrix(QQ,old['basis_height_gram']);assert len(basis)==12 and H0.det()==189
    node=-3*b[12]/(2*a[8]);tangent=QQ(3*node).sqrt()
    def infinity(f,weight):
        if not f:return K(0)
        n,d=f.numerator(),f.denominator();valuation=d.degree()-n.degree()+weight
        return None if valuation<0 else K(0) if valuation>0 else n.leading_coefficient()/d.leading_coefficient()
    def profile(point):
        if not point:return (QQ(0),0,0)
        x,y=point.xy();n,d=x.numerator(),x.denominator();twice=max(d.degree(),n.degree()-4);assert twice>=0 and twice%2==0
        comp=0
        if infinity(x,4)==node and infinity(y,6)==0:
            dx=x-node*T**4
            vx=10**6 if not dx else 4+dx.denominator().degree()-dx.numerator().degree()
            vy=10**6 if not y else 6+y.denominator().degree()-y.numerator().degree()
            if min(vx,vy)>=2:comp=2
            else:
                slope=infinity(y/dx,2);assert slope in (tangent,-tangent);comp=1 if slope==tangent else 3
        return (QQ(4+twice)-QQ(comp*(4-comp))/4,int(twice//2),comp)
    points=basis+[P];profiles=[profile(Q) for Q in points];H=matrix(QQ,13)
    for i,Q in enumerate(points):
        H[i,i]=profiles[i][0];assert profile(4*Q)[0]==16*H[i,i]
        for j in range(i):
            plus=profile(Q+points[j])[0];minus=profile(Q-points[j])[0]
            assert plus+minus==2*(profiles[i][0]+profiles[j][0]);H[i,j]=H[j,i]=(plus-profiles[i][0]-profiles[j][0])/2
    assert H[:12,:12]==H0 and H.is_positive_definite()
    anti=P-sigma;trace=P+sigma;anti_profile=profile(anti)
    word=H0.solve_right(2*vector(QQ,H[:12,12].column(0)));assert all(c in ZZ for c in word)
    assert trace==sum((ZZ(c)*Q for c,Q in zip(word,basis)),E(0))
    # Any nonzero geometric section has height >=3 on20I1+I4.
    # Thus anti of height4 is primitive. For any geometric R, R-bar(R)
    # is an integer multiple of anti; subtract that multiple of P to get
    # a rational section in the previously certified full old basis.
    assert anti_profile[0]==4 and profiles[-1]==(3,0,2)
    G=matrix(QQ,18);G[0,1]=G[1,0]=1;G[1,1]=-2
    for i in range(2,5):G[i,i]=-2
    for i in (2,3):G[i,i+1]=G[i+1,i]=1
    for j,(hh,oo,c) in enumerate(profiles):
        idx=5+j;G[idx,idx]=-2;G[0,idx]=G[idx,0]=1;G[1,idx]=G[idx,1]=oo
        for k in (1,2,3):G[1+k,idx]=G[idx,1+k]=int(c==k)
        for i in range(j):
            c0=profiles[i][2];correction=min(c,c0)-QQ(c*c0)/4
            value=2+oo+profiles[i][1]-H[i,j]-correction
            G[5+i,idx]=G[idx,5+i]=value
    assert all(c in ZZ for c in G.list()) and G.det()==-756 and H.det()==189
    sigmaprofile=profile(sigma);assert sigmaprofile==profiles[-1]
    rhs=list(G[:,17].column(0));rhs[-1]=2-profile(P+sigma)[0]/2+3-1
    # Direct formula: <P,sigma P>=(h(P+sigma P)-6)/2;
    # intersection=2-<P,sigma P>-1.
    rhs[-1]=1-(profile(P+sigma)[0]-6)/2
    action=matrix.identity(QQ,18);action[17]=G.solve_right(vector(QQ,rhs))
    assert all(c in ZZ for c in action.list()) and action**2==matrix.identity(QQ,18) and action*G*action.transpose()==G
    assert 18-(action-matrix.identity(QQ,18)).rank()==17
    anti_x=anti[0];anti_y=anti[1]/r
    assert all(c in QQ for c in anti_x.numerator().list()+anti_x.denominator().list()+anti_y.numerator().list()+anti_y.denominator().list())
    return {'status':'PASS','constant_field':'Q(sqrt(-3))','generic_rational_MW_rank':12,'generic_geometric_MW_rank':13,
        'full_geometric_NS_Gram':[list(map(int,row)) for row in G.rows()],'geometric_NS_determinant':-756,
        'geometric_NS_discriminant_group':list(map(int,matrix(ZZ,G).elementary_divisors())),
        'Galois_action':[list(map(int,row)) for row in action.rows()],'Galois_fixed_NS_rank':17,
        'full_geometric_MW_height_Gram':[list(map(str,row)) for row in H.rows()],'height_determinant':'189',
        'section_height':'3','anti_height':'4','trace_word':list(map(int,word)),
        'new_section':[str(Qx),str(Qy)],'anti_x':str(anti_x),'anti_y_div_sqrt_minus3':str(anti_y),
        'scope':'First Kihara parent only. Exact nonrational section, Galois trace and primitive anti-direction give a full geometric MW basis from the previously proved full rational basis and exact geometric rank13. The resulting full geometric NS determinant is-756, with fixed rank17. No rational rank gain, new parent, parameter scan, or universal parent-family field identity.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path);p.add_argument('--check',type=Path);a=p.parse_args();r=verify(a.input)
    if a.output:
        with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    if a.check:assert r==json.loads(a.check.read_text())
    print('PASS first Kihara full geometric MW13, NS determinant-756, exact constant field Q(sqrt(-3))',flush=True)
