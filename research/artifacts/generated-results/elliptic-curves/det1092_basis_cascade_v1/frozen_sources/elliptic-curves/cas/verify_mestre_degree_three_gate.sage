#!/usr/bin/env sage-python
"""Standalone rational-curve, divisor and exhaustive frame-root proof."""
import argparse,json,math
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector

def enumerate_exact(gram,bound):
    """All lattice vectors by rational LDL and exact recursive intervals."""
    n=len(gram);u=[[Q(int(i==j)) for j in range(n)] for i in range(n)];d=[]
    for i in range(n):
        di=Q(gram[i][i])-sum(d[k]*u[k][i]**2 for k in range(i));assert di>0;d.append(di)
        for j in range(i+1,n):u[i][j]=(Q(gram[i][j])-sum(d[k]*u[k][i]*u[k][j] for k in range(i)))/di
    v=[0]*n;out=set();nodes=0
    def descend(i,left):
        nonlocal nodes
        nodes+=1
        if i<0:
            if any(v):out.add(tuple(v))
            return
        shift=sum(u[i][j]*v[j] for j in range(i+1,n));radius=left/d[i]
        # An integer upper bound on sqrt(radius); over-enumeration is harmless.
        b=math.isqrt(radius.numerator//radius.denominator)+1
        lo=math.ceil(-shift-b);hi=math.floor(-shift+b)
        for a in range(lo,hi+1):
            term=d[i]*(a+shift)**2
            if term<=left:v[i]=a;descend(i-1,left-term)
        v[i]=0
    descend(n-1,Q(bound))
    return out,nodes

def main(path):
    bundle=json.loads(path.read_text());ns=bundle['ns'];h=bundle['heights'];frames=bundle['frames'];pencil=bundle['pencil']
    oldG=matrix(ZZ,ns['rational_NS_Gram']);H=matrix(QQ,h['seed_height_gram']);oldinv=oldG.inverse()
    assert oldG.det()==-468 and H.det()==QQ(117)/4
    R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();Z=PolynomialRing(QQ,'z');z=Z.gen();L=Z.fraction_field()
    A=R(h['curve_A_coefficients']);B=R(h['curve_B_coefficients'])
    oldseed=[(K(x),K(y)) for x,y in ns['section_points']]
    pole=[K(c) for c in pencil['pole_point']];den=R(pencil['q']);offset=R(pencil['c'])
    O=vector(ZZ,[int(i==1) for i in range(18)])
    def section_class(word):
        w=vector(ZZ,word)
        if not any(w):return O
        comp=[sum(w[i]*ns['section_profiles'][i]['components'][j] for i in range(11))%n for j,n in enumerate((2,2,4))]
        corr=QQ(comp[0]+comp[1])/2+QQ(comp[2]*(4-comp[2]))/4
        oo=(w*H*w-4+corr)/2;assert oo in ZZ
        cross=[]
        for i,profile in enumerate(ns['section_profiles']):
            other=profile['components']
            correction=sum(QQ(comp[j]*other[j])/2 for j in (0,1))+min(comp[2],other[2])-QQ(comp[2]*other[2])/4
            cross.append(2+oo+profile['zero_section_intersection']-(w*H)[i]-correction)
        rhs=[1,oo,*comp[:2],*[int(comp[2]==j) for j in (1,2,3)],*cross]
        v=oldinv*vector(QQ,rhs);assert all(c in ZZ for c in v) and v*oldG*v==-2
        return v
    options=[(i,j,a,b,c) for i,j in product(range(3),repeat=2) for a,b,c in product(range(3),repeat=3) if a+b+c<=2]
    def identify(oo,cross):
        result=[]
        for comp in options:
            v=oldinv*vector(QQ,[2,oo,*comp,*cross])
            if all(c in ZZ for c in v) and v*oldG*v==-2:result.append(v)
        assert len(result)==1
        return result[0]
    bases={}
    for b in bundle['base_bisections']:
        t=L(b['old_T']);x=L(b['old_x']);y=L(b['old_y']);deck=L(b['deck'])
        E=EllipticCurve(L,[A(t),B(t)]);P=E([x,y]);seed=[E([a(t),b(t)]) for a,b in oldseed]
        assert max(t.numerator().degree(),t.denominator().degree())==2 and t(deck)==t and deck!=z
        assert (den(t)*(y+pole[1](t))/(x-pole[0](t))+offset(t))/den(t)**2==z
        # The inverse equals the alternate fibration parameter: this is a
        # section of that smooth K3 fibration, hence an embedded (-2)-curve.
        def oi(xx):
            dx=xx.denominator();dt=t.denominator()
            val=dx.degree()-dx.gcd(dt**4).degree()+max(0,xx.numerator().degree()-dx.degree()-4*max(0,t.numerator().degree()-dt.degree()))
            assert val%2==0;return val//2
        v=identify(oi(x),[oi((P-Q)[0]) for Q in seed]);bases[b['index']]=(v,t,P,seed,oi)
    G=matrix(ZZ,frames['geometric_NS_gram']);action=matrix(ZZ,frames['Galois_action'])
    glue=ns['eligible_index_two_geometric_glues'][0]
    assert G==matrix(ZZ,glue['geometric_Gram']) and action==matrix(ZZ,glue['Galois_action'])
    change=matrix.identity(QQ,19);change[glue['replace_basis_index']]=vector(QQ,glue['fixed_numerator']+[1])/2
    results=[]
    for i,row in enumerate(frames['rows']):
        base,t,P,seed,oi=bases[row['source_index']];w=vector(ZZ,row['subtract_section_word']);S=section_class(w)
        assert base*oldG*S==2
        cross=[base*oldG*section_class(w+vector(ZZ,[int(k==j) for k in range(11)])) for j in range(11)]
        v=identify(2,cross);assert list(v)==row['bisection_class']
        f=O+v;assert list(f)==row['fibre_class'] and f*oldG*f==0 and f*oldG*O==0
        gf=vector(QQ,list(f)+[0])*change.inverse();gs=vector(ZZ,row['geometric_zero_section_class'])
        assert list(gf)==row['geometric_fibre_class'] and gs*G*gs==-2 and gf*G*gs==1
        assert gf*action==gf and gs*action==gs
        # Identify the claimed degree-one curve in the explicit old component list.
        old_section=vector(QQ,gs)*change
        witness=bundle['visible_curve_classes'][row['visible_section_indices'][0]]
        assert list(old_section)==witness+[0]
        F=matrix(ZZ,row['frame_basis']);Hf=matrix(ZZ,row['frame_gram']);U=matrix(ZZ,row['reduction_change']);small=matrix(ZZ,row['reduced_frame_gram'])
        assert not any(F*G*gf) and not any(F*G*gs)
        whole=matrix(ZZ,[gf,gs+gf,*F.rows()]);assert abs(whole.det())==1
        assert Hf==-F*G*F.transpose() and Hf.det()==468 and abs(U.det())==1 and small==U.transpose()*Hf*U
        exact,nodes=enumerate_exact(row['reduced_frame_gram'],2)
        retained={tuple(v) for v in row['root_vectors_up_to_sign']};retained|={tuple(-c for c in v) for v in retained}
        assert exact==retained and len(exact)==row['geometric_root_count']
        roots=matrix(ZZ,sorted(exact))*U.transpose()*F
        rg=roots.rank();rq=rg-(roots*(action-matrix.identity(ZZ,19))).rank()
        assert rg==row['geometric_root_rank'] and rq==row['rational_root_rank']
        assert (16-rq,17-rg)==(row['generic_Q_MW_rank'],row['generic_geometric_MW_rank'])
        if i==0:
            actual=P-sum((a*Q for a,Q in zip(w,seed)),P.curve()(0));given=frames['selected_bisection']
            assert str(t)==given['old_T'] and actual==P.curve()([L(given['old_x']),L(given['old_y'])])
            assert oi(actual[0])==2 and [oi((actual-Q)[0]) for Q in seed]==cross
        results.append([int(16-rq),int(17-rg),len(exact),nodes])
    print('PASS standalone exact bisections, primitive degree3 pencils and exhaustive roots',results,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);args=p.parse_args();main(args.input)
