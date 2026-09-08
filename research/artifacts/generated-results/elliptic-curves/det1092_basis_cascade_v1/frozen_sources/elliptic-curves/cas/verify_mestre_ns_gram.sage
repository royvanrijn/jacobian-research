#!/usr/bin/env sage-python
"""Independent NS intersection and unique-glue replay; only Sage and JSON."""
import argparse, json
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, vector, block_diagonal_matrix


def main(bundle, certificate):
    data=json.loads(bundle.read_text()); ns=json.loads(certificate.read_text())
    assert ns['status']=='PASS' and len(ns['rows'])==6
    common=None
    for source, row in zip(data['rows'],ns['rows']):
        assert source['outer_u']==row['outer_u']
        h=source['generic_heights'];R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field()
        A=R(h['curve_A_coefficients']);B=R(h['curve_B_coefficients'])
        E=EllipticCurve(K,[A,B]);points=[E([K(x),K(y)]) for x,y in row['section_points']]
        raw=[E([K(x),K(y)]) for x,y in h['covariant_points']]
        assert len(points)==11 and points[0]==raw[0]
        assert all(2*points[j]==raw[i]-raw[0] for j,i in enumerate(h['seed_indices']) if j)
        delta=-16*(4*A**3+27*B**2)
        assert delta.degree()==20 and delta.gcd(A).degree()==0
        double=next(f for f,m in delta.squarefree_decomposition() if m==2)
        bases=sorted(double.roots(QQ,multiplicities=False));assert len(bases)==2
        node=-3*B[12]/(2*A[8]);tangent=(3*node).sqrt();assert tangent in QQ
        assert A[7]==B[11]==0 and A[6]*node+B[10]==0
        def intersection_zero(P):
            x=P[0];n,d=x.numerator(),x.denominator()
            twice=max(d.degree(),n.degree()-4);assert twice>=0 and twice%2==0
            return twice//2
        def narrow_height(P):
            if not P:return QQ(0)
            Q=4*P;assert Q
            return QQ(4+2*intersection_zero(Q))/16
        def at_infinity(f,weight):
            # Reciprocal rational functions, evaluated at the local parameter 0.
            numerator=R(list(reversed(f.numerator().list())))
            denominator=R(list(reversed(f.denominator().list())))
            return K(T**(weight+f.denominator().degree()-f.numerator().degree())*numerator/denominator)
        G=matrix(ZZ,row['rational_NS_Gram']);assert G.dimensions()==(18,18) and G==G.transpose()
        assert G[:2,:2]==matrix(ZZ,[[0,1],[1,-2]])
        root=matrix.diagonal(ZZ,[-2]*5);root[2,3]=root[3,2]=root[3,4]=root[4,3]=1
        assert G[2:7,2:7]==root and not any(G[i,j] for i in range(2) for j in range(2,7))
        heights=[narrow_height(P) for P in points];H=matrix(QQ,11,11)
        for i,P in enumerate(points):
            H[i,i]=heights[i]
            for j in range(i):H[i,j]=H[j,i]=(narrow_height(P+points[j])-heights[i]-heights[j])/2
            assert G[7+i,7+i]==-2 and G[0,7+i]==1 and G[1,7+i]==intersection_zero(P)
            components=[]
            for base in bases:
                c=0
                if P[0].denominator()(base):
                    c=int(P[0](base)==-3*B(base)/(2*A(base)) and P[1](base)==0)
                components.append(c)
            X=at_infinity(P[0],4);Y=at_infinity(P[1],6);c=0
            if X.denominator()(0) and Y.denominator()(0) and X(0)==node and Y(0)==0:
                dx=X-node
                nx=dx.numerator().valuation() if dx else 100000
                ny=Y.numerator().valuation() if Y else 100000
                assert min(nx,ny)>=1
                if min(nx,ny)>=2:c=2
                else:
                    ratio=Y/dx;assert ratio.denominator()(0)
                    if ratio(0)==tangent:c=1
                    else:assert ratio(0)==-tangent;c=3
            components.append(c)
            assert components==row['section_profiles'][i]['components']
            expected=components[:2]+[int(c==k) for k in (1,2,3)]
            assert list(G[2:7,7+i].column(0))==expected
        assert H==matrix(QQ,h['seed_height_gram']) and H.det()==QQ(117)/4
        trivial=G[:7,:7];cross=G[7:,:7]
        assert G[7:,7:]-cross*trivial.inverse()*cross.transpose()==-H
        assert G.det()==-468 and list(G.elementary_divisors())==[1]*16+[6,78]
        split=matrix(ZZ,row['U_split_change']);assert abs(split.det())==1
        split=split*G*split.transpose();assert split[:2,:2]==matrix(ZZ,[[0,1],[1,0]])
        assert not any(split[i,j] for i in range(2) for j in range(2,18)) and (-split[2:,2:]).is_positive_definite()
        kernel=matrix(GF(2),G).left_kernel();assert kernel.dimension()==2
        eligible=[]
        for vv in kernel:
            v=vector(ZZ,[int(c) for c in vv])
            if any(v) and (v*G*v-4)%8==0:eligible.append(v)
        assert len(eligible)==1 and len(row['eligible_index_two_geometric_glues'])==1
        v=eligible[0];glue=row['eligible_index_two_geometric_glues'][0]
        assert list(v)==glue['fixed_numerator'] and all(c%2==0 for c in G*v)
        C=matrix.identity(QQ,19);j=glue['replace_basis_index'];assert v[j]==1
        C[j]=vector(QQ,list(v)+[1])/2
        extended=C*block_diagonal_matrix(G,matrix(ZZ,[[-4]]))*C.transpose()
        assert extended==matrix(ZZ,glue['geometric_Gram']) and extended.det()==468
        assert all(extended[i,i]%2==0 for i in range(19))
        assert list(extended.elementary_divisors())==[1]*17+[3,156]
        action=matrix(ZZ,glue['Galois_action'])
        assert action==C*matrix.diagonal(QQ,[1]*18+[-1])*C.inverse()
        assert action*extended*action.transpose()==extended and action**2==matrix.identity(QQ,19)
        assert 19-(action-matrix.identity(QQ,19)).rank()==18
        if common is None:common=(G,extended,action)
        else:assert common==(G,extended,action)
        print('PASS u'+row['outer_u']+' intersections, 66 heights, unique glue and Galois action',flush=True)
    print('PASS6 COMMON ARITHMETIC AND GEOMETRIC NS MATRICES',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--bundle',type=Path,default=Path('mestre_468_replay_bundle_v1.json'))
    p.add_argument('--certificate',type=Path,default=Path('mestre_rational_ns_gram_v2.json'))
    a=p.parse_args();main(a.bundle,a.certificate)
