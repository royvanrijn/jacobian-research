#!/usr/bin/env sage-python
"""Standalone exact replay of the six determinant468 parent certificates.

Requires only Sage and the adjacent JSON bundle; no repository imports.
Standard geometric lemmas and the Artin-Tate argument are written in the
canonical proof note. This script checks their arithmetic witnesses.
"""
import argparse,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,lcm,gcd

def finite_rank(model,points,primes,ell):
    E=EllipticCurve(QQ,model);P=[E([QQ(x),QQ(y)]) for x,y in points];projective=[]
    for point in P:
        den=lcm([c.denominator() for c in point]);v=[ZZ(c*den) for c in point];g=gcd(v);projective.append([c//g for c in v])
    rows=[];records=[]
    for prime in primes:
        F=GF(prime);e=EllipticCurve(F,[F(c) for c in model])
        if not e.discriminant():raise ArithmeticError('good finite specialization required')
        key=lambda point:tuple(int(c) for c in point)
        elements=e.points();multiples={key(ell*P):ell*P for P in elements};mask={key(P):0 for P in multiples.values()};reps=[e(0)]
        while len(mask)<len(elements):
            P=next(P for P in elements if key(P) not in mask);old=list(reps);size=len(old)
            for digit in range(1,ell):
                for i,R in enumerate(old):
                    rep=R+digit*P;reps.append(rep)
                    for T in multiples.values():
                        k=key(rep+T)
                        if k in mask:raise ArithmeticError('quotient cosets overlap')
                        mask[k]=i+size*digit
        dimension=ZZ(len(reps)).valuation(ell)
        if ell**dimension!=len(reps) or dimension>2:raise ArithmeticError('elliptic quotient dimension differs')
        reduced=[e([F(c) for c in P]) for P in projective]
        for j in range(dimension):rows.append([(mask[key(P)]//ell**j)%ell for P in reduced])
        records.append({'prime':prime,'group_order':len(elements),'ell_multiple_subgroup_order':len(multiples),'quotient_dimension':int(dimension)})
    rank=int(matrix(GF(ell),rows).rank())
    if rank!=11:raise ArithmeticError('selected generic seed is not injective in finite ell quotients')
    return {'modulus':ell,'rank':rank,'groups':records}

def main(path):
    data=json.loads(path.read_text());certificates=data['certificates'];verified=0
    for row in data['rows']:
        u=row['outer_u'];h=row['generic_heights'];seed=row['specialized_seed']
        R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();A=R(h['curve_A_coefficients']);B=R(h['curve_B_coefficients']);delta=-16*(4*A**3+27*B**2)
        assert A.degree()==8 and B.degree()==12 and delta.degree()==20 and delta.gcd(A).degree()==0
        assert sorted((int(f.degree()),int(m)) for f,m in delta.squarefree_decomposition())==[(2,2),(16,1)]
        assert not any(A[i] for i in range(1,9,2)) and not any(B[i] for i in range(1,13,2))
        a=R([A[2*i] for i in range(5)]);b=R([B[2*i] for i in range(7)]);dd=-16*(4*a**3+27*b**2)
        assert sorted((int(f.degree()),int(m)) for f,m in dd.squarefree_decomposition())==[(1,2),(8,1)] and dd(0)!=0
        double=next(f for f,m in delta.squarefree_decomposition() if m==2)
        assert len(double.roots(QQ,multiplicities=False))==2
        E=EllipticCurve(K,[A,B]);P=[E([K(x),K(y)]) for x,y in h['covariant_points']]
        def height(point):
            if point.is_zero():return QQ(0)
            Q=4*point;assert not Q.is_zero();x=Q[0];n,d=x.numerator(),x.denominator()
            twice=max(int(d.degree()),int(n.degree())-4);assert twice>=0 and twice%2==0
            return QQ(4+twice)/16
        H=[height(point) for point in P];G=matrix(QQ,14,14)
        for i in range(14):
            G[i,i]=H[i]
            for j in range(i):G[i,j]=G[j,i]=(height(P[i]+P[j])-H[i]-H[j])/2
        assert G==matrix(QQ,h['covariant_height_gram']) and G.rank()==11
        for word in h['generic_covariant_relations']:
            assert sum((ZZ(c)*point for c,point in zip(word,P)),E(0))==E(0)
        C=matrix(QQ,h['divisor_change_matrix']);indices=h['seed_indices'];C0=C.matrix_from_rows_and_columns(indices,indices)
        H=C*G*C.transpose();basis_height=H.matrix_from_rows_and_columns(indices,indices)
        assert basis_height==matrix(QQ,h['seed_height_gram']) and basis_height.is_positive_definite() and basis_height.det()==QQ(117)/4
        scale=(A(1)/QQ(seed['curve'][3])).sqrt().sqrt();assert scale in QQ and B(1)==scale**6*QQ(seed['curve'][4])
        for point,original in zip(P,seed['covariant_images']):
            assert point[0](1)==scale**2*QQ(original[0]) and point[1](1)==scale**3*QQ(original[1])
        specialized=EllipticCurve(QQ,seed['curve']);images=[specialized([QQ(c) for c in P]) for P in seed['covariant_images']]
        cloud=[specialized([QQ(c) for c in P]) for P in seed['divisor_cloud']]
        assert cloud[0]==images[0]
        assert all(2*cloud[i]==images[i]-images[0] for i in range(1,14))
        assert [seed['divisor_cloud'][i] for i in indices]==seed['points']
        sigma=next(r for r in certificates['mestre_base_involution_v1.json']['rows'] if r['outer_u']==u)
        M=matrix(QQ,sigma['action_matrix']);raw_basis=[P[i] for i in indices]
        for j,relation in enumerate(sigma['exact_group_relations']):
            target=E([c(-T) for c in raw_basis[j].xy()])
            assert relation['multiplier']*target==sum((ZZ(c)*point for c,point in zip(relation['coefficients'],raw_basis)),E(0))
            assert list(M.column(j))==[QQ(c)/relation['multiplier'] for c in relation['coefficients']]
        identity=matrix.identity(QQ,11);assert M*M==identity and 11-(M-identity).rank()==5 and 11-(M+identity).rank()==6
        info=next(r for r in certificates['mestre_parent_picard_and_saturation_v1.json']['rows'] if r['outer_u']==u)
        for finite in info['finite_saturation_injections']:
            actual=finite_rank(seed['curve'],seed['points'],[g['prime'] for g in finite['groups']],finite['modulus'])
            assert actual==finite
        classes=[];arithmetic_bounds=[]
        for counts in row['finite_surface_counts']:
            p=counts['prime'];totals=[]
            for record in counts['counts']:
                degree=record['extension_degree'];q=p**degree
                if degree==1:F=GF(p);z=None
                else:
                    S=PolynomialRing(GF(p),'z');F=GF(q,'z',modulus=S.gen()**2-record['quadratic_nonsquare']);z=F.gen()
                Qring=PolynomialRing(F,'t');tvar=Qring.gen();aa=Qring(A);bb=Qring(B);disc=-16*(4*aa**3+27*bb**2)
                assert aa.degree()==8 and bb.degree()==12 and disc.degree()==20 and disc.gcd(aa).degree()==0
                assert sorted((int(f.degree()),int(m)) for f,m in disc.squarefree_decomposition())==[(2,2),(16,1)]
                repairs=[];total=0;smooth=0
                for i,expected in enumerate(record['weierstrass_fibre_counts']):
                    t=None if i==q else (F(i) if degree==1 else F(i%p)+F(i//p)*z)
                    av,bv=(aa[8],bb[12]) if t is None else (aa(t),bb(t))
                    if 4*av**3+27*bv**2:
                        actual=int(EllipticCurve(F,[av,bv]).cardinality());correction=0;smooth+=1
                    else:
                        actual=1
                        for x in F:
                            value=x**3+av*x+bv
                            actual+=1 if value==0 else (2 if value.is_square() else 0)
                        if t is None:
                            cubic=tvar**3+av*tvar+bv;node=cubic.gcd(cubic.derivative()).monic()
                            assert node.degree()==1 and (3*(-node[0])).is_square();correction=3*q
                        else:
                            f=disc;mult=0
                            while f(t)==0:f=f//(tvar-t);mult+=1
                            assert mult in (1,2);correction=(mult-1)*q
                    assert actual==expected
                    total+=actual+correction
                    if correction:repairs.append({'base_index':i,'correction':correction})
                assert total==record['surface_point_count'] and repairs==record['resolution_corrections'] and smooth==record['smooth_fibres_independently_counted']
                totals.append(total)
            trace=totals[0]-1-p*p-18*p;second=totals[1]-1-p**4-18*p*p;mid=QQ(trace*trace-second)/2;assert mid.denominator()==1
            sign=1 if mid else -1
            if not mid:assert p<abs(trace)<=2*p
            S=PolynomialRing(QQ,'z');z=S.gen();f=z**4-trace*z**3+mid*z*z-sign*p*p*trace*z+sign*p**4;ones=0;algebraic=0
            for s in (1,-1):
                while f(s*p)==0:
                    f=f//(z-s*p);algebraic+=1
                    if s==1:ones+=1
            assert algebraic==2 and f.degree()==2 and f[0]==p*p
            pair_trace=-f[1];assert pair_trace not in (-2*p,-p,0,p,2*p) and abs(pair_trace)<2*p
            classes.append(pair_trace*pair_trace-4*p*p);arithmetic_bounds.append(18+ones)
        assert min(arithmetic_bounds)==18 and len(classes)==2 and not QQ(classes[0]/classes[1]).is_square()
        assert info['geometric_NS_rank']==19 and info['arithmetic_NS_rank']==18 and info['arithmetic_generic_MW_rank']==11 and info['geometric_generic_MW_rank']==12
        assert [i for i in range(1,22) if 468%(i*i)==0]==[1,2,3,6]
        # The finite injections exclude every nontrivial index. The quotient
        # rank6 lattice has determinant1/4. Its fixed rank5 determinant1/2
        # determines a narrow anti generator of height2 and index2 glue.
        action=C0.transpose().inverse()*M*C0.transpose();assert all(c.denominator()==1 for c in action.list())
        fixed=matrix(ZZ,action-identity).right_kernel().basis_matrix()
        fixed_det=(fixed*basis_height*fixed.transpose()).det()/2**5;assert fixed_det==QQ(1)/2
        eligible=[(i,QQ(i*i)/4/fixed_det) for i in (1,2) if QQ(i*i)/4/fixed_det>=1];assert eligible==[(2,QQ(2))]
        geo=next(r for r in certificates['mestre_geometric_discriminant_v2.json']['rows'] if r['outer_u']==u)
        assert geo['geometric_NS_discriminant']==468*4//(2*2)==468 and geo['arithmetic_NS_discriminant']==-468
        verified+=1;print('PASS u'+u,'714-height-chain / counts / rank / saturation / discriminant witnesses',flush=True)
    assert verified==6;print('PASS6 INDEPENDENT DET468 PARENT REPLAYS',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=Path(__file__).with_name('mestre_468_replay_bundle_v1.json'));a=p.parse_args();main(a.input)
