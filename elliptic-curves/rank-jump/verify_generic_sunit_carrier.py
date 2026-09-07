#!/usr/bin/env python3
"""Independent cyclic-ring, regularized-symbol and generic-factor certificates."""
import argparse
from pathlib import Path
import retrospective as r
import generic_sunit_carrier as run
from verify_half_ideal_artin import jacobi

OUTPUT=r.OUT/'rank_jump_generic_sunit_carrier_verification_v1.json'


def compute():
    from sage.all import QQ,ZZ,GF,NumberField,PolynomialRing,matrix,vector,pari
    source=r.read(run.half.INPUT);out=r.read(run.OUTPUT);old=r.read(run.old_artin.OUTPUT)
    for path,sha in out['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    verified=[]
    for data in out['rows']:
        assert data['status']=='PASS'
        c=next(z for z in source['cases'] if z['case_index']==data['case_index'])
        oldrow=next(z for z in old['rows'] if z['case_index']==data['case_index'])
        factor=r.read(run.rem.INPUT)['cases'][data['case_index']]['factor']
        ps=[p for p,e in factor['factors']];R=PolynomialRing(QQ,'z');f=R(c['integral_cubic_ascending'])
        assert factor['factorization_complete'] and __import__('math').prod(p**e for p,e in factor['factors'])==abs(16*int(f.discriminant()))
        K=NumberField(f,'theta');theta=K.gen();basis=[K(R(z)) for z in c['maximal_order_basis']]
        B=matrix(QQ,[list(z) for z in basis]).transpose();Binv=B.inverse()
        assert B.det()**2*f.discriminant()==QQ(c['field_discriminant'])
        nf=pari.nfinit([pari(f),ps]);assert str(nf.disc())==c['field_discriminant']
        primes=[(p,j,P) for p in ps for j,P in enumerate(pari.idealprimedec(nf,p))]
        pts=c['points'];n=len(pts);m=data['generic_dimension'];k=data['strict_dimension']
        bits=matrix(GF(2),n,n);symbols_checked=0;regularized=0
        for j,column in enumerate(data['columns']):
            assert column['point_position']==j
            H=matrix(QQ,column['coprime_hnf']);N=ZZ(column['norm'])
            assert H.det()==N and N%2 and all(N%p for p in ps)
            assert H[0,0]==N and H[1,1]==H[2,2]==1 and all(H[a,b]==0 for a in range(3) for b in range(a))
            residues=vector(QQ,[1,-H[0,1],-H[0,2]])
            for a in range(3):
                for b in range(3):
                    cc=Binv*vector(QQ,list(basis[a]*basis[b]));assert all(x in ZZ for x in cc)
                    assert ZZ(residues.dot_product(cc)-residues[a]*residues[b])%N==0
            root=ZZ(residues.dot_product(Binv*vector(QQ,list(theta))))%N
            assert root==ZZ(column['root_residue']) and f(root)%N==0 and ZZ(f.derivative()(root)).gcd(N)==1
            product=pari(H)
            for rec,(p,h,P) in zip(column['removed_S_factors'],primes):
                assert rec[:2]==[p,h]
                if rec[2]:product=pari.idealmul(nf,product,pari.idealpow(nf,P,rec[2]))
            assert matrix(QQ,product)==matrix(QQ,pts[j]['gcd_ideal_hnf'])
            # Recheck the localized half-ideal identity, allowing all S valuations.
            gamma=K(pts[j]['gamma_coordinates']);I2=pari.idealpow(nf,pari(H),2)
            quotient=pari.idealdiv(nf,I2,pari(R(list(gamma))))
            for p,h,P in primes:
                e=int(pari.idealval(nf,quotient,P))
                if e:quotient=pari.idealmul(nf,quotient,pari.idealpow(nf,P,-e))
            assert pari.idealhnf(nf,quotient)==pari.idealhnf(nf,1)
            for i,(P,rec) in enumerate(zip(pts,column['point_evaluations'])):
                a,b,d=[ZZ(P[z]) for z in ('a','b','d')]
                gamma=a-d*d*theta;other=a*a+a*d*d*theta+d**4*(theta**2+K(f[1]))
                assert gamma*other==b*b  # Exact identity underlying every nonunit regularization.
                value=ZZ((a-d*d*root)%N);u=ZZ(rec['unit_modulus']);v=ZZ(rec['derivative_modulus'])
                assert value==ZZ(rec['residue']) and u*v==N and u.gcd(v)==1 and value.gcd(u)==1
                leftover=v
                while leftover>1:
                    g=leftover.gcd(value);assert g>1;leftover//=g
                derivative=3*a*a+ZZ(f[1])*d**4
                assert derivative==ZZ(rec['derivative_numerator']) and derivative.gcd(v)==1 and d.gcd(v)==1
                sign=jacobi(int(value),int(u))*jacobi(int(derivative),int(v))
                assert sign==rec['symbol'] and sign in (-1,1)
                bits[i,j]=int(sign==-1);symbols_checked+=1;regularized+=int(v>1)
        masks=data['strict_masks'];C=matrix(GF(2),[[(mask>>i)&1 for i in range(n)] for mask in masks])
        evaluations=C*bits
        assert evaluations==matrix(GF(2),[c['strict_character_bits'] for c in data['columns']]).transpose()
        known=evaluations*C.transpose();assert known==matrix(GF(2),oldrow['matrix_rows'])
        generic=evaluations.matrix_from_columns(range(m))
        assert generic.rank()==k==data['generic_detected_rank']
        assert known.rank()==data['strict_detected_rank']
        assert generic.augment(known).rank()==k and data['strict_dimensions_outside_generic_plus_S_units_lower_bound']==0
        # Explicit elementary factor generated entirely by generic half ideals.
        pivots=list(generic.pivots());assert len(pivots)==k
        square=generic.matrix_from_columns(pivots);inverse=square.inverse()
        dual=[]
        for i in range(k):
            v=vector(GF(2),m)
            for j,pivot in enumerate(pivots):v[pivot]=inverse[j,i]
            assert generic*v==vector(GF(2),[int(j==i) for j in range(k)])
            dual.append(sum(int(v[j])<<j for j in range(m)))
        verified.append({'case_index':data['case_index'],'id':data['id'],'point_symbols_verified':symbols_checked,
            'nonunit_symbols_regularized':regularized,'old_strict_Artin_entries_reproduced':k*k,
            'generic_pivot_positions':pivots,'dual_generic_half_ideal_masks':dual,
            'elementary_S_class_factor_from_generic_half_ideals_dimension':k,
            'generic_strict_dimension':oldrow['generic_strict_dimension'],
            'additional_strict_character_dimensions_lower_bound':k-oldrow['generic_strict_dimension'],
            'strict_dimensions_outside_generic_plus_S_units_detected':0})
    return {'schema':'rank-jump.generic-sunit-carrier-verification.v1','status':'PASS','rows':verified,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),run.OUTPUT,run.half.INPUT,run.old_artin.OUTPUT,Path(__file__).with_name('verify_half_ideal_artin.py'))},
        'boundary':'Retrospective oracle characters certify elementary factors generated by generic half ideals. Neither class identities modulo this factor nor an equation-only independence proof are supplied.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['build','check']);args=parser.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS',[(x['id'],x['elementary_S_class_factor_from_generic_half_ideals_dimension']) for x in result['rows']])
