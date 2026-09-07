#!/usr/bin/env python3
"""Verify the twist bounds, transported CT block and explicit torsor models."""
import argparse
from pathlib import Path
from fractions import Fraction as Q
import retrospective as r
import fixed_cubic_minus_reference as source
import minus_reference_carriers as covers
import verify_bounded_gain_reference as upstream
from verify_unpointed_governing_norm import Algebra

OUTPUT=r.OUT/'rank_jump_fixed_cubic_minus_reference_verification_v1.json'


def determinant(M):
    A=[list(map(Q,row)) for row in M];d=Q(1)
    for i in range(len(A)):
        j=next((j for j in range(i,len(A)) if A[j][i]),None)
        if j is None:return Q(0)
        if i!=j:A[i],A[j]=A[j],A[i];d=-d
        p=A[i][i];d*=p
        for j in range(i+1,len(A)):
            s=A[j][i]/p
            for k in range(i,len(A)):A[j][k]-=s*A[i][k]
    return d


def compute():
    from sage.all import QQ,ZZ,pari,PolynomialRing,AA,GF,matrix
    # This replays all six original strict ideals and every repaired Artin entry.
    assert upstream.compute()==r.read(upstream.OUTPUT)
    data=r.read(source.OUTPUT);report=r.read(source.REPORT);carrier=r.read(covers.OUTPUT)
    for obj in (data,report,carrier):
        for p,sha in obj['bindings'].items():assert r.digest((r.ROOT/p).read_bytes())==sha
    R=PolynomialRing(QQ,'z');fs=[R(row['cubic_ascending']) for row in data['arms']]
    assert fs[1]==-fs[0](-R.gen()) and fs[1].discriminant()==fs[0].discriminant()
    old=r.read(source.completed.OUTPUT)['stages']['local'];primes=old['S_finite'];Ss=[]
    for arm,f in zip(data['arms'],fs):
        alg=Algebra(list(map(str,f.list())));delta=Q(str(f.discriminant()));nf=pari.nfinit([pari(f),primes])
        E=pari.ellinit([0,0,0,pari(f[1]),pari(f[0])]);S=[2];ell=1
        for row in arm['local']:
            p=row['prime'];cv=int(pari.elllocalred(E,p)[0]);assert cv==row['conductor_valuation']
            if p!=2 and cv:S.append(p)
            dec=pari.idealprimedec(nf,p);dim=len(dec)-1+int(p==2)
            assert dim==row['local_point_dimension'] and [[int(P[2]),int(P[3])] for P in dec]==row['splitting']
            if p==2 or cv:ell+=dim
        assert S==arm['S_finite'] and ell==arm['local_point_product_dimension']==11
        assert arm['Selmer_boundary_upper_bound']==ell-1==10
        beta=alg.elt(arm['derivative_coefficients']);assert beta==alg.elt([-delta*alg.f[1],0,-3*delta])
        assert alg.norm(beta)==delta**4==Q(arm['derivative_norm'])
        roots=f.roots(AA,multiplicities=False);assert len(roots)==3
        assert [int(-f.discriminant()*f.derivative()(x)<0) for x in roots]==[1,0,1]
        assert str(nf.disc())==arm['field_discriminant'] and arm['omitted_good_prime_valuations']==[]
        Ss.append(S)
    assert Ss[0]==Ss[1]==primes
    M=data['twist_CT_matrix'];assert M==old['minus_twist_CT_matrix']
    assert matrix(GF(2),M).rank()==6==data['twist_Sha2_dimension_unconditional_lower_bound']
    assert data['transported_block_rational_dimension_on_twist']==0
    indices=[v for pair in carrier['symplectic_pairs_in_strict_basis'] for v in pair]
    T=matrix(GF(2),[[v>>i&1 for i in range(6)] for v in indices]);assert T.rank()==6
    assert T*matrix(GF(2),M)*T.transpose()==matrix(GF(2),carrier['standard_CT_matrix'])
    # Independent rational-algebra reconstruction of all quadric coefficients.
    alg=Algebra(list(map(str,fs[1].list())));theta=alg.elt([0,1,0]);coefficient_checks=0
    for row,rec in zip(carrier['rows'],old['class_records']):
        beta=alg.elt([Q(c)*(-1)**j for j,c in enumerate(rec['beta_ascending'])])
        assert beta==alg.elt(row['beta_minus_ascending'])
        assert alg.norm(beta)>0
        assert ZZ(alg.norm(beta).numerator).is_square() and ZZ(alg.norm(beta).denominator).is_square()
        Hs=[[[Q(c) for c in rr] for rr in H] for H in row['quadric_matrices']]
        for index,H in enumerate(Hs):
            for i in range(4):
                for j in range(4):
                    expected=alg.mul(beta,alg.power(theta,i+j))[index+1] if i<3 and j<3 else Q(int(index==0 and i==j==3))
                    assert H[i][j]==expected;coefficient_checks+=1
        pencil=list(map(Q,row['pencil_ascending']));assert len(pencil)==4 and pencil[3]
        for t in range(5):assert determinant([[Hs[0][i][j]+t*Hs[1][i][j] for j in range(4)] for i in range(4)])==sum(c*t**i for i,c in enumerate(pencil))
        d,c,b,a=pencil;disc=b*b*c*c-4*a*c**3-4*b**3*d-27*a*a*d*d+18*a*b*c*d
        assert disc==Q(row['pencil_discriminant']) and disc
        assert row['genus']==1 and row['projective_degree']==4
    result=report['result'];assert result['shared_strict_Selmer_dimension_conditional']==12
    assert result['twist_rank_conditional_upper_bound']==12+10-6==16
    assert result['rank_drop_conditional_lower_bound']==22-16==6
    assert result['twist_strict_rational_dimension_conditional_upper_bound']==12-6==6
    assert result['twist_Selmer_dimension_conditional_interval']==[12,22]
    files=(Path(__file__),source.OUTPUT,source.REPORT,covers.OUTPUT,upstream.OUTPUT,
           Path(upstream.__file__),Path(__file__).with_name('verify_unpointed_governing_norm.py'))
    return {'schema':'rank-jump.fixed-cubic-minus-reference-verification.v1','status':'PASS',
        'upstream_strict_Artin_replay':True,'local_arms':2,'common_boundary_upper_bound':10,
        'unconditional_CT_rank':6,'quadric_coefficient_checks':coefficient_checks,
        'independent_pencil_determinants':30,'smooth_genus_one_covers':6,
        'conditional_twist_rank_upper_bound':16,
        'carrier_minimality':'Mathematical argument in proof note; numerical checks certify its six independent torsor-class hypothesis.',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print({k:v for k,v in result.items() if k!='bindings'},flush=True)
