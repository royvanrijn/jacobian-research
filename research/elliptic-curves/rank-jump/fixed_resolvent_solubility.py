#!/usr/bin/env python3
"""Field incidence versus exact-resolvent rational solubility."""
import argparse
from pathlib import Path
import retrospective as r
import bounded_gain_reference_completion as completed
import fixed_cubic_minus_reference as switch

PROTOCOL=Path(__file__).with_name('FIXED_RESOLVENT_SOLUBILITY_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_fixed_resolvent_solubility_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_fixed_resolvent_solubility_v1.json'


def bindings(paths):return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def export():
    old=r.read(completed.OUTPUT)['stages']['local'];ct=r.read(switch.OUTPUT)
    assert old['strict_generic_dimension']==ct['twist_CT_restriction_rank']==6
    r.write_new(INPUT,{'schema':'rank-jump.fixed-resolvent-solubility-inputs.v1',
        'cubic_ascending':ct['arms'][0]['cubic_ascending'],
        'field_discriminant':old['field_discriminant'],
        'classes':[{'generic_mask':x['generic_mask'],'alpha':x['beta_ascending']} for x in old['class_records']],
        'twist_CT_matrix':ct['twist_CT_matrix'],
        'bindings':bindings([Path(__file__),PROTOCOL,completed.OUTPUT,switch.OUTPUT]),
        'boundary':'Projection of equation and already certified generic strict classes only; no exceptional coordinates, classes or outcome labels.'})


def universal():
    from sage.all import QQ,PolynomialRing,matrix
    R=PolynomialRing(QQ,names=('A','B','a','b','c','p','q','w','s','v','lam'));A,B,a,b,c,p,q,w,s,v,lam=R.gens()
    P=PolynomialRing(R,'X');X=P.gen()
    C=matrix(R,[[0,0,-B],[1,0,-A],[0,1,0]]);I=matrix.identity(R,3)
    H=a*I+b*C+c*C**2;T=H.trace();S=(T*T-(H*H).trace())/2;N=H.det()
    change=matrix(R,[[1,a,(H*H)[0,0]],[0,b,(H*H)[1,0]],[0,c,(H*H)[2,0]]]).det()
    assert change==b**3+A*b*c*c+B*c**3
    minimal=X**3-T*X**2+S*X-N
    assert minimal.discriminant()==(-4*A**3-27*B**2)*change**2
    quartic=X**4+p*X**2+q*X+w
    resolvent=X**3-p*X**2-4*w*X+4*p*w-q*q
    assert quartic.discriminant()==resolvent.discriminant()
    tracefree=resolvent(X+p/3)
    assert tracefree==X**3+(-p*p/3-4*w)*X-2*p**3/27+8*p*w/3-q*q
    # A representative alpha=s+lam*theta has this norm. Its quartic's trace-free resolvent is lam*theta.
    normlinear=s**3+lam*lam*A*s-lam**3*B
    assert (s*I+lam*C).det()==normlinear
    qq=X**4-3*s*X**2/2-v*X-3*s*s/16-lam*lam*A/4
    rr=X**3-qq[2]*X**2-4*qq[0]*X+4*qq[2]*qq[0]-qq[1]**2
    rem=rr(X+qq[2]/3)-(X**3+lam*lam*A*X+lam**3*B)
    assert rem==normlinear-v*v
    return {'general_quartic':'Z^4 - T/2*Z^2 - n*Z + (T^2-4*S)/16',
        'T':str(T),'S':str(S),'norm_alpha':str(N),'change_of_generator_determinant':str(change),
        'resolvent_identity':'R_Q(Z)=minpoly_alpha(Z+T/2)',
        'tracefree_resolvent_generator':'alpha - Trace(alpha)/3',
        'polynomial_discriminant_ratio_to_original_cubic':'Norm(b-c*theta)^2',
        'scaled_generator_norm':str(normlinear),
        'scaled_generator_quartic':str(qq),'remaining_resolvent_coefficient_condition':str(rem),
        'fixed_oriented_resolvent_gate':'alpha*z^2=s+lam*theta iff alpha is rational on E^(-lam), via v^2=s^3+lam^2*A*s-lam^3*B; lam nonzero',
        'scope':'Exact oriented trace-free generator. Arbitrary isomorphic cubic-field generators are not restricted by this gate.'}


def compute():
    from sage.all import QQ,GF,PolynomialRing,matrix
    data=r.read(INPUT);P=PolynomialRing(QQ,'X');X=P.gen();f=P(data['cubic_ascending']);A=f[1];B=f[0]
    assert f.degree()==3 and f[2]==0 and f.is_irreducible() and not f.discriminant().is_square()
    C=matrix(QQ,[[0,0,-B],[1,0,-A],[0,1,0]]);I=matrix.identity(QQ,3)
    M=matrix(GF(2),data['twist_CT_matrix']);assert M.rank()==6 and M==M.transpose() and all(M[i,i]==0 for i in range(6))
    exclusions=[]
    for word in range(1,64):
        witness=next(j for j in range(6) if sum(((word>>i)&1)*int(M[i,j]) for i in range(6))%2)
        exclusions.append({'word':word,'CT_witness_basis_index':witness,'CT_value':1})
    rows=[]
    for i,rec in enumerate(data['classes']):
        a,b,c=map(QQ,rec['alpha']);H=a*I+b*C+c*C**2
        T=H.trace();S=(T*T-(H*H).trace())/2;N=H.det();assert N>0 and N.is_square();n=N.sqrt()
        change=b**3+A*b*c*c+B*c**3;assert change!=0 and c!=0
        q=X**4-T*X**2/2-n*X+(T*T-4*S)/16
        m=X**3-T*X**2+S*X-N
        resolvent=X**3-q[2]*X**2-4*q[0]*X+4*q[2]*q[0]-q[1]**2
        assert resolvent==m(X+T/2)
        tracefree=resolvent(X+q[2]/3);assert tracefree==m(X+T/3)
        assert q.discriminant()==m.discriminant()==f.discriminant()*change**2
        assert m(H)==matrix(QQ,3,3)
        witnesses={}
        for p in r.primes(503):
            try:qp=q.change_ring(GF(p));fp=f.change_ring(GF(p))
            except (ValueError,ZeroDivisionError):continue
            if not qp.discriminant() or not fp.discriminant():continue
            if qp.is_irreducible() and 'quartic_irreducible' not in witnesses:witnesses['quartic_irreducible']=p
            if fp.is_irreducible() and 'cubic_irreducible' not in witnesses:witnesses['cubic_irreducible']=p
            if len(witnesses)==2:break
        assert len(witnesses)==2
        rows.append({'index':i,'generic_mask':rec['generic_mask'],'alpha':rec['alpha'],
            'quartic_ascending':list(map(str,q.list())),
            'tracefree_resolvent_ascending':list(map(str,tracefree.list())),
            'tracefree_generator_ascending':[str(a-T/3),str(b),str(c)],
            'norm_alpha':str(N),'norm_square_root':str(n),
            'change_of_generator_determinant':str(change),'prime_witnesses':witnesses,
            'field_discriminant':data['field_discriminant'],
            'field_discriminant_evidence':'Deduced from the retained strict unramified class certificate and the norm-square S4 correspondence; no new maximal-order computation.',
            'Galois_group':'S4','fixed_plus_resolvent_presentation':'IMPOSSIBLE by twist CT',
            'fixed_minus_resolvent_presentation':'EXISTS because this is an original generic rational Kummer class; no new point coordinates computed.'})
    return {'schema':'rank-jump.fixed-resolvent-solubility.v1','status':'PASS',
        'universal':universal(),'quartics':rows,'twist_CT_rank':6,'all_nonzero_block_exclusions':exclusions,
        'unramified_S4_classes_excluded_by_fixed_plus_resolvent':63,
        'exceptional_classes_constructed':0,'bindings':bindings([Path(__file__),PROTOCOL,INPUT,Path(r.__file__)]),
        'boundary':'A verified false-negative mechanism for fixed-polynomial resolvent constructors. These are existing generic strict classes; no new specialization class or rank predictor is supplied.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','build','check']);args=p.parse_args()
    if args.mode=='export':export()
    else:
        result=compute()
        if args.mode=='build':r.write_new(OUTPUT,result)
        else:assert result==r.read(OUTPUT)
        print('PASS six explicit unramified S4 quartics; all 63 nonzero classes excluded by fixed oriented resolvent')
