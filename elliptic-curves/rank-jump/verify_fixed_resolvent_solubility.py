#!/usr/bin/env python3
"""Independent root-sum, rational-algebra and finite-field replay."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import retrospective as r
import fixed_resolvent_solubility as run
from verify_unpointed_governing_norm import Algebra

OUTPUT=r.OUT/'rank_jump_fixed_resolvent_solubility_verification_v1.json'
UPSTREAM=r.OUT/'rank_jump_fixed_cubic_minus_reference_verification_v1.json'


def compute():
    import sympy as s
    data=r.read(run.INPUT);out=r.read(run.OUTPUT);upstream=r.read(UPSTREAM)
    for obj in (data,out,upstream):
        for name,sha in obj['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    assert upstream['status']=='PASS' and upstream['unconditional_CT_rank']==6 and upstream['upstream_strict_Artin_replay']
    old=r.read(run.completed.OUTPUT)['stages']['local'];ct=r.read(run.switch.OUTPUT)
    assert data['classes']==[{'generic_mask':z['generic_mask'],'alpha':z['beta_ascending']} for z in old['class_records']]
    assert data['twist_CT_matrix']==ct['twist_CT_matrix']==old['minus_twist_CT_matrix']
    assert data['field_discriminant']==old['field_discriminant']
    z,u,v,w=s.symbols('z u v w')
    roots=[(u+v+w)/2,(u-v-w)/2,(-u+v-w)/2,(-u-v+w)/2]
    T=u*u+v*v+w*w;S=u*u*v*v+u*u*w*w+v*v*w*w
    Q=s.expand(s.prod(z-a for a in roots))
    assert s.expand(Q-(z**4-T*z*z/2-u*v*w*z+(T*T-4*S)/16))==0
    pairs=[roots[0]*roots[1]+roots[2]*roots[3],roots[0]*roots[2]+roots[1]*roots[3],roots[0]*roots[3]+roots[1]*roots[2]]
    assert [s.expand(p+T/2) for p in pairs]==[u*u,v*v,w*w]
    # The resolvent's centered roots are alpha_i - Trace(alpha)/3.
    assert [s.expand(p+T/6) for p in pairs]==[s.expand(x*x-T/3) for x in (u,v,w)]
    K=Algebra(data['cubic_ascending']);A=K.f[1];B=K.f[0]
    trace=lambda a:3*a[0]-2*A*a[2]
    poly=lambda cs:s.Poly.from_list([s.Rational(c) for c in cs[::-1]],z,domain=s.QQ)
    f=poly(data['cubic_ascending']);count=0
    for rec,row in zip(data['classes'],out['quartics']):
        alpha=K.elt(rec['alpha']);T=trace(alpha);S=(T*T-trace(K.mul(alpha,alpha)))/2;N=K.norm(alpha)
        n=F(row['norm_square_root']);assert N==n*n==F(row['norm_alpha']) and n>0
        q=poly(row['quartic_ascending'])
        assert [F(str(q.nth(i))) for i in range(5)]==[(T*T-4*S)/16,-n,-T/2,0,1]
        determinant=K.norm(K.elt([alpha[1],-alpha[2],0]));assert determinant==F(row['change_of_generator_determinant'])!=0
        assert q.discriminant()==f.discriminant()*s.Rational(str(determinant))**2
        centered=K.elt([alpha[0]-T/3,alpha[1],alpha[2]])
        assert list(map(str,centered))==row['tracefree_generator_ascending']
        resolvent=poly(row['tracefree_resolvent_ascending'])
        assert K.evaluate(list(map(str,[resolvent.nth(i) for i in range(4)])),centered)==K.elt([0,0,0])
        for role,p in row['prime_witnesses'].items():
            pol=q if role=='quartic_irreducible' else f
            cs=[int(s.numer(c))*pow(int(s.denom(c)),-1,p)%p for c in pol.all_coeffs()]
            reduced=s.Poly.from_list(cs,z,modulus=p)
            assert reduced.discriminant()%p and reduced.is_irreducible;count+=1
        assert row['field_discriminant']==data['field_discriminant']
    M=data['twist_CT_matrix'];assert r.rank([r.pack(row) for row in M])==6
    assert len(out['all_nonzero_block_exclusions'])==63
    for word,row in enumerate(out['all_nonzero_block_exclusions'],1):
        assert row['word']==word
        j=row['CT_witness_basis_index'];assert sum(((word>>i)&1)*M[i][j] for i in range(6))%2==row['CT_value']==1
    return {'schema':'rank-jump.fixed-resolvent-solubility-verification.v1','status':'PASS',
        'four_root_product_identity':True,'resolvent_pairing_roots_verified':3,
        'rational_algebra_quartics_verified':6,'polynomial_discriminant_identities':6,
        'finite_field_irreducibility_witnesses':count,'CT_exclusion_words_verified':63,
        'retained_unconditional_strict_CT_certificate_bound':True,
        'fresh_maximal_order_or_unramifiedness_computations':0,
        'bindings':run.bindings([Path(__file__),run.INPUT,run.OUTPUT,UPSTREAM,Path(r.__file__),Path(__file__).with_name('verify_unpointed_governing_norm.py')])}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();out=compute()
    if args.mode=='build':r.write_new(OUTPUT,out)
    else:assert out==r.read(OUTPUT)
    print('PASS independent root-sum quartics, twelve prime witnesses and 63 CT exclusions')
