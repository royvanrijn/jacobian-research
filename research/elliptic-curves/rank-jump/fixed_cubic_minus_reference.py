#!/usr/bin/env python3
"""A fixed-incidence, simultaneous-solubility control from the masked reference."""
import argparse
from pathlib import Path
from math import prod
import subprocess
import sys
import retrospective as r
import bounded_gain_reference as ref
import bounded_gain_reference_completion as completed

PROTOCOL=Path(__file__).with_name('FIXED_CUBIC_MINUS_REFERENCE_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_fixed_cubic_minus_reference_v1.json'
REPORT=r.OUT/'rank_jump_fixed_cubic_minus_reference_comparison_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-fixed-cubic-minus-reference-v1'
THEOREM=Path(__file__).with_name('INDEPENDENT_SCALAR_CUP_AND_TWIST_BLOCKS.md')


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
        (Path(__file__),PROTOCOL,ref.INPUT,completed.OUTPUT,Path(ref.__file__),THEOREM)}


def compute():
    from sage.all import QQ,ZZ,pari,EllipticCurve,AA,GF,matrix,PolynomialRing
    pari.allocatemem(64000000,268435456,silent=True);ref.configure()
    f,pts,scale=ref.base.model_data(ref.TOKEN);R=f.parent();z=R.gen();g=-f(-z)
    assert g==z**3+f[1]*z-f[0] and g.discriminant()==f.discriminant()
    prior=r.read(completed.OUTPUT)['stages'];fac=prior['factor'];local=prior['local']
    primes=[p for p,e in fac['factors']]
    assert prod(p**e for p,e in fac['factors'])==abs(16*f.discriminant())
    assert all(ZZ(p).is_prime(proof=True) for p in primes)
    arms=[]
    for name,h in [('original',f),('minus_one',g)]:
        E=EllipticCurve([0,0,0,h[1],h[0]]);nf=pari.nfinit([pari(h),primes]);th=pari.Mod('z',pari(h))
        finite=[];S=[2];ell=0
        for p in primes:
            local_data=E.local_data(p,proof=True);cv=int(local_data.conductor_valuation())
            if p!=2 and cv:S.append(p)
            dec=pari.idealprimedec(nf,p)
            point_dim=len(dec)-1+int(p==2)
            if p==2 or cv:ell+=point_dim
            finite.append({'prime':p,'conductor_valuation':cv,'local_point_dimension':point_dim,
                           'splitting':[[int(P[2]),int(P[3])] for P in dec]})
        roots=h.roots(AA,multiplicities=False);assert len(roots)==3;ell+=1
        delta=h.discriminant();beta=-pari(delta)*pari(h.derivative())(th)
        assert pari.nfeltnorm(nf,beta)==delta**4
        signs=[int(-delta*h.derivative()(x)<0) for x in roots];assert signs==[1,0,1]
        omitted=[]
        for p in sorted(set(primes)-set(S)):
            vals=[int(pari.idealval(nf,beta,P)) for P in pari.idealprimedec(nf,p)]
            assert all(v%2==0 for v in vals);omitted.append({'prime':p,'valuations':vals})
        arms.append({'name':name,'model':['0','0','0',str(h[1]),str(h[0])],
            'cubic_ascending':list(map(str,h.list())),'field_discriminant':str(nf.disc()),
            'S_finite':S,'local':finite,'local_point_product_dimension':ell,
            'derivative_coefficients':[str(pari.lift(beta).polcoef(i)) for i in range(3)],
            'derivative_norm':str(delta**4),'real_derivative_signs':signs,
            'omitted_good_prime_valuations':omitted,'Selmer_boundary_upper_bound':ell-1})
    assert arms[0]['S_finite']==arms[1]['S_finite']==local['S_finite']
    assert arms[0]['field_discriminant']==arms[1]['field_discriminant']==local['field_discriminant']
    A=local['Artin_matrix'];M=[[A[i][j]^A[j][i] for j in range(len(A))] for i in range(len(A))]
    assert M==local['minus_twist_CT_matrix'];rank=int(matrix(GF(2),M).rank());k=len(M)
    assert rank==k==local['strict_generic_dimension']==6
    # Check the transported representatives directly in the twist cubic.
    nf=pari.nfinit([pari(g),primes]);th=pari.Mod('z',pari(g));power_checks=0;norm_checks=0
    for rec in local['class_records']:
        coefficients=[QQ(c)*(-1)**i for i,c in enumerate(rec['beta_ascending'])]
        beta=pari(R(coefficients))(th);N=ZZ(pari.nfeltnorm(nf,beta));assert N>0 and N.is_square();norm_checks+=1
        for p in arms[1]['S_finite']:
            for P in pari.idealprimedec(nf,p):
                assert pari.nfislocalpower(nf,P,beta,2)==1;power_checks+=1
        for root in g.roots(AA,multiplicities=False):assert R(coefficients)(root)>0
    return {'schema':'rank-jump.fixed-cubic-minus-reference.v1','status':'PASS','bindings':bindings(),
        'arms':arms,'field_identification':'theta_original = -theta_minus',
        'same_strict_class_group':True,'transported_strict_dimension':k,
        'twist_CT_matrix':M,'twist_CT_restriction_rank':rank,
        'transported_local_power_checks':power_checks,'transported_norm_square_checks':norm_checks,
        'twist_Sha2_dimension_unconditional_lower_bound':rank,
        'transported_block_rational_dimension_on_original':k,
        'transported_block_rational_dimension_on_twist':0,
        'twist_strict_rational_dimension_upper_formula':f'c_S - {rank}',
        'twist_rank_upper_formula':f'c_S + {arms[1]["Selmer_boundary_upper_bound"]} - {rank}',
        'boundary':'Unconditional witnessed Sha[2] subspace and bounds in c_S. Original strict classes are generic and rational; their transported nondegenerate CT block intersects the twist rational Kummer image trivially. No total c_S or numerical whole-twist rank upper bound computed here.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);path=WORK/'worker.json'
    if not path.exists():
        error=None
        with (WORK/'worker.log').open('x') as log:
            try:
                p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker'],stdout=log,stderr=log,timeout=30)
                if p.returncode:error='Worker failure'
            except subprocess.TimeoutExpired:error='Bounded worker timeout'
        if error:r.write_new(path,{'status':'UNKNOWN','reason':error,'bindings':bindings()})
    result=r.read(path);assert result['bindings']==bindings();r.write_new(OUTPUT,result)
    print(result['status'],result.get('twist_Sha2_dimension_unconditional_lower_bound'),flush=True)


def report():
    d=r.read(OUTPUT);label=r.read(ref.REPORT)['result'];assert d['status']=='PASS' and d['same_strict_class_group']
    assert label['localized_class_dimension_conditional_interval']==[12,12]
    c=12;b=d['arms'][1]['Selmer_boundary_upper_bound'];rho=d['twist_CT_restriction_rank']
    result={'assumption':label['assumption'],'shared_strict_Selmer_dimension_conditional':c,
        'original_rank_exact_conditional':label['label_rank_exact_conditional'],
        'twist_Selmer_dimension_conditional_interval':[c,c+b],
        'twist_rank_conditional_upper_bound':c+b-rho,
        'original_strict_rational_dimension_exact_conditional':c,
        'twist_strict_rational_dimension_conditional_upper_bound':c-rho,
        'rank_drop_conditional_lower_bound':label['label_rank_exact_conditional']-(c+b-rho),
        'twist_Sha2_dimension_lower_bound_unconditional':rho,
        'rank_of_twist_exact':'UNKNOWN','new_point_searches':0,
        'family_parameter_transfer':'UNKNOWN; this is a quadratic twist, not a K3 base specialization'}
    files=(Path(__file__),PROTOCOL,OUTPUT,ref.REPORT,ref.LABEL)
    r.write_new(REPORT,{'schema':'rank-jump.fixed-cubic-minus-reference-comparison.v1','result':result,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files},
        'boundary':'Conditional class/rank label joined after masked arithmetic. The six-dimensional Sha[2] lower bound is unconditional; the twist rank upper bound and full strict-incidence comparison use the stated GRH assumption.'})
    print(result,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['worker','capture','report']);args=p.parse_args()
    if args.mode=='worker':r.write_new(WORK/'worker.json',compute())
    else:globals()[args.mode]()
