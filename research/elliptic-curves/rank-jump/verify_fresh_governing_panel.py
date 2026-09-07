#!/usr/bin/env python3
"""Independent integer/modular and local-power replay of the masked production panel."""
import argparse
from fractions import Fraction
from pathlib import Path
import sys
import retrospective as r
import fresh_governing_panel as base
import fresh_governing_completion as completion
import fresh_governing_octics as octics
import verify_explicit_governing_octic as portable
import explicit_governing_octic as group

OUTPUT=r.OUT/'rank_jump_fresh_governing_panel_verification_v1.json'


def jacobi(a,n):
    assert n>0 and n%2
    a%=n;sign=1
    while a:
        while not a%2:
            a//=2
            if n%8 in (3,5):sign=-sign
        a,n=n,a
        if a%4==n%4==3:sign=-sign
        a%=n
    return sign if n==1 else 0


def verify():
    from sage.all import QQ,ZZ,pari,matrix,GF,AA,PolynomialRing
    sys.set_int_max_str_digits(20000)
    inputs=r.read(base.INPUT);assert set(inputs)=={'schema','cases'}
    records=[r.read(p) for p in (base.OUTPUT,completion.OUTPUT,octics.OUTPUT)]
    for data in records:
        for name,digest in data['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==digest
    completed={x['token']:x for x in records[1]['rows']};octic_rows={x['token']:x for x in records[2]['rows']}
    rows=[];modular_count=0;local_power_count=0;ideal_count=0;jacobi_count=0
    for old in records[0]['rows']:
        token=old['token'];f,pts,scale=base.model_data(token);model=['0','0','0',str(f[1]),str(f[0])]
        local=old['local'];patch=completed.get(token)
        used=patch if patch and 'local' in patch else local
        assert used.get('local')
        vectors=[]
        for place in used['local']:
            sigs=place['signatures']
            vectors.extend([[s[j] for s in sigs] for j in range(len(sigs[0]))])
        kernel=matrix(GF(2),vectors).right_kernel();k=int(kernel.dimension())
        masks=used.get('generic_strict_masks',used.get('tested_kernel_masks'))
        assert len(masks)==k
        for mask in masks:assert matrix(GF(2),[[mask>>i&1 for i in range(len(pts))]]).row(0) in kernel
        blocks=[(p,tuple(roots)) for p,roots in used['independence_blocks']]
        sigs=[r.point_signature(model,list(map(str,P)),blocks) for P in pts]
        assert sigs==used['independence_signatures'] and r.rank(sigs)==len(pts)
        # Reconstruct strict representatives from generic points and check local
        # squareness by PARI's local-power API, separately from its character logs.
        tested=used.get('S_finite',used.get('tested_bad_primes'))
        maximal_primes=[p for p,e in old['factor']['factors']] if old['factor']['status']=='PASS' else tested
        nf=pari.nfinit([pari(f),maximal_primes]);theta=pari.Mod('z',pari(f))
        gammas=[]
        for x,y in pts:
            d=ZZ(x.denominator()).sqrt();gammas.append(pari(ZZ(x*d*d))-pari(d*d)*theta)
        betas=[]
        for j,mask in enumerate(masks):
            beta=pari.Mod(1,pari(f))
            for i,gamma in enumerate(gammas):
                if mask>>i&1:beta*=gamma
            betas.append(beta)
            for p in tested:
                for P in pari.idealprimedec(nf,p):
                    assert pari.nfislocalpower(nf,P,beta,2)==1;local_power_count+=1
            for root in f.roots(AA,multiplicities=False):
                assert f.parent()(list(QQ(pari.lift(beta).polcoef(i)) for i in range(3)))(root)>0
            if 'class_records' in local:
                rec=local['class_records'][j]
                assert [str(pari.lift(beta).polcoef(i)) for i in range(3)]==rec['beta_ascending']
                I=pari(matrix(QQ,rec['half_ideal_hnf']))
                assert pari.idealpow(nf,I,2)==pari.idealhnf(nf,beta);ideal_count+=1
        matrix_source=patch if patch and 'Artin_matrix' in patch else local
        M=matrix_source.get('minus_twist_CT_matrix');A=matrix_source.get('Artin_matrix')
        if M is not None:
            assert M==[[A[i][j]^A[j][i] for j in range(len(A))] for i in range(len(A))]
            rank=int(matrix(GF(2),M).rank()) if M else 0
            assert rank==matrix_source['minus_twist_CT_rank'] and rank%2==0
            for j,c in enumerate(local.get('artin_columns',[])):
                for i,e in enumerate(c['evaluations']):
                    if 'artin_bit' in e:
                        assert int(jacobi(int(e['residue']),int(c['norm']))==-1)==A[i][j];jacobi_count+=1
            if patch:
                for repair in patch.get('repairs',[]):
                    value=jacobi(int(repair['residue']),int(repair['cofactor']));assert value==repair['jacobi'];jacobi_count+=1
                    bit=int(value==-1)
                    for pl in repair['local']:
                        for c in pl['contributions']:bit^=(c['exponent']%2)*c['independent_generic_character_bit']
                    assert bit==A[repair['row']][repair['column']]
        # Governing octics: independent Sylvester determinant and finite-field
        # polynomial gcd/Tonelli--Shanks implementation in ordinary Python.
        o=octic_rows[token];assert o['status']=='PASS'
        H=list(map(int,o['integral_octic_ascending']));assert portable.discriminant(H)==int(o['integral_octic_discriminant'])
        x,y=map(Fraction,map(str,pts[0]));u,v=map(Fraction,map(str,pts[1]));c=u-x
        expect=[c**6,0,-4*c**3*(v-y),0,6*c*c*(x+u),0,-4*(y+v),0,1]
        assert list(map(Fraction,o['rational_octic_ascending']))==expect
        transform={'cubic_ascending':list(map(str,f.list())),'integral_octic_ascending':o['integral_octic_ascending'],'scaled_points':[list(map(str,P)) for P in pts[:2]]}
        for e in o['inert_prime_table']:
            portable.prime_replay(transform,{'prime':e['prime'],'octic_factor_degrees':e['factor_degrees'],'radical_norm_mod_p':e['radical_norm'],'psi':e['psi'],'independent_radical_psi':e['psi']});modular_count+=1
        rows.append({'token':token,'generic_independence_rank':len(pts),'verified_tested_kernel_dimension':k,
                     'full_kernel_certified':local['status'] in ('PASS','PARTIAL') or k==0,
                     'CT_switch_rank':matrix_source.get('minus_twist_CT_rank',patch.get('minus_twist_CT_rank') if patch else None),
                     'octic_and_all_inert_primes_verified':True})
        print('verified',token,k,flush=True)
    files=(Path(__file__),base.INPUT,base.OUTPUT,completion.OUTPUT,octics.OUTPUT,Path(base.__file__),Path(r.__file__),Path(portable.__file__),Path(group.__file__))
    return {'schema':'rank-jump.fresh-governing-panel-verification.v1','status':'PASS',
      'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files},'rows':rows,
      'independent_octic_modular_replays':modular_count,'independent_local_power_checks':local_power_count,
      'ideal_square_certificates':ideal_count,'independent_Jacobi_checks':jacobi_count,
      'universal_group_certificate':group.finite_certificate(),
      'boundary':'Certifies the inherited baseline and explicitly partial kernels. Does not fill additional Selmer or full CT data.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();out=verify()
    if a.mode=='build':r.write_new(OUTPUT,out)
    else:assert r.read(OUTPUT)==out
    print('PASS masked governing panel verification')
