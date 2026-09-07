#!/usr/bin/env python3
"""Hilbert-dual replay and certified small quadratic S-class-group controls."""
import argparse
from math import prod
from pathlib import Path
import retrospective as r
import resolvent_genus_component as source
import standard_s3_class_module as module

OUTPUT=r.OUT/'rank_jump_resolvent_genus_component_verification_v1.json'
COMPARISON=r.OUT/'rank_jump_resolvent_genus_component_comparison_v1.json'


def compute():
    from sage.all import ZZ,QQ,GF,matrix,pari,PolynomialRing
    from sage.env import SAGE_VERSION
    data=r.read(source.INPUT);out=r.read(source.OUTPUT);finite=r.read(module.OUTPUT)
    for obj in (data,out,finite):
        for p,sha in obj['bindings'].items():assert r.digest((r.ROOT/p).read_bytes())==sha,p
    assert finite==module.compute()
    def rank(words,n):return int(matrix(GF(2),[[v>>i&1 for i in range(n)] for v in words]).rank()) if words else 0
    def hilbert(a,b,p):
        return -1 if a<0 and b<0 else 1 if p=='infinity' else int(pari.hilbert(a,b,p))
    # Keep the archimedean case separate from the finite sign of a,b.
    def symbol(a,b,p):return (-1 if a<0 and b<0 else 1) if p=='infinity' else int(pari.hilbert(a,b,p))
    primes_checked=set();local_checks=0
    def verify_genus(row):
        nonlocal local_checks
        D=int(row['quadratic_discriminant']);atoms=list(map(int,row['prime_discriminants']));n=len(atoms)
        assert prod(atoms)==D and D%4 in (0,1)
        odd=[a for a in atoms if a not in (-4,8,-8)]
        assert all(abs(a)%2==1 and a%4==1 for a in odd)
        assert len([a for a in atoms if a in (-4,8,-8)])<=1
        for a in odd:
            assert abs(a) not in [1]
            if abs(a) not in primes_checked:assert ZZ(abs(a)).is_prime(proof=True);primes_checked.add(abs(a))
        assert len(set(map(abs,odd)))==len(odd)
        constraints=[]
        for loc in row['places']:
            p=loc['place']
            if p=='infinity':basis=[-1]
            elif p==2:basis=[2,-1,5]
            else:
                u=next(z for z in range(2,258) if pow(z,(p-1)//2,p)==p-1);basis=[p,u]
            width=len(basis);characters=[]
            # All Hilbert-dual functionals annihilating the class of D.
            for mask in range(1,1<<width):
                c=prod(v for i,v in enumerate(basis) if mask>>i&1)
                if symbol(D,c,p)!=1:continue
                characters.append(sum(int(symbol(a,c,p)==-1)<<i for i,a in enumerate(atoms)))
            assert rank(characters,n)==rank(loc['constraint_rows'],n)
            assert rank(characters+loc['constraint_rows'],n)==rank(characters,n)
            constraints.extend(characters);local_checks+=len(characters)*n
        k=n-rank(constraints,n);assert k-1==row['S_genus_dimension']
        full=row['full_kernel_basis_masks'];assert len(full)==k and rank(full,n)==k
        assert all((v&c).bit_count()%2==0 for v in full for c in constraints)
        total=(1<<n)-1;assert row['trivial_total_product_mask']==total
        selected=row['S_genus_basis_masks'];assert rank([total]+selected,n)==len(selected)+1==k
        assert all((v&c).bit_count()%2==0 for v in selected for c in constraints)
        assert row['S_genus_generators']==[str(prod(a for i,a in enumerate(atoms) if v>>i&1)) for v in selected]
        sign_constraint=[sum(int(a<0)<<i for i,a in enumerate(atoms))] if D>0 else []
        assert row['ordinary_genus_dimension']==n-1-rank(sign_constraint,n)
        assert row['finite_unramified_genus_dimension']==n-1
        return {'token':row['token'],'quadratic_discriminant':str(D),'ordinary_genus_dimension':row['ordinary_genus_dimension'],
            'S_genus_dimension':k-1,'standard_component_contribution':0}
    rows=[]
    for x,row in zip(data['rows'],out['rows']):
        assert x['token']==row['token']
        if x['status']=='UNKNOWN':assert x==row;rows.append(row);continue
        c=list(map(int,x['cubic_ascending']));delta=-4*c[1]**3-27*c[0]**2
        assert c[2:]==[0,1] and 16*delta==int(x['elliptic_discriminant'])
        assert prod(p**e for p,e in x['discriminant_factors'])==abs(16*delta)
        for p,e in x['discriminant_factors']:
            if p not in primes_checked:assert ZZ(p).is_prime(proof=True);primes_checked.add(p)
        p=x['cubic_irreducibility_prime'];assert ZZ(p).is_prime(proof=True)
        assert all((z**3+c[1]*z+c[0])%p for z in range(p))
        d=(-1 if delta<0 else 1)*prod(p for p,e in x['discriminant_factors'] if e%2)
        D=d if d%4==1 else 4*d;assert d!=1 and str(D)==row['quadratic_discriminant']
        assert [v['place'] for v in row['places']]==x['S_finite']+['infinity']
        rows.append({'status':'PASS',**verify_genus(row)})
    controls=[];P=PolynomialRing(QQ,'z');z=P.gen()
    for row in out['small_controls']:
        g=verify_genus(row);D=int(row['quadratic_discriminant']);d=D if D%4==1 else D//4
        bnf=pari.bnfinit(pari(z*z-d),1);assert pari.bnfcertify(bnf)==1 and int(bnf.disc())==D
        cyc=list(map(int,bnf.bnf_get_cyc()));even=[i for i,v in enumerate(cyc) if v%2==0];coords=[]
        for loc in row['places']:
            if loc['place']=='infinity':continue
            for ideal in pari.idealprimedec(bnf,loc['place']):
                v=list(map(int,pari.bnfisprincipal(bnf,ideal,0)));coords.append(sum((v[j]%2)<<i for i,j in enumerate(even)))
        direct=len(even)-rank(coords,len(even));assert direct==g['S_genus_dimension'] and len(even)==g['ordinary_genus_dimension']
        controls.append({**g,'ordinary_class_group_invariants':cyc,'S_prime_class_coordinates':coords,'bnfcertify':True,
            'direct_S_class_two_dimension':direct})
    return {'schema':'rank-jump.resolvent-genus-component-verification.v1','status':'PASS','rows':rows,'small_controls':controls,
        'local_Hilbert_bit_checks':local_checks,'certified_input_primes':len(primes_checked),'software':{'sage':SAGE_VERSION,'pari':str(pari('version()'))},
        'bindings':source.bindings([Path(__file__),source.INPUT,source.OUTPUT,module.OUTPUT,Path(module.__file__)]),
        'boundary':'Exact quadratic-resolvent S-class 2-dimensions by genus theory and independent local Hilbert constraints; four small certified class-group regressions. No research-field class group or additional elliptic class is constructed.'}


def compare(result):
    cap_path=r.OUT/'rank_jump_strict_cover_creation_capacity_v1.json';caps=r.read(cap_path)
    rows=[];bytoken={x['token']:x for x in caps['rows']}
    for x in result['rows']:
        if x['token'] not in bytoken:continue
        c=bytoken[x['token']];row={'token':x['token'],'id':c['id'],'family':c['family'],
            'retained_rank_lower_bound':c['retained_rank_lower_bound'],'generic_dimension':c['generic_dimension'],'status':x['status']}
        if x['status']=='PASS':
            g=x['S_genus_dimension'];k=c['generic_strict_dimension'];lower=c['total_strict_rational_dimension_lower_bound']
            extra=c['strict_rational_classes_outside_global_pool_lower_bound']
            row.update({'quadratic_resolvent_S_class_two_dimension':g,'cubic_S_class_two_dimension':'UNKNOWN',
                'sextic_S_class_module':'F2^'+str(g)+' + V^c_S',
                'sextic_S_class_two_dimension_formula':'2*c_S + '+str(g),
                'rank_label_forced_sextic_S_class_two_dimension_lower_bound':2*lower+g,
                'rank_label_forced_new_standard_copies_beyond_global_pool':extra,
                'rank_label_forced_new_standard_vector_dimensions':2*extra,
                'resolvent_genus_contribution_to_elliptic_strict_classes':0})
        else:row['reason']=x['reason']
        rows.append(row)
    return {'schema':'rank-jump.resolvent-genus-component-comparison.v1','rows':rows,
        'bindings':source.bindings([Path(__file__),OUTPUT,cap_path]),
        'boundary':'Exact genus component joined with previous rank-derived strict necessities. Standard multiplicities remain UNKNOWN; no rank predictor is produced.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result);r.write_new(COMPARISON,compare(result))
    else:assert result==r.read(OUTPUT);assert compare(result)==r.read(COMPARISON)
    print('PASS quadratic genus/local audit;',len(result['small_controls']),'certified small class groups;',result['local_Hilbert_bit_checks'],'Hilbert bits')
