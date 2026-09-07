#!/usr/bin/env python3
"""Independently verify the relation-root obstruction by lattice membership."""
import argparse
from pathlib import Path
from fractions import Fraction as Q
import retrospective as r
import relation_root_class as original
import complete_relation_root_class as run
import early_relation_pool as pool
import retained_norm_inherited_hit as hit
import bounded_gain_reference as ref
from verify_unpointed_governing_norm import Algebra

OUTPUT=r.OUT/'rank_jump_relation_root_class_verification_v1.json'


def compute():
    from sage.all import QQ,ZZ,PolynomialRing,pari,matrix,vector
    d=r.read(original.INPUT);out=r.read(run.OUTPUT)
    for obj in (d,out):
        for name,sha in obj['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    source=r.read(pool.INPUT);h=r.read(hit.OUTPUT)
    rel=next(x for x in source['relations'] if [x['m'],x['n']]==h['address'])
    assert d['probe_primes']==sorted({source['columns'][i]['p'] for i,e in rel['ideal_factorization']}-set(d['S_finite']))[:32]
    assert d['root_ascending']==h['square_root_ascending'] and d['generic_product_mask']==h['generic_product_mask']
    ref.configure();fm,pts,scale=ref.base.model_data(ref.TOKEN)
    b,a=map(QQ,source['affine_map_old_root_from_masked'])
    assert d['generic_classes_ascending']==[[str(a*x+b),'-1','0'] for x,y in pts]
    A=Algebra(d['cubic_ascending']);w=A.elt(d['root_ascending']);alpha=A.elt(d['alpha_ascending'])
    gamma=list(map(A.elt,d['generic_classes_ascending']))
    parent=A.mul(A.elt([A.norm(alpha)]),alpha)
    for i,g in enumerate(gamma):
        if d['generic_product_mask']>>i&1:parent=A.mul(parent,g)
    assert A.mul(w,w)==parent
    N=A.norm(w);assert str(N)==out['norm_root']
    beta=A.mul(A.elt([N]),w)
    assert beta==A.elt(out['projection_ascending']) and A.norm(beta)==N**4
    R=PolynomialRing(QQ,'z');f=R(d['cubic_ascending']);nf=pari.nfinit([pari(f),d['S_finite']])
    basis=list(map(A.elt,out['maximal_order_basis']));B=matrix(QQ,basis).transpose()
    assert basis==[A.elt([str(pari.lift(z).polcoef(i)) for i in range(3)]) for z in nf.nf_get_zk()]
    assert B.det()**2*f.discriminant()==QQ(d['field_discriminant'])
    # Multiply ideal lattices independently in the certified integral basis.
    def multiply_ideals(H,J):
        columns=[]
        for i in range(3):
            for j in range(3):
                left=A.elt(map(str,B*H.column(i)));right=A.elt(map(str,B*J.column(j)))
                columns.append(B.inverse()*vector(QQ,A.mul(left,right)))
        # Sage row HNF provides a canonical basis for the column span.
        return matrix(ZZ,columns).hermite_form(include_zero_rows=False).transpose()
    def same_lattice(H,J):
        return all(x in ZZ for x in (H.inverse()*J).list()) and all(x in ZZ for x in (J.inverse()*H).list())
    checks=0;witnesses=[];rows=[]
    for p in d['probe_primes']:
        assert ZZ(p).is_prime(proof=True) and p not in d['S_finite'] and f.discriminant()%p
        recorded=[x for x in out['rows'] if x['p']==p];dec=list(pari.idealprimedec(nf,p))
        assert len(recorded)==len(dec)
        for row,P in zip(recorded,dec):
            assert row['prime_hnf']==str(pari.idealhnf(nf,P)) and row['e']==int(P[2])==1 and row['f']==int(P[3])
            H=matrix(ZZ,pari.idealhnf(nf,P));assert abs(H.det())==p**row['f']
            powers=[matrix(ZZ,3,3,1)]
            def power(k):
                while len(powers)<=k:powers.append(multiply_ideals(powers[-1],H))
                return powers[k]
            for label,z in [('root',w),('projection',beta)]:
                cert=row[label+'_membership'];D=ZZ(cert['denominator']);v=row[label+'_valuation'];k=cert['integral_valuation']
                coords=B.inverse()*vector(QQ,z)*D
                assert list(map(str,coords))==cert['integral_coordinates'] and all(c in ZZ for c in coords)
                assert k==v+D.valuation(p)
                Hk=matrix(QQ,cert['ideal_power_hnf']);Hnext=matrix(QQ,cert['next_ideal_power_hnf'])
                assert same_lattice(Hk,power(k)) and same_lattice(Hnext,power(k+1))
                assert abs(Hk.det())==p**(row['f']*k) and abs(Hnext.det())==p**(row['f']*(k+1))
                assert all(c in ZZ for c in Hk.inverse()*coords)
                assert any(c not in ZZ for c in Hnext.inverse()*coords)
                checks+=1
            # The retained witness primes are units on every generic class.
            assert row['generic_valuations']==[0]*16
            for g in gamma:
                coords=B.inverse()*vector(QQ,g)
                assert all(QQ(x).denominator()%p for x in coords)
                assert any(QQ(x).denominator()%p==0 for x in H.inverse()*coords)
            vp=QQ(str(N)).valuation(p)
            assert row['projection_valuation']==row['root_valuation']+vp
            defect=(2*row['root_valuation']+2*vp)%4
            if defect:witnesses.append(out['rows'].index(row))
            rows.append({'p':p,'residue_degree':row['f'],'projection_valuation':row['projection_valuation'],
                'parent_norm_projection_valuation_mod4':int(defect)})
    assert witnesses==out['odd_valuation_witness_rows'] and witnesses
    assert out['candidate_status']=='EXCLUDED_AT_GOOD_PRIME' and out['generic_corrections_excluded']==65536
    return {'schema':'rank-jump.relation-root-class-verification.v1','status':'PASS','rows':rows,
        'independent_valuation_checks':checks,'generic_unit_checks':16*len(rows),
        'odd_good_prime_ideal_witnesses':len(witnesses),'rational_witness_primes':len({x['p'] for x in rows if x['parent_norm_projection_valuation_mod4']}),
        'method':'Pure rational root and norm identities; independently multiplied integral ideal lattices and membership/nonmembership. PARI is used only for the maximal order and prime ideal identification, not valuations or ideal powers.',
        'bindings':original.bindings([Path(__file__),original.INPUT,run.OUTPUT,pool.INPUT,hit.OUTPUT,ref.INPUT,
            Path(__file__).with_name('verify_unpointed_governing_norm.py')])}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS',result['independent_valuation_checks'],'lattice valuations;',result['odd_good_prime_ideal_witnesses'],'odd witnesses')
