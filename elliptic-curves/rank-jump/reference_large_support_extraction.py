#!/usr/bin/env python3
"""Positive class extraction from exact, potentially incomplete principal relations."""
import argparse
from pathlib import Path
from math import prod, gcd
import re
import subprocess
import sys
import retrospective as r
import seeded_reference_class as original
import early_relation_pool as early

PROTOCOL=Path(__file__).with_name('REFERENCE_LARGE_SUPPORT_EXTRACTION_PROTOCOL.json')
RELWORK=r.ROOT/'artifacts/local/rank-jump-seeded-reference-relations-memory-v1'
WORK=r.ROOT/'artifacts/local/rank-jump-reference-large-support-extraction-v1'
OUTPUT=r.OUT/'rank_jump_reference_large_support_extraction_v1.json'
ACCEPTED=r.ROOT/'artifacts/local/rank-jump-reference-small-representative-class-v1/monitor_011_elements.json'

def indices(mask):
    while mask:
        bit=mask&-mask;yield bit.bit_length()-1;mask^=bit

def compute():
    from sage.all import AA,QQ,ZZ,PolynomialRing,pari
    sys.path.insert(0,str(r.ROOT/'elliptic-curves/cas'))
    from research_runtime.local_kummer import LocalSquareclasses
    pari.allocatemem(64000000,268435456,silent=True)
    ref=r.read(original.REFERENCE);pool=r.read(original.POOL);accepted=r.read(ACCEPTED)
    R=PolynomialRing(QQ,'z');f=R(ref['cubic_ascending']);nf=pari.nfinit([pari(f),ref['S_finite']])
    assert str(nf.disc())==ref['field_discriminant']
    assert list(map(str,nf.nf_get_zk()))==accepted['basis_GP']
    S=set(ref['S_finite']);disc=abs(ZZ(f.discriminant()))
    for p in S:
        while disc%p==0:disc//=p
    assert disc==1
    columns=list(pool['columns']);lookup={(c['p'],c['hnf']):i for i,c in enumerate(columns)}
    blocks={};dec={}
    for i,c in enumerate(columns):blocks.setdefault(c['p'],[]).append(i)
    atoms=[];keys={}
    def primitive_key(a):
        v=pari.nfalgtobasis(nf,a);v*=pari.denominator(v);v/=pari.content(v)
        w=tuple(map(int,v));sign=1 if next(c for c in w if c)>0 else -1
        return tuple(sign*c for c in w)
    for source in pool['relations']:
        a=pari.Mod(pari(R(source['alpha_ascending'])),pari(f));key=primitive_key(a)
        assert key not in keys
        N=ZZ(pari.nfeltnorm(nf,a));fac=source['ideal_factorization']
        assert prod(columns[i]['p']**(columns[i]['f']*e) for i,e in fac)==abs(N)
        keys[key]=len(atoms);atoms.append({'alpha_ascending':source['alpha_ascending'],
            'norm':str(N),'ideal_factorization':fac,'source':'eligible_early_pool'})
    def get_block(p):
        if p not in dec:dec[p]=list(pari.idealprimedec(nf,p))
        for P in dec[p]:
            h=str(pari.idealhnf(nf,P));key=(p,h)
            if key not in lookup:
                i=len(columns);lookup[key]=i;blocks.setdefault(p,[]).append(i)
                columns.append({'p':p,'hnf':h,'e':int(P[2]),'f':int(P[3])})
        return dec[p]
    added=0
    for encoded in accepted['elements_GP']:
        a=pari.nfbasistoalg(nf,pari(encoded));key=primitive_key(a)
        if key in keys:continue
        factor=pari.idealfactor(nf,a);fac=[];N=ZZ(pari.nfeltnorm(nf,a))
        for j in range(factor.nrows()):
            P,e=factor[j,0],int(factor[j,1]);p=int(P[0]);get_block(p)
            fac.append([lookup[(p,str(pari.idealhnf(nf,P)))],e])
        assert all(e>=0 for i,e in fac)
        assert prod(columns[i]['p']**(columns[i]['f']*e) for i,e in fac)==abs(N)
        assert pari.idealhnf(nf,pari.idealfactorback(nf,factor))==pari.idealhnf(nf,a)
        keys[key]=len(atoms);atoms.append({'alpha_ascending':[str(pari.lift(a).polcoef(i)) for i in range(3)],
            'norm':str(N),'ideal_factorization':fac,'source':'seeded_relation_worker'});added+=1
    masks=[]
    for atom in atoms:
        vals=dict(atom['ideal_factorization']);mask=0
        for p in {columns[i]['p'] for i in vals}-S:
            # Complete prime blocks are essential: norm projection adds v_p(N)
            # to every prime above p, including primes absent from (alpha).
            if p not in blocks:raise ArithmeticError('missing complete prime block')
            assert sum(columns[i]['e']*columns[i]['f'] for i in blocks[p])==3
            bit=sum(columns[i]['f']*vals.get(i,0) for i in blocks[p])%2
            for i in blocks[p]:mask^=((vals.get(i,0)+bit*columns[i]['e'])%2)<<i
        masks.append(mask)
    outside_rank,kernel=early.elimination(masks)
    used=sorted({i for m in kernel for i in indices(m)})
    local=[LocalSquareclasses(nf,p) for p in sorted(S)];roots=f.roots(AA,multiplicities=False)
    assert len(roots)==3
    polys=[];labels=[]
    for j,g in enumerate(ref['generic_classes']):
        beta=R(g['beta_ascending']);N=ZZ(pari.nfeltnorm(nf,pari.Mod(pari(beta),pari(f))))
        assert str(N)==g['norm'] and N.is_square() and beta.degree()==1
        d=-beta[1];assert d>0 and d.is_square() and d.denominator()==1
        polys.append(beta);labels.append({'kind':'generic','index':j})
    for i in used:
        polys.append(ZZ(atoms[i]['norm'])*R(atoms[i]['alpha_ascending']))
        labels.append({'kind':'projected_atom','index':i})
    locbits=[]
    for beta in polys:
        a=pari.Mod(pari(beta),pari(f));bits=[int(b) for L in local for b in L.signature(a)]
        signs=[int(beta(x)<0) for x in roots];assert sum(signs)%2==0
        locbits.append(r.pack(bits+signs))
    finite=[0]*len(polys);proof_primes=[];count=0
    for p0 in __import__('sage.all',fromlist=['prime_range']).prime_range(50001,55001):
        p=int(p0)
        if p in S or f.discriminant()%p==0:continue
        if any(c.denominator()%p==0 for b in polys for c in b):continue
        ff=f.change_ring(__import__('sage.all',fromlist=['GF']).GF(p))
        rr=ff.roots(multiplicities=False)
        if len(rr)!=3:continue
        vals=[[int(b.change_ring(ff.base_ring())(x)) for x in rr] for b in polys]
        if any(v==0 for row in vals for v in row):continue
        for j,row in enumerate(vals):
            for k,v in enumerate(row):finite[j]|=int(pow(v,(p-1)//2,p)==p-1)<<(count+k)
        count+=3;proof_primes.append(p)
    assert r.rank(finite[:16])==16
    amap={i:16+j for j,i in enumerate(used)}
    candidates=[1<<j for j in range(16)]
    candidates += [sum(1<<amap[i] for i in indices(m)) for m in kernel]
    def xor_values(mask,values):
        v=0
        for j in indices(mask):v^=values[j]
        return v
    local_rows=[xor_values(m,locbits) for m in candidates]
    local_rank,strict_coeff=early.elimination(local_rows)
    strict=[]
    for coef in strict_coeff:
        m=xor_values(coef,candidates);assert xor_values(m,locbits)==0
        strict.append(m)
    chars=[xor_values(m,finite) for m in strict]
    generic_basis=r.basis(finite[:16]);extra=[j for j,c in enumerate(chars) if r.reduce(c,generic_basis)]
    generic_local_rank,generic_strict=early.elimination(locbits[:16]);assert len(generic_strict)==6
    # These six controls are strict rational classes, constructed from G only.
    control_chars=[xor_values(m,finite) for m in generic_strict];assert r.rank(control_chars)==6
    strict_records=[{'factor_labels':[labels[j] for j in indices(m)],'proof_character':str(c),
                     'certified_outside_generic_span':j in extra} for j,(m,c) in enumerate(zip(strict,chars))]
    return {'schema':'rank-jump.reference-large-support-extraction.v1','status':'PASS',
        'bindings':original.bindings([Path(__file__),PROTOCOL,original.REFERENCE,original.POOL,ACCEPTED,
            RELWORK/'reference.log',RELWORK/'launch.json',Path(early.__file__),
            r.ROOT/'elliptic-curves/cas/research_runtime/local_kummer.py']),
        'equation':{'cubic_ascending':list(map(str,f.list())),'S_finite':sorted(S),'discriminant':str(nf.disc())},
        'accepted_worker_occurrences':len(accepted['elements_GP']),'new_distinct_principal_elements':added,
        'total_distinct_principal_elements':len(atoms),'outside_S_parity_rank':outside_rank,
        'outside_S_kernel_supports':[list(indices(m)) for m in kernel],
        'used_atoms':[{'index':i,**atoms[i]} for i in used],
        'new_atoms':[{'index':i,**a} for i,a in enumerate(atoms) if a['source']=='seeded_relation_worker'],
        'columns':columns if added else [],'local_S_order':sorted(S),
        'factor_labels':labels,'factor_local_signatures':[str(v) for v in locbits],
        'proof_primes':proof_primes,'factor_proof_characters':[str(v) for v in finite],
        'generic_character_rank':16,'generic_strict_dimension':6,
        'generic_strict_masks':[str(m) for m in generic_strict],
        'strict_coefficient_kernel_dimension':len(strict),'strict_character_rank':r.rank(chars),
        'strict_candidates':strict_records,'additional_independence_witnesses':extra,
        'positive_endpoint':'CANDIDATE_PENDING_INDEPENDENT_VERIFICATION' if extra else 'NOT_REACHED',
        'rational_or_Sha':'UNKNOWN for any additional class',
        'boundary':'Coefficient dependencies and fixed characters do not prove global class-group completeness. No new class is certified unless independence and all strict conditions replay independently.'}

def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    with (WORK/'worker.log').open('x') as log:
        try:
            p=subprocess.run([sys.executable,__file__,'worker'],stdout=log,stderr=log,
                timeout=r.read(PROTOCOL)['bounds']['worker_seconds'])
            error='worker failure' if p.returncode else None
        except subprocess.TimeoutExpired:error='bounded timeout'
    if error and not OUTPUT.exists():r.write_new(OUTPUT,{'status':'UNKNOWN','reason':error})
    result=r.read(OUTPUT);print({k:v for k,v in result.items() if k in ['status','positive_endpoint',
        'new_distinct_principal_elements','outside_S_parity_rank','generic_strict_dimension','additional_independence_witnesses']},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);a=p.parse_args()
    if a.mode=='worker':r.write_new(OUTPUT,compute())
    else:capture()
