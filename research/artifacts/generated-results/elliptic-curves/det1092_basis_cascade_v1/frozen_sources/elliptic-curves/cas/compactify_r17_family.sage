#!/usr/bin/env sage-python
"""Exact R17 base compactification using bounded auxiliary bad-prime analysis.

Strip small factors from gcd(disc(A),disc(B),Res(A,B)); extract its largest
exact power before bounded factorization. Feed only proved primes to PARI's
prime-local auxiliary minimization. No curve search or global minimality claim.
"""
import argparse
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys
from sage.all import PolynomialRing,QQ,ZZ,pari,Matrix,gcd,lcm,prime_range

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
import audit_r17_constant_scaling as scaling
from research_runtime.store import checkpoint
helpers=SourceFileLoader('compact_base_helpers',str(CAS/'reduce_r17_family_base.sage')).load_module()
INPUT=helpers.INPUT
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/r17-compactification-v1'


def sources():
    return {**helpers.sources(),str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__).resolve())}


def prepare(directory):
    if (directory/'protocol.json').exists():raise FileExistsError('protocol already frozen')
    checkpoint(directory/'protocol.json',{'schema':'elliptic-curves.r17-compactification.v1','sources':sources(),
        'families':['103b2','11952','074d9','07ca9','08234','08f72'],
        'prime_trial_bound':10000,'maximum_power_exponent':132,'composite_factor_bit_gate':128,'prime_root_bit_gate':512,
        'precision_bits':8192,'worker_wall_seconds':120,'worker_rss_bytes':1073741824,'pari_stack_bytes':256000000,
        'maximum_concurrent_workers':2,
        'gate':'Run 103b2 first. Run other five only after exact base identities prove at least a 25-percent reduction of maximum coefficient bits.',
        'method':'Exact gcd of primitive A/B discriminants and resultant, trial division below 10000, largest exact perfect power through 132, proved factorization within the bit gate; prime-local auxiliary minimization and Cremona-Stoll reduction propose a base map; exact weighted identities decide validity.',
        'scope':'Six possible family coordinate computations; no parameter census, point search, global minimality or rank increase. Failed arithmetic is censored.'})


def primitive(f):
    den=lcm(q.denominator() for q in f);g=f*den
    return g/gcd(ZZ(q) for q in g)


def run(directory,family):
    protocol=cert.read(directory/'protocol.json')
    if protocol['sources']!=sources() or family not in protocol['families']:raise ArithmeticError('source mismatch')
    output=directory/(family+'.json')
    if output.exists():raise FileExistsError('preserve prior attempt')
    row=next(r for r in cert.read(INPUT)['rows'] if r['family']==family)
    R=PolynomialRing(QQ,'t');A,B=(R(list(map(QQ,row[k]))) for k in ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
    P,Q=primitive(A),primitive(B)
    invariants=[ZZ(P.discriminant()),ZZ(Q.discriminant()),ZZ(P.resultant(Q))]
    if any(x==0 for x in invariants):raise ArithmeticError('degenerate auxiliary invariant')
    common=abs(gcd(invariants));remainder=common;small=[]
    for p in prime_range(protocol['prime_trial_bound']):
        e=remainder.valuation(p)
        if e:small.append([int(p),int(e)]);remainder//=p**e
    exponent,root=1,remainder
    if remainder!=1:
        for k in range(2,protocol['maximum_power_exponent']+1):
            candidate,exact=remainder.nth_root(k,truncate_mode=True)
            if exact:exponent,root=k,candidate
    diagnostic={'common_invariant_gcd':str(common),'small_prime_factors':small,'remaining_cofactor':str(remainder),
        'largest_power_exponent':exponent,'power_root':str(root),'power_root_bits':int(root.nbits())}
    checkpoint(output,{'status':'RUNNING_BAD_PRIME_ANALYSIS','family':family,'protocol_sha256':cert.hashed(directory/'protocol.json'),'diagnostic':diagnostic})
    if root==1:factors=[]
    elif root.nbits()<=protocol['prime_root_bit_gate'] and root.is_prime(proof=True):factors=[(root,1)]
    elif root.nbits()<=protocol['composite_factor_bit_gate']:factors=list(root.factor(proof=True))
    else:raise ArithmeticError('remaining composite root exceeds the declared factor gate')
    product=ZZ(1)
    for p,e in factors:
        if not p.is_prime(proof=True):raise ArithmeticError('unproved prime')
        product*=p**e
    if product!=root or root**exponent!=remainder:raise ArithmeticError('factor reconstruction failed')
    restored=root**exponent
    for p,e in small:restored*=ZZ(p)**e
    if restored!=common:raise ArithmeticError('common invariant factorization mismatch')
    primes=sorted(set([ZZ(p) for p,e in small]+[p for p,e in factors]))
    diagnostic['root_factorization']=[[str(p),int(e)] for p,e in factors]
    checkpoint(output,{'status':'RUNNING_LOCAL_MINIMIZATION','family':family,'diagnostic':diagnostic,'primes':list(map(str,primes)),
        'protocol_sha256':cert.hashed(directory/'protocol.json')})
    pari.allocatemem(protocol['pari_stack_bytes'],silent=True);pari.set_real_precision_bits(protocol['precision_bits'])
    reducer=pari('(P,L)->{my(m,n);my(Q=hyperellminimalmodel(P,&m,L));my(S=hyperellred(Q,&n));[Q,m,S,n]}')
    minimal,m,reduced,n=reducer(pari(P),pari(primes))
    matrices=[Matrix(QQ,2,2,[z[1][i,j] for i,j in ((0,0),(0,1),(1,0),(1,1))]) for z in (m,n)]
    M=matrices[0]*matrices[1];a,b,c,d=M.list();det=M.det()
    if not det:raise ArithmeticError('singular base map')
    AA,BB=helpers.homogeneous(A,8,a,b,c,d),helpers.homogeneous(B,12,a,b,c,d)
    u,content=scaling.scale_for([cert.F(str(x)) for x in AA],[cert.F(str(x)) for x in BB])
    new_A,new_B=AA/QQ(str(u))**4,BB/QQ(str(u))**6
    if helpers.homogeneous(new_A,8,d,-b,-c,a)*QQ(str(u))**4!=det**8*A:raise ArithmeticError('inverse A failed')
    if helpers.homogeneous(new_B,12,d,-b,-c,a)*QQ(str(u))**6!=det**12*B:raise ArithmeticError('inverse B failed')
    bits=scaling.bits([cert.F(str(x)) for x in new_A.list()+new_B.list()])
    result={'status':'PASS_EXACT_BASE_CHANGE','family':family,'protocol_sha256':cert.hashed(directory/'protocol.json'),
        'diagnostic':diagnostic,'primes':list(map(str,primes)),'before_bits':row['after_bits'],'after_bits':bits,
        'base_matrix_a_b_c_d':list(map(str,M.list())),'total_scale_from_literal_source':str(u*cert.F(row['scale_u'])),
        'A_coefficients_low_to_high':list(map(str,new_A.list())),'B_coefficients_low_to_high':list(map(str,new_B.list())),
        'auxiliary_minimal_model':str(minimal),'auxiliary_minimization_transform':str(m),
        'auxiliary_reduced_model':str(reduced),'auxiliary_reduction_transform':str(n),'constant_scaling':content,
        'source':row['source'],'source_sha256':row['source_sha256'],
        'claim_boundary':'Exact base and coefficient identities. Generic-section transport must precede any new search. No global minimality or rank increase claim.'}
    checkpoint(output,result);print('COMPACT FAMILY',family,row['after_bits'],'->',bits,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run']);p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--family',default='103b2');a=p.parse_args()
    prepare(a.directory) if a.stage=='prepare' else run(a.directory,a.family)
