#!/usr/bin/env sage-python
"""Use bounded prime-local auxiliary minimization to propose rational base maps.

PARI's explicit prime list avoids global discriminant factorization. The
auxiliary hyperelliptic model only proposes a map; both elliptic coefficient
identities and their eventual section transport are separate exact gates.
"""
import argparse
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys
from sage.all import PolynomialRing, QQ, ZZ, pari, Matrix

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
import audit_r17_constant_scaling as scaling
import minimize_r17_small_prime_base as local
from research_runtime.store import checkpoint
previous=SourceFileLoader('base_reduction_helpers',str(CAS/'reduce_r17_family_base.sage')).load_module()
INPUT=previous.INPUT
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/r17-prime-local-base-v1'


def sources():
    return {**previous.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in
            (Path(__file__).resolve(),Path(local.__file__).resolve())}}


def prepare(directory):
    if (directory/'protocol.json').exists():raise FileExistsError('protocol already frozen')
    checkpoint(directory/'protocol.json',{'schema':'elliptic-curves.r17-prime-local-base.v1','sources':sources(),
        'families':['103b2','11952','074d9','07ca9','08234','08f72'],'precision_bits':8192,
        'method':'hyperellminimalmodel of primitive A at explicitly listed primes <=997 with common A/B multiplicities >=4/6, then Cremona-Stoll hyperellred. Compose base maps, transform both coefficients, remove exact weighted constant content.',
        'gate':'Run 103b2 first; proceed to the other five only if largest coefficient size falls at least 25 percent.',
        'worker_wall_seconds':120,'worker_rss_bytes':1073741824,'pari_stack_bytes':256000000,
        'maximum_concurrent_workers':2,
        'scope':'At most six bounded auxiliary minimizations. No global factorization, parameter sweep, point search, minimality of elliptic surface, or new-rank claim.'})


def run(directory,family):
    protocol=cert.read(directory/'protocol.json')
    if protocol['sources']!=sources() or family not in protocol['families']:raise ArithmeticError('source mismatch')
    output=directory/(family+'.json')
    if output.exists():raise FileExistsError('preserve prior attempt')
    row=next(r for r in cert.read(INPUT)['rows'] if r['family']==family)
    R=PolynomialRing(QQ,'t');t=R.gen()
    A,B=(R(list(map(QQ,row[k]))) for k in ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
    aa,bb=[[int(QQ(x)) for x in row[k]] for k in ('A_coefficients_low_to_high','B_coefficients_low_to_high')]
    if any(QQ(x) not in ZZ for x in row['A_coefficients_low_to_high']+row['B_coefficients_low_to_high']):raise ArithmeticError('nonintegral input')
    primes=[p for p in range(2,998) if local._is_prime(p) and local.roots(aa,bb,p)]
    primitive=A/QQ(scaling.content(list(map(cert.F,row['A_coefficients_low_to_high']))))
    if primitive.gcd(primitive.derivative()).degree()!=0:raise ArithmeticError('auxiliary repeated root')
    checkpoint(output,{'status':'RUNNING','family':family,'primes':primes,'protocol_sha256':cert.hashed(directory/'protocol.json')})
    pari.allocatemem(protocol['pari_stack_bytes'],silent=True);pari.set_real_precision_bits(protocol['precision_bits'])
    reducer=pari('(P,L)->{my(m,n);my(Q=hyperellminimalmodel(P,&m,L));my(S=hyperellred(Q,&n));[Q,m,S,n]}')
    minimized,m,reduced,n=reducer(pari(primitive),pari(primes))
    matrices=[Matrix(QQ,2,2,[z[1][i,j] for i,j in ((0,0),(0,1),(1,0),(1,1))]) for z in (m,n)]
    M=matrices[0]*matrices[1];a,b,c,d=M.list();det=M.det()
    if not det:raise ArithmeticError('singular base map')
    AA,BB=previous.homogeneous(A,8,a,b,c,d),previous.homogeneous(B,12,a,b,c,d)
    u,diag=scaling.scale_for([cert.F(str(x)) for x in AA],[cert.F(str(x)) for x in BB])
    new_A,new_B=AA/QQ(str(u))**4,BB/QQ(str(u))**6
    if previous.homogeneous(new_A,8,d,-b,-c,a)*QQ(str(u))**4 != det**8*A:raise ArithmeticError('inverse A identity failed')
    if previous.homogeneous(new_B,12,d,-b,-c,a)*QQ(str(u))**6 != det**12*B:raise ArithmeticError('inverse B identity failed')
    bits=scaling.bits([cert.F(str(x)) for x in new_A.list()+new_B.list()])
    result={'status':'PASS_EXACT_BASE_CHANGE','family':family,'primes':primes,'protocol_sha256':cert.hashed(directory/'protocol.json'),
        'before_bits':row['after_bits'],'after_bits':bits,'base_matrix_a_b_c_d':list(map(str,M.list())),
        'total_scale_from_literal_source':str(u*cert.F(row['scale_u'])),'additional_scale_u':str(u),
        'A_coefficients_low_to_high':list(map(str,new_A.list())),'B_coefficients_low_to_high':list(map(str,new_B.list())),
        'auxiliary_minimal_model':str(minimized),'auxiliary_minimization_transform':str(m),
        'auxiliary_reduced_model':str(reduced),'auxiliary_reduction_transform':str(n),
        'constant_scaling':diag,'source':row['source'],'source_sha256':row['source_sha256'],
        'claim_boundary':'Exact base and coefficient identities only. Section transport must be checked before a prospective search. No globally minimal family assertion.'}
    checkpoint(output,result);print('PRIME LOCAL BASE',family,row['after_bits'],'->',bits,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run']);p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--family',default='103b2');a=p.parse_args()
    prepare(a.directory) if a.stage=='prepare' else run(a.directory,a.family)
