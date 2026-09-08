#!/usr/bin/env sage-python
"""Bounded local-only PARI image construction, with exact proxy margins.

Reuse the existing local squareclass engine, but never its cubic-specific
point_kummer_dimension property. Independent finite-ring replay is required.
"""
import argparse,hashlib,json,signal,sys
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,pari,floor
from sage.version import version as sage_version
CAS=Path(__file__).resolve().parent;sys.path.insert(0,str(CAS))
from research_runtime.local_kummer import LocalSquareclasses
ROOT=CAS.parents[1];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_rr_dyadic_images_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def construct(i):
    sources=[ART/'det1092_rr_good_local_images_v1'/('case-%02d.json'%j) for j in range(10)]
    prior=ART/'det1092_rr_real_dyadic_panel_replay_v1.json'
    protocol={'classification':'local2-adic image construction; independent replay required',
              'inputs':{str(p.relative_to(ROOT)):sha(p) for p in sources+[prior,CAS/'research_runtime/local_kummer.py']},
              'limits':{'cases':10,'seconds_per_case':25,'proxy_precision':512,
                        'generic_divisors_per_case':17,'local_maximal_order_primes':[2],
                        'nfinit_calls_per_case':1,'global_class_groups':0,'global_unit_groups':0,
                        'full_discriminant_factorizations':0,'point_searches':0,'remote_arithmetic':0,'pilot_changes':0},
              'selector':'Use all17 inherited generic divisors, then determine independence after quotienting the three rational squareclasses -1,2,5. No marked point.',
              'script_sha256':sha(Path(__file__))}
    OUT.mkdir(parents=True,exist_ok=True);retain(OUT/'protocol.json',protocol)
    source=json.loads(sources[i].read_text());dimensions=json.loads(prior.read_text())['cases'][i]
    R=PolynomialRing(QQ,'x');x=R.gen();q=R(source['q']).monic();N=512;modulus=ZZ(2)**N
    lift=lambda a: ZZ(QQ(a).numerator()*QQ(a).denominator().inverse_mod(modulus)%modulus)
    s=min(floor(a.valuation(2)/(6-j)) for j,a in enumerate(q.list()[:-1]) if a)
    f=R(q(2**s*x)/2**(6*s));D=int(f.discriminant().valuation(2));assert N>2*D
    proxy=R([lift(a) for a in f]);assert proxy.discriminant().valuation(2)==D
    primes=json.loads((ART/'det1092_rr_good_local_images_v1/protocol.json').read_text())['primes']
    irreducibility_prime=next((p for p in primes if PolynomialRing(GF(p),'x')(proxy.list()).is_irreducible()),None)
    assert irreducibility_prime is not None,'NO_BOUNDED_PROXY_IRREDUCIBILITY_WITNESS'
    x0=QQ(source['base_x'])
    def normalized(P):
        v=min(a.valuation(2) for a in P if a);integral=R(P/2**v)
        small=R([lift(a) for a in integral]);norm=proxy.resultant(small)
        vn=int(norm.valuation(2));assert vn>=0 and N-D>vn+2
        return small,{'input_coefficient_min_v2':int(v),'normalized_residue_polynomial':list(map(str,small.list())),
                      'norm_v2':vn,'root_matching_lower_bound':N-D,
                      'relative_error_v2_strictly_greater_than_two':True}
    b,b_record=normalized(R([x0,-2**s]));a_records=[];as_poly=[]
    for j,row in enumerate(source['generic_divisors']):
        g=R(row);a,record=normalized(R(g(2**s*x)))
        record.update({'generic_divisor_index':j,'degree':int(g.degree())});a_records.append(record);as_poly.append(a)
    prefix={'classification':'verified proxy arithmetic margins; local images require independent replay',
            'case_index':i,'source':str(sources[i].relative_to(ROOT)),'source_sha256':sha(sources[i]),
            'T_equals_two_power_times_x':int(s),'proxy_polynomial':list(map(str,proxy.list())),
            'coefficient_precision':N,'integral_monic_discriminant_v2':D,
            'proxy_irreducibility_prime':irreducibility_prime,'anchor':b_record,'generic_evaluations':a_records,
            'dimensions':dimensions,'limits':protocol['limits']}
    retain(OUT/('case-%02d-input.json'%i),prefix)
    # The explicit [2] argument excludes global discriminant factorization.
    nf=pari.nfinit(pari([pari(proxy),pari([2])]))
    basis=[R(a) for a in nf.nf_get_zk()]
    engine=LocalSquareclasses(nf,2)
    theta=pari(x).Mod(pari(proxy));base=pari(b)(theta)
    values=[pari(-1),pari(2),pari(5)]+[(-1)**r['degree']*pari(a)(theta)/base**r['degree'] for a,r in zip(as_poly,a_records)]
    signatures=[list(engine.signature(a)) for a in values]
    M=matrix(GF(2),signatures).transpose();scalar_rank=M[:,:3].rank()
    rank=int(M.rank()-scalar_rank);target=dimensions['Q2_fake_Kummer_dimension'];assert rank<=target
    basis_indices=[];B=M[:,:3];current=B.rank()
    for j in range(17):
        trial=B.augment(M[:,3+j])
        if trial.rank()>current:basis_indices.append(j);B=trial;current+=1
    assert len(basis_indices)==rank
    prime_data=[]
    for prime,uniformizer,bid in engine.data:
        prime_data.append({'e':int(prime[2]),'f':int(prime[3]),
                           'uniformizer':list(map(str,R(uniformizer.lift()).list())),
                           'unit_cyclic_orders':[str(a) for a in bid.bid_get_cyc()]})
    result={**prefix,'status':'PASS_CONSTRUCTED_COMPLETE_DYADIC_IMAGE' if rank==target else 'INCOMPLETE_INHERITED_DYADIC_SPAN',
            'two_maximal_order_basis':[list(map(str,b.list())) for b in basis],
            'prime_data':prime_data,'signatures':signatures,'rational_scalar_image_rank':int(scalar_rank),
            'generic_fake_image_rank':rank,'independent_generic_divisor_indices':basis_indices,
            'generic_true_image_rank_including_D0':rank+dimensions['Q2_true_Kummer_dimension']-target,
            'full_global_Selmer_group':'NOT_COMPUTED','software':{'sage':sage_version,'pari':str(pari.version())},
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [sources[i],prior,OUT/'protocol.json',Path(__file__)]}}
    retain(OUT/('case-%02d.json'%i),result)
    print('case',i,result['status'],'fake rank',rank,'target',target,'prime e,f',[(r['e'],r['f']) for r in prime_data],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--case',type=int,required=True,choices=range(10));args=ap.parse_args()
    signal.alarm(25);construct(args.case)
