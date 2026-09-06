#!/usr/bin/env sage-python
"""Full finite factor base and exact principal relations for the fixed512-box."""
import argparse
import json
from math import prod
from pathlib import Path
import time
from sage.all import QQ, pari
import extend_small_conductor_norm_batch as batch
from audit_bnf_free_s_class_quotient import bounds, packed_rank
from research_runtime.store import checkpoint

ROOT, ART, D, cert = batch.ROOT, batch.ART, batch.D, batch.cert
OUTPUT = ART / 'small_conductor_norm_batch_relations_v1.json'


def expected(check=False):
    protocol = batch.protocol()
    frozen = cert.read(D/'audit_protocol.json')
    for name,digest in frozen['sources'].items():
        if cert.hashed(ROOT/name)!=digest:
            raise ArithmeticError('frozen audit input changed')
    summary = cert.read(batch.OUT)
    scalar = cert.read(D/'scalar_replay.json')
    if scalar['status']!='PASS' or scalar['source_sha256']!=cert.hashed(batch.OUT):
        raise ArithmeticError('scalar batch replay required')
    for chunk in summary['chunks']:
        if cert.hashed(ROOT/chunk['path'])!=chunk['sha256']:
            raise ArithmeticError('retained batch chunk changed')
    data = batch.forms.expected()
    if json.loads(json.dumps(data)) != cert.read(batch.forms.OUT):
        raise ArithmeticError('norm identity changed')
    proof = cert.read(batch.forms.target.PROOF)
    polynomial = pari.Pol(list(map(int,data['original_monic_cubic_descending'])))
    primes = [int(p) for p,e in proof['discriminant_factorization']]
    pari.addprimes(primes)
    nf = pari.nfinit([polynomial,primes])
    if len(pari.nfcertify(nf)) or str(nf.disc())!=data['field_discriminant']:
        raise ArithmeticError('certified maximal order differs')
    class FieldBoundInput:
        def degree(self): return 3
        def signature(self): return (3,0)
        def discriminant(self): return int(nf.disc())
    minkowski, bach = map(int,bounds(FieldBoundInput()))
    B = protocol['smooth_bound']
    if B < bach:
        raise ArithmeticError('factor base below Bach upper endpoint')
    columns, blocks, canonical = [], [], []
    started = time.monotonic()
    for number,p in enumerate(batch.prior.primes_to(B)):
        ideals = pari.idealprimedec(nf,p)
        product = pari.idealhnf(nf,1)
        sparse = []
        total_degree = 0
        for ideal in ideals:
            e,f = int(ideal[2]),int(ideal[3])
            total_degree += e*f
            hnf = str(pari.idealhnf(nf,ideal))
            if int(pari.idealnorm(nf,ideal))!=p**f:
                raise ArithmeticError('prime ideal norm differs')
            if e%2: sparse.append(len(columns))
            columns.append({'p':p,'e':e,'f':f,'hnf':hnf})
            product = pari.idealmul(nf,product,pari.idealpow(nf,ideal,e))
        if total_degree!=3 or pari.idealhnf(nf,product)!=pari.idealhnf(nf,p):
            raise ArithmeticError('complete rational-prime factorization differs')
        if not sparse:
            raise ArithmeticError('odd degree requires an odd exponent')
        canonical.append(sum(1<<i for i in sparse))
        blocks.append({'p':p,'parity_columns':sparse})
        if number%2048==0:
            print('FACTOR BASE',number+1,'RATIONAL PRIMES',len(columns),'IDEALS',flush=True)
            if not check:
                checkpoint(D/'audit_progress.json',{'stage':'factor_base','last_prime':p,'ideal_columns':len(columns),'wall_seconds':time.monotonic()-started})
    canonical_rank = packed_rank(canonical)
    if canonical_rank!=len(blocks):
        raise ArithmeticError('disjoint rational-prime rows unexpectedly dependent')
    lookup = {(r['p'],r['hnf']):i for i,r in enumerate(columns)}
    x = pari.Mod('x',polynomial)
    a = int(data['integral_norm_generator']['fixed_a'])
    w = sum(pari(QQ(c))*x**i for i,c in enumerate(data['integral_norm_generator']['w_power_basis']))
    M = data['sl2_matrix']
    relations, actual = [], []
    for index,record in enumerate(summary['smooth_records']):
        m,n,value = record['m'],record['n'],int(record['value'])
        factors = record['factorization']
        if prod(p**e for p,e in factors)!=abs(value) or any(p>B or not pari.isprime(p) for p,e in factors):
            raise ArithmeticError('smooth factors differ')
        beta = a*(M[0][0]*m+M[0][1]*n)+(M[1][0]*m+M[1][1]*n)*w
        if pari.nfeltnorm(nf,beta)!=a*a*value:
            raise ArithmeticError('exact norm transport differs')
        factor = pari.idealfactor(nf,beta)
        if pari.idealhnf(nf,pari.idealfactorback(nf,factor))!=pari.idealhnf(nf,beta):
            raise ArithmeticError('principal ideal product differs')
        row,mask = [],0
        for j in range(factor.nrows()):
            ideal,e = factor[j,0],int(factor[j,1])
            if e<0:
                raise ArithmeticError('nonintegral norm generator')
            column = lookup[(int(ideal[0]),str(pari.idealhnf(nf,ideal)))]
            row.append([column,e])
            mask ^= (e%2)<<column
        actual.append(mask)
        relations.append({'m':m,'n':n,'beta_power_basis':[str(pari.polcoef(pari.lift(beta),i)) for i in range(3)],'ideal_factorization':row})
        if index%128==0:
            print('EXACT RELATIONS',index+1,flush=True)
    snapshots = []
    for box in [64,128,256,512]:
        selected = [mask for r,mask in zip(relations,actual) if abs(r['m'])<=box and r['n']<=box]
        gain = packed_rank(canonical+selected)-canonical_rank
        snapshots.append({'box':box,'smooth_relations':len(selected),'additional_relation_rank':gain,
                          'quotient_dimension':len(columns)-canonical_rank-gain})
    if snapshots[0]['smooth_relations']!=18 or snapshots[0]['additional_relation_rank']!=18:
        raise ArithmeticError('prior18 relation regression differs')
    dimension = snapshots[-1]['quotient_dimension']
    rank_bound = dimension+7
    rank_bound -= rank_bound%2
    sources = [Path(__file__).resolve(),batch.OUT,D/'scalar_replay.json',D/'audit_protocol.json',batch.forms.OUT,
               ROOT/'elliptic-curves/cas/audit_bnf_free_s_class_quotient.py']
    result = {'schema':'elliptic-curves.small-conductor-norm-batch-relations.v1','status':'PASS',
              'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in sources},
              'pari_version':str(pari.version()),'field_discriminant':str(nf.disc()),
              'factor_base_rational_bound':B,'bach_bound_upper_endpoint':bach,
              'minkowski_bound_upper_endpoint':str(minkowski),
              'generation_hypothesis':'GRH sufficient for Bach ideal-class generation; see Klagsbrun-Sherman-Weigandt section4.2. All primes of ideal norm<=B occur because all ideals above every rational prime<=B are included.',
              'reference':'https://arxiv.org/html/1606.07178#S4.SS2',
              'columns':columns,'canonical_rational_relations':blocks,
              'canonical_rank':canonical_rank,'relations':relations,'box_snapshots':snapshots,
              'conditional_on_grh_class_two_rank_upper_bound':dimension,
              'conditional_on_grh_curve_rank_upper_bound':rank_bound,
              'additional_independent_rows_needed_for_target':max(0,dimension-16),
              'unconditional_generation_certified':False,'unconditional_rank_upper_bound':None,
              'exact_rank':None,
              'claim_boundary':'Exact finite supported relation matrix. Its upper bound for Cl(K)/2 and hence curve rank uses GRH generation; no unconditional upper bound or extra point is proved.'}
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=expected(args.check)
    if args.check:
        if result!=cert.read(OUTPUT):
            raise ArithmeticError('factor-base relation replay differs')
    else:
        if OUTPUT.exists():
            raise FileExistsError('preserve relation certificate')
        cert.write(OUTPUT,result)
    print('MATRIX AUDIT PASS',len(result['columns']),'COLUMNS',result['box_snapshots'],flush=True)
