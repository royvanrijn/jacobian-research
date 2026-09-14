#!/usr/bin/env sage -python
"""Post-evaluation prime decomposition; forbidden as a prospective feature input."""
from collections import Counter
from fractions import Fraction as F
import gzip
import json
from math import prod, log2
from pathlib import Path
import time

from cypari2 import Pari
from sympy import isprime
from finite_cancellation_corpus import OUT, ART, canonical, digest, write

def main():
    start=time.process_time();pari=Pari();counts=Counter();residuals={};case_inputs={}
    for p in sorted((OUT/'cases').glob('*.json.gz')):
        b=p.read_bytes();case_inputs[p.name]=digest(b);d=json.loads(gzip.decompress(b))
        if d['status']!='PASS_EXACT_ACCESSIBILITY':continue
        for row in d['observations']:
            counts[row['g']]+=1;residuals[int(row['unprocessed_g_cofactor'])]=None
    # The recorded largest residual is below2^64. This cheap retrospective
    # decomposition does not authorize factoring a model discriminant.
    assert max(residuals)<2**64
    for n in sorted(residuals):
        fs=pari(n).factor();factors=[(int(fs[i,0]),int(fs[i,1])) for i in range(fs.nrows())]
        assert prod(p**v for p,v in factors)==n and all(isprime(p) and bool(pari(p).isprime()) for p,v in factors)
        residuals[n]=factors
    primes=(2,3,5,7,11,13,17,19,23,29,31);decompositions={};prime_occurrences=Counter()
    for text,count in counts.items():
        g=int(text);rest=g;factors=[]
        for p in primes:
            v=0
            while rest%p==0:rest//=p;v+=1
            if v:factors.append((p,v))
        factors+=residuals[rest]
        assert prod(p**v for p,v in factors)==g
        decompositions[text]={'factors':factors,'occurrences':count}
        for p,v in factors:prime_occurrences[p]+=count
    with gzip.GzipFile(str(OUT/'g_factorizations.json.gz'),'wb',mtime=0) as f:f.write(canonical(decompositions))
    development=[]
    for rank in [27,29,30]:
        folder=ART/'height_model_next_direction_v1'/f'curve302-M{rank}'
        r=json.loads((folder/'report.json').read_text());models=[]
        for row in r['models']:
            g=int(row['actual_cancellation_gcd']);fs=pari(g).factor();factors=[(int(fs[i,0]),int(fs[i,1])) for i in range(fs.nrows())]
            assert prod(p**v for p,v in factors)==g and all(isprime(p) for p,v in factors)
            C=F(row['finite_slack_factor'])*g;S=F(row['actual_real_norm']);H=int(row['parameter_height']);Hx=int(r['x_height'])
            assert F(H**4,Hx)==F(g)/S
            models.append({'model':row['model'],'H':H,'g':str(g),'g_factors':factors,'g_over_S':str(F(g)/S),
                'log2_g':log2(g),'log2_C_over_g':log2(C.numerator)-log2(C.denominator)-log2(g)})
        development.append({'case':f'Curve302-M{rank}','models':models,'source_sha256':digest((folder/'report.json').read_bytes())})
    write(OUT/'prime_decomposition.json',{'status':'PASS_COMPLETE_OBSERVED_G_FACTORIZATION','distinct_g':len(counts),
        'distinct_unprocessed_cofactors':len(residuals),'largest_unprocessed_cofactor':str(max(residuals)),
        'observations':sum(counts.values()),'prime_occurrence_counts':dict(sorted(prime_occurrences.items())),
        'factorizations_sha256':digest((OUT/'g_factorizations.json.gz').read_bytes()),'case_files_sha256':case_inputs,
        'development_curve302_original_models':development,'cpu_seconds':time.process_time()-start,
        'source_sha256':digest(Path(__file__).read_bytes()),
        'boundary':'Retrospective decomposition only, after predictor evaluation. All observed g reconstructed exactly; remaining factors below2^64 checked by PARI and SymPy. No discriminant or uniform-bound factorization; scheduler cannot read this artifact.'})
    print(json.dumps({'status':'PASS_COMPLETE_OBSERVED_G_FACTORIZATION','distinct_g':len(counts),'residuals':len(residuals),'cpu_seconds':time.process_time()-start}))

if __name__=='__main__':main()
