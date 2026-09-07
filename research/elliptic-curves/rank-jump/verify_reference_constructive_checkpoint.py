#!/usr/bin/env python3
"""Independent local-power and matrix checks on the constructive checkpoint."""
import argparse
from pathlib import Path
import retrospective as r
import seeded_reference_class as seed

SOURCE=r.OUT/'rank_jump_reference_targeted_class_extraction_v1.json'
OUTPUT=r.OUT/'rank_jump_reference_constructive_checkpoint_verification_v1.json'

def compute():
    from sage.all import QQ,ZZ,AA,PolynomialRing,pari,GF,matrix
    data=r.read(SOURCE)
    for name,sha in data['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha
    ref=r.read(seed.REFERENCE);R=PolynomialRing(QQ,'z');f=R(ref['cubic_ascending'])
    nf=pari.nfinit([pari(f),ref['S_finite']]);roots=f.roots(AA,multiplicities=False)
    factors={('generic',j):R(g['beta_ascending']) for j,g in enumerate(ref['generic_classes'])}
    for atom in data['used_atoms']:
        factors[('projected_atom',atom['index'])]=ZZ(atom['norm'])*R(atom['alpha_ascending'])
    primes=[P for p in ref['S_finite'] for P in pari.idealprimedec(nf,p)]
    checks=0
    for row in data['strict_candidates']:
        beta=pari.Mod(1,pari(f))
        for label in row['factor_labels']:beta*=pari.Mod(pari(factors[(label['kind'],label['index'])]),pari(f))
        assert ZZ(pari.nfeltnorm(nf,beta)).is_square()
        # Distinct interface from the producer's idealstar/ideallog signatures.
        for P in primes:assert pari.nfislocalpower(nf,P,beta,2)==1;checks+=1
        poly=R([QQ(pari.lift(beta).polcoef(i)) for i in range(3)])
        assert all(poly(x)>0 for x in roots)
    def rows(values):
        values=list(map(int,values));width=max(v.bit_length() for v in values)
        return matrix(GF(2),[[(v>>j)&1 for j in range(width)] for v in values])
    assert rows(data['factor_proof_characters'][:16]).rank()==16
    assert rows([z['proof_character'] for z in data['strict_candidates']]).rank()==6
    assert rows(data['factor_local_signatures'][:16]).rank()==10
    assert data['generic_strict_dimension']==6 and not data['additional_independence_witnesses']
    return {'schema':'rank-jump.reference-constructive-checkpoint-verification.v1','status':'PASS',
        'strict_candidate_products':len(data['strict_candidates']),'individual_local_power_checks':checks,
        'generic_character_rank':16,'generic_local_rank':10,'detected_strict_character_rank':6,
        'additional_class_certified':False,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),SOURCE,seed.REFERENCE]},
        'boundary':'Independent local-power and linear-algebra verification of retained candidate products. No full class-group computation or extra class asserted.'}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['build','check']);a=parser.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS independent local-power checks; positive endpoint remains unmet')
