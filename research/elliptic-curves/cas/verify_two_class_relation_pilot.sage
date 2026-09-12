#!/usr/bin/env sage-python
"""Replay frozen candidate coverage and ideal identities without class groups."""
import argparse
from collections import Counter
from hashlib import sha256
from math import gcd, prod
from pathlib import Path

from sage.all import GF, QQ, ZZ, PolynomialRing, matrix, pari, prime_range
import run_two_class_relation_pilot as pilot
import two_class_relation_core as core


def verify_field(out, row, plan):
    key = row['curve_key']
    dest = out/key
    inp = pilot.read(out/'inputs'/f'{key}.json')
    supervisor = pilot.read(dest/'supervisor.json')
    if supervisor['status'] != 'PASS_BOUNDED_WORKER':
        return {'curve_key':key,'status':supervisor['status'],'global_class_2rank_estimate':None}
    setup, result = pilot.read(dest/'setup.json'), pilot.read(dest/'result.json')
    R = PolynomialRing(QQ, 'x')
    f = pari(R(list(map(QQ, inp['integral_monic_cubic_ascending']))))
    hints = list(map(ZZ, inp['certified_bad_rational_primes']))
    assert all(p.is_prime(proof=True) for p in hints)
    nf = pari.nfinit([f, sorted(set([ZZ(2)]+hints))])
    assert pari.nfcertify(nf) == []
    assert str(nf.disc()) == inp['field_discriminant']
    assert list(map(int, nf.nf_get_sign())) == inp['field_signature']
    assert [str(v) for v in nf.nf_get_zk()] == setup['form']['integral_basis_polynomials']
    beta = nf.nfbasistoalg(pari(list(map(QQ, setup['form']['generator_in_integral_basis']))).Col())
    c = setup['form']['reduced_cubic_ascending']
    assert sum(c[i]*beta**i for i in range(4)) == 0
    assert str(R(c).discriminant()) == inp['field_discriminant']
    assert all(QQ(v).denominator()==1 for v in nf.nfalgtobasis(c[3]*beta))
    denominator_primes = [ZZ(p) for p,e in setup['denominator_factorization']]
    assert all(p.is_prime(proof=True) for p in denominator_primes)
    assert prod(ZZ(p)**e for p,e in setup['denominator_factorization']) == c[3]
    support = sorted(set(map(int, prime_range(inp['protocol']['factor_base_bound']+1))) | set(map(int,denominator_primes)))
    assert support == setup['rational_support']
    fb, ideals = setup['factor_base'], []
    for p in support:
        for P in sorted(nf.idealprimedec(p),key=lambda P:str(nf.idealhnf(P))):
            record=fb[len(ideals)]
            assert record['column']==len(ideals) and record['rational_prime']==p
            assert (record['ramification_index'],record['residue_degree'])==(int(P[2]),int(P[3]))
            assert record['norm']==str(ZZ(p)**int(P[3]))
            assert record['hnf']==[[str(v) for v in r] for r in nf.idealhnf(P).sage().rows()]
            ideals.append(P)
    assert len(ideals)==len(fb)
    budget=inp['protocol']['candidate_budget']
    relations=pilot.read(dest/f'relations-{budget:06}.json')['rows']
    canonical=setup['canonical_count']
    assert canonical==len(support) and len(relations)==canonical+result['noncanonical_relation_count']
    masks=[]
    sparse_entries={}
    smooth_by_attempt={}
    for j, relation in enumerate(relations):
        gamma=nf.nfbasistoalg(pari(list(map(QQ,relation['generator_in_integral_basis']))).Col())
        if j<canonical:
            p=support[j]
            assert relation['kind']=='canonical_rational_prime' and relation['rational_prime']==p and gamma==p
            assert relation['signed_ideal_exponents']==[[i, fb[i]['ramification_index']] for i in range(len(fb)) if fb[i]['rational_prime']==p]
        else:
            assert relation['kind']=='smooth_binary_cubic'
            assert gamma==relation['a']+relation['b']*beta
            assert relation['norm']==str(QQ(nf.nfeltnorm(gamma)))
            assert relation['attempt'] not in smooth_by_attempt
            smooth_by_attempt[relation['attempt']]=relation
        ideal=nf.idealhnf(1)
        exponents=relation['signed_ideal_exponents']
        assert len({i for i,e in exponents})==len(exponents)
        mask=0
        for i,e in exponents:
            assert e and int(nf.idealval(gamma,ideals[i]))==e
            ideal=nf.idealmul(ideal,nf.idealpow(ideals[i],e))
            if e%2:
                mask ^= 1 << i
                sparse_entries[j,i]=1
        assert nf.idealhnf(ideal)==nf.idealhnf(gamma)
        assert hex(mask)==relation['parity_hex']
        masks.append(mask)
    # Independent Sage GF(2) matrix implementation, not packed-row insertion.
    M=matrix(GF(2),len(relations),len(fb),sparse_entries,sparse=False)
    assert M.rank()==result['final_matrix']['rank_mod2']
    assert len(fb)-M.rank()==result['final_matrix']['deficiency']
    assert result['final_matrix']==core.matrix_stats(masks,len(fb),canonical)
    trial_hash=sha256()
    norm_bits,residual_bits=Counter(),Counter()
    smooth_attempts=[]
    product=prod(support)
    for attempted,(a,b) in enumerate(core.primitive_pairs(budget),1):
        numerator=sum(c[i]*a**i*(-b)**(3-i) for i in range(4))
        remainder=abs(numerator)
        while True:
            common=gcd(remainder,product)
            if common==1:
                break
            remainder//=common
        trial_hash.update(f'{a},{b},{numerator},{remainder}\n'.encode())
        norm_bits[str(abs(numerator).bit_length())]+=1
        residual_bits[str(remainder.bit_length())]+=1
        if remainder==1:
            smooth_attempts.append(attempted)
            r=smooth_by_attempt[attempted]
            assert (r['a'],r['b'])==(a,b)
            assert QQ(r['norm'])==QQ(numerator)/c[3]
        if attempted in inp['protocol']['candidate_prefixes']:
            saved=pilot.read(dest/f'prefix-{attempted:06}.json')
            assert saved['trial_sha256']==trial_hash.hexdigest()
            assert saved['attempted']==attempted and saved['smooth_relations']==len(smooth_attempts)
            count=canonical+len(smooth_attempts)
            assert pilot.read(dest/f'relations-{attempted:06}.json')['rows']==relations[:count]
            assert all(saved[k]==v for k,v in core.matrix_stats(masks[:count],len(fb),canonical).items())
    assert set(smooth_attempts)==set(smooth_by_attempt)
    assert trial_hash.hexdigest()==result['trial_sha256']
    assert dict(norm_bits)==result['norm_numerator_bit_histogram']
    assert dict(residual_bits)==result['residual_bit_histogram']
    assert result['global_class_2rank_estimate'] is None and result['global_class_2rank_upper_bound'] is None
    bindings={p.name:pilot.sha(p) for p in sorted(dest.iterdir()) if p.is_file()}
    return {'curve_key':key,'t':row['t'],'status':'PASS_EXACT_BOUNDED_RELATION_REPLAY',
            'cohort':row['cohort'],'selection_mode':row['selection_mode'],
            'candidates_replayed':budget,'ideal_relations_replayed':len(relations),
            'noncanonical_relations':len(smooth_attempts),'matrix':result['final_matrix'],
            'bindings':bindings,'global_class_2rank_estimate':None}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=pilot.OUT)
    args=ap.parse_args()
    out=args.output.resolve()
    plan=pilot.check(out,complete=True)
    rows=[]
    for row in plan['rows']:
        got=verify_field(out,row,plan)
        rows.append(got)
        print(f"MOD2_REPLAY|t={row['t']}|{got['status']}",flush=True)
    pilot.save(out/'REPLAY.json',{'schema':'elliptic-curves.two-class-pilot-replay.v1',
        'status':'PASS_BOUNDED_EVIDENCE_REPLAY','rows':rows,
        'checker_sha256':pilot.sha(Path(__file__)),'plan_sha256':pilot.sha(out/'plan.json'),
        'boundary':'Exact candidate replay, principal-ideal identities and independent GF(2) matrix rank. Shared Sage/PARI; no class-group generation or full g result.'})


if __name__=='__main__':
    main()
