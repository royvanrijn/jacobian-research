#!/usr/bin/env python3
"""Modular irreducibility and the canonical discriminant-cover gate."""
import argparse
from pathlib import Path
import retrospective as r
import parameter_cover_capacity as source
import bad_prime_support as bad
from retrospective_secant_pencils import rational_square

PROTOCOL=Path(__file__).with_name('PARAMETER_COVER_VERIFICATION_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_parameter_cover_capacity_verification_v1.json'


def possible_degrees(degrees):
    sums={0}
    for d in degrees:sums|={s+d for s in list(sums)}
    return sums-{0,sum(degrees)}


def compute():
    from sage.all import QQ,GF,PolynomialRing
    inp=r.read(source.INPUT);data=r.read(source.OUTPUT)
    for record in [inp['source_hashes'],data['bindings']]:
        for name,digest in record.items():assert r.digest((r.ROOT/name).read_bytes())==digest
    R=PolynomialRing(QQ,'u');u=R.gen()
    A=R(list(map(QQ,inp['A'])));B=R(list(map(QQ,inp['B'])))
    D=-64*A**3-432*B**2
    assert list(map(str,D))==data['discriminant_coefficients_ascending']
    candidates=set(range(1,24));modular=[]
    for p in r.primes(r.read(PROTOCOL)['limits']['maximum_modular_prime']):
        fp=D.change_ring(GF(p))
        if fp.degree()!=24 or fp.gcd(fp.derivative())!=1:
            modular.append({'prime':p,'status':'SKIP_BAD_REDUCTION'})
            continue
        factors=fp.factor();degrees=[int(f.degree()) for f,e in factors for _ in range(e)]
        assert all(e==1 for f,e in factors)
        candidates &= possible_degrees(degrees)
        modular.append({'prime':p,'factor_degrees':degrees,'remaining_proper_factor_degrees':sorted(candidates)})
        if not candidates:break
    assert data['smallest_singular_closed_point_degree']==24
    for q,cert in zip(inp['covers'],data['retained_covers'],strict=True):
        d=R(list(map(QQ,q['q'])));b=int(d.degree())
        assert d.gcd(d.derivative())==1 and d.gcd(D)==1
        # Local Euler numbers: 24 nodal fibres plus b I0* fibres.
        euler=24+6*b;chi=euler//12;roots=4*b
        assert euler%12==0
        assert cert['twist_chi']==chi and cert['twist_root_rank']==roots
        assert cert['twist_geometric_MW_rank_upper_bound']==10*chi-2-roots==18+b
    params=[x['parameter'] for x in r.read(source.branch_blocks.INPUT)['fibres']]
    params += [x['parameter'] for x in bad.cases()[4:6]]
    assert len(params)==len(set(params))==7
    rows=[]
    for parameter in params:
        value=D(QQ(parameter));assert value
        square=bool(value.is_square());assert square==rational_square(str(value))
        rows.append({'parameter':parameter,'discriminant_value':str(value),'square':square,
            'canonical_discriminant_cover_has_rational_lift':square})
    return {'schema':'rank-jump.parameter-cover-capacity-verification.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,source.INPUT,source.OUTPUT,source.branch_blocks.INPUT,r.INPUT)},
        'modular_irreducibility_certificate':modular,
        'independent_irreducibility_status':'PASS' if not candidates else 'UNKNOWN',
        'remaining_proper_factor_degrees':sorted(candidates),
        'cover_counts':{str(b):sum(x['branch_points']==b for x in data['retained_covers']) for b in [2,4]},
        'retained_discriminant_cover_tests':rows,
        'discriminant_twist':source.capacity(0,24,24),
        'boundary':'No branch point of a rational divisor of degree <24 can be a singular R17 fibre. Capacity is an upper bound; a failed discriminant-cover lift does not exclude other covers or high rank.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();data=compute()
    if a.mode=='check':assert r.read(OUTPUT)==data;print('PASS independent parameter-cover verification')
    else:r.write_new(OUTPUT,data)
    print(data['modular_irreducibility_certificate'])
    print([(x['parameter'],x['square']) for x in data['retained_discriminant_cover_tests']])
