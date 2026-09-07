#!/usr/bin/env python3
"""Replay the small certified compact-unit control, never the timed-out fields."""
import argparse
from pathlib import Path
import re
import sys
import retrospective as r
import compact_sunit_basis as run

OUTPUT=r.OUT/'rank_jump_compact_sunit_basis_verification_v1.json'


def compute():
    from sage.all import QQ,ZZ,PolynomialRing,pari,GF,matrix
    data=r.read(run.OUTPUT)
    for name,sha in data['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha
    inputs=r.read(run.INPUT)
    for name,sha in inputs['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha
    control=data['rows'][0];assert control['token']=='control' and control['terminal']['status']=='CERTIFIED'
    bnf=pari(control['stages']['bnf']['compact_bnf']);nf=bnf[6];units=pari(control['stages']['units']['compact_units'])
    R=PolynomialRing(QQ,'z');f=R([-11,0,0,1]);assert str(nf.disc())=='-3267' and pari.bnfcertify(bnf)==1
    assert list(map(int,bnf.bnf_get_cyc()))==[2]
    ps=[2,3,11];primes=[P for p in ps for P in pari.idealprimedec(nf,p)]
    assert len(primes)==4 and len(units[0])==6
    sys.path.insert(0,str(r.ROOT/'elliptic-curves/cas'))
    from research_runtime.local_kummer import LocalSquareclasses
    local=[LocalSquareclasses(nf,p) for p in ps]
    rows=[];local_checks=0
    for i,factors in enumerate(units[0]):
        actual=pari.Mod(1,pari(f));local_bits=None;atoms=[]
        for j in range(factors.nrows()):
            base=factors[j,0];exponent=int(factors[j,1])
            a=pari.nfbasistoalg(nf,base) if base.type()=='t_COL' else pari.Mod(base,pari(f))
            norm=pari.nfeltnorm(nf,a);projected=norm*a
            bits=[int(b) for c in local for b in c.signature(projected)]
            if local_bits is None:local_bits=[0]*len(bits)
            for k,b in enumerate(bits):local_bits[k]^=(exponent%2)*b
            actual*=a**exponent
            atoms.append({'coefficients':[str(pari.lift(a).polcoef(k)) for k in range(3)],
                'exponent':exponent,'projected_local_bits':bits})
        # Exact S-unit verification is cheap on this small control only.
        ideal=pari.idealhnf(nf,actual)
        valuations=[]
        for P in primes:
            e=int(pari.idealval(nf,ideal,P));valuations.append(e)
            if e:ideal=pari.idealmul(nf,ideal,pari.idealpow(nf,P,-e))
        assert pari.idealhnf(nf,ideal)==pari.idealhnf(nf,1)
        projected=actual*pari.nfeltnorm(nf,actual)
        direct=[int(b) for c in local for b in c.signature(projected)]
        assert direct==local_bits
        # A separate local-power interface verifies every local zero/nonzero status.
        for chars in local:
            signature=chars.signature(projected)
            squares=[int(pari.nfislocalpower(nf,P,projected,2)) for P in chars.primes]
            assert (not any(signature))==all(squares);local_checks+=len(squares)
        rows.append({'generator':i,'factors':atoms,'S_valuations':valuations,
            'norm_square_projected_local_bits':local_bits})
    rank=matrix(GF(2),[z['norm_square_projected_local_bits'] for z in rows]).rank()
    e_S=1+1-1+4-3
    assert rank==e_S==2
    # Thus the projected units span the norm-square S-unit group and its strict kernel is zero.
    class_coords=[list(map(int,pari.bnfisprincipal(bnf,P,0))) for P in primes]
    s_rank=matrix(GF(2),class_coords).rank();assert s_rank==1
    assert 1-s_rank==0  # The S-class group has trivial two-primary part.
    failures=[]
    for row,source in zip(data['rows'][1:],inputs['cases']):
        assert row['token']==source['token'] and row['terminal']=={'reason':'30-second timeout','status':'UNKNOWN'}
        assert set(row['stages'])=={'setup'}
        assert row['stages']['setup']['field_discriminant']==source['field_discriminant']
        assert r.digest(row['log'].encode())==row['log_sha256']
        requests=[list(map(int,z)) for z in re.findall(r'Look for (\d+) relations in (\d+) ideals',row['log'])]
        assert requests and len(set(x[0] for x in requests))==1
        bases=[list(map(int,z)) for z in re.findall(r'KCZ = (\d+), KC = (\d+), n = (\d+)',row['log'])]
        assert bases
        failures.append({'token':row['token'],'factor_base_initial_columns':bases[0][1],
            'logged_relation_requests':requests,'constant_requested_relation_count':requests[0][0],
            'returned_bnf':False,'returned_S_unit_basis':False,'status':'UNKNOWN'})
    return {'schema':'rank-jump.compact-sunit-basis-verification.v1','status':'PASS',
        'control':{'cubic_ascending':[-11,0,0,1],'S_finite':ps,'ordinary_class_group':[2],
            'ordinary_S_unit_generators':rows,'norm_square_S_unit_dimension':2,'projected_local_rank':2,
            'strict_S_unit_dimension':0,'S_class_two_rank':0,'local_power_checks':local_checks},
        'research_fields':failures,
        'bindings':run.binding([Path(__file__),run.INPUT,run.OUTPUT,r.ROOT/'elliptic-curves/cas/research_runtime/local_kummer.py']),
        'boundary':'Only the small compact-unit control is certified. No full research-field computation is replayed; research class groups, S-units and extra strict classes remain UNKNOWN.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['build','check']);args=parser.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS compact S-unit control; research fields UNKNOWN')
