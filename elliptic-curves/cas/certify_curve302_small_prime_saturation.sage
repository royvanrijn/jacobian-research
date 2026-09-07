#!/usr/bin/env sage-python
"""Exact displayed-group saturation test at3,5,7 from finite quotients.

Good primes<=3000,120 seconds,one worker; stop when all three ranks reach31.
Failure to reach full rank is UNKNOWN, never nonsaturation evidence.
"""
import argparse
from hashlib import sha256
import itertools
import json
from pathlib import Path
import signal
import sys
from sage.all import EllipticCurve,GF,QQ,matrix,prime_range,prod

ROOT=Path(__file__).resolve().parents[2];CAS=Path(__file__).resolve().parent;sys.path.insert(0,str(CAS))
from icarm_curve302 import GENERAL_WEIERSTRASS_COEFFICIENTS,POINTS
OUT=ROOT/'artifacts/generated-results/elliptic-curves/curve302_small_prime_saturation_v1.json'
REPLAY=OUT.with_name('curve302_small_prime_saturation_replay_v1.json')
ELLS=[3,5,7]


def key(P):return None if P.is_zero() else (int(P[0]),int(P[1]))
def model():
    E0=EllipticCurve(QQ,list(map(QQ,GENERAL_WEIERSTRASS_COEFFICIENTS)));E=E0.short_weierstrass_model();iso=E0.isomorphism_to(E)
    return E,[iso(E0(list(map(QQ,P)))) for P in POINTS]


def coordinates(Ep,gens,orders):
    multiples=[]
    for g,n in zip(gens,orders):
        values=[];P=Ep(0)
        for _ in range(n):values.append(P);P+=g
        assert P.is_zero();multiples.append(values)
    table={}
    for cs in itertools.product(*(range(n) for n in orders)):
        P=sum((multiples[i][c] for i,c in enumerate(cs)),Ep(0));k=key(P);assert k not in table;table[k]=cs
    assert len(table)==int(Ep.cardinality());return table


def reduce_public(Ep,public,p):
    F=GF(p);out=[]
    for P in public:
        if P.is_zero() or QQ(P[0]).denominator()%p==0:out.append(None)
        else:out.append(key(Ep(F(P[0]),F(P[1]))))
    return out


def build():
    E,public=model();rows={ell:[] for ell in ELLS};ranks={ell:0 for ell in ELLS};records=[];examined=[]
    for p0 in prime_range(5,3001):
        p=int(p0);F=GF(p)
        if F(E.discriminant())==0:continue
        Ep=EllipticCurve(F,list(map(F,E.a_invariants())));order=int(Ep.cardinality());examined.append(p)
        wanted=[ell for ell in ELLS if ranks[ell]<31 and order%ell==0]
        if not wanted:continue
        group=Ep.abelian_group();gens=[g.element() for g in group.gens()];orders=[int(g.order()) for g in gens]
        table=coordinates(Ep,gens,orders);labels=[table[k] for k in reduce_public(Ep,public,p)];blocks={}
        for ell in wanted:
            block=[[int(c[i])%ell for c in labels] for i,n in enumerate(orders) if n%ell==0]
            rank=int(matrix(GF(ell),rows[ell]+block).rank())
            if rank>ranks[ell]:rows[ell]+=block;ranks[ell]=rank;blocks[str(ell)]=block
        if blocks:
            records.append({'prime':p,'order':order,'generator_points':[key(g) for g in gens],'generator_orders':orders,'blocks':blocks})
            print('RANKS',p,ranks,flush=True)
            (ROOT/'artifacts/local/elliptic-curves/curve302-wide-packets/saturation_checkpoint.json').write_text(json.dumps({'records':records,'ranks':ranks},sort_keys=True)+'\n')
        if all(rank==31 for rank in ranks.values()):break
    return {'schema':'curve302.small-prime-saturation.v1','status':'PASS_SATURATED_AT_3_5_7' if all(v==31 for v in ranks.values()) else 'INCOMPLETE_FINITE_QUOTIENT_RANKS',
        'input_sha256':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),CAS/'icarm_curve302.py']},
        'limits':{'saturation_primes':ELLS,'good_prime_bound':3000,'seconds':120,'workers':1},'examined_good_primes':examined,'records':records,
        'combined_ranks':{str(k):v for k,v in ranks.items()},'matrix_rows':{str(k):v for k,v in rows.items()},
        'torsion_authority':'ECR31','boundary':'Full finite-quotient rank proves D intersection ell E(Q)=ell D; trivial rational torsion then proves ell-saturation, and induction handles powers. This does not prove saturation at every prime or an exact rank upper bound. A deficient finite matrix is unresolved, not a nonsaturation certificate.'}


def replay():
    d=json.loads(OUT.read_text())
    for p,h in d['input_sha256'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    E,public=model();rows={ell:[] for ell in ELLS}
    for rec in d['records']:
        p=rec['prime'];assert p in prime_range(5,3001);F=GF(p);assert F(E.discriminant())
        Ep=EllipticCurve(F,list(map(F,E.a_invariants())));assert int(Ep.cardinality())==rec['order']
        gens=[Ep(*pt) for pt in rec['generator_points']];orders=rec['generator_orders'];assert int(prod(orders))==rec['order']
        table=coordinates(Ep,gens,orders);labels=[table[k] for k in reduce_public(Ep,public,p)]
        for ell,expected in rec['blocks'].items():
            ell=int(ell);block=[[int(c[i])%ell for c in labels] for i,n in enumerate(orders) if n%ell==0];assert block==expected;rows[ell]+=block
    ranks={str(ell):int(matrix(GF(ell),rows[ell]).rank()) for ell in ELLS}
    assert ranks==d['combined_ranks'] and {str(k):v for k,v in rows.items()}==d['matrix_rows']
    out={'schema':'curve302.small-prime-saturation.replay.v1','status':'PASS_FINITE_QUOTIENT_REPLAY','input_sha256':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),OUT]},'combined_ranks':ranks,
         'boundary':'Reconstructs every recorded finite group from its generators and enumerates all quotient labels; shared group-law implementation. Global saturation additionally uses trivial rational torsion from ECR31.'}
    REPLAY.write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--replay',action='store_true');args=p.parse_args();signal.alarm(120)
    if args.replay:d=replay()
    else:d=build();OUT.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
    print(d['status'],flush=True)
